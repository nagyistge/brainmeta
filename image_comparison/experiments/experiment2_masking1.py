#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 2 Masking:

# Here I want to do a simple experiment to show how pearson correlations change with 
# three masking strategies:

# brain mask [everything in brain mask]
# pairwise inclusion [voxels in either mask]
# pairwise deletion [voxels in both mask] 

# 1: Define 144 unthresholded Z score images
# 2: Threshold at 1.96, 2.58, 3.02 (equivalent to p=0.05, 0.005, 0.001 levels)
# 3: Calculate pearsonr similarity for each masking strategy
# 4: Plot differences in scores comparing strategies

import os
import time
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
import similarity_metrics as SM
import image_transformations as IT
from scipy.stats import pearsonr
from nilearn.masking import apply_mask

image_id = sys.argv[1]
threshold = float(sys.argv[2])
absolute_value = sys.argv[3]

# Make sure we have boolean
if absolute_value == "True": absolute_value = True
else: absolute_value = False

basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
outdirectory = "%s/masking_scores" %(basedir)
indirectory = "%s/mr" %(basedir)
tmpdirectory = "%s/tmp" %(basedir)
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)
brain_mask = nib.load(standard)

# Input file
input_file = "%s/openfmri_labels.tsv" %(basedir)
input_delim = "\t"

# Read in input file
inputs = pandas.read_csv(input_file,sep=input_delim)

# Remove our query image
inputs = inputs[inputs.ID!=int(image_id)]
image_ids = inputs.ID.tolist()

# Read in all images - these similarities will be gold standard
image_path = "%s/000%s.nii.gz" %(indirectory,image_id)
mr1 = nib.load(image_path)
mrs = [nib.load("%s/000%s.nii.gz" %(indirectory,ii)) for ii in inputs.ID]
  
# We also want the maps thresholded
thresholded = []
for mr in mrs:
  if absolute_value: thresholded.append(IT.threshold_abs(mr,thresholds=[threshold])[threshold])
  else: thresholded.append(IT.threshold_pos(mr,thresholds=[threshold])[threshold])
  
# Calculate "gold standard" list of pearson scores
# These are unthresholded maps vs. unthresholded maps
# No absolute value is taken!
pearsons_gs = []
print "Calculating gold standard for %s..." %(mr1)
for mr2 in mrs:
  pdmask = IT.get_pairwise_deletion_mask(mr1,mr2,brain_mask)
  data = apply_mask([mr1,mr2],pdmask)  
  pearsons_gs.append(pearsonr(data[0],data[1])[0])

# Now calculate same image vs thresholded maps:

# PAIRWISE DELETION, PAIRWISE INCLUSION, BRAINMASK
pearsons_pd = []
pearsons_pi = []
pearsons_bm = []

# Edge cases:
# No surviving voxels for pdmask: we append 0
# Thresholded image is 0: all correlations 0
# Fewer than 3 remaining voxels to compare: we append 0

# We also will save sizes of each
sizes = pandas.DataFrame(columns=["pd","pi","bm"])
print "Calculating mask varieties [PD,PI,BM] vs thresholded..."
idx = 0
size_ids = []
for mr2 in thresholded:
  # If the image is empty thresholded, we must append zeros by default
  if len(np.unique(mr2.get_data()))==1:
    pearsons_bm.append(0)
    pearsons_pi.append(0)
    pearsons_pd.append(0)
  else:
    pdmask = IT.get_pairwise_deletion_mask(mr1,mr2,brain_mask)
    pimask = IT.get_pairwise_inclusion_mask(mr1,mr2,brain_mask)
    # Calculate correlation if there is overlap, otherwise it is 0
    if len(np.unique(pdmask.get_data())) == 2:
      datapd = apply_mask([mr1,mr2],pdmask) 
      # We need at least 3 values
      if np.shape(datapd)[1] > 2: pearsons_pd.append(pearsonr(datapd[0],datapd[1])[0])
      else: pearsons_pd.append(0)
    else: pearsons_pd.append(0)
    # Calculate correlation if there is overlap, otherwise it is 0
    if len(np.unique(pimask.get_data())) == 2:
      datapi = apply_mask([mr1,mr2],pimask)  
      # We need at least 3 values
      if np.shape(datapi)[1] > 2: pearsons_pi.append(pearsonr(datapi[0],datapi[1])[0])
      else: pearsons_pi.append(0)
    else: pearsons_pi.append(0)
    databm = apply_mask([mr1,mr2],brain_mask)  
    pearsons_bm.append(pearsonr(databm[0],databm[1])[0])
    sizes.loc[idx] = [len(datapd[0]),len(datapi[0]),len(databm[0])]
    size_ids.append(image_ids[idx])
  idx+=1
# Save all data to output dictionary
output = {"ids":inputs.ID.tolist(),"pearson_gs":pearsons_gs,"mr_vs_thresh_pearson_pd":pearsons_pd,
          "mr_vs_thresh_pearson_pi":pearsons_pi,"mr_vs_thresh_pearson_bm":pearsons_bm,
          "sizes":sizes,"size_ids":size_ids}

pickle.dump(output,open("%s/thresh_%s_%s/000%s_masking_scores_thresh_%s.pkl" %(outdirectory,threshold,absolute_value,image_id,threshold),"wb"))
