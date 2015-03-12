#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 3 How to handle missing data:

# Here I want to do a simple experiment to show how pearson and spearman correlations change with 
# three strategies for handling missing data.

# brain mask [everything in brain mask]
# union (pairwise inclusion) [voxels in either mask]
# intersection (pairwise deletion) [voxels in both mask] 

# 1: Define 144 unthresholded Z score images
# 2: Threshold at various levels (eg 1.96, 2.58, 3.02 equivalent to p=0.05, 0.005, 0.001)
# 3: Calculate pearsonr, spearmanr, similarity for each masking strategy

import os
import time
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
import similarity_metrics as SM
import image_transformations as IT
from scipy.stats import pearsonr, spearmanr
from nilearn.masking import apply_mask

basedir = sys.argv[1]
image_id = sys.argv[2]
filename1 = sys.argv[3]
threshold = float(sys.argv[4])
direction = sys.argv[5]
output_file = sys.argv[6]

# Make sure we have boolean
if direction == "posneg": absolute_value = True
else: absolute_value = False

# standard brain mask that will be used for "brain mask" and to eliminate out of brain voxels for all masks
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)
brain_mask = nib.load(standard)

# Input file
input_file = "%s/doc/hcp_groupmaps.tsv" %(basedir)
inputs = pandas.read_csv(input_file,sep="\t")

# Remove our query image
#inputs = inputs[inputs.ID!=int(image_id)]
image_ids = inputs.uid.tolist()
filenames = inputs.files.tolist()

# Read in main image
mr1 = nib.load(filename1)
  
# We will threshold maps as we go, because it takes up too much memory to store 1720 maps
# 3/12/2015: confirmed that thresholding at 0.0 is equivalent to original image
  
# PAIRWISE DELETION, PAIRWISE INCLUSION, BRAINMASK, for each of spearman and pearson
pearsons_pd = []
pearsons_pi = []
pearsons_bm = []
spearman_pd = []
spearman_pi = []
spearman_bm = []

# Edge cases:
# No surviving voxels for pdmask: we append 0
# Thresholded image is 0: all correlations 0
# Fewer than 3 remaining voxels to compare: we append 0

# We also will save sizes of each - this can confirm if the comparison was possible, period
# (to distinguish values of 0 == no correlation vs 0 == no comparison possible)
sizes = pandas.DataFrame(columns=["pd","pi","bm"])

print "Calculating mask varieties [PD,PI,BM] vs thresholded..."
idx = 0
size_ids = []

for filename2 in filenames:
  # Load the image and calculate the thresholded version
  mr2 = nib.load(filename2)
  if absolute_value: mrthresh = IT.threshold_abs(mr2,thresholds=[threshold])[threshold]
  else: mrthresh = IT.threshold_pos(mr2,thresholds=[threshold])[threshold]
  # 3/12/2015: confirmed that first returns +/- values, second returns only positive  

  # If the image is empty thresholded, we must append zeros by default
  if len(np.unique(mrthresh.get_data()))==1:
    pearsons_bm.append(0)
    pearsons_pi.append(0)
    pearsons_pd.append(0)
    spearman_pd.append(0)
    spearman_pd.append(0)
    spearman_pd.append(0)
  else:
    # Generate a union (pi) and intersection (pd) mask
    pdmask = IT.get_pairwise_deletion_mask(mr1,mrthresh,brain_mask)
    pimask = IT.get_pairwise_inclusion_mask(mr1,mrthresh,brain_mask)

    # PAIRWISE DELETION (intersection)
    # Calculate correlation if there is overlap
    if len(np.unique(pdmask.get_data())) == 2:
      datapd = apply_mask([mr1,mrthresh],pdmask) 
      
      # We need at least 3 values
      if np.shape(datapd)[1] > 2: 
        pearsons_pd.append(pearsonr(datapd[0],datapd[1])[0])
        spearman_pd.append(spearmanr(datapd[0],datapd[1])[0])
      else: 
        pearsons_pd.append(0)
        spearman_pd.append(0)
    
    # Otherwise (no overlap) it is zero
    else: 
      pearsons_pd.append(0)
      spearman_pd.append(0)
    
    # PAIRWISE INCLUSION (union)
    # Calculate correlation if there is overlap
    if len(np.unique(pimask.get_data())) == 2:
      datapi = apply_mask([mr1,mrthresh],pdmask) 
      
      # We need at least 3 values
      if np.shape(datapi)[1] > 2: 
        pearsons_pi.append(pearsonr(datapi[0],datapi[1])[0])
        spearman_pi.append(spearmanr(datapi[0],datapi[1])[0])
      else: 
        pearsons_pi.append(0)
        spearman_pi.append(0)
    
    # Otherwise (no overlap) it is zero
    else: 
      pearsons_pi.append(0)
      spearman_pi.append(0)

    # BRAIN MASK
    databm = apply_mask([mr1,mrthresh],brain_mask)  
    # We need at least 3 values
    if np.shape(databm)[1] > 2: 
      pearsons_bm.append(pearsonr(databm[0],databm[1])[0])
      spearman_bm.append(spearmanr(databm[0],databm[1])[0])
    else: 
      pearsons_bm.append(0)
      spearman_bm.append(0)    

    # Save sizes of all masks
    sizes.loc[idx] = [len(datapd[0]),len(datapi[0]),len(databm[0])]
    size_ids.append(image_ids[idx])
  idx+=1

# Save all data to output dictionary
output = {"uid":inputs.uid.tolist(),
          "mr_vs_thresh_pearson_pd":pearsons_pd,
          "mr_vs_thresh_pearson_pi":pearsons_pi,
          "mr_vs_thresh_pearson_bm":pearsons_bm,
          "mr_vs_thresh_spearman_pd":spearman_pd,
          "mr_vs_thresh_spearman_pi":spearman_pi,
          "mr_vs_thresh_spearman_bm":spearman_bm,
          "sizes":sizes,"size_ids":size_ids}

pickle.dump(output,open(outfile,"wb"))
