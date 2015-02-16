#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 2 Masking:

# Here I want to do a simple experiment to show how pearson correlations change with 
# three masking strategies:

# brain mask [everything in brain mask]
# pairwise inclusion [voxels in either mask]
# pairwise deletion [voxels in both mask] 

# 1: Define 144 unthresholded Z score images
# 2: Threshold at 1.96 (equivalent to p=0.05 level)
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

basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
outdirectory = "%s/masking_scores" %(basedir)
indirectory = "%s/mr" %(basedir)
tmpdirectory = "%s/tmp" %(basedir)
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)
brain_mask = nib.load(standard)

# Input file
input_file = "%s/openfmri_labels.tsv" %(basedir)
input_delim = "\t"

# Here is the threshold for p=0.05
threshold = 1.96

# Read in input file
inputs = pandas.read_csv(input_file,sep=input_delim)

# Read in all images - these similarities will be gold standard
image_path = "%s/000%s.nii.gz" %(indirectory,image_id)
mr1 = nib.load(image_path)
mrs = [nib.load("%s/000%s.nii.gz" %(indirectory,ii)) for ii in inputs.ID]

# We also want the maps thresholded
thresholded = []
for mr in mrs:
  thresholded.append(IT.threshold_abs(mr,thresholds=[threshold])[threshold])
  
# Calculate "gold standard" list of pearson scores
# These are unthresholded maps vs. unthresholded maps
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

# We also will save sizes of each
sizes = pandas.DataFrame(columns=["pd","pi","bm"])
idx = 0
print "Calculating mask varieties [PD,PI,BM] vs thresholded..."
for mr2 in thresholded:
  pdmask = IT.get_pairwise_deletion_mask(mr1,mr2,brain_mask)
  pimask = IT.get_pairwise_inclusion_mask(mr1,mr2,brain_mask)
  datapd = apply_mask([mr1,mr2],pdmask)  
  datapi = apply_mask([mr1,mr2],pimask)  
  databm = apply_mask([mr1,mr2],brain_mask)  
  pearsons_pd.append(pearsonr(datapd[0],datapd[1])[0])
  pearsons_pi.append(pearsonr(datapi[0],datapi[1])[0])
  pearsons_bm.append(pearsonr(databm[0],databm[1])[0])
  sizes.loc[idx] = [len(datapd[0]),len(datapi[0]),len(databm[0])]
  idx+=1

# Save all data to output dictionary
output = {"ids":inputs.ID.tolist(),"pearson_gs":pearsons,"mr_vs_thresh_pearson_pd":pearsons_pd,
          "mr_vs_thresh_pearson_pi":pearsons_pi,"mr_vs_thresh_pearson_bm":pearsons_bm,"sizes":sizes}

pickle.dump(output,open("%s/000%s_masking_scores.pkl" %(outdirectory,image_id),"wb"))
