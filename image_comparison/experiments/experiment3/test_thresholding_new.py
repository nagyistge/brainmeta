#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 3 How to handle missing data:

# Here I want to do a simple experiment to show how pearson and spearman correlations change with 
# three strategies for handling missing data.

# complete case analysis [intersection, voxels in both mask]
# single-value imputation [union, voxels in either mask]
# multiple-imputation [considered gold standard]

# 1: Define 470 unthresholded Z score images
# 2: Threshold at various levels (eg 1.96, 2.58, 3.02 equivalent to p=0.05, 0.005, 0.001)
# 3: Compare unthresholded to thresholded
# 4: Calculate pearsonr, spearmanr, similarity for each masking strategy

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
groups_pkl = sys.argv[7]

# Make sure we have boolean
if direction == "posneg": absolute_value = True
else: absolute_value = False

# standard brain mask that will be used for "brain mask" and to eliminate out of brain voxels for all masks
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)
brain_mask = nib.load(standard)

# Input file
input_file = "%s/doc/hcp_groupmaps_filter.tsv" %(basedir)
inputs = pandas.read_csv(input_file,sep="\t")

# Groups file tells us the dof for each group
groups = pickle.load(open(groups_pkl,"rb"))

# Remove our query image
#inputs = inputs[inputs.ID!=int(image_id)]
image_ids = inputs.uid.tolist()
filenames = inputs.files.tolist()

# Read in main image
mr1 = nib.load(filename1)

# Convert to Z score map
group1 = inputs.groups[inputs.files==filename1].tolist()[0]
df = len(groups[group1]) - 2
mr1 = IT.t_to_z(mr1,df)
  
# We will threshold maps as we go, because it takes up too much memory to store 1720 maps
# 3/12/2015: confirmed that thresholding at 0.0 is equivalent to original image
  
# Complete case analysis, single value imputation, multiple-imputation
cca_pearson = []
svi_pearson = []
cca_spearman = []
svi_spearman = []

# Edge cases:
# No surviving voxels: we append nan
# Thresholded image is 0: all correlations nan
# Fewer than 3 remaining voxels to compare: we append nan
# Suggested by Russ 3/13/2015, and I agree

# We also will save sizes of each, and a log to keep track of nan behavior
sizes = pandas.DataFrame(columns=["cca","svi"])
nanlog_cca = []
nanlog_svi = []

print "Calculating mask varieties [CCA,SVI,MI] vs thresholded..."
idx = 0
size_ids = []

for filename2 in filenames:
  # Load the image, convert to Z, and calculate the thresholded version
  mr2 = nib.load(filename2)
  group2 = inputs.groups[inputs.files==filename2].tolist()[0]
  dof = len(groups[group2]) - 2
  mr2 = IT.t_to_z(mr2,dof)
  if absolute_value: mrthresh = IT.threshold_abs(mr2,thresholds=[threshold])[threshold]
  else: mrthresh = IT.threshold_pos(mr2,thresholds=[threshold])[threshold]
  # 3/12/2015: confirmed that first returns +/- values, second returns only positive  

  # If the image is empty thresholded, we append NaN
  if len(np.unique(mrthresh.get_data()))==1:
    cca_pearson.append(np.nan)
    svi_pearson.append(np.nan)
    cca_spearman.append(np.nan)
    svi_spearman.append(np.nan)
    sizes.loc[idx] = [0,0]
    nanlog_cca.append("nan_mrthresh_empty")
    nanlog_svi.append("nan_mrthresh_empty")
    size_ids.append(image_ids[idx])
  else:
    # Generate a union (svi) and intersection (cca) mask
    ccamask = IT.get_pairwise_deletion_mask(mr1,mrthresh,brain_mask)
    svimask = IT.get_pairwise_inclusion_mask(mr1,mrthresh,brain_mask)

    # COMPLETE CASE ANALYSIS (OLD PAIRWISE DELETION, intersection)
    # Calculate correlation if there is overlap
    if len(np.unique(ccamask.get_data())) == 2:
      datacca = apply_mask([mr1,mrthresh],ccamask) 
      
      # We need at least 3 values
      if np.shape(datacca)[1] > 2: 
        cca_pearson.append(pearsonr(datacca[0],datacca[1])[0])
        cca_spearman.append(spearmanr(datacca[0],datacca[1])[0])
        nanlog_cca.append("success")
      else: 
        cca_pearson.append(np.nan)
        cca_spearman.append(np.nan)
        nanlog_cca.append("nan_fewer_3_values")
    
    # Otherwise (no overlap) it is nan
    else: 
      cca_pearson.append(np.nan)
      cca_spearman.append(np.nan)
      nanlog_cca.append("nan_no_overlap")
    
    # SINGLE VALUE IMPUTATION (old pairwise inclusion, union)
    # Calculate correlation if there is overlap
    if len(np.unique(svimask.get_data())) == 2:
      datasvi = apply_mask([mr1,mrthresh],svimask) 
      
      # We need at least 3 values
      if np.shape(datasvi)[1] > 2: 
        svi_pearson.append(pearsonr(datasvi[0],datasvi[1])[0])
        svi_spearman.append(spearmanr(datasvi[0],datasvi[1])[0])
        nanlog_svi.append("success")
      else: 
        svi_pearson.append(np.nan)
        svi_spearman.append(np.nan)
        nanlog_svi.append("nan_fewer_3_values")
    
    # Otherwise (no overlap) it is nan
    else: 
      svi_pearson.append(np.nan)
      svi_spearman.append(np.nan)
      nanlog_svi.append("nan_no_overlap")

    # Save sizes of all masks
    sizes.loc[idx] = [len(datacca[0]),len(datasvi[0])]
    size_ids.append(image_ids[idx])
  idx+=1

# Save all data to output dictionary
output = {"uid":inputs.uid.tolist(),
          "cca_pearson":cca_pearson,
          "svi_pearson":svi_pearson,
          "cca_spearman":cca_spearman,
          "svi_spearman":svi_spearman,
          "sizes":sizes,"size_ids":size_ids,
          "nanlog_cca":nanlog_cca,
          "nanlog_svi":nanlog_svi,
          "mr_df":df}

pickle.dump(output,open(output_file,"wb"))
