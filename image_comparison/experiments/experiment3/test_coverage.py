#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 3 Handling Missing Values: Assessing Ranking and Distribtions of Scores as we Increase Coverage (using Random Sampling)

# Algorithm

# We will sample a certain percentage of the regions after the imputation / masking strategy.

# For each unthresholded brain map, Xi
#  For a level of coverage (represented by a random sampling percentage)
#      and calculate pearson r,spearman r, for each sample vs. Xi
 
import os
import pandas
import sys
import pickle
import random
import numpy as np
import nibabel as nib
from glob import glob
from scipy.stats import pearsonr, spearmanr
import image_transformations as IT
from nilearn.masking import apply_mask

percent_sample = float(sys.argv[1]) # eg, 0.2
input_file = sys.argv[2]       # mr image input file
outfile = sys.argv[3]          # output pickle file
standard = sys.argv[4]         # standard image

# Read in our brain map paths
inputs = pandas.read_csv(input_file,sep="\t")
# mrs = [nib.load(i) for i in inputs.files]

# Unzip functiom
unzip = lambda l:tuple(zip(*l))

# Get our standard brain mask
brain_mask = nib.load(standard)

# We will save a distribution of correlations for each image
# each of these will be a list of lists (distributions for 5000 permutations)
pd_pearsons = pandas.DataFrame(columns=inputs.uid)
pi_pearsons = pandas.DataFrame(columns=inputs.uid)
bm_pearsons = pandas.DataFrame(columns=inputs.uid)
pd_spearmans = pandas.DataFrame(columns=inputs.uid)
pi_spearmans = pandas.DataFrame(columns=inputs.uid)
bm_spearmans = pandas.DataFrame(columns=inputs.uid)

# And a distribution of the sizes for each
pd_sizes = pandas.DataFrame(columns=inputs.uid)
pi_sizes = pandas.DataFrame(columns=inputs.uid)
bm_sizes = pandas.DataFrame(columns=inputs.uid)

files = inputs.files.tolist()

# For each unthresholded brain map, Xi
for row1 in inputs.iterrows():
  filename1 = row1[1].files
  uid1 = row1[1].uid
  mr1 = nib.load(filename1)
  pd_pearson_single = []
  pi_pearson_single = []
  bm_pearson_single = []
  pd_spearman_single = []
  pi_spearman_single = []
  bm_spearman_single = []
  # We also need to save the distributions of total sizes
  roi_sizes = []
  bm_size_single = []
  pd_size_single = []
  pi_size_single = [] 
  image_ids = []
  # For each other image, randomly sample the percentage of voxels
  for row2 in inputs.iterrows():
    filename2 = row2[1].files
    uid2 = row2[1].uid
    # Don't compare images to themselves
    if uid1 != uid2:
      image_ids.append(uid2)
      print "Calculating for %s vs %s" %(uid1,uid2)
      mr2 = nib.load(filename2)
      # Prepare pd and pi masks
      pdmask = IT.get_pairwise_deletion_mask(mr1,mr2,brain_mask)
      pimask = IT.get_pairwise_inclusion_mask(mr1,mr2,brain_mask)      
      # Make random voxel mask
      empty_nii = np.zeros(brain_mask.shape)
      x,y,z = np.where(brain_mask.get_data()==1)
      idx = zip(x,y,z)
      np.random.shuffle(idx) # mix it up! shake it up!
      if percent_sample != 1.0: # if threshold is 1, we take all voxels 
        number_voxels = int(np.floor(percent_sample * len(idx)))
        random_idx = idx[0:number_voxels]
        random_idx = unzip(random_idx)
        empty_nii[random_idx] = 1
      else: empty_nii[brain_mask.get_data()==1] = 1
      # For each of brain mask, pdmask, and pimask, combine with random selection
      pdmask_sub = np.logical_and(pdmask, empty_nii).astype(int)
      pimask_sub = np.logical_and(pimask, empty_nii).astype(int)
      brain_mask_sub = np.logical_and(brain_mask, empty_nii).astype(int)  
      # How many total voxels?
      pi_size_single.append(len(np.where(pimask_sub[pimask_sub==1])[0]))
      pd_size_single.append(len(np.where(pdmask_sub[pdmask_sub==1])[0]))
      bm_size_single.append(len(np.where(brain_mask_sub[brain_mask_sub==1])[0]))
      # Make into nifti images
      pdmask_sub = nib.nifti1.Nifti1Image(pdmask_sub,affine=brain_mask.get_affine(),header=brain_mask.get_header())
      pimask_sub = nib.nifti1.Nifti1Image(pimask_sub,affine=brain_mask.get_affine(),header=brain_mask.get_header())
      brain_mask_sub = nib.nifti1.Nifti1Image(brain_mask_sub,affine=brain_mask.get_affine(),header=brain_mask.get_header())
      # PAIRWISE DELETION (INTERSECTION)
      # Calculate correlation if there is overlap, otherwise it is nan
      if len(np.unique(pdmask_sub.get_data())) == 2:
        datapd = apply_mask([mr1,mr2],pdmask_sub) 
        # We need at least 3 values
        if np.shape(datapd)[1] > 2: 
          pd_spearman_single.append(spearmanr(datapd[0],datapd[1])[0])
          pd_pearson_single.append(spearmanr(datapd[0],datapd[1])[0])
        else: 
          pd_pearson_single.append(np.nan)
          pd_spearman_single.append(np.nan)
      else: 
        pd_pearson_single.append(np.nan)
        pd_spearman_single.append(np.nan)
      # PAIRWISE INCLUSION (UNION)
      # Calculate correlation if there is overlap, otherwise it is 0
      if len(np.unique(pimask_sub.get_data())) == 2:
        datapi = apply_mask([mr1,mr2],pimask_sub)  
        # We need at least 3 values
        if np.shape(datapi)[1] > 2: 
          pi_spearman_single.append(spearmanr(datapi[0],datapi[1])[0])
          pi_pearson_single.append(pearsonr(datapi[0],datapi[1])[0])
        else: 
          pd_pearson_single.append(np.nan)
          pd_spearman_single.append(np.nan)
      else: 
        pd_pearson_single.append(np.nan)
        pd_spearman_single.append(np.nan)

      # BRAIN MASK
      # There will always be overlap, but why not check anyway
      if len(np.unique(brain_mask_sub.get_data())) == 2:    
        databm = apply_mask([mr1,mr2],brain_mask_sub)  
        # We will have three, but why not check and be consistent :)
        if np.shape(databm)[1] > 2: 
          bm_pearson_single.append(pearsonr(databm[0],databm[1])[0])
          bm_spearman_single.append(spearmanr(databm[0],databm[1])[0])
        else: 
          bm_pearson_single.append(np.nan)
          bm_spearman_single.append(np.nan)
      else: 
        bm_pearson_single.append(np.nan)
        bm_spearman_single.append(np.nan)


  # Add each to our data frames
  pd_pearsons.loc[uid1,image_ids] = pd_pearson_single
  pi_pearsons.loc[uid1,image_ids] = pi_pearson_single
  bm_pearsons.loc[uid1,image_ids] = bm_pearson_single
  pd_spearmans.loc[uid1,image_ids] = pd_spearman_single
  pi_spearmans.loc[uid1,image_ids] = pi_spearman_single
  bm_spearmans.loc[uid1,image_ids] = bm_spearman_single
  pd_sizes.loc[uid1,image_ids] = pd_size_single
  pi_sizes.loc[uid1,image_ids] = pi_size_single
  bm_sizes.loc[uid1,image_ids] = bm_size_single

# Save all data to output file
output = {"ids":inputs.uid.tolist(),"pd_pearsons":pd_pearsons,"pi_pearsons":pi_pearsons,
          "bm_pearsons":bm_pearsons,"pd_spearmans":pd_spearmans,"pi_spearmans":pi_spearmans,
          "bm_spearmans":bm_spearmans,"pd_sizes":pd_sizes,"pi_sizes":pi_sizes,"bm_sizes":bm_sizes,
          "percent_sample":percent_sample}
pickle.dump(output,open(outfile,"wb"))
