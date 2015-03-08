#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 2 Masking: Assessing Variance in Scores As we Increase Coverage (using Random Sampling)

# Algorithm

# We will sample a certain percentage of the regions after the imputation / masking strategy.

# For each unthresholded brain map, Xi
#  For a level of coverage (represented by a random sampling percentage)
#      and calculate pearson r for each sample vs. Xi
 
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
gs_file = sys.argv[5]          # gold standard file (will be created if doesn't exist)

# Read in our brain map paths
inputs = pandas.read_csv(input_file,sep="\t")
mrs = [nib.load(i) for i in inputs.FILE]

# Get our standard brain mask
brain_mask = nib.load(standard)

# We will save a distribution of correlations for each image
# each of these will be a list of lists (distributions for 5000 permutations)
pd_correlations = pandas.DataFrame(columns=inputs.ID)
pi_correlations = pandas.DataFrame(columns=inputs.ID)
bm_correlations = pandas.DataFrame(columns=inputs.ID)

# And a distribution of the sizes for each
pd_sizes = pandas.DataFrame(columns=inputs.ID)
pi_sizes = pandas.DataFrame(columns=inputs.ID)
bm_sizes = pandas.DataFrame(columns=inputs.ID)

# For the gold standard correlation, we only get one value for each image
# Only do this if the file does not exist yet
if not os.path.exists(gs_file):
  # The gold standard comparison data frame will be calculated once
  gs_correlations = pandas.DataFrame(columns=inputs.ID) # each image has a list, length N images (144)
  for m in range(0,len(mrs)):
    gs_single = np.zeros(len(mrs))
    mr1 = mrs[m]  
    print "Calculating %s of %s" %(m,len(mrs))
    for mm in range(0,len(mrs)):
      mr2 = mrs[mm]
      data = apply_mask([mr1,mr2],brain_mask) 
      gs_single[mm] = spearmanr(data[0],data[1])[0]
    gs_correlations.loc[m] = gs_single
  gs_correlations.index = inputs.ID
  gs_correlations.to_csv(gs_file,sep="\t")


# For each unthresholded brain map, Xi
for m in range(0,len(mrs)):
  # Here is the current image ID
  image1_id = inputs.ID[m]
  mr1 = mrs[m] # This is the unthresholded image
  pd_single = []
  pi_single = []
  bm_single = []
  # We also need to save the distributions of total sizes
  roi_sizes = []
  bm_size_single = []
  pd_size_single = []
  pi_size_single = [] 
  image_ids = []
  # For each other image, randomly sample the percentage of voxels
  for mm in range(0,len(mrs)):
    image2_id = inputs.ID[mm]
    # Don't compare images to themselves
    if image1_id != image2_id:
      image_ids.append(image2_id)
      print "Calculating for %s vs %s" %(image1_id,image2_id)
      mr2 = mrs[mm] # this is an unthresholded image
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
        random_idx = self.unzip(random_idx)
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

      # PAIRWISE DELETION
      # Calculate correlation if there is overlap, otherwise it is 0
      if len(np.unique(pdmask_sub.get_data())) == 2:
        datapd = apply_mask([mr1,mr2],pdmask_sub) 
        # We need at least 3 values
        if np.shape(datapd)[1] > 2: pd_single.append(spearmanr(datapd[0],datapd[1])[0])
        else: pd_single.append(0)
      else: pd_single.append(0)
      # PAIRWISE INCLUSION
      # Calculate correlation if there is overlap, otherwise it is 0
      if len(np.unique(pimask_sub.get_data())) == 2:
        datapi = apply_mask([mr1,mr2],pimask_sub)  
        # We need at least 3 values
        if np.shape(datapi)[1] > 2: pi_single.append(spearmanr(datapi[0],datapi[1])[0])
        else: pi_single.append(0)
      else: pi_single.append(0)
      databm = apply_mask([mr1,mr2],brain_mask_sub)  
      bm_single.append(spearmanr(databm[0],databm[1])[0])
  # Add each to our data frames
  pd_correlations.loc[image1_id,image_ids] = pd_single
  pi_correlations.loc[image1_id,image_ids] = pi_single
  bm_correlations.loc[image1_id,image_ids] = bm_single
  pd_sizes.loc[image1_id,image_ids] = pd_size_single
  pi_sizes.loc[image1_id,image_ids] = pi_size_single
  bm_sizes.loc[image1_id,image_ids] = bm_size_single

# Save all data to output file
output = {"ids":inputs.ID.tolist(),"pd_correlations":pd_correlations,
          "pi_correlations":pi_correlations,"bm_correlations":bm_correlations,
          "pd_sizes":pd_sizes,"pi_sizes":pi_sizes,"bm_sizes":bm_sizes,"num_regions":num_regions}
pickle.dump(output,open(outfile,"wb"))
