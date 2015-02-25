#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 2 Masking: Assessing Variance in Scores As we Increase Coverage (using ROIs)

# Algorithm

# We will sample N ROIs at a time. I.e., if we have 150 regions, randomly sample 5 each time (and then you can systematically vary the number of regions, which is a proxy for increasing the amount of coverage).

# We will get a distribution of r's for each of N images, for each set of ROIs, for each approach (PD vs. PI).

# For each unthresholded brain map, Xi
#  For a level of coverage (represented by a number of rois, N)
#    Introduce noise to the process [which means] for each other map in X:
#      randomly sample N of the regions
#      and calculate pearson r for each sample vs. Xi
 
import os
import pandas
import sys
import pickle
import random
import numpy as np
import nibabel as nib
from glob import glob
from scipy.stats import pearsonr
import image_transformations as IT
from nilearn.masking import apply_mask

# Atlas to use for regions
atlas_nifti = sys.argv[1]
num_regions = int(sys.argv[2]) # eg, 2
atlas_file = sys.argv[3]       # file with list of atlas regions
input_file = sys.argv[4]       # mr image input file
outfile = sys.argv[5]          # output pickle file
standard = sys.argv[6]         # standard image
gs_file = sys.argv[7]          # gold standard file (will be created if doesn't exist)

# Read in the atlas image and file
atlas = nib.load(atlas_nifti)
atlas_labels = pandas.read_csv(atlas_file,sep="\t",header=None)

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
      gs_single[mm] = pearsonr(data[0],data[1])[0]
    gs_correlations.loc[m] = gs_single
  gs_correlations.index = inputs.ID
  gs_correlations.to_csv(gs_file,sep="\t")


# For each unthresholded brain map, Xi
#  For a level of coverage (represented by a number of rois, N)
#    Introduce noise to the process [which means] for each other map in X:
#      randomly sample N of the regions
#      and calculate pearson r for each sample vs. Xi


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
  # For each other image, randomly sample N regions
  for mm in range(0,len(mrs)):
    image2_id = inputs.ID[mm]
    # Don't compare images to themselves
    if image1_id != image2_id:
      image_ids.append(image2_id)
      print "Calculating for %s vs %s" %(image1_id,image2_id)
      mr2 = mrs[mm] # this is an unthresholded image
      # Select a random set of N masks
      roi_labels = random.sample(xrange(len(atlas_labels[2])), num_regions)
      roi_labels = atlas_labels[2][roi_labels].tolist()
      # Create binary mask
      roi = np.zeros(atlas.shape)
      for roi_label in roi_labels:
        roi[atlas.get_data() == int(roi_label)] = 1
      # How many total voxels?
      region_size = len(roi[roi==1])
      roi_sizes.append(region_size)
      # Mask the image with the region
      masked_mr = IT.get_masked_images(mr2,roi)[0]
      pdmask = IT.get_pairwise_deletion_mask(mr1,masked_mr,brain_mask)
      pimask = IT.get_pairwise_inclusion_mask(mr1,masked_mr,brain_mask)      
      # Save the sizes of the masks
      pi_size_single.append(len(np.where(pimask.get_data()==1)[0]))
      pd_size_single.append(len(np.where(pdmask.get_data()==1)[0]))
      bm_size_single.append(len(np.where(brain_mask.get_data()==1)[0])) # always the same
      # PAIRWISE DELETION
      # Calculate correlation if there is overlap, otherwise it is 0
      if len(np.unique(pdmask.get_data())) == 2:
        datapd = apply_mask([mr1,masked_mr],pdmask) 
        # We need at least 3 values
        if np.shape(datapd)[1] > 2: pd_single.append(pearsonr(datapd[0],datapd[1])[0])
        else: pd_single.append(0)
      else: pd_single.append(0)
      # PAIRWISE INCLUSION
      # Calculate correlation if there is overlap, otherwise it is 0
      if len(np.unique(pimask.get_data())) == 2:
        datapi = apply_mask([mr1,masked_mr],pimask)  
        # We need at least 3 values
        if np.shape(datapi)[1] > 2: pi_single.append(pearsonr(datapi[0],datapi[1])[0])
        else: pi_single.append(0)
      else: pi_single.append(0)
      databm = apply_mask([mr1,masked_mr],brain_mask)  
      bm_single.append(pearsonr(databm[0],databm[1])[0])
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
          "pd_sizes":pd_sizes,"pi_sizes":pi_sizes,"bm_sizes":bm_sizes,"num_regions",num_regions}
pickle.dump(output,open(outfile,"wb"))
