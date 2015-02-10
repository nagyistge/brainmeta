#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 1 Analysis:
# Individual images compared to transformations of themselves.

# INPUT: a single brainmap image
# OUTPUT: similarity metrics for the image across different transformations (compared to original)

# 1. generate multiple transformations
# 2. calculate similarity metrics
# 3. save results


# This is the worker script to run for a single author (on SLURM cluster)
# Usage : authorSynth_cluster.py uuid "author" email outdirectory
# This will look up papers on pubmed, and cross list with
# neurosynth database. You should run this with run_authorSynth_cluster.py

import sys
import pandas
import numpy as np
import nibabel as nib
import similarity_metrics as SM
import image_transformations as IT

# Get arguments
image_id = sys.argv[1]
indirectory = sys.argv[2]
output_file = sys.argv[3]
standard_mask = sys.argv[4]
input_file = sys.argv[5]
input_delim = sys.argv[6]

print "Processing image %s" %(image_path)

# Load other image paths
inputs = pandas.read_csv(inputfile,sep=input_delim)
image_path = "%s/000%s.nii.gz" %(indirectory,image_id)
original = nib.load(image_path)
mask = nib.load(standard_mask)

# Produce thresholdings (these are all Z score maps, this includes original image (thresh 0))
thresholded1 = IT.threshold_abs(original)
thresholds1 = np.sort(thresholded1.keys())
image1_labels = ["%s_thr_%s" %(image_id,thresh) for thresh in thresholds1]

# We will have a matrix of image threshold combinations (rows) by similarity metrics (columns)
similarity_metrics = pandas.DataFrame()

# Extract a column (list of similarity metrics) for each image vs original (index 0 == original)
for t in range(0,len(thresholds1)):
  thresh = thresholds1[t]
  image1 = thresholded1[thresh]
  label1 = image1_labels[t]
    
  # Do a comparison for each pairwise set at each threshold
  for i in inputs["IMAGE_ID"]:

    image2_path = "%s/000%s.nii.gz" %(indirectory,i)
    image2 = nib.load(image2_path)
  
    # Only proceed if image dimensions are equal, and in same space
    if ((image1.shape == image2.shape) and np.all(image1.get_affine() == image2.get_affine())):
      thresholded2 = IT.threshold_abs(image2)
      thresholds2 = np.sort(thresholded2.keys())
      image2_labels = ["%s_thr_%s" %(i,th) for th in thresholds2]

      for tt in range(0,len(thresholds2)):
        thresh2 = thresholds2[t]
        image2 = thresholded2[thresh2]
        label2 = image2_labels[tt]
        single_metrics,pairwise_metrics = SM.run_all(image1=image1,image2=image2,image1_label=label1,
                                            image2_label=image2_label,brain_mask=mask)    
        # Here we need to format into a data frame to save, and add single_metrics MUST TEST THIS!
        similarity_metrics["%s_%s" %(label1,label2)] = pairwise_metrics
        

    else:
      print "ERROR: mask %s and image %s are not the same shape! Exiting." %(standard_mask,image_path)

  # Save the similarity metrics to file
  similarity_metrics.to_csv(output_file,sep="\t")

