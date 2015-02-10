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
import numpy as np
import nibabel as nib
import similarity_metrics as SM
import image_transformations as IT

# Get arguments
image_id = sys.argv[1]
image_path = sys.argv[2]
output_file = sys.argv[3]
standard_mask = sys.argv[4]

print "Processing image %s" %(image_path)

original = nib.load(image_path)
mask = nib.load(standard_mask)

# Only proceed if image dimensions are equal, and in same space
if ((mask.shape == standard.shape) and np.all(mask.get_affine() == original.get_affine())):

  # Produce thresholdings (these are all Z score maps, this includes original image (thresh 0))
  thresholded = IT.threshold_abs(original)
  thresholds = np.sort(thresholded.keys())
  image_labels = ["%s_thresh_%s" %(image_id,thresh) for thresh in thresholds]

  # We will have a matrix of similarity scores (rows) by transformations (columns) compared to gold standard
  similarity_metrics = pandas.DataFrame()

  # Extract a column (list of similarity metrics) for each image vs original (index 0 == original)
  for t in range(0,len(thresholds)):
    thresh = thresholds[t]
    image = thresholded[thresh]
    label = image_labels[t]
    similarity_metrics[thresh] = SM.run_all(image1=image,image2=original,image1_label=label,
                                            image2_label=image_id,brain_mask=mask)    

  # Save the similarity metrics to file
  similarity_metrics.to_csv(output_file,sep="\t")

else:
  print "ERROR: mask %s and image %s are not the same shape! Exiting." %(standard_mask,image_path)
