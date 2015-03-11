#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 1 Scoring:
# Scoring similarity metrics based on "gold standard" whole brain (unthresholded) metric

# 1. read in and format data
# 2. calculate gold standard [ranked list of order of similarities for unthresholded maps, pairwise deletion]
# 3. generate new rankings for all transformations for each metric
# 4. calculate scores for each metric, save to file

# IN PROGRESS

import sys
import os
import pandas
import pickle
import numpy as np
from glob import glob
from scoring_functions import run_scoring

# Read in inputs file
input_file = ""
inputs = pandas.read_csv(input_file,sep="\t")
input_folder = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis",
output_folder = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/compiled_scores",
metrics = ["covariance","correlation_coefficient","correlation_ratio","correlation_ratio_norm",
           "mutual_information_norm","mutual_information","cosine","euclidean","minkowski",
           "cityblock","seuclidean","sqeuclidean","chebyshev","canberra","braycurtis",
           "skl_linear_kernel","skl_l1","skl_l2","skl_sigmoid_kernel","skl_polynomial_kernel",
           "skl_rbf_kernel"]

for i in inputs.ID:
  run_scoring(image_id=i,input_folder=input_folder,output_folder=output_folder,
              num_images=144,num_thresh=16,metrics=metrics)


