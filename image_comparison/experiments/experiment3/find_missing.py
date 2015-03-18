#!/usr/bin/env python

import os
import pickle
import pandas
# This script will find missing files

# Input file with map paths
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3"
outdirectory = "%s/scores" %(basedir)
input_file = "%s/doc/hcp_groupmaps.tsv" %(basedir)
inputs = pandas.read_csv(input_file,sep="\t")
direction = "posneg"
groups = "%s/doc/hcp_10groups460_alltasks.pkl" %(basedir)

# We will run for a set of thresholds
thresholds = [0.0,0.5,1.0,1.5,1.65,1.7,1.75,1.8,1.85,1.9,1.96,2.0,2.58,3.02,3.5,4.0]

# Find missing result files
missing = []
missing_outfile = []
missing_thresh = []
for thresh in thresholds:
  for groupmap in inputs.iterrows():
    image_id = groupmap[1].uid
    filename = groupmap[1].files
    topdir = "%s/thresh_%s_%s" %(outdirectory,thresh,direction)
    outfile = "%s/%s_scores_thresh_%s.pkl" %(topdir,image_id,thresh)
    if not os.path.exists(outfile):
      missing.append(groupmap)
      missing_thresh.append(thresh)
      missing_outfile.append(outfile)
      print "%s is missing!" %(outfile)


      

