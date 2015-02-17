#/usr/bin/python Experiment 1:

# This script will check for missing runs, and create a list of images and thresholds to resubmit!

import os
import pandas
import re
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
outdirectory = "%s/similarity_scores" %(basedir)
from glob import glob

# Input file
input_file = "%s/openfmri_labels.tsv" %(basedir)
input_delim = "\t"
inputs = pandas.read_csv(input_file,sep=input_delim)

# We are going to do these one at a time!
# thresholds = [0.0,0.5,1.0,1.5,1.65,1.7,1.75,1.8,1.85,1.9,1.96,2,2.58,3,3.5,4.0]
thresholds = [4.0]

# Generate a set of all that we should have
should_have = set([ "%s_%s" %(thresh,ii) for thresh in thresholds for ii in inputs.ID])

missing = list()
for i in inputs.ID:
  subdir = "%s/%s/" %(outdirectory,i)
  files = glob("%s*.tsv" %(subdir))
  files = set([x.replace("%s000%s_thr_" %(subdir,i),"").replace("_pairwise_metrics.tsv","") for x in files])
  tmp = list(should_have.difference(files))
  print "%s is missing %s." %(i,len(tmp))
  if len(tmp) > 0:
    missing.append(i)
