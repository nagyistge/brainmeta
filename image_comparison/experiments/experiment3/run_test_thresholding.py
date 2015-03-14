#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import time
import pandas
import pickle

# Input file with map paths
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3"
outdirectory = "%s/scores" %(basedir)
input_file = "%s/doc/hcp_groupmaps.tsv" %(basedir)
inputs = pandas.read_csv(input_file,sep="\t")

# Our groups will tells us the number of degrees of freedom for each map
groups = "%s/doc/hcp_10groups460_alltasks.pkl" %(basedir)

# We will run for a set of thresholds
thresholds = [0.0,0.5,1.0,1.5,1.65,1.7,1.75,1.8,1.85,1.9,1.96,2.0,2.58,3.02,3.5,4.0]
thresholds = [0.0,0.5,1.0,1.5,1.65]

# Take absolute value (include negative values): True
# Do not take absolute value (include only positive values): False
absolute_value = True
if absolute_value == True:
  direction = "posneg"
else:
  direction = "pos"

# Prepare and submit a job for each
for thresh in thresholds:
  for groupmap in inputs.iterrows():
    image_id = groupmap[1].uid
    filename = groupmap[1].files
    # if the output directory doesn't exist, make it
    topdir = "%s/thresh_%s_%s" %(outdirectory,thresh,direction)
    if not os.path.exists(topdir): os.mkdir(topdir)
    outfile = "%s/%s_scores_thresh_%s.pkl" %(topdir,image_id,thresh)
    if not os.path.exists(outfile):
      filey = ".job/%s_masking_%s.job" %(image_id,thresh)
      filey = open(filey,"w")
      filey.writelines("#!/bin/bash\n")
      filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
      filey.writelines("#SBATCH --output=.out/%s_%s.out\n" %(image_id,thresh))
      filey.writelines("#SBATCH --error=.out/%s_%s.err\n" %(image_id,thresh))
      filey.writelines("#SBATCH --time=2-00:00\n")
      filey.writelines("#SBATCH --mem=64000\n")
      filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/test_thresholding.py %s %s %s %s %s %s %s" %(basedir,image_id,filename,thresh,direction,outfile,groups))
      filey.close()
      os.system("sbatch -p russpold " + ".job/%s_masking_%s.job" %(image_id,thresh))
