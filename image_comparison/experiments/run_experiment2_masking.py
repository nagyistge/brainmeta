#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import time
import pandas

# Input file
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
outdirectory = "%s/masking_scores" %(basedir)
input_file = "%s/openfmri_labels.tsv" %(basedir)
input_delim = "\t"
inputs = pandas.read_csv(input_file,sep=input_delim)

# We will run for a set of thresholds
#thresholds = [1,1.65,1.96,2.58,3.02]
thresholds = [0.5,1.0,1.5,1.65,1.7,1.75,1.8,1.85,1.9,1.96,2.0,2.58,3.02,3.5,4.0]

# Do not limit to only positive values
absolute_value = False

# Prepare and submit a job for each
for thresh in thresholds:
  for image_id in inputs.ID:
    # if the output directory doesn't exist, make it
    topdir = "%s/thresh_%s_%s" %(outdirectory,thresh,absolute_value)
    if not os.path.exists(topdir): os.mkdir(topdir)
    outfile = "%s/000%s_masking_scores_thresh_%s.pkl" %(topdir,image_id,thresh)
    if not os.path.exists(outfile):
      filey = ".job/%s_masking_%s.job" %(image_id,thresh)
      filey = open(filey,"w")
      filey.writelines("#!/bin/bash\n")
      filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
      filey.writelines("#SBATCH --output=.out/%s_%s.out\n" %(image_id,thresh))
      filey.writelines("#SBATCH --error=.out/%s_%s.err\n" %(image_id,thresh))
      filey.writelines("#SBATCH --time=2-00:00\n")
      filey.writelines("#SBATCH --mem=64000\n")
      filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment2_masking1.py %s %s %s" %(image_id,thresh,absolute_value))
      filey.close()
      os.system("sbatch -p russpold " + ".job/%s_masking_%s.job" %(image_id,thresh))
