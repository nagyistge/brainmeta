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

# Prepare and submit a job for each
for i in inputs.ID:
  image_id = i
  time.sleep(1)
  outfile = "%s/000%s_masking_scores.pkl" %(outdirectory,image_id)
  if not os.path.exists(outfile):
    filey = ".job/%s_masking.job" %(image_id)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
    filey.writelines("#SBATCH --output=.out/%s.out\n" %(image_id))
    filey.writelines("#SBATCH --error=.out/%s.err\n" %(image_id))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment2_masking.py %s" %(image_id)
    filey.close()
    os.system("sbatch -p russpold " + ".job/%s_masking.job" %(image_id))
