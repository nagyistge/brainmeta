#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import time
import pandas
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
#basedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON"
outdirectory = "%s/similarity_scores" %(basedir)
indirectory = "%s/mr" %(basedir)
tmpdirectory = "%s/tmp" %(basedir)
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

# Input file
input_file = "%s/openfmri_labels.tsv" %(basedir)
input_delim = "\t"

# Here is the threshold to run
thresholds = 4.0

# Read in input file
inputs = pandas.read_csv(input_file,sep=input_delim)
# IMAGE_ID should correspond to integer ID

# Prepare and submit a job for each
for image_id in missing:
  time.sleep(1)
  output_directory = "%s/%s" %(outdirectory,image_id)
  if not os.path.exists(output_directory):
    os.mkdir(output_directory)
  filey = ".job/%s.job" %(image_id)
  filey = open(filey,"w")
  filey.writelines("#!/bin/bash\n")
  filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
  filey.writelines("#SBATCH --output=.out/%s.out\n" %(image_id))
  filey.writelines("#SBATCH --error=.out/%s.err\n" %(image_id))
  filey.writelines("#SBATCH --time=2-00:00\n")
  filey.writelines("#SBATCH --mem=64000\n")
  filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment1_analysis1.py %s %s %s %s %s %s %s" %(image_id,indirectory,tmpdirectory,output_directory,standard,input_file,thresholds))
  filey.close()
  os.system("sbatch -p russpold " + ".job/%s.job" %(image_id))
