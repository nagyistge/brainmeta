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

# Read in input file
inputs = pandas.read_csv(input_file,sep=input_delim)
# IMAGE_ID should correspond to integer ID

empty = [115,116,117,118,119,120,121,122,123,124,125,126,127,129,130,131,132,142,143,144,145,146,147,148,149,150,151,152,153,154,155,173,174,175,177,178,179,180,181,299,300,303,305,306,307,308,309,311,436,437,439,440,441,442,443,446,449,451,452,454,456,457,459,460,461,463,464,466,472,526,528,529,530,531,532,533,534,535,536,538,539,540,542,543,544]

# Prepare and submit a job for each
for i in empty:
  image_id = i
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
  filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment1_analysis1.py %s %s %s %s %s %s" %(image_id,indirectory,tmpdirectory,output_directory,standard,input_file))
  filey.close()
  os.system("sbatch -p russpold " + ".job/%s.job" %(image_id))
