#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
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

# Prepare and submit a job for ach
for i in inputs["ID"]:
  image_id = i
  output_directory = "%s/%s" %(outdirectory,image_id)
  if not os.path.exists(output_directory):
    os.mkdir(output_directory)
  if not os.path.isfile(output_metrics):
    filey = ".job/%s.job" %(image_id)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
    filey.writelines("#SBATCH --output=.out/%s.out\n" %(image_id))
    filey.writelines("#SBATCH --error=.out/%s.err\n" %(image_id))
    filey.writelines("#SBATCH --time=1-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    # Usage : experiment1_analysis1.py image_id base_path output_metrics, single_metrics,standard input_file input_delim
    filey.writelines("/home/vsochat/python-lapack-blas/bin/python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment1_analysis1.py %s %s %s %s %s %s" %(image_id,indirectory,tmpdirectory,output_directory,standard,input_file))
    filey.close()
    os.system("sbatch " + ".job/%s.job" %(image_id))
