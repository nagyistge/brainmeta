#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import pandas
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
#basedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/mr"
outdirectory = "%s/similarity_scores" %(basedir)
indirectory = "%s/mr" %(basedir)
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

# Input file
input_file = "%s/openfmri_labels.tsv" %(basedir)
input_delim = "\t"

# Read in input file
inputs = pandas.read_csv(input_file,sep=input_delim)
# IMAGE_ID should correspond to integer ID

# Prepare and submit a job for each
for i in inputs["ID"]:
  image_id = i
  single_metrics = "%s/000%s.pkl" %(outdirectory,i)
  output_metrics = "%s/000%s.tsv" %(outdirectory,i)
  if not os.path.isfile(output_file):
    filey = ".job/" + uuids[i] + ".job"
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=" + uuids[i] + "\n")
    filey.writelines("#SBATCH --output=.out/" + uuids[i] + ".out\n")
    filey.writelines("#SBATCH --error=.out/" + uuids[i] + ".err\n")
    filey.writelines("#SBATCH --time=1-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    # Usage : experiment1_analysis1.py image_id base_path output_file standard input_file input_delim
    filey.writelines("/home/vsochat/python-lapack-blas/bin/python /home/vsochat/SCRIPT/python/brainmeta/experiment1/experiment1_analysis1.py %s %s %s" %(image_id,indirectory,output_metrics,single_metrics,standard,input_file,input_delim))
    filey.close()
    os.system("sbatch " + ".job/" + image_id + ".job")
