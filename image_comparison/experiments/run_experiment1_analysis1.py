#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
#TODO: decide if should run on TACC or Sherlock!

import os
import pandas
#basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
basedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/mr"
outdirectory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1/similarity_scores"
indirectory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1/data/mrs"
standard = "/scratch/users/vsochat/DATA/BRAINMETA/data/MNI152_T1_2mm_brain_mask.nii.gz"

# Read in input file
inputs = pandas.read_csv("%s/openfmri_labels.tsv" %(basedir),sep="\t")

# Prepare and submit a job for each
for i in inputs["IMAGE_ID"]:
  mr_file = "%s/000%s.nii.gz" %(indirectory,i)
  image_id = i
  output_file = "%s/000%s.pkl" %(outdirectory,i)
  if not os.path.isfile(output_file):
    filey = ".job/" + uuids[i] + ".job"
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=" + uuids[i] + "\n")
    filey.writelines("#SBATCH --output=.out/" + uuids[i] + ".out\n")
    filey.writelines("#SBATCH --error=.out/" + uuids[i] + ".err\n")
    filey.writelines("#SBATCH --time=1-00:00\n")
    filey.writelines("#SBATCH --mem=12000\n")
    # Usage : experiment1_analysis1.py image_id image_path output_file standard
    filey.writelines("/home/vsochat/python-lapack-blas/bin/python /home/vsochat/SCRIPT/python/brainmeta/experiment1/experiment1_analysis1.py %s %s %s" %(image_id,mr_file,output_file,standard))
    filey.close()
    os.system("sbatch " + ".job/" + uuids[i] + ".job")
