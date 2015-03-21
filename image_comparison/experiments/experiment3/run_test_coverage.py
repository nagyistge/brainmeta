#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster

import os
import time
import pandas
import numpy as np
import nibabel as nib
import pickle

# Input file with map paths
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3"
outdirectory = "%s/coverage" %(basedir)
input_file = "%s/doc/hcp_groupmaps_filter.tsv" %(basedir)
inputs = pandas.read_csv(input_file,sep="\t")
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

# Our groups will tells us the number of degrees of freedom for each map
groups = "%s/doc/hcp_10groups460_alltasks.pkl" %(basedir)

# Range of sampling percentages from 0.05 to 1.0
percentages = np.divide(range(1,101),100.0)

for percent_sample in percentages:
  for groupmap in inputs.iterrows():
    image_id = groupmap[1].uid
    filename = groupmap[1].files
    # if the output directory doesn't exist, make it
    topdir = "%s/coverage_%s" %(outdirectory,percent_sample)
    if not os.path.exists(topdir): os.mkdir(topdir)
    outfile = "%s/%s_coverage_%s.pkl" %(topdir,image_id,percent_sample)
    if not os.path.exists(outfile):
      filey = ".job/coverage_%s_%s.job" %(image_id,percent_sample)
      filey = open(filey,"w")
      filey.writelines("#!/bin/bash\n")
      filey.writelines("#SBATCH --job-name=coverage_%s_%s\n" %(image_id,percent_sample))
      filey.writelines("#SBATCH --output=.out/coverage_%s_%s.out\n" %(image_id,percent_sample))
      filey.writelines("#SBATCH --error=.out/coverage_%s_%s.err\n" %(image_id,percent_sample))
      filey.writelines("#SBATCH --time=2-00:00\n")
      filey.writelines("#SBATCH --mem=64000\n")
      filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/test_coverage.py %s %s %s %s %s" %(percent_sample,input_file,outfile,standard,image_id))
      filey.close()
      os.system("sbatch -p russpold " + ".job/coverage_%s_%s.job" %(image_id,percent_sample))
