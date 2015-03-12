#!/usr/bin/python

# This script will submit a job to run randomise for each of our 4D maps corresponding to a GROUP_TASK_CONTRAST for experiment 3

from glob import glob
import pandas
import fnmatch
import filecmp
import sys
import os

search_directory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3/group_maps"
doc_directory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc"

# Find all 4D volumes
volumes = []
for root, dirnames, filenames in os.walk(search_directory):
  for filename in fnmatch.filter(filenames, '*4D.nii.gz'):
      print "Found %s/%s" %(root,filename)
      volumes.append(os.path.join(root, filename))

volume_df = pandas.DataFrame(volumes)
volume_df.columns = ["image_path"]
volume_df.to_csv("%s/hcp_10groups_4Dlist.txt" %(doc_directory))

for volume in volumes:
  contrast_directory = os.path.split(volume)[0]
  output_nii = os.path.split(volume)[1].replace("_4D.nii.gz",".nii.gz")
  output_nii_path = "%s/%s.nii" %(contrast_directory,output_nii)
  if not os.path.exists(output_nii_path):
    filey = ".job/randomise_%s.job" %(output_nii)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=randomise_%s\n" %(output_nii))
    filey.writelines("#SBATCH --output=.out/randomise_%s.out\n" %(output_nii))
    filey.writelines("#SBATCH --error=.out/randomise_%s.err\n" %(output_nii))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("module load fsl\n")
    filey.writelines("randomise -i %s -o %s -1 -T\n" %(volume,output_nii_path))
    filey.close()
    os.system("sbatch -p russpold .job/randomise_%s.job" %(output_nii))
