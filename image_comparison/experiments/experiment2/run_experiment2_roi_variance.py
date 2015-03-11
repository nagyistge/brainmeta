#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster

import os
import time
import pandas
import numpy as np
import nibabel as nib

# Input file
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment1"
atlas_nifti = "/scratch/users/vsochat/DATA/ATLAS/AAL/ROI_MNI_V4.nii"
atlas_file = "/scratch/users/vsochat/DATA/ATLAS/AAL/ROI_MNI_V4.txt"
outdirectory = "%s/atlas_permutations" %(basedir)
input_file = "%s/openfmri_labels.tsv" %(basedir)
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

# Read in the atlas file
atlas_labels = pandas.read_csv(atlas_file,sep="\t",header=None)

# Read in the atlas file
atlas = nib.load(atlas_nifti)

# Sanity check - 116 regions in both
len(np.unique(atlas.get_data()))-1 == len(np.unique(atlas_labels[2]))
number_regions = len(np.unique(atlas_labels[2]))

# This will be a gold standard correlation data file
gs_file = "%s/gs_comparisons.tsv" %(basedir)

# We will sample some number of regions each time, from 1 through 116, to represent different mask sizes
# This will simulate different levels of coverage!
for num_regions in range(1,number_regions):
  outfile = "%s/region_%s.pkl" %(outdirectory,num_regions)
  if not os.path.exists(outfile):
    filey = ".job/coverage_%s.job" %(num_regions)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=coverage_%s\n" %(num_regions))
    filey.writelines("#SBATCH --output=.out/coverage_%s.out\n" %(num_regions))
    filey.writelines("#SBATCH --error=.out/coverage_%s.err\n" %(num_regions))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment2_atlas_variance.py %s %s %s %s %s %s %s" %(atlas_nifti,num_regions,atlas_file,input_file,outfile,standard,gs_file))
    filey.close()
    os.system("sbatch -p russpold " + ".job/coverage_%s.job" %(num_regions))
