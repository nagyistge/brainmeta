#!/usr/bin/python

import pandas
import sys
import os

subject_path_pkl = "/scratch/users/vsochat/DATA/BRAINMETA/experiment5/doc/HCP_465_motor_cope_input.pkl"
motor_ss_maps = pandas.read_pickle(subject_path_pkl)
data_directory = "/scratch/PI/russpold/work/HCP/data_matrices/tfMRI_MOTOR"
contrasts = list(motor_ss_maps.columns)

# We will create a brain mask
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3"
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

# Boolean to calculate original matrix or not
do_data_prep = False

for c in range(0,len(contrasts)):
    con = contrasts[c]
    print "Processing %s" %(con)
    data_outfile_prefix = "%s/%s" %(data_directory,con)
    # Write job to file
    filey = ".job/data_matrix_%s.job" %(con)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=data_matrix_%s\n" %(con))
    filey.writelines("#SBATCH --output=.out/data_matrix_%s.out\n" %(con))
    filey.writelines("#SBATCH --error=.out/data_matrix_%s.err\n" %(con))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment5/2.cluster_subjects.py %s %s %s %s %s\n" %(con,data_outfile_prefix,standard,subject_path_pkl,do_data_prep))
    filey.close()
    os.system("sbatch -p russpold .job/data_matrix_%s.job" %(con))
