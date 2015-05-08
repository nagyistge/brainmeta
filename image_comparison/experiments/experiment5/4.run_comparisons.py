#!/usr/bin/python

import pandas
import sys
import os

subject_path_pkl = "/scratch/users/vsochat/DATA/BRAINMETA/experiment5/doc/HCP_465_motor_cope_input_similar.pkl"
motor_ss_maps = pandas.read_pickle(subject_path_pkl)
data_directory = "/scratch/PI/russpold/work/HCP/data_matrices/tfMRI_MOTOR"
output_directory = "%s/comparisons" %(data_directory)
contrasts = list(motor_ss_maps.columns)
metrics = pandas.read_csv("/home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment5/4.comparison_metrics.txt")[0]

# We will create a brain mask
basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3"
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

for m in metrics:
    for c in range(0,len(contrasts)):
        con = contrasts[c]
        print "Processing contrast %s for metric %s" %(con,m)
        data_outfile_prefix = "%s/%s" %(output_directory,con)
        raw_data_path = "%s/%s_brainmask.pkl" %(data_directory,con)
        # Write job to file
        filey = ".job/comp_%s_%s.job" %(con,m)
        filey = open(filey,"w")
        filey.writelines("#!/bin/bash\n")
        filey.writelines("#SBATCH --job-name=comp_%s_%s\n" %(con,m))
        filey.writelines("#SBATCH --output=.out/comp_%s_%s.out\n" %(con,m))
        filey.writelines("#SBATCH --error=.out/comp_%s.err\n" %(con,m))
        filey.writelines("#SBATCH --time=2-00:00\n")
        filey.writelines("#SBATCH --mem=64000\n")
        filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment5/4.comparisons/4.%s.py %s %s %s %s %s %s %s\n" %(m,con,data_outfile_prefix,standard,
                                                         subject_path_pkl,raw_data_path,m))
        filey.close()
        os.system("sbatch -p russpold .job/comp_%s_%s.job" %(con,m))
