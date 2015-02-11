#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import pandas
basedir = "/corral-repl/utexas/poldracklab/data/brainmeta/image_comparison/experiment1"
outdirectory = "%s/similarity_scores" %(basedir)
indirectory = "%s/mr" %(basedir)
standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

# Input file
input_file = "%s/openfmri_labels.tsv" %(basedir)

# Read in input file
inputs = pandas.read_csv(input_file,sep="\t")
# IMAGE_ID should correspond to integer ID

# Prepare a job script and submit
jobname = "exp1_test"
filey = ".job/%s.job" %(jobname)
filey = open(filey,"w")
for i in inputs["ID"]:
  image_id = i
  single_metrics = "%s/000%s.pkl" %(outdirectory,i)
  output_metrics = "%s/000%s.tsv" %(outdirectory,i)
  filey.writelines("/home1/02092/vsochat/SOFTWARE/python-venv/bin/python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment1_analysis1.py %s %s %s %s %s %s\n" %(image_id,indirectory,output_metrics,single_metrics,standard,input_file))

filey.close()
os.system("launch -s .job/%s.job -r 04:00:00 -p 1728 -e 1way -n exp1_test -j Analysis_Lonestar -m vsochat@stanford.edu" %(jobname))
