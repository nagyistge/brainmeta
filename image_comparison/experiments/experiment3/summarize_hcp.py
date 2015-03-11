#!/usr/bin/python

# This script will summarize the tasks that we have for HCP, and output how many subjects per, etc.

import os
import sys
import nibabel
import pandas
import numpy as np
from glob import glob

# We will extract timeseries for each 4D file from out task data
inputs_nii = glob("/scratch/PI/russpold/work/HCP/group_maps/nii/*4D.nii.gz")

# Directory with subject lists
copes_directory = "/scratch/PI/russpold/work/HCP/group_maps/copes"
copes_files = glob("%s/*" %(copes_directory))

# How many copes?
len(copes_files)
# 86

# Contrasts file corresponding to the copes files
contrasts_file = pandas.read_csv("/scratch/PI/russpold/work/HCP/group_maps/doc/hcp_contrasts.tsv",sep="\t")

# This is raw data for all tasks, resting, all subjects. We can get all subject ids from here
timeseries_file = pandas.read_csv("/scratch/PI/russpold/work/HCP/group_maps/doc/hcp_timeseries_paths.tsv",sep="\t")

# We will save a data frame to keep a lot of who has data for each task!
# Get subids from the complete task data
subids = np.unique([x.split("/")[7] for x in timeseries_file.timeseries_files])
# 501

df = pandas.DataFrame(columns=subids)
# Here we will save list of task contrasts
task_contrasts = []
count=1
# We now want to find a set of subjects that have data for all the tasks
for cope_file in copes_files:
  task_contrast = os.path.split(cope_file)[1].replace("_copes.txt","")
  task_contrasts.append(task_contrast)
  cope_file = pandas.read_csv(cope_file,sep="\t",header=None)
  subids = [x.split("/")[7] for x in cope_file[0]]
  df.loc[count,subids] = 1
  count = count+1

# Add the task_contrasts
df.index = task_contrasts
df.to_csv("/scratch/PI/russpold/work/HCP/group_maps/doc/hcp_ss_tasklog.tsv",sep="\t")

# For which subjects do we have data for all tasks?
subids = df.sum(axis=0)==df.shape[0]
subids = subids.index[subids==True]
subids = subids.tolist()
len(subids)
hcp_subjects = pandas.DataFrame()
hcp_subjects["id"] = subids
hcp_subjects.to_csv("/scratch/PI/russpold/work/HCP/group_maps/doc/hcp_465_with_all_tasks.tsv",sep="\t")
# 465
