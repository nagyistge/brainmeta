#!/usr/bin/python

import fnmatch
import pandas
import filecmp
import sys
import os

# This script will produce a file with a list of actual file path inputs for run_test_thresholding.py

search_directory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3"
doc_directory = "%s/doc" %(search_directory)

# First we will find all the group map volumes
volumes = []
for root, dirnames, filenames in os.walk(search_directory):
  for filename in fnmatch.filter(filenames, '*nii_tstat1.nii.gz'):
      print "Found %s/%s" %(root,filename)
      volumes.append(os.path.join(root, filename))

# Make some labels!
tasks = [os.path.split(v)[1].split(".")[0].split("_")[1] for v in volumes]
contrasts = [os.path.split(v)[1].split(".")[0].split("_")[2] for v in volumes]
groups = [os.path.split(v)[1].split(".")[0].split("_")[0] for v in volumes]
uid = ["%s_%s_%s" %(groups[v],tasks[v],contrasts[v]) for v in range(0,len(volumes))]
tasks_contrasts = ["%s_%s" %(tasks[v],contrasts[v]) for v in range(0,len(volumes))]
groups_tasks = ["%s_%s" %(groups[v],tasks[v]) for v in range(0,len(volumes))]
groups_contrasts = ["%s_%s" %(groups[v],contrasts[v]) for v in range(0,len(volumes))]

df = pandas.DataFrame()
df["uid"] = uid
df["files"] = volumes
df["tasks"] = tasks
df["groups"] = groups
df["contrasts"] = contrasts
df["tasks_contrasts"] = tasks_contrasts
df["groups_tasks"] = groups_tasks
df["groups_contrasts"] = groups_contrasts

# Sort by uid
df.to_csv("%s/hcp_groupmaps.tsv" %(doc_directory) ,sep="\t",index=False)



