#!/usr/bin/python

import pandas
import nibabel
import numpy
import sys
import os
import re

search_directory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3/group_maps"
doc_directory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc"
volumes = pandas.read_csv("%s/hcp_10groups_unrelated_4Dlist.txt" %(doc_directory))
image_paths = volums.image_path.tolist()
expression = re.compile("TASK04_CON70")
idx = []

for i in range(0,len(image_paths)):
  if expression.search(image_paths[i]):
      idx.append(i)
image_paths = [image_paths[x] for x in idx]
sizes = []

for volume in image_paths:
  nii = nibabel.load(volume)
  sizes.append(nii.shape[3])

# [50, 51, 46, 51, 51, 47, 52, 43, 52, 49]
