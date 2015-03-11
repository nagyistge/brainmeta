#!/usr/bin/python

# This script will generate a group map from a set of subjects for a particular contrast / task.

import sys
import nibabel
import os
import pandas

group_id = sys.argv[1]
contrast_name = sys.argv[2]
output_nii = sys.argv[3]
contrast_map = sys.argv[4]
contrast_list = sys.argv[5]
subs = sys.argv[6]

# Parse the list of subjects
subjects = subs.split(",")

# Read in the contrast map nifti, and list of subjects
nii = nibabel.load(contrast_map)

# Extract subject ids
nii_subjects = pandas.read_csv(contrast_list,header=None)
nii_subjects = nii_subjects[0].tolist()
expression = re.compile("/[0-9]+/")
subject_ids = []
for nii_sub in nii_subjects:
  match = expression.search(nii_sub)
  subject_ids.append(nii_sub[match.start()+1:match.end()-1])

# Check #1: parsing the subject ids from the copes list
if nii.shape[3] != len(subject_ids):
  print "ERROR: did not extract subject IDS - expecting all numerical in path! Exiting."
  sys.exit(32)

# Now we need to find the timepoints with our subjects (0 index)
indices = [subject_ids.index(x) for x in subjects]

# Check #2: that we have all subjects in the 4D file
if len(indices) != len(subjects):
  print "ERROR: did not find all subjects in 4D group map. Exiting!"
  sys.exit(32)

# Extract only those timepoints, and make a new image!
data = nii.get_data()
subset = data[:,:,:,indices]

# Write a new nifti image!
affine = nii.get_affine()
header = nii.get_header()
nii_group = nibabel.Nifti1Image(subset,affine=affine,header=header) 
nibabel.save(nii_group,file=output_nii)
