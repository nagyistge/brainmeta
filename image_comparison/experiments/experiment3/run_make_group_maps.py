#!/usr/bin/python

# This script will first define 10 groups per task (N=46 per group) for a total of 860 groups for the 465 individuals in HCP.  Group membership will be constant across tasks/contrasts. Then we will derive the maps.

from glob import glob
import numpy as np
import pickle
import pandas
import os

##### PREP ############################################################

# Here is the list of hcp subjects with all tasks
subjects = pandas.read_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_465_with_all_tasks.tsv",sep="\t")
subjects = subjects.id.tolist()
random.shuffle(subjects)

# Break into an even number, and save remaining 5 to add after
remaining_five = subjects[0:5]
all_subjects = subjects
subjects = subjects[5:len(subjects)]

groups = dict()
index = 0
for g in range(1,11):
  groups["GRP%02d" %(g)] = subjects[index:index+46]
  # If we are in first five groups
  if g in range(1,6):
    sub = remaining_five.pop()
    groups["GRP%02d" %(g)].append(sub)
  index = index + 46

# Check that we have everyone, and that there aren't any repeated
check = []
for grp,subs in groups.iteritems():
  check = check + subs 

check.sort()
all_subjects.sort()
check == all_subjects # True

# Save groups to file
pickle.dump(groups,open("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_10groups_alltasks.pkl","wb"))

# Read in the file with contrasts. We now need to make unique IDs for each!
contrasts = pandas.read_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_contrasts.tsv",sep="\t")
tasks = np.unique(contrasts.task)

task_ids = contrasts.task.tolist()
contrast_ids = contrasts.index.tolist() + 1
contrast_ids = ["CON%02d" %(c+1) for c in contrast_ids]

count = 1
for task in tasks:
  task_ids = [t.replace(task,"TASK%02d" %(count)) for t in task_ids]
  count +=1

contrasts["task_id"] = task_ids
contrasts["contrast_ids"] = contrast_ids

# Now generate what will be the full file names / identifiers
full_names = []
for con in contrasts.iterrows():
  name = "%s_%s" %(con[1].task_id,con[1].contrast_ids)
  full_names.append(name)

contrasts["id"] = full_names
contrasts.to_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_contrasts_id.tsv",sep="\t")

##### RUN! ############################################################

# Variables to run jobs
contrasts = pandas.read_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_contrasts_id.tsv",sep="\t")
output_directory = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3/group_maps"
group_maps_directory = "/scratch/PI/russpold/work/HCP/group_maps/nii"
groups = pickle.load(open("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_10groups_alltasks.pkl","rb"))
group_list_directory = "/scratch/PI/russpold/work/HCP/group_maps/copes"

# Now for each group, for each contrast, we will read in the image file, select the subjects, and run randomise to generate the group maps
for grp,subs in groups.iteritems():
  for con in contrasts.iterrows():
    contrast_name = con[1].id
    contrast_directory = "%s/%s" %(output_directory,contrast_name)
    if not os.path.exists(contrast_directory): os.mkdir(contrast_directory)
    subs.sort()
    subs = ["%s" %(s) for s in subs]
    subs = ",".join(subs)
    contrast_subjects = "%s/%s_%s_copes.txt" %(group_list_directory,con[1].task,con[1].contrasts)
    contrast_map = "%s/%s_%s_copes_4D.nii.gz" %(group_maps_directory,con[1].task,con[1].contrasts)
    output_nii = "%s/%s_%s_4D.nii.gz" %(output_directory,grp,contrast_name)    
    # Write job to file
    filey = ".job/mk_%.job" %(contrast_name)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=mk_%s\n" %(contrast_name))
    filey.writelines("#SBATCH --output=.out/mk_%s.out\n" %(contrast_name))
    filey.writelines("#SBATCH --error=.out/mk_%s.err\n" %(contrast_name))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python /scratch/users/vsochat/SCRIPT/python/brainmeta/experiments/experiment3/make_group_maps.py %s %s %s %s\n" %(grp,contrast_name,output_nii,contrast_map,contrast_list,subs))
    filey.close()
    os.system("sbatch -p russpold .job/mk_%s.job" %(contrast_name))
