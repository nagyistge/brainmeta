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
