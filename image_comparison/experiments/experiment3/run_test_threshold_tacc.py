#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import pandas
from clusterhcp.database import get_hcp_paths
from clusterhcp.stats import select_unrelated_samples

# This is the top directory of the HCP data
top_directory = "/corral-tacc/tacc/HCP/"
basedir = "/scratch/02092/vsochat/DATA/BRAINMETA/IMAGE_COMPARISON/experiments/experiment3"
outdirectory = "%s/permutations" %(basedir)

# This is an optional dictionary to replace disk names {"lookup":"replacement"}
disks = {"Disk1of5":"disk1",
        "Disk2of5":"disk2",
        "Disk3of5":"disk3",
        "Disk4of5":"disk4",
        "Disk5of5":"disk5"}

# This is the size of the groups that we want to generate
size = 46

# Read in contrast
input_file = "%s/doc/hcp_contrasts_id_filter.csv" %(basedir)

# Read in input file
contrasts = pandas.read_csv(input_file,sep="\t")


# We will submit a file with all jobs
jobname = "experiment3_maps"
filey = ".job/%s.job" %(jobname)
filey = open(filey,"w")


# STEP 1: First we will make groups A and B maps for 500 runs x 96 contrasts
nruns = 500
for i in range(0,nruns):
    # get unrelated samples, default is two groups
    A,B = select_unrelated_samples(paths,size=size)
    groupA = ",".join(A)
    groupB = ",".join(B)

    # Top level of output directory is for the iteration
    output_directory = "%s/%s" %(outdirectory,i)
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        os.mkdir("%s/maps" %(output_directory))
        os.mkdir("%s/comparisons" %(output_directory))
    maps_directory = "%s/maps" %(output_directory)

    # Now generate images for each contrast, group
    for con in contrasts.iterrows():
        task = con[1]["task"]
        contrast = con[1]["contrasts"]
        map_id = con[1]["id"]
        paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast, disks=disks)
        filey.writelines("/work/02092/vsochat/SOFTWARE/python-venv/bin/python /work/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_maps_tacc.py %s %s %s %s\n" %(groupA,groupB,maps_directory,map_id))

# Close the file
filey.close()

# Submit the job
os.system("launch -s .job/%s.job -r 04:00:00 -p 1728 -e 1way -n exp3_maps -j Analysis_Lonestar -m vsochat@stanford.edu" %(jobname))



# STEP 2: Now run analysis over each iteration folder

# We will run over a set of thresholds
thresholds = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
thresholds = ",".join(thresholds)

standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

jobname = "experiment3_analysis"
filey = ".job/%s.job" %(jobname)
filey = open(filey,"w")

for i in range(0,nruns):    
    output_directory = "%s/%s" %(outdirectory,i)
    maps_directory = "%s/maps" %(output_directory)
    # There will be one of group A for each of group B, unrelated
    image_pairs = pandas.DataFrame()
    image_pairs["groupA"] = glob("%s/*_groupA_tstat1.nii.gz")
    image_pairs["groupB"] = glob("%s/*_groupB_tstat1.nii.gz")
    for run in image_pairs.iterrows():
        #TODO: test to make sure the groups are ordered correctly (matching rows)
        groupA_path = run[1]["groupA"]
        groupB_path = run[1]["groupB"]
        contrast_task = groupA_path.split("/")[-1].replace("_groupA_tstat1.nii.gz","")
        dofA = int(open("%s/maps/%s_N.txt" %(output_directory,contrast_task),"r").readlines()[0].strip("\n"))
        dofB = int(open("%s/maps/%s_N.txt" %(output_directory,contrast_task),"r").readlines()[0].strip("\n"))
        output_pkl = "%s/comparisons/%s.pkl" %(output_directory,contrast_task)
        filey.writelines("/work/02092/vsochat/SOFTWARE/python-venv/bin/python /work/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/run_test_threshold_tacc.py %s %s %s %s %s %s %s %s\n" %(groupA_path,groupB_path,thresholds,standard,output_pkl,dofA,dofB,contrast_task))


# Close the file
filey.close()

# Submit the job
os.system("launch -s .job/%s.job -r 04:00:00 -p 1728 -e 1way -n exp3_analysis -j Analysis_Lonestar -m vsochat@stanford.edu" %(jobname))
