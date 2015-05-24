#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import pickle
import time
import pandas
import numpy
from glob import glob
from clusterhcp.database import get_hcp_paths
from clusterhcp.stats import select_unrelated_samples

# This is the top directory of the HCP data
top_directory = "/scratch/projects/UT/poldracklab/data/HCP"
basedir = "/work/02092/vsochat/wrangler/DATA/BRAINMETA/IMAGE_COMPARISON/experiments/experiment3"
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

# STEP 1: First we will make groups A and B maps for 500 runs x 96 contrasts
nruns = 500
number_jobs = int(os.popen('squeue -u vsochat | wc -l').read().replace("\n"))-1
for i in range(0,nruns):
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
        # get unrelated samples, default is two groups
        A,B = select_unrelated_samples(paths,size=size)
        groupA = ",".join(A)
        groupB = ",".join(B)
        output_file_check = "%s/%s_groupB_tstat1.nii.gz" %(maps_directory,map_id)
        if not os.path.exists(output_file_check):
            filey = ".job/%s_%s.job" %(i,map_id)
            filey = open(filey,"w")
            filey.writelines("#!/bin/bash\n")
            filey.writelines("#SBATCH --job-name=%s_%s\n" %(i,map_id))
            filey.writelines("#SBATCH --output=.out/%s_%s.out\n" %(i,map_id))
            filey.writelines("#SBATCH --error=.out/%s_%s.err\n" %(i,map_id))
            filey.writelines("#SBATCH --time=48:00\n")
            filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_maps_tacc.py %s %s %s %s\n" %(groupA,groupB,maps_directory,map_id)) 
            filey.close()
            while number_jobs >= 50:
                time.sleep(15)
            os.system("sbatch -p normal -A TG-CCR130001 -n 1 .job/%s_%s.job" %(i,map_id))
            number_jobs = int(os.popen('squeue -u vsochat | wc -l').read().replace("\n",""))-1

# Or use launch:
nruns = 500
counter = 0
jobnum = 1
launch_file = "exp3_%s.job" %(jobnum)
filey = open(launch_file,"wb")
for i in range(10,nruns):
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
        # get unrelated samples, default is two groups
        A,B = select_unrelated_samples(paths,size=size)
        groupA = ",".join(A)
        groupB = ",".join(B)
        output_file_check = "%s/%s_groupB_tstat1.nii.gz" %(maps_directory,map_id)
        if not os.path.exists(output_file_check):
            if counter < 4096:
                filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_maps_tacc.py %s %s %s %s\n" %(groupA,groupB,maps_directory,map_id)) 
                counter = counter + 1
            else:
                filey.close()
                os.system("launch -s %s -r 48:00 -e 1way -n %s -j Analysis_Lonestar -m vsochat@stanford.edu" %(launch_file, launch_file.replace(".job","")))
                jobnum = jobnum + 1
                launch_file = "exp3_%s.job" %(jobnum)
                print "Writing to new %s" %(launch_file)
                filey = open(launch_file,"wb")
                filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_maps_tacc.py %s %s %s %s\n" %(groupA,groupB,maps_directory,map_id)) 
                counter = 1


# STEP 2: Now run analysis over each iteration folder

# We will run over a set of thresholds
thresholds = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
thresholds = ",".join([str(x) for x in thresholds])

standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

# We will again use launcher, and submit 4095 jobs at once
counter = 0
jobnum = 1
launch_file = "sim3_%s.job" %(jobnum)
filey = open(launch_file,"wb")

for i in range(0,nruns):    
    output_directory = "%s/%s" %(outdirectory,i)
    maps_directory = "%s/maps" %(output_directory)
    # There will be one of group A for each of group B, unrelated
    image_pairs = pandas.DataFrame()
    image_pairs["groupA"] = numpy.sort(glob("%s/*_groupA_tstat1.nii.gz" %(maps_directory)))
    image_pairs["groupB"] = numpy.sort(glob("%s/*_groupB_tstat1.nii.gz" %(maps_directory)))
    for run in image_pairs.iterrows():
        mapnameA = run[1].groupA.replace("_groupA_tstat1.nii.gz","")
        mapnameB = run[1].groupB.replace("_groupB_tstat1.nii.gz","")
        if mapnameA == mapnameB:
            groupA_path = run[1]["groupA"]
            groupB_path = run[1]["groupB"]
            contrast_task = groupA_path.split("/")[-1].replace("_groupA_tstat1.nii.gz","")
            dofA = int(open("%s/maps/%s_groupA_N.txt" %(output_directory,contrast_task),"r").readlines()[0].strip("\n"))
            dofB = int(open("%s/maps/%s_groupB_N.txt" %(output_directory,contrast_task),"r").readlines()[0].strip("\n"))
            output_pkl = "%s/comparisons/%s.pkl" %(output_directory,contrast_task)
            if not os.path.exists(output_pkl):
                if counter < 4096:
                    filey.writelines("python /work/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/run_test_threshold_tacc.py %s %s %s %s %s %s %s %s\n" %(groupA_path,groupB_path,thresholds,standard,output_pkl,dofA,dofB,contrast_task))        
                    counter = counter + 1
                else:
                    filey.close()
                    os.system("launch -s %s -r 48:00 -e 1way -n %s -j Analysis_Lonestar -m vsochat@stanford.edu" %(launch_file, launch_file.replace(".job","")))
                    jobnum = jobnum + 1
                    launch_file = "sim3_%s.job" %(jobnum)
                    print "Writing to new %s" %(launch_file)
                    filey = open(launch_file,"wb")
                    filey.writelines("python /work/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/run_test_threshold_tacc.py %s %s %s %s %s %s %s %s\n" %(groupA_path,groupB_path,thresholds,standard,output_pkl,dofA,dofB,contrast_task))        
                    counter = 1
        else:
            print "Error, mismatch for run %s" %(run)


