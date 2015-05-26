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



### STEP 1: First we will make groups A and B maps for 500 runs x 96 contrasts ############################

# Launch single jobs...

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
for i in range(0,nruns):
    # Top level of output directory is for the iteration
    output_directory = "%s/%s" %(outdirectory,i)
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        os.mkdir("%s/maps" %(output_directory))
        os.mkdir("%s/comparisons" %(output_directory))
    maps_directory = "%s/maps" %(output_directory)
    # Get a group of subjects for one of the tasks/contrasts
    task = 'tfMRI_WM'
    contrast = '0BK_BODY'
    paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast, disks=disks) # default is to return 465 subjects with all tasks
    # get unrelated samples, default is two groups
    A,B = select_unrelated_samples(paths,size=size)
    subjectsA = [int(s.split("/")[8]) for s in A]
    subjectsB = [int(s.split("/")[8]) for s in B]
    for con in contrasts.iterrows():     
        task = con[1]["task"]
        contrast = con[1]["contrasts"]
        map_id = con[1]["id"]
        paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast, disks=disks)
        # Filter down to A and B
        A = [s for s in A if int(s.split("/")[8]) in subjectsA]
        B = [s for s in B if int(s.split("/")[8]) in subjectsB]
        groupA = ",".join(A)
        groupB = ",".join(B)
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


### STEP 2: Find missing maps (2047 out of ~23,000 were not generated) ###########################################################

counter = 1
jobnum = 1
launch_file = "map3_%s.job" %(jobnum)
filey = open(launch_file,"wb")

for i in range(271,nruns):    
    output_directory = "%s/%s" %(outdirectory,i)
    maps_directory = "%s/maps" %(output_directory)
    # There will be one of group A for each of group B, unrelated
    groupA_files = numpy.sort(glob("%s/*_groupA_tstat1.nii.gz" %(maps_directory)))
    groupB_files = numpy.sort(glob("%s/*_groupB_tstat1.nii.gz" %(maps_directory)))
    if len(groupA_files) != len(groupB_files):
            groupA_id = [x.split("/")[-1].replace("_groupA_tstat1.nii.gz","") for x in groupA_files]
            groupB_id = [x.split("/")[-1].replace("_groupB_tstat1.nii.gz","") for x in groupB_files]
            missingA = [x for x in groupA_id if x not in groupB_id]
            missingB = [x for x in groupB_id if x not in groupA_id]
            missing = numpy.unique(missingA + missingB)
            print "Iteration %s is missing %s maps" %(i,len(missing))
            missing_df = contrasts[contrasts.id.isin(missing)]
            # Load subject groups from directory
            copelistA = glob("%s/*_groupA_cope_inputs.txt" %(maps_directory))[0]
            copelistB = glob("%s/*_groupB_cope_inputs.txt" %(maps_directory))[0]
            groupA = open(copelistA,"r").readlines()
            groupA = [g.strip("\n") for g in groupA]
            groupB = open(copelistB,"r").readlines()
            groupB = [g.strip("\n") for g in groupB]
            groupA = ",".join(groupA)
            groupB = ",".join(groupB)
            for con in missing_df.iterrows():
                map_id = con[1].id
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


### STEP 3: Calculate similarities for maps that have all images generated ############################


# We will run over a set of thresholds
thresholds = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
thresholds = ",".join([str(x) for x in thresholds])

standard = "%s/standard/MNI152_T1_brain_mask_exp3.nii.gz" %(basedir)

# We will again use launcher, and submit 4095 jobs at once
counter = 1
jobnum = 1
launch_file = "sim3_%s.job" %(jobnum)
filey = open(launch_file,"wb")

for i in range(0,nruns):    
    output_directory = "%s/%s" %(outdirectory,i)
    maps_directory = "%s/maps" %(output_directory)
    # There will be one of group A for each of group B, unrelated
    image_pairs = pandas.DataFrame()
    groupA_files = numpy.sort(glob("%s/*_groupA_tstat1.nii.gz" %(maps_directory)))
    groupB_files = numpy.sort(glob("%s/*_groupB_tstat1.nii.gz" %(maps_directory)))
    if len(groupA_files) == len(groupB_files):
	    image_pairs["groupA"] = groupA_files
	    image_pairs["groupB"] = groupB_files
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
		            filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/test_thresholding_tacc.py %s %s %s %s %s %s %s %s\n" %(groupA_path,groupB_path,thresholds,standard,output_pkl,dofA,dofB,contrast_task))        
		            counter = counter + 1
		        else:
		            filey.close()
		            os.system("launch -s %s -r 48:00 -e 1way -n %s -j Analysis_Lonestar -m vsochat@stanford.edu" %(launch_file, launch_file.replace(".job","")))
		            jobnum = jobnum + 1
		            launch_file = "sim3_%s.job" %(jobnum)
		            print "Writing to new %s" %(launch_file)
		            filey = open(launch_file,"wb")
		            filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/test_thresholding_tacc.py %s %s %s %s %s %s %s %s\n" %(groupA_path,groupB_path,thresholds,standard,output_pkl,dofA,dofB,contrast_task))        
		            counter = 1
		else:
		    print "Error, mismatch for run %s %s" %(run[1].groupA,run[1].groupB)
    else:
        print "Error: iteration %s is missing maps, fix!" %(i)

