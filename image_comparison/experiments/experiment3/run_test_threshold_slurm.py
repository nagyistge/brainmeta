#!/usr/bin/python

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import pickle
import time
import pandas
import numpy
import nibabel
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


### STEP 0: Try loading each file
for con in contrasts.iterrows():     
    task = con[1]["task"]
    contrast = con[1]["contrasts"]
    map_id = con[1]["id"]
    print "Testing %s" %(map_id)
    paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast, disks=disks)
    for p in paths:
        try:
            nibabel.load(p)
        except:
            print "ERROR: problem with image %s" %(p)


### STEP 1: First we will define groups A and B for 500 runs ############################
nruns=500
# Generate the subject groups first, and then run them over iterations
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
    outputA = "\n".join([str(x) for x in subjectsA])
    outputB = "\n".join([str(x) for x in subjectsB])
    fileyA = open("%s/copesA.txt" %(output_directory),"wb")
    fileyB = open("%s/copesB.txt" %(output_directory),"wb") 
    fileyA.writelines(outputA)
    fileyB.writelines(outputB)
    fileyA.close()
    fileyB.close()
    

### STEP 2: Now generate grid runs to be run ACROSS runs ############################
nruns = 500
counter = 1  
jobnum = 1
launch_file = ".job/experiment3_%s.job" %(jobnum)
filey = open(launch_file,"wb")    
for con in contrasts.iterrows(): 
    print "Processing contrast %s" %(con[1].id)  
    task = con[1]["task"]
    contrast = con[1]["contrasts"]
    map_id = con[1]["id"]
    paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast, disks=disks)
    A = [s for s in paths if s.split("/")[8] in subjectsA]
    B = [s for s in paths if s.split("/")[8] in subjectsB]
    groupA = ",".join(A)
    groupB = ",".join(B)
    for i in range(0,nruns):
        # Top level of output directory is for the iteration
        output_directory = "%s/%s" %(outdirectory,i)
        maps_directory = "%s/maps" %(output_directory)
        # Load the groups
        subjectsA = open("%s/copesA.txt" %(output_directory),"r").readlines()
        subjectsB = open("%s/copesB.txt" %(output_directory),"r").readlines()
        subjectsA = [x.replace("\n","") for x in subjectsA]
        subjectsB = [x.replace("\n","") for x in subjectsB]
        if counter < 4096:
            filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_maps_tacc.py %s %s %s %s\n" %(subjectsA,subjectsB,maps_directory,map_id)) 
            counter = counter + 1
        else:
            filey.close()
            os.system("launch -s %s -r 20:00 -e 1way -j Analysis_Lonestar -m vsochat@stanford.edu" %(launch_file))
            jobnum = jobnum + 1
            launch_file = ".job/experiment3_%s.job" %(jobnum)
            filey = open(launch_file,"wb")    
            filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_maps_tacc.py %s %s %s %s\n" %(subjectsA,subjectsB,maps_directory,map_id)) 
            counter = 1 


### STEP 3: Find any missing maps ###########################################################
counter = 1  
jobnum = 1
launch_file = ".job/miss3_%s.job" %(jobnum)
filey = open(launch_file,"wb")    
for con in contrasts.iterrows(): 
    print "Processing contrast %s" %(con[1].id)  
    task = con[1]["task"]
    contrast = con[1]["contrasts"]
    map_id = con[1]["id"]
    for i in range(0,nruns):
        # Top level of output directory is for the iteration
        output_directory = "%s/%s" %(outdirectory,i)
        maps_directory = "%s/maps" %(output_directory)
        outfileA = "%s/%s_groupA_tstat1.nii.gz" %(maps_directory,map_id)
        outfileB = "%s/%s_groupA_tstat1.nii.gz" %(maps_directory,map_id)
        if not os.path.exists(outfileA):
            print "Missing %s" %(outfileA)        
            subjectsA = open("%s/copesA.txt" %(output_directory),"r").readlines()
            subjectsA = [x.replace("\n","") for x in subjectsA]
            if counter < 4096:
                filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_map_single_tacc.py %s %s %s %s\n" %(subjectsA,"A",maps_directory,map_id)) 
                counter = counter + 1
            else:
                filey.close()
                os.system("launch -s %s -r 20:00 -e 1way -j Analysis_Lonestar -m vsochat@stanford.edu" %(launch_file))
                jobnum = jobnum + 1
                launch_file = ".job/miss3_%s.job" %(jobnum)
                filey = open(launch_file,"wb")    
                filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_map_single_tacc.py %s %s %s %s\n" %(subjectsA,"A",maps_directory,map_id)) 
                counter = 1 
        if not os.path.exists(outfileB):
            print "Missing %s" %(outfileB)        
            subjectsB = open("%s/copesB.txt" %(output_directory),"r").readlines()
            subjectsB = [x.replace("\n","") for x in subjectsB]
            if counter < 4096:
                filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_map_single_tacc.py %s %s %s %s\n" %(subjectsB,"B",maps_directory,map_id)) 
                counter = counter + 1
            else:
                filey.close()
                os.system("launch -s %s -r 20:00 -e 1way -j Analysis_Lonestar -m vsochat@stanford.edu" %(launch_file))
                jobnum = jobnum + 1
                launch_file = ".job/miss3_%s.job" %(jobnum)
                filey = open(launch_file,"wb")    
                filey.writelines("python /home1/02092/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_map_single_tacc.py %s %s %s %s\n" %(subjectsB,"B",maps_directory,map_id)) 
                counter = 1 
    
    

### STEP 4: Which runs are done? ############################


done = []
for i in range(0,nruns):
   num = glob("/work/02092/vsochat/wrangler/DATA/BRAINMETA/IMAGE_COMPARISON/experiments/experiment3/permutations/%s/maps/*_tstat1.nii.gz" %(i))
   if len(num)==94: done.append(i)

redo = [x for x in range(0,500) if x not in done]
for r in redo:
  os.system("rm -rf /work/02092/vsochat/wrangler/DATA/BRAINMETA/IMAGE_COMPARISON/experiments/experiment3/permutations/%s" %(r))


### STEP 5: Calculate similarities for maps that have all images generated ############################


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


