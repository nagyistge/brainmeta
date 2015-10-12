#!/usr/bin/python
from glob import glob
import pickle
import os

# For each concept, (N=140):
#    Select a number G from 1...[total "in" group] as the size of the set to investigate
#        For each image in (entire) "in"set:
#            For some number of iterations:
#                Randomly select G other images for "in" set, calculate P(concept|image)
#                Take mean score of iterations as P(concept|image)


base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON/v2"
data = "%s/data" %base        # mostly images
likelihood_pickles = glob("%s/likelihood/*.pkl" %(data))
scores_folder = "%s/group_size_vary_scores" %(data)     # output folder for group size scores

if not os.path.exists(scores_folder):
    os.mkdir(scores_folder)

# Generate an output file for each image
for i in range(0,len(likelihood_pickles)):
    node = likelihood_pickles[i]
    group = pickle.load(open(node,"rb"))
    # We will vary the size of the "in" set from 1 to the number of "in" images.
    in_group = group["in"]
    for j in range(1,len(in_group)-1):
        all_images = group["in"] + group["out"]
        for image in all_images:
            image_id = os.path.split(image)[1].replace(".nii.gz","")
            # Output convention is [node]_size_[size]_[image_id].pkl
            run_id = "%s_size_%s_%s" %(group["nid"],in_group_index,image_id)
            output_pkl = "%s/%s.pkl" %(scores_folder,run_id)
            if not os.path.exists(output_pkl):
                filey = ".jobs/ri_%s.job" %(run_id)
                filey = open(filey,"w")
                filey.writelines("#!/bin/bash\n")
                filey.writelines("#SBATCH --job-name=%s\n" %(run_id))
                filey.writelines("#SBATCH --output=.out/%s.out\n" %(run_id))
                filey.writelines("#SBATCH --error=.out/%s.err\n" %(run_id))
                filey.writelines("#SBATCH --time=2-00:00\n")
                filey.writelines("#SBATCH --mem=64000\n")
                filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/ontological_comparison/cluster/classification-framework/4.explore_group_sizes.py %s %s %s %s" %(image, node, output_pkl,j))
                filey.close()
                os.system("sbatch -p russpold " + ".jobs/ri_%s.job" %(image_id))
