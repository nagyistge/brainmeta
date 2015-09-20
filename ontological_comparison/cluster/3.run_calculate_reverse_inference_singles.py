#!/usr/bin/python
from glob import glob
import os

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON"
data = "%s/data" %base        # mostly images
priors_folder = "%s/priors" %(data)
tables_folder = "%s/priors/tables" %(data)         # folder with priors tables
scores_folder = "%s/individual_scores" %(data)     # output folder for individual scores
images_folder = "%s/resampled_z" %data

if not os.path.exists(scores_folder):
    os.mkdir(scores_folder)

images = glob("%s/*"%images_folder)

# Generate an output file for each image
for i in range(0,len(images)):
    input_image = images[i]
    image_id = os.path.split(images[i])[-1].replace(".nii.gz","")
    filey = ".jobs/inf_%s.job" %(image_id)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
    filey.writelines("#SBATCH --output=.out/%s.out\n" %(image_id))
    filey.writelines("#SBATCH --error=.out/%s.err\n" %(image_id))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/ontological_comparison/cluster/3.calculate_reverse_inference_singles.py %s %s %s %s" %(input_image, tables_folder, scores_folder, priors_folder))
    filey.close()
    os.system("sbatch -p russpold " + ".jobs/inf_%s.job" %(image_id)) 
