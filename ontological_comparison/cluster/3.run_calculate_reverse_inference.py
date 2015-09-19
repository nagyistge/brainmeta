#!/usr/bin/python
from glob import glob

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON"
data = "%s/data" %base        # mostly images
priors_pickles = glob("%s/priors/*.pkl" %(data))
tables_folder = "%s/priors/tables" %(data) # output folder for priors tables
scores_folder = "%s/scores" %(data)        # output folder for scores

if not os.path.exists(tables_folder):
    os.mkdir(tables_folder)

if not os.path.exists(scores_folder):
    os.mkdir(scores_folder)

for pkl in priors_pickles:
    contrast_id = os.path.split(pkl)[-1].split("_")[-1].replace(".pkl","")
    filey = ".job/revinf_%s.job" %(contrast_id)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
    filey.writelines("#SBATCH --output=.out/%s_%s.out\n" %(image_id,thresh))
    filey.writelines("#SBATCH --error=.out/%s_%s.err\n" %(image_id,thresh))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/ontological_comparison/cluster/3.calculate_reverse_inference.py %s %s %s" %(pkl, tables_folder, scores_folder))
    filey.close()
    os.system("sbatch -p russpold " + ".job/revinf_%s.job" %(contrast_id)) 
