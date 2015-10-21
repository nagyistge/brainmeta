#!/usr/bin/python
from glob import glob
import pickle
import os

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON/v2"
data = "%s/data" %base        # mostly images
likelihood_pickles = glob("%s/likelihood/*.pkl" %(data))
scores_folder = "%s/individual_scores" %(data)     # output folder for individual scores
tables_folder = "%s/likelihood/tables" %(data)

if not os.path.exists(scores_folder):
    os.mkdir(scores_folder)

# Generate an output file for each image
for i in range(0,len(likelihood_pickles)):
    node = likelihood_pickles[i]
    group = pickle.load(open(node,"rb"))
    all_images = group["in"] + group["out"]
    for image in all_images:
        image_id = os.path.split(image)[1].replace(".nii.gz","")
        output_pkl = "%s/%s_%s.pkl" %(scores_folder,group["nid"],image_id)
        if not os.path.exists(output_pkl):
            filey = ".jobs/ri_%s.job" %(image_id)
            filey = open(filey,"w")
            filey.writelines("#!/bin/bash\n")
            filey.writelines("#SBATCH --job-name=%s\n" %(image_id))
            filey.writelines("#SBATCH --output=.out/%s.out\n" %(image_id))
            filey.writelines("#SBATCH --error=.out/%s.err\n" %(image_id))
            filey.writelines("#SBATCH --time=2-00:00\n")
            filey.writelines("#SBATCH --mem=64000\n")
            filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/ontological_comparison/cluster/classification-framework/3.calculate_reverse_inference.py %s %s %s %s" %(image, node, output_pkl, tables_folder))
            filey.close()
            os.system("sbatch -p russpold " + ".jobs/ri_%s.job" %(image_id))


# This will produce data for a LOO cross validation procedure:
# For each image node (defined with the likelihood group pickles above)
#     For each image: select him to leave out
#     With the remaining images, calculate likelihood tables, prior, and, RI for the query image
#     Condition A: With actual labels, do reverse inference procedure for correct/real tags. For each query image:
# We want to save a pickle/data object with:
# result["id"]: "trm*"
# result["prior"]: {"in": 0.98, "out": 0.02}
# result["query_RI"]: panda data frame with:
#               reverse_inference_score  bayes_factor
#  queryimage1
#  queryimage2
#  queryimage3

# Then for analysis:
# For each image node (defined with the likelihood group pickles above)
#     For each image: select him to leave out
#     Condition A: load scores of image for "correct/real" tags
#     Condition B: randomly select incorrect tags, load scores, compare
#     For each of the above, calculate a CE score, and compare distributions.
#     We would want "correct" tags to have higher scores
# This means that we need, for each node: to save a reverse inference score for ALL query images in the databases (against the node) and to save the priors values (not including the image) and then a reverse inference score.
