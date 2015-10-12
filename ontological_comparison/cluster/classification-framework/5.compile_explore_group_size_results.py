#!/usr/bin/python
from glob import glob
import numpy
import pickle
import pandas
import os

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON/v2"
data = "%s/data" %base        # mostly images
scores_folder = "%s/group_size_vary_scores" %(data)     # output folder for group size scores
likelihood_pickles = glob("%s/likelihood/*.pkl" %(data))
scores = glob("%s/*.pkl" %scores_folder)

nodes = []

# First make a list of all the nodes, images
for i in range(0,len(likelihood_pickles)):
    node = likelihood_pickles[i]
    group = pickle.load(open(node,"rb"))
    all_images = group["in"] + group["out"]
    nodes.append(group["nid"])

# Parse image IDS
image_ids = [os.path.split(x)[1].replace(".nii.gz","") for x in all_images]

# Let's make a flat data frame this time
ri_scores = pandas.DataFrame(columns=["image_id","node","ri_score","bayes_factor","in_count","out_count"])

for s in range(0,len(scores)):
    print "Parsing %s of %s" %(s,len(scores))
    r = pickle.load(open(scores[s],"rb"))
    result_id = scores[s].replace(".pkl","")
    ri_scores.loc[result_id] = [r["image_id"],r["nid"],r["ri_query"],r["bayes_factor"],r["in_count"],r["out_count"]]
    
# Save tables with result to file
ri_scores.to_csv("%s/ri_explore_size_results.tsv" %data,sep="\t")
