#!/usr/bin/python
from glob import glob
import numpy
import pickle
import pandas
import os

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON/v2"
data = "%s/data" %base        # mostly images
scores_folder = "%s/individual_scores" %(data)   # output folder for individual scores
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

ri_score = pandas.DataFrame(columns=nodes,index=image_ids)
ri_bayes = pandas.DataFrame(columns=nodes,index=image_ids)
counts_in = pandas.DataFrame(columns=nodes,index=image_ids)
counts_out = pandas.DataFrame(columns=nodes,index=image_ids)

for s in range(0,len(scores)):
    print "Parsing %s of %s" %(s,len(scores))
    # Read in each score table, we will save to one master data frame
    result = pickle.load(open(scores[s],"rb"))
    image_id = result["image_id"]
    node = result["nid"]
    ri_score.loc[image_id,node] = result["ri_query"] 
    ri_bayes.loc[image_id,node] = result["bayes_factor"]
    counts_in.loc[image_id,node] = result["in_count"]
    counts_out.loc[image_id,node] = result["out_count"]
    
# Save tables with result to file
ri_score.to_csv("%s/reverse_inference_scores.tsv" %data,sep="\t")
ri_bayes.to_csv("%s/reverse_inference_bayes.tsv" %data,sep="\t")
counts_in.to_csv("%s/reverse_inference_counts_in.tsv" %data,sep="\t")
counts_out.to_csv("%s/reverse_inference_counts_out.tsv" %data,sep="\t")
