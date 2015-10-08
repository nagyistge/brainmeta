#!/usr/bin/python
from glob import glob
import numpy
import pickle
import pandas
import os

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON"
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

ri_ranges = pandas.DataFrame(columns=nodes,index=image_ids)
ri_binary = pandas.DataFrame(columns=nodes,index=image_ids)
ri_priors_in = pandas.DataFrame(columns=nodes,index=image_ids)
ri_priors_out = pandas.DataFrame(columns=nodes,index=image_ids)
bayes_in_ranges = pandas.DataFrame(columns=nodes,index=image_ids)
bayes_out_ranges = pandas.DataFrame(columns=nodes,index=image_ids)
bayes_in_bin = pandas.DataFrame(columns=nodes,index=image_ids)
bayes_out_bin = pandas.DataFrame(columns=nodes,index=image_ids)

for s in range(0,len(scores)):
    print "Parsing %s of %s" %(s,len(scores))
    # Read in each score table, we will save to one master data frame
    result = pickle.load(open(scores[s],"rb"))
    range_table = result["ranges_table"]    
    image_id = result["image_id"]
    node = result["nid"]
    ri_ranges.loc[image_id,node] = result["ri_ranges"] 
    ri_binary.loc[image_id,node] = result["ri_binary"]
    ri_priors_in.loc[image_id,node] = result["prior_in"]
    ri_priors_out.loc[image_id,node] = result["prior_out"]
    bayes_in_ranges.loc[image_id,node] = result["bayes_factor_ranges_in"]
    bayes_out_ranges.loc[image_id,node] = result["bayes_factor_ranges_out"]
    bayes_in_bin.loc[image_id,node] = result["bayes_factor_bin_in"]
    bayes_out_bin.loc[image_id,node] = result["bayes_factor_ranges_out"]

# Save tables with result to file
ri_ranges.to_csv("%s/reverse_inference_scores_ranges.tsv" %data,sep="\t")
ri_binary.to_csv("%s/reverse_inference_scores_binary.tsv" %data,sep="\t")
ri_priors_in.to_csv("%s/reverse_inference_priors_in.tsv" %data,sep="\t")
ri_priors_out.to_csv("%s/reverse_inference_priors_out.tsv" %data,sep="\t")
bayes_in_ranges.to_csv("%s/reverse_inference_bayes_in_ranges" %data,sep="\t")
bayes_out_ranges.to_csv("%s/reverse_inference_bayes_out_ranges.tsv" %data,sep="\t")
bayes_in_bin.to_csv("%s/reverse_inference_bayes_in_binary.tsv" %data,sep="\t")
bayes_out_bin.to_csv("%s/reverse_inference_bayes_out_binary.tsv" %data,sep="\t")
