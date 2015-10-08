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

# Read in results files
# Save tables with result to file
ri_ranges = pandas.read_csv("%s/reverse_inference_scores_ranges.tsv" %data,sep="\t")
ri_binary = pandas.read_csv("%s/reverse_inference_scores_binary.tsv" %data,sep="\t")
ri_priors_in = pandas.read_csv("%s/reverse_inference_priors_in.tsv" %data,sep="\t")
ri_priors_out = pandas.read_csv("%s/reverse_inference_priors_out.tsv" %data,sep="\t")
bayes_in_ranges = pandas.read_csv("%s/reverse_inference_bayes_in_ranges" %data,sep="\t")
bayes_out_ranges = pandas.read_csv("%s/reverse_inference_bayes_out_ranges.tsv" %data,sep="\t")
bayes_in_bin = pandas.read_csv("%s/reverse_inference_bayes_in_binary.tsv" %data,sep="\t")
bayes_out_bin = pandas.read_csv("%s/reverse_inference_bayes_out_binary.tsv" %data,sep="\t")

count_bin = pandas.DataFrame(0,columns=["for","against"],index=bayes_in_bin.columns)
count_range = pandas.DataFrame(0,columns=["for","against"],index=bayes_in_bin.columns)

# For each likelihood pickle, read in the "in" and "out" groups, and take a count for when evidence for > evidence against
for i in range(0,len(likelihood_pickles)):
    print "Parsing %s of %s" %(i,len(likelihood_pickles))
    node = likelihood_pickles[i]
    group = pickle.load(open(node,"rb"))
    nid = group["nid"]
    # Look at bayes for range and bin given "in" group
    for image in group["in"]:
        image_id = os.path.split(image)[1].replace(".nii.gz","")
        # If we have data yet!
        if not numpy.isnan(ri_ranges.loc[image_id,nid]): 
            if bayes_in_ranges.loc[image_id,nid] > bayes_out_ranges.loc[image_id,nid]:
                count_range.loc[nid,"for"] = count_range.loc[nid,"for"] + 1
            else:
                count_range.loc[nid,"against"] = count_range.loc[nid,"against"] + 1
            if bayes_in_bin.loc[image_id,nid] > bayes_out_bin.loc[image_id,nid]:
                count_bin.loc[nid,"for"] = count_bin.loc[nid,"for"] + 1
            else:
                count_bin.loc[nid,"against"] = count_bin.loc[nid,"against"] + 1
    # Look at bayes for range and bin given "out" group
    for image in group["out"]:
        image_id = os.path.split(image)[1].replace(".nii.gz","")
        # If we have data yet!
        if not numpy.isnan(ri_ranges.loc[image_id,nid]): 
            if bayes_in_ranges.loc[image_id,nid] < bayes_out_ranges.loc[image_id,nid]:
                count_range.loc[nid,"for"] = count_range.loc[nid,"for"] + 1
            else:
                count_range.loc[nid,"against"] = count_range.loc[nid,"against"] + 1
            if bayes_in_bin.loc[image_id,nid] < bayes_out_bin.loc[image_id,nid]:
                count_bin.loc[nid,"for"] = count_bin.loc[nid,"for"] + 1
            else:
                count_bin.loc[nid,"against"] = count_bin.loc[nid,"against"] + 1

