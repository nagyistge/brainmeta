#!/usr/bin/python

# This first script will calculate threshold specific probabilities for an entire node,
# and then one score for all the mean image at the node, and all unique images below it
# Another script will be used to calculate reverse inference scores for all images at all nodes

from pybraincompare.ontology.inference import save_priors_df, calculate_reverse_inference, calculate_reverse_inference_threshes
from pybraincompare.mr.datasets import get_standard_mask
from pybraincompare.compare.mrutils import get_images_df
from glob import glob
import pickle
import pandas
import os
import re
import sys

input_image = sys.argv[1]    
tables_folder = sys.argv[2] # output folder for priors tables
scores_folder = sys.argv[3] # folder to output individual results to
priors_folder = sys.argv[4] # folder with all priors

# Load in the image
standard_mask = get_standard_mask()
image_df = get_images_df(input_image,standard_mask)
single_scores = pandas.DataFrame(columns=["ri_ranges","ri_binary","name","count_in","count_out"])

# We get the concepts and group sizes from the original priors pickle folder
priors_pickles = glob("%s/*.pkl" %priors_folder)

for p in range(0,len(priors_pickles)):
    print "Calculating RI for priors %s of %s" %(p,len(priors_pickles))
    group = pickle.load(open(priors_pickles[p],"rb"))
    nid = group["nid"]

    # Find the calculated priors tables
    priors_tables = glob("%s/*%s*" %(tables_folder,nid))
    priors_tables_ranges_out = [x for x in priors_tables if re.search("out_ranges",x)][0]
    priors_tables_ranges_in = [x for x in priors_tables if re.search("in_ranges",x)][0]
    priors_tables_binary = [x for x in priors_tables if x not in priors_tables_ranges_in + priors_tables_ranges_out]
    priors_tables_binary_out = [x for x in priors_tables_binary if re.search("df_out",x)][0]
    priors_tables_binary_in = [x for x in priors_tables_binary if re.search("df_in",x)][0]
    
    print "Priors tables ranges out is %s" % priors_tables_ranges_out
    print "Priors tables ranges in is %s" % priors_tables_ranges_in
    print "Priors tables binary out is %s" % priors_tables_binary_out
    print "Priors tables binary in is %s" % priors_tables_binary_in

    # Read in priors tables
    priors_in_ranges = pandas.read_pickle(priors_tables_ranges_in)
    priors_out_ranges = pandas.read_pickle(priors_tables_ranges_out)
    priors_in_bin = pandas.read_pickle(priors_tables_binary_in)
    priors_out_bin = pandas.read_pickle(priors_tables_binary_out)
    in_count = len(group["in"])
    out_count = len(group["out"])
    range_table = group["range_table"]

    # Use priors table, calculate a reverse inference score for the image
    ri_single_ranges = calculate_reverse_inference(image_df,priors_in_ranges,priors_out_ranges,in_count,out_count,range_table)
    ri_single_binary = calculate_reverse_inference(image_df,priors_in_bin,priors_out_bin,in_count,out_count)
    single_scores.loc[nid] = [ri_single_ranges,ri_single_binary,group["name"],in_count,out_count]

# Save all to file
image_id = int(os.path.split(input_image)[-1].replace(".nii.gz",""))
result = {}
result["single_scores"] = single_scores
result["image_id"] = image_id
result["image_file"] = input_image
pickle.dump(result,open("%s/%s_riscores.pkl" %(scores_folder,image_id),"wb"))
