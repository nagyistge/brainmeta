#!/usr/bin/python

# This first script will calculate threshold specific probabilities for an entire node,
# and then one score for all the mean image at the node, and all unique images below it
# Another script will be used to calculate reverse inference scores for all images at all nodes

from pybraincompare.ontology.inference import save_likelihood_df, calculate_reverse_inference, calculate_reverse_inference_threshes
from pybraincompare.mr.datasets import get_standard_mask
from pybraincompare.compare.mrutils import get_images_df
import pickle
import pandas
import sys

likelihood_pkl = sys.argv[1]    # input pickle with groups
tables_folder = sys.argv[2]     # output folder for priors tables

group = pickle.load(open(likelihood_pkl,"rb"))
standard_mask = get_standard_mask()

###### 3.2 Calculate likelihood tables
# P(node mental process|activation) = P(activation|mental process) * P(mental process)
# divided by
# P(activation|mental process) * P(mental process) + P(A|~mental process) * P(~mental process)
# P(activation|mental process): my voxelwise prior map

# Save likelihood tables to file
range_table = group["range_table"]
likelihood_df_files = save_likelihood_df(group["nid"],group["in"],group["out"],standard_mask,tables_folder,range_table)
