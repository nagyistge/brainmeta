#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 2 Masking:

# This script will compile results from piecewise analyses
import os
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
from glob import glob

input_folders = glob("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores_all/*")

# We will save all results to these lists
pearsons_gs = []
pearsons_pd = []
pearsons_pi = []
pearsons_bm = []
pd_sizes = []
pi_sizes = []
bm_sizes = []

for input_folder in input_folders:
  inputs = glob("%s/*.pkl" %input_folder)
  thresh = os.path.split(input_folder)[1].split("_")[1]
  only_pos = os.path.split(input_folder)[1].split("_")[2]
  input_ids =  [int(i.replace("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores_all/thresh_%s_%s/000" %(thresh,only_pos),"").replace("_masking_scores_thresh_%s.pkl" %(thresh),"")) for i in inputs]
  input_ids.sort()
  column_labels = input_ids + ["thresh","pos_only"]
  # Load the first to get image ids and order
  pearson_gs = pandas.DataFrame(columns=column_labels)
  pearson_pd = pandas.DataFrame(columns=column_labels)
  pearson_pi = pandas.DataFrame(columns=column_labels)
  pearson_bm = pandas.DataFrame(columns=column_labels)
  # We will also save matrices of size differences
  pd_size = pandas.DataFrame(columns=column_labels)
  pi_size = pandas.DataFrame(columns=column_labels)
  bm_size = pandas.DataFrame(columns=column_labels)
  # Each row corresponds to one image, there will be nan for images vs. themselves
  for i in inputs:
    print "Processing %s" %i
    input_id = int(i.replace("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores_all/thresh_%s_%s/000" %(thresh,only_pos),"").replace("_masking_scores_thresh_%s.pkl" %(thresh),""))
    tmp = pickle.load(open(i,"rb"))
    pearson_gs.loc[input_id,tmp["ids"]] = tmp["pearson_gs"]
    pearson_pi.loc[input_id,tmp["ids"]] = tmp["mr_vs_thresh_pearson_pi"]
    pearson_pd.loc[input_id,tmp["ids"]] = tmp["mr_vs_thresh_pearson_pd"]
    pearson_bm.loc[input_id,tmp["ids"]] = tmp["mr_vs_thresh_pearson_bm"]
    # Save all mask sizes differences, again will be nan for image vs itself.
    pd_size.loc[input_id,tmp["sizes"].index.tolist()] = tmp["sizes"]["pd"].tolist()
    pi_size.loc[input_id,tmp["sizes"].index.tolist()] = tmp["sizes"]["pi"].tolist()
    bm_size.loc[input_id,tmp["sizes"].index.tolist()] = tmp["sizes"]["bm"].tolist()
  # Add thresholding and directions
  pd_size["thresh"] = thresh; pd_size["pos_only"] = only_pos
  pi_size["thresh"] = thresh; pi_size["pos_only"] = only_pos
  bm_size["thresh"] = thresh; bm_size["pos_only"] = only_pos
  pearson_gs["thresh"] = thresh; pearson_gs["pos_only"] = only_pos
  pearson_pd["thresh"] = thresh; pearson_pd["pos_only"] = only_pos
  pearson_bm["thresh"] = thresh; pearson_bm["pos_only"] = only_pos
  pearson_pi["thresh"] = thresh; pearson_pi["pos_only"] = only_pos
  # Append to lists
  pearsons_gs.append(pearsons_gs)
  pearsons_pd.append(pearsons_pd)
  pearsons_pi.append(pearsons_pi)
  pearsons_bm.append(pearsons_bm)
  pd_sizes.append(pd_size)
  pi_sizes.append(pi_size)
  bm_sizes.append(bm_size)

# Finally, merge them into one!

# Save all data matrices to file
pearsons_gs.to_csv("144_masking_gs_%s.tsv" %(thresh),sep="\t")
pearsons_pd.to_csv("144_masking_pd_%s.tsv" %(thresh),sep="\t")
pearsons_pi.to_csv("144_masking_pi_%s.tsv" %(thresh),sep="\t")
pearsons_bm.to_csv("144_masking_bm_%s.tsv" %(thresh),sep="\t")
pd_vs_bm.to_csv("144_pd_vs_bm_sizediff_%s.tsv" %(thresh),sep="\t")
pi_vs_bm.to_csv("144_pi_vs_bm_sizediff_%s.tsv" %(thresh),sep="\t")
pd_vs_pi.to_csv("144_pd_vs_pi_sizediff_%s.tsv" %(thresh),sep="\t")
