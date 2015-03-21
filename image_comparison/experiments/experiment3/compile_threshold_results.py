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

input_folders = np.sort(glob("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/scores/*")).tolist()

# We will save all results to these lists
pearsons_pd = []  # pearsonr scores
pearsons_pi = []
pearsons_bm = []
spearmans_pd = [] # spearmanr scores
spearmans_pi = []
spearmans_bm = []
pd_sizes = []     # size of final masks
pi_sizes = []
bm_sizes = []
nanlog_pd = []    # "success","nan_fewer_3_values","nan_no_overlap"
nanlog_pi = []
nanlog_bm = []

for input_folder in input_folders:
  os.chdir(input_folder)
  inputs = glob("%s/*.pkl" %input_folder)
  thresh = os.path.split(input_folder)[1].split("_")[1]
  direction = os.path.split(input_folder)[1].split("_")[2]
  input_ids =  [i.replace("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/scores/thresh_%s_%s/" %(thresh,direction),"").replace("_scores_thresh_%s.pkl" %(thresh),"") for i in inputs]
  #groups = [x.split("_")[0] for x in input_ids]
  #tasks = [x.split("_")[1] for x in input_ids]
  #contrasts = [x.split("_")[2] for x in input_ids]
  column_labels = input_ids + ["thresh","direction","dof_mr1"]
  # Load the first to get image ids and order
  pearson_pd = pandas.DataFrame(columns=column_labels)
  pearson_pi = pandas.DataFrame(columns=column_labels)
  pearson_bm = pandas.DataFrame(columns=column_labels)
  spearman_pd = pandas.DataFrame(columns=column_labels)
  spearman_pi = pandas.DataFrame(columns=column_labels)
  spearman_bm = pandas.DataFrame(columns=column_labels)
  # We will also save matrices of size differences
  pd_size = pandas.DataFrame(columns=column_labels)
  pi_size = pandas.DataFrame(columns=column_labels)
  bm_size = pandas.DataFrame(columns=column_labels)
  # The nan logs keep track of if the comparison was successful, or if we appended nan
  nanlog_pd_single = pandas.DataFrame(columns=column_labels)
  nanlog_pi_single = pandas.DataFrame(columns=column_labels)
  nanlog_bm_single = pandas.DataFrame(columns=column_labels)
  # Each row corresponds to one image, there will be nan for images vs. themselves
  for i in inputs:
    print "Processing %s" %i
    input_id =  i.replace("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/scores/thresh_%s_%s/" %(thresh,direction),"").replace("_scores_thresh_%s.pkl" %(thresh),"")
    tmp = pickle.load(open(i,"rb"))
    pearson_pi.loc[input_id,tmp["uid"]] = tmp["mr_vs_thresh_pearson_pi"]
    pearson_pd.loc[input_id,tmp["uid"]] = tmp["mr_vs_thresh_pearson_pd"]
    pearson_bm.loc[input_id,tmp["uid"]] = tmp["mr_vs_thresh_pearson_bm"]
    spearman_pi.loc[input_id,tmp["uid"]] = tmp["mr_vs_thresh_spearman_pi"]
    spearman_pd.loc[input_id,tmp["uid"]] = tmp["mr_vs_thresh_spearman_pd"]
    spearman_bm.loc[input_id,tmp["uid"]] = tmp["mr_vs_thresh_spearman_bm"]
    # Save all mask sizes differences, again will be nan for image vs itself.
    pd_size.loc[input_id,tmp["size_ids"]] = tmp["sizes"]["pd"].tolist()
    pi_size.loc[input_id,tmp["size_ids"]] = tmp["sizes"]["pi"].tolist()
    bm_size.loc[input_id,tmp["size_ids"]] = tmp["sizes"]["bm"].tolist()
    # Nanlog - did we append a nan and why?
    nanlog_pd_single.loc[input_id,tmp["uid"]] = tmp["nanlog_pd"]
    nanlog_pi_single.loc[input_id,tmp["uid"]] = tmp["nanlog_pi"]
    nanlog_bm_single.loc[input_id,tmp["uid"]] = tmp["nanlog_bm"]
  # Add thresholding, directions, mr1 degrees of freedom
  pd_size["thresh"] = thresh; pd_size["direction"] = direction
  pi_size["thresh"] = thresh; pi_size["direction"] = direction
  bm_size["thresh"] = thresh; bm_size["direction"] = direction
  pearson_pd["thresh"] = thresh; pearson_pd["direction"] = direction
  pearson_bm["thresh"] = thresh; pearson_bm["direction"] = direction
  pearson_pi["thresh"] = thresh; pearson_pi["direction"] = direction
  spearman_pd["thresh"] = thresh; spearman_pd["direction"] = direction
  spearman_bm["thresh"] = thresh; spearman_bm["direction"] = direction
  spearman_pi["thresh"] = thresh; spearman_pi["direction"] = direction
  pearson_pd["mr_df"] = tmp["mr_df"] 
  pearson_bm["mr_df"] = tmp["mr_df"]
  pearson_pi["mr_df"] = tmp["mr_df"]
  spearman_pd["mr_df"] = tmp["mr_df"]
  spearman_bm["mr_df"] = tmp["mr_df"]
  spearman_pi["mr_df"] = tmp["mr_df"]
  # Append to lists
  pearsons_pd.append(pearson_pd)
  pearsons_pi.append(pearson_pi)
  pearsons_bm.append(pearson_bm)
  spearmans_pd.append(pearson_pd)
  spearmans_pi.append(pearson_pi)
  spearmans_bm.append(pearson_bm)
  pd_sizes.append(pd_size)
  pi_sizes.append(pi_size)
  bm_sizes.append(bm_size)
  nanlog_pd.append(nanlog_pd_single)
  nanlog_pi.append(nanlog_pi_single)
  nanlog_bm.append(nanlog_bm_single)

# Finally, merge them into one!
pbm = pandas.tools.merge.concat(pearsons_bm,axis=0)
ppi = pandas.tools.merge.concat(pearsons_pi,axis=0)
ppd = pandas.tools.merge.concat(pearsons_pd,axis=0)
sbm = pandas.tools.merge.concat(spearmans_bm,axis=0)
spi = pandas.tools.merge.concat(spearmans_pi,axis=0)
spd = pandas.tools.merge.concat(spearmans_pd,axis=0)
bmsizes = pandas.tools.merge.concat(bm_sizes,axis=0)
pdsizes = pandas.tools.merge.concat(pd_sizes,axis=0)
pisizes = pandas.tools.merge.concat(pi_sizes,axis=0)
nanlogpd = pandas.tools.merge.concat(nanlogpd,axis=0)
nanlogpi = pandas.tools.merge.concat(nanlogpi,axis=0)
nanlogbm = pandas.tools.merge.concat(nanlogbm,axis=0)

# Save all data matrices to file
os.chdir("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/scores/")
ppd.to_csv("860_pearson_pd.tsv",sep="\t")
ppi.to_csv("860_pearson_pi.tsv",sep="\t")
pbm.to_csv("860_pearson_bm.tsv",sep="\t")
spd.to_csv("860_spearman_pd.tsv",sep="\t")
spi.to_csv("860_spearman_pi.tsv",sep="\t")
sbm.to_csv("860_spearman_bm.tsv",sep="\t")
bmsizes.to_csv("860_bm_sizes.tsv",sep="\t")
pdsizes.to_csv("860_pd_sizes.tsv",sep="\t")
pisizes.to_csv("860_pi_sizes.tsv",sep="\t")
nanlogpd.to_csv("860_nanlog_pd.tsv",sep="\t")
nanlogpi.to_csv("860_nanlog_pi.tsv",sep="\t")
nanlogbm.to_csv("860_nanlog_bm.tsv",sep="\t")
