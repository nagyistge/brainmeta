#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 2 Handling Missing Values:

# This script will compile results from piecewise analyses, likely needs big memory node
import os
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
from glob import glob

top_directory = "/work/02092/vsochat/wrangler/DATA/BRAINMETA/IMAGE_COMPARISON/experiments/experiment3"
inputs_directory = "%s/permutations" %(top_directory)
output_directory = "%s/results" %(top_directory)
nruns = 500

for r in range(0,nruns):
    inputs = np.sort(glob("%s/%s/comparisons*.pkl" %(inputs_directory,r))).tolist()
    input_ids = ["%s_%s" %(r,x.replace(".pkl","")) for x in inputs]
    # We will save all results to these lists
    pearsons_pi = []; pearsons_pd = []   # pearsonr scores
    spearmans_pi = []; spearmans_pd = [] # spearmanr scores
    pd_sizes = []; pi_sizes = []         # size of final masks
    nanlog_pd = []; nanlog_pi = []       # "success","nan_fewer_3_values","nan_no_overlap"
    for i in inputs:
        tmp = pickle.load(open(i,"rb"))
        input_id = "%s_%s" %(r,i.replace(".pkl",""))
        column_labels = input_ids + ["thresh","direction","dof_mr1"]
        # Load the first to get image ids and order
        pearson_pd = pandas.DataFrame(columns=column_labels)
        pearson_pi = pandas.DataFrame(columns=column_labels)
        spearman_pd = pandas.DataFrame(columns=column_labels)
        spearman_pi = pandas.DataFrame(columns=column_labels)
        # We will also save matrices of size differences
        pd_size = pandas.DataFrame(columns=column_labels)
        pi_size = pandas.DataFrame(columns=column_labels)
        # The nan logs keep track of if the comparison was successful, or if we appended nan
        nanlog_pd_single = pandas.DataFrame(columns=column_labels)
        nanlog_pi_single = pandas.DataFrame(columns=column_labels)
        # Each row corresponds to one image, there will be nan for images vs. themselves
        uids = [x for x in tmp["uid"] if x in column_labels]
        size_ids = [x for x in tmp["size_ids"] if x in column_labels]
        uid_index = [tmp["uid"].index(x) for x in uids]
        size_index = [tmp["size_ids"].index(x) for x in size_ids]
        pearson_pi.loc[input_id,uids] = [tmp["svi_pearson"][x] for x in uid_index]
        pearson_pd.loc[input_id,uids] = [tmp["cca_pearson"][x] for x in uid_index]
        spearman_pi.loc[input_id,uids] = [tmp["svi_spearman"][x] for x in uid_index]
        spearman_pd.loc[input_id,uids] = [tmp["cca_spearman"][x] for x in uid_index]
        # Save all mask sizes differences, again will be nan for image vs itself.
        pd_size.loc[input_id,size_ids] = [tmp["sizes"]["cca"].tolist()[x] for x in size_index]
        pi_size.loc[input_id,size_ids] = [tmp["sizes"]["svi"].tolist()[x] for x in size_index]
        # Nanlog - did we append a nan and why?
        nanlog_pd_single.loc[input_id,uids] = [tmp["nanlog_cca"][x] for x in uid_index]
        nanlog_pi_single.loc[input_id,uids] = [tmp["nanlog_svi"][x] for x in uid_index]
        # Add thresholding, directions, mr1 degrees of freedom
  

  pd_size["thresh"] = thresh; pd_size["direction"] = direction
  pi_size["thresh"] = thresh; pi_size["direction"] = direction
  #bm_size["thresh"] = thresh; bm_size["direction"] = direction
  nanlog_pd_single["thresh"] = thresh; nanlog_pd_single["direction"] = direction
  nanlog_pi_single["thresh"] = thresh; nanlog_pi_single["direction"] = direction
  #nanlog_bm_single["thresh"] = thresh; nanlog_bm_single["direction"] = direction
  pearson_pd["thresh"] = thresh; pearson_pd["direction"] = direction
  #pearson_bm["thresh"] = thresh; pearson_bm["direction"] = direction
  pearson_pi["thresh"] = thresh; pearson_pi["direction"] = direction
  spearman_pd["thresh"] = thresh; spearman_pd["direction"] = direction
  #spearman_bm["thresh"] = thresh; spearman_bm["direction"] = direction
  spearman_pi["thresh"] = thresh; spearman_pi["direction"] = direction
  pearson_pd["mr_df"] = tmp["mr_df"] 
  #pearson_bm["mr_df"] = tmp["mr_df"]
  pearson_pi["mr_df"] = tmp["mr_df"]
  spearman_pd["mr_df"] = tmp["mr_df"]
  #spearman_bm["mr_df"] = tmp["mr_df"]
  spearman_pi["mr_df"] = tmp["mr_df"]
  # Append to lists
  pearsons_pd.append(pearson_pd)
  pearsons_pi.append(pearson_pi)
  #pearsons_bm.append(pearson_bm)
  spearmans_pd.append(spearman_pd)
  spearmans_pi.append(spearman_pi)
  #spearmans_bm.append(spearman_bm)
  pd_sizes.append(pd_size)
  pi_sizes.append(pi_size)
  #bm_sizes.append(bm_size)
  nanlog_pd.append(nanlog_pd_single)
  nanlog_pi.append(nanlog_pi_single)
  #nanlog_bm.append(nanlog_bm_single)

# Finally, merge them into one!
#pbm = pandas.tools.merge.concat(pearsons_bm,axis=0)
ppi = pandas.tools.merge.concat(pearsons_pi,axis=0)
ppd = pandas.tools.merge.concat(pearsons_pd,axis=0)
#sbm = pandas.tools.merge.concat(spearmans_bm,axis=0)
spi = pandas.tools.merge.concat(spearmans_pi,axis=0)
spd = pandas.tools.merge.concat(spearmans_pd,axis=0)
#bmsize = pandas.tools.merge.concat(bm_sizes,axis=0)
pdsize = pandas.tools.merge.concat(pd_sizes,axis=0)
pisize = pandas.tools.merge.concat(pi_sizes,axis=0)
nanlogpd = pandas.tools.merge.concat(nanlog_pd,axis=0)
nanlogpi = pandas.tools.merge.concat(nanlog_pi,axis=0)
#nanlogbm = pandas.tools.merge.concat(nanlog_bm,axis=0)

# Save all data matrices to file
os.chdir(output_directory)
ppd.to_csv("860_pearson_cca.tsv",sep="\t")
ppi.to_csv("860_pearson_svi.tsv",sep="\t")
#pbm.to_csv("860_pearson_bm_v2.tsv",sep="\t")
spd.to_csv("860_spearman_cca.tsv",sep="\t")
spi.to_csv("860_spearman_svi.tsv",sep="\t")
#sbm.to_csv("860_spearman_bm_v2.tsv",sep="\t")
#bmsize.to_csv("860_bm_sizes_v2.tsv",sep="\t")
pdsize.to_csv("860_cca_sizes.tsv",sep="\t")
pisize.to_csv("860_svi_sizes.tsv",sep="\t")
nanlogpd.to_csv("860_nanlog_cca.tsv",sep="\t")
nanlogpi.to_csv("860_nanlog_svi.tsv",sep="\t")
nanlogbm.to_csv("860_nanlog_bm_v2.tsv",sep="\t")
