#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 2 Masking: Assessing how variance changes with different levels of coverage (region based)

# This script will compile results from piecewise analyses
import os
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
from glob import glob

input_files = glob("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/atlas_scores_all/region*.pkl")

gs_file = glob("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/atlas_scores_all/gs*.pkl")
gs = pickle.load(open(gs_file[0],"rb"))

# We will save the masking strategy, number of rois, and size to the data frame
column_labels = gs.columns.tolist() + ["strategy","n_roi"]

# We will save all results to these lists
pearsons_pd = pandas.DataFrame(columns=column_labels)
pearsons_pi = pandas.DataFrame(columns=column_labels)
pearsons_bm = pandas.DataFrame(columns=column_labels)
pd_sizes = pandas.DataFrame(columns=column_labels)
pi_sizes = pandas.DataFrame(columns=column_labels)
bm_sizes = pandas.DataFrame(columns=column_labels)

for i in range(0,len(input_files)):
  input_file = input_files[i]
  "Processing %s of %s" %(i,len(input_files))
  result = pickle.load(open(input_file,"rb"))
  result["pd_sizes"]["strategy"] = "PD"
  result["pi_sizes"]["strategy"] = "PI"
  result["bm_sizes"]["strategy"] = "BM"
  result["pi_correlations"]["strategy"] = "PI"
  result["pd_correlations"]["strategy"] = "PD"
  result["bm_correlations"]["strategy"] = "BM"
  result["pd_sizes"]["n_roi"] = result["num_regions"]
  result["pi_sizes"]["n_roi"] = result["num_regions"]
  result["bm_sizes"]["n_roi"] = result["num_regions"]
  result["pi_correlations"]["n_roi"] = result["num_regions"]
  result["pd_correlations"]["n_roi"] = result["num_regions"]
  result["bm_correlations"]["n_roi"] = result["num_regions"]
  pearsons_pd = pearsons_pd.append(result["pd_correlations"])
  pearsons_pi = pearsons_pi.append(result["pi_correlations"])
  pearsons_bm = pearsons_bm.append(result["bm_correlations"])
  pd_sizes = pd_sizes.append(result["pd_sizes"])
  pi_sizes = pi_sizes.append(result["pi_sizes"])
  bm_sizes = bm_sizes.append(result["bm_sizes"])

pearsons = pandas.DataFrame()
pearsons = pearsons.append(pearsons_pd) 
pearsons = pearsons.append(pearsons_pi) 
pearsons = pearsons.append(pearsons_bm) 

sizes = pandas.DataFrame()
sizes = sizes.append(pd_sizes)
sizes = sizes.append(pi_sizes)
sizes = sizes.append(bm_sizes)

# Save all data matrices to file
os.chdir("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/atlas_scores_all/")
pearsons.to_csv("pearsons_all.tsv",sep="\t")
gs.to_csv("pearsons_gs.tsv",sep="\t")
sizes.to_csv("pearsons_sizes.tsv",sep="\t")
