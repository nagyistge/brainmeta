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

# For each threshold
thresholds = [1.96,2.58,3.02]

for thresh in thresholds:
  inputs = glob("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores/thresh_%s/*.pkl" %(thresh))

  # Load the first to get image ids and order
  tmp = pickle.load(open(inputs[0],"rb"))
  pearsons_gs = pandas.DataFrame(columns=tmp["ids"])
  pearsons_pd = pandas.DataFrame(columns=tmp["ids"])
  pearsons_pi = pandas.DataFrame(columns=tmp["ids"])
  pearsons_bm = pandas.DataFrame(columns=tmp["ids"])

  # We will also save matrices of size differences
  pd_vs_bm = pandas.DataFrame(columns=tmp["ids"])
  pi_vs_bm = pandas.DataFrame(columns=tmp["ids"])
  pd_vs_pi = pandas.DataFrame(columns=tmp["ids"])

  for i in inputs:
    print "Processing %s" %i
    input_id = int(i.replace("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores/thresh_%s/000" %(thresh),"").replace("_masking_scores_thresh_%s.pkl" %(thresh),""))
    tmp = pickle.load(open(i,"rb"))
    pearsons_gs.loc[input_id] = tmp["pearson_gs"]
    pearsons_pi.loc[input_id] = tmp["mr_vs_thresh_pearson_pi"]
    pearsons_pd.loc[input_id] = tmp["mr_vs_thresh_pearson_pd"]
    pearsons_bm.loc[input_id] = tmp["mr_vs_thresh_pearson_bm"]
    # Calculate size differences
    pd_vs_bm.loc[input_id] = abs(tmp["sizes"]["pd"] - tmp["sizes"]["bm"]).tolist()
    pd_vs_pi.loc[input_id] = abs(tmp["sizes"]["pd"] - tmp["sizes"]["pi"]).tolist()
    pi_vs_bm.loc[input_id] = abs(tmp["sizes"]["pi"] - tmp["sizes"]["bm"]).tolist()

  # Save all data matrices to file
  pearsons_gs.to_csv("144_masking_gs_%s.tsv" %(thresh),sep="\t")
  pearsons_pd.to_csv("144_masking_pd_%s.tsv" %(thresh),sep="\t")
  pearsons_pi.to_csv("144_masking_pi_%s.tsv" %(thresh),sep="\t")
  pearsons_bm.to_csv("144_masking_bm_%s.tsv" %(thresh),sep="\t")
  pd_vs_bm.to_csv("144_pd_vs_bm_sizediff_%s.tsv" %(thresh),sep="\t")
  pi_vs_bm.to_csv("144_pi_vs_bm_sizediff_%s.tsv" %(thresh),sep="\t")
  pd_vs_pi.to_csv("144_pd_vs_pi_sizediff_%s.tsv" %(thresh),sep="\t")
