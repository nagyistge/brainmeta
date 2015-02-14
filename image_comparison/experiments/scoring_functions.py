#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 1 Scoring:
# Scoring similarity metrics based on "gold standard" whole brain (unthresholded) metric
# This module also includes functions for better formatting data.

# IN PROGRESS

import sys
import os
import pandas
import pickle
import numpy as np
from glob import glob

# Run all scoring, including saving compiled dataframe, calculating gold standard, ranking
def run_scoring(image_id,input_folder,output_folder,num_images,num_thresh,metrics):

  # If we haven't compiled the metrics yet
  single_output_folder = "%s/%s" %(output_folder,image_id)

  if not os.path.exists(single_output_folder):
    save_compiled_df(image_folder=single_output_folder,output_folder=output_folder)

  thresholds = glob(single_output_folder)
  for thresh_file in thresholds:
    # For each result, calculate gold standard rankings 
    metric_df = pandas.read_csv(thresh_file,sep="\t")
    # If we have all expected results
    if check_result_count(num_images,num_thresh,df):
      # Calculate the gold standard rankings - the ordering for the images against unthresholded versions
      unthresholded_comparisons = []
      for row in metric_df.iterrows():
        if row[1].ID.split("_")[2] == "0.0":
          print "Found %s" % row[0]
          unthresholded_comparisons.append(row[1].ID)
          if len(unthresholded_comparisons != num_images: 
            print "ERROR %s does not have %s gold standard images!" %(image_id)
          else:
            gs_df = metric_df[metric_df.ID.isin(unthresholded_comparisons)]
            gs_rankings = get_gold_standard_rankings(gs_df,metrics) 

    else:
      print "ERROR %s does not have all results!" %(thresh_file)



# Save compiled df for each threshold of an image across all metrics and "other image" transforms
def save_compiled_df(image_folder,output_folder):
  results = glob("%s/*.tsv" %(image_folder))
  # Get the column headers for the similarity metrics
  column_names = pandas.read_csv(results[0],sep="\t").columns.tolist()
  # These will be new column names to add at the end
  column_names_new = pandas.read_csv(results[0],sep="\t").columns.tolist()
  column_names_new[0] = "ID"    
  # We will save a compiled tsv for each image_id_thresh [this is what we score]
  metric_df = pandas.DataFrame(columns=column_names)
  thresh = os.path.split(results[0])[1].split("_")[2]
  for res in results:
    # Split filepath into parts
    # ['000134', 'thr', '0.0', '541', 'pairwise', 'metrics.tsv']
    parts = os.path.split(res)[1].split("_")
    image1_id = int(parts[0])
    current_thresh = parts[2]
    # If we are working with a new threshold
    if current_thresh != thresh:
      metric_df.columns = column_names_new
      metric_df.to_csv("%s/%s_%s.tsv" %(output_folder,image1_id,thresh),sep="\t",index=False)
      metric_df = pandas.DataFrame(columns=column_names)
      thresh = current_thresh
    # Append similarity scores to current data frame
    tmp_df = pandas.read_csv(res,sep="\t")
    metric_df = metric_df.append(tmp_df)
  
# Check that we have #images by #thresholds results in a df
def check_result_count(num_images,num_thresh,df):
  if df.shape[0] != num_images*num_thresh: return False
  else: return True

# Calculate a gold standard ordering
def get_gold_standard_rankings(metric_df,metrics):
  df = metric_df[metrics]
  #TODO: need to visualize the distribution of metric scores for each so we know how to rank





  

