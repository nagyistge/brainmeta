#!/usr/bin/python
from glob import glob
import pickle
import pandas
import os

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON"
data = "%s/data" %base        # mostly images
priors_folder = "%s/priors" %(data)
scores_folder = "%s/scores" %(data)                  # output folder for node scores
indscores_folder = "%s/individual_scores" %(data)   # output folder for individual scores
scores = glob("%s/*.pkl" %scores_folder)

#images = glob("%s/*"%images_folder)

# Read in each score table, we will save to one master data frame
result = pickle.load(open(scores[0],"rb"))
column_names = [r for r in result.keys() if r not in ["range_table","single_scores","out","in","meta","scores_ranges"]]
scores_df = pandas.DataFrame(columns=column_names + ["count_in","count_out","threshold"] + result["scores_ranges"].index.tolist())
for score_file in scores:
    result = pickle.load(open(score_file,"rb")) 
    for column_name in column_names:
        if isinstance(result[column_name],pandas.core.series.Series):
            score = float(result[column_name])
        # Threshold is stored with "score_binary"
        if column_name == "scores_binary": 
            scores_df.loc[result["nid"],"threshold"] = result[column_name].index[0]
        else:
            score = result[column_name]
        scores_df.loc[result["nid"],column_name] = score
    scores_df.loc[result["nid"],"count_in"] =  len(result["in"])
    scores_df.loc[result["nid"],"count_out"] =  len(result["out"])
    scores_df.loc[result["nid"],result["scores_ranges"].index.tolist()] = result["scores_ranges"].tolist()

# Save table with result to file
scores_df.to_csv("%s/reverse_inference_scores_compiled.tsv" %data,sep="\t")
