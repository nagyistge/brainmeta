#!/usr/bin/python
from glob import glob
import numpy
import pickle
import pandas
import os

base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON"
data = "%s/data" %base        # mostly images
priors_folder = "%s/priors" %(data)
scores_folder = "%s/scores" %(data)                  # output folder for node scores
indscores_folder = "%s/individual_scores" %(data)   # output folder for individual scores
scores = glob("%s/*.pkl" %scores_folder)

# Read in each score table, we will save to one master data frame
result = pickle.load(open(scores[0],"rb"))
column_names = [r for r in result.keys() if r not in ["range_table","single_scores","out","in","meta","scores_ranges"]]
scores_df = pandas.DataFrame(columns=column_names + ["count_in","count_out","threshold"] + result["scores_ranges"].index.tolist())

# We want to get a list of all unique images
images_list = []

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
    images_list = images_list + result["single_scores"].index.tolist()
    images_list = numpy.unique(images_list).tolist()

# Now get individual scores
image_ids = [int(os.path.split(i)[-1].replace(".nii.gz","")) for i in images_list]
individual_scores_binary = pandas.DataFrame(columns=image_ids)
individual_scores_ranges = pandas.DataFrame(columns=image_ids)
for score_file in scores:
    result = pickle.load(open(score_file,"rb"))
    images_at_node = result["single_scores"].index.tolist() 
    images_at_node = [int(os.path.split(i)[-1].replace(".nii.gz","")) for i in images_at_node]
    individual_scores_binary.loc[result["nid"],images_at_node] = result["single_scores"]["ri_binary"].tolist()
    individual_scores_ranges.loc[result["nid"],images_at_node] = result["single_scores"]["ri_ranges"].tolist()

# Save tables with result to file
scores_df.to_csv("%s/reverse_inference_scores_compiled.tsv" %data,sep="\t")
individual_scores_ranges.to_csv("%s/reverse_inference_scores_ind_ranges.tsv" %data,sep="\t")
individual_scores_binary.to_csv("%s/reverse_inference_scores_ind_binary.tsv" %data,sep="\t")

