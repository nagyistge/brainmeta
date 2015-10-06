#!/usr/bin/python

# This prepares data for a LOO cross validation procedure:
# For each image node (defined with the likelihood group pickles above)
#     For each image: select him to leave out
#     With the remaining images, calculate likelihood tables, prior, and, RI for the query image
#     Condition A: With actual labels, do reverse inference procedure for correct/real tags. For each query image:
# We want to save a pickle/data object with:
# result["id"]: "trm*"
# result["prior"]: {"in": 0.98, "out": 0.02}
# result["query_RI"]: panda data frame with:
#               reverse_inference_score  bayes_factor
#  queryimage1
#  queryimage2
#  queryimage3

from pybraincompare.ontology.inference import get_likelihood_df, calculate_reverse_inference, calculate_reverse_inference_threshes
from pybraincompare.mr.datasets import get_standard_mask
from pybraincompare.compare.mrutils import get_images_df
import pickle
import pandas
import sys
import os

image = sys.argv[1]             # Full path of input image
node = sys.argv[2]       # pickle with groups for concept node
output_pkl = sys.argv[3]        # Path to save output file pickle

group = pickle.load(open(node,"rb"))
standard_mask = get_standard_mask()

# Remove image from the in and out groups
in_group = [x for x in group["in"] if x != image]
out_group = [x for x in group["out"] if x != image]

# We will save to a result object
result = dict()

# Calculate likelihood tables
# P(node mental process|activation) = P(activation|mental process) * P(mental process)
# divided by
# P(activation|mental process) * P(mental process) + P(A|~mental process) * P(~mental process)
# P(activation|mental process): my voxelwise prior map
range_table = group["range_table"]
likelihood_dfs = get_likelihood_df(group["nid"],in_group,out_group,standard_mask,range_table)

# Calculate reverse inference (threshold specific) scores for the concept node
# This is the p(activation | mental process)
like_in_ranges = likelihood_dfs["in_ranges"]
like_out_ranges = likelihood_dfs["out_ranges"]
like_in_bin = likelihood_dfs["in_bin"]
like_out_bin = likelihood_dfs["out_bin"]
in_count = len(in_group)
out_count = len(out_group)

# The prior for each group is N/total
prior_in = in_count / float(in_count+out_count)
prior_out = out_count / float(in_count+out_count)
result["prior_in"] = prior_in
result["prior_out"] = prior_out

# Now calculate the reverse inference scores for the query image
# This is a reverse inference score, the p(cognitive process | activation in range [x1..x2])
# p(cognitive process | an activation "level")
image_df = get_images_df(file_paths=image,mask=standard_mask)
image_id = os.path.split(image)[1].replace(".nii.gz","")    

# REVERSE INFERENCE TO CLASSIFY IMAGES - USING THRESHOLDS
# When we provide the range_table: we use the likelihood tables as "lookups" to create a vector of probabilities
# (one per voxel) matched to the appropriate probability [voxel,threshold] in the priors lookup tables. 
# We can calculate a score using the "in" likelihood table (the images labeled with the concept)
# and the "out" likelihood table (everything else). When we generate this score for the query image and compare
# it to the prior, we can determine if the image adds new evidence for the concept (or not)
ri_ranges = calculate_reverse_inference(image_df,like_in_ranges,like_out_ranges,in_count,out_count,range_table)
result["ri_ranges"] = ri_ranges
result["bayes_factor_ranges_in"] = ri_ranges / 0.5
result["bayes_factor_ranges_out"] = ri_ranges / 0.5

# REVERSE INFERENCE TO CLASSIFY IMAGES - USING BINARY ACTIVATION THRESHOLD
# When we don't provide range table, threshold is extracted from priors_* column name
ri_bin = calculate_reverse_inference(image_df,like_in_bin,like_out_bin,in_count,out_count)
result["ri_binary"] = ri_bin
result["bayes_factor_bin_in"] = ri_bin / 0.5
result["bayes_factor_bin_out"] = ri_bin / 0.5

# Save rest of varibles to result object
result["in_count"] = in_count
result["out_count"] = out_count
result["ranges_table"] = range_table
result["image"] = image
result["nid"] = group["nid"]
result["concept_node"] = node
result["image_id"] = image_id

# Save result to file
pickle.dump(result,open(output_pkl,"wb"))
