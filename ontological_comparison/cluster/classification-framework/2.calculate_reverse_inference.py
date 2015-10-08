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

from pybraincompare.ontology.inference import calculate_reverse_inference_distance
from pybraincompare.mr.datasets import get_standard_mask
from pybraincompare.compare.mrutils import get_images_df
import pickle
import pandas
import sys
import os

image = sys.argv[1]             # Full path of input image
node = sys.argv[2]              # pickle with groups for concept node
output_pkl = sys.argv[3]        # Path to save output file pickle

group = pickle.load(open(node,"rb"))
standard_mask = get_standard_mask()
image_id = os.path.split(image)[1].replace(".nii.gz","")

# Remove image from the in and out groups
in_group = [x for x in group["in"] if x != image]
out_group = [x for x in group["out"] if x != image]

# We will save to a result object
result = dict()

# Calculate reverse inference (posterior) for query image
# P(node mental process|activation) = P(activation|mental process) * P(mental process)
# divided by
# P(activation|mental process) * P(mental process) + P(A|~mental process) * P(~mental process)
# P(activation|mental process): my voxelwise prior map

# This is a reverse inference score, the p(cognitive process | query)
ri = calculate_reverse_inference_distance(image,in_group,out_group,standard_mask)
result["ri_query"] = ri
result["bayes_factor"] = ri / 0.5

# Save rest of varibles to result object
result["in_count"] = len(group["in"])
result["out_count"] = len(group["out"])
result["image"] = image
result["nid"] = group["nid"]
result["concept_node"] = node
result["image_id"] = image_id

# Save result to file
pickle.dump(result,open(output_pkl,"wb"))
