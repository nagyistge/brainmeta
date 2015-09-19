#!/usr/bin/python

# This first script will calculate threshold specific probabilities for an entire node,
# and then one score for all the mean image at the node, and all unique images below it
# Another script will be used to calculate reverse inference scores for all images at all nodes

from pybraincompare.ontology.inference import save_priors_df, calculate_reverse_inference, calculate_reverse_inference_threshes
from pybraincompare.mr.datasets import get_standard_mask
from pybraincompare.compare.mrutils import get_images_df
import pickle
import pandas

priors_pkl = sys.argv[1]    # input pickle with groups
tables_folder = sys.argv[2] # output folder for priors tables
scores_folder = sys.argv[3] # folder to output results to
                            # "scores" subdirectory will be 

group = pickle.load(open(priors_pkl,"rb"))
standard_mask = get_standard_mask()

###### 3.2 Calculate priors tables, and then reverse inference 
# P(node mental process|activation) = P(activation|mental process) * P(mental process)
# divided by
# P(activation|mental process) * P(mental process) + P(A|~mental process) * P(~mental process)
# P(activation|mental process): my voxelwise prior map

# Save priors tables to file
range_table = group["range_table"]
priors_df_files = save_priors_df(group["nid"],group["in"],group["out"],standard_mask,tables_folder,range_table)

# Calculate reverse inference (threshold specific) scores for the concept node
priors_in_ranges = pandas.read_pickle(priors_df_files["in_ranges"])
priors_out_ranges = pandas.read_pickle(priors_df_files["out_ranges"])
priors_in_bin = pandas.read_pickle(priors_df_files["in_bin"])
priors_out_bin = pandas.read_pickle(priors_df_files["out_bin"])
in_count = len(group["in"])
out_count = len(group["out"])
scores_ranges = calculate_reverse_inference_threshes(priors_in_ranges,priors_out_ranges,in_count,out_count)
scores_bin = calculate_reverse_inference_threshes(priors_in_bin,priors_out_bin,in_count,out_count)

# This is a reverse inference score, the p(cognitive process | activation in range [x1..x2])
# p(cognitive process | an activation "level") based on an entire group of images with the tag

# Add to group data structure
group["scores_ranges"] = scores_ranges
group["scores_binary"] = scores_bin

# Now let's pretend we have a "query" image to calculate a score for, the default takes the mean of a group of images
mrin = get_images_df(file_paths=group["in"],mask=standard_mask)
mrout = get_images_df(file_paths=group["out"],mask=standard_mask)
    
# REVERSE INFERENCE TO CLASSIFY IMAGES - USING THRESHOLDS
# When we provide the range_table: we use the priors tables as "lookups" to create a vector of probabilities
# (one per voxel) matched to the appropriate probability [voxel,threshold] in the priors lookup tables. 
# We can calculate a score using the "in" priors table (the images labeled with the concept)
# and the "out" priors table (everything else), and calculate a bayes factor to determine if they are different.
ri_in_range = calculate_reverse_inference(mrin,priors_in_ranges,priors_out_ranges,in_count,out_count,range_table)
ri_out_range = calculate_reverse_inference(mrout,priors_in_ranges,priors_out_ranges,in_count,out_count,range_table)
group["ri_in_ranges"] = ri_in_range
group["ri_out_ranges"] = ri_out_range
group["bayes_factor_ranges"] = ri_in_range/ri_out_range


# REVERSE INFERENCE TO CLASSIFY IMAGES - USING BINARY ACTIVATION THRESHOLD
# When we don't provide range table, threshold is extracted from priors_* column name
ri_in_bin = calculate_reverse_inference(mrin,priors_in_bin,priors_out_bin,in_count,out_count)
ri_out_bin = calculate_reverse_inference(mrout,priors_in_bin,priors_out_bin,in_count,out_count)
group["ri_in_binary"] = ri_in_bin
group["ri_out_binary"] = ri_out_bin
group["bayes_factor_binary_%s" %(priors_in_bin.columns[0])] = ri_in_bin/ri_out_bin

# For each image, calculate a reverse inference score
single_scores = pandas.DataFrame(columns=["ri_ranges","ri_binary"])
for image in group["in"] + group["out"]:
    image_df = get_images_df(file_paths = image,mask=standard_mask)
    ri_single_ranges = calculate_reverse_inference(image_df,priors_in_ranges,priors_out_ranges,in_count,out_count,range_table)
    ri_single_binary = calculate_reverse_inference(image_df,priors_in_bin,priors_out_bin,in_count,out_count)
    single_scores.loc[image] = [ri_single_ranges,ri_single_binary]

group["single_scores"] = single_scores

# Save all to file
pickle.dump(group,open("%s/%s_riscores.pkl" %(scores_folder,group["nid"]),"wb"))
