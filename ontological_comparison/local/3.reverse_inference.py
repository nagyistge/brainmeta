#!/usr/bin/python
### REVERSE INFERENCE

import pandas

# Read in images metadata
outfolder = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference"
images = pandas.read_csv("%s/contrast_defined_images.tsv" %outfolder,sep="\t")


## STEP 1: GENERATE TRIPLES DATA STRUCTURE
from cognitiveatlas.datastructure import concept_node_triples, get_concept_categories

'''
  id    parent  name
  1 none BASE                   # there is always a base node
  2 1   MEMORY                  # high level concept groups
  3 1   PERCEPTION              
  4 2   WORKING MEMORY          # concepts
  5 2   LONG TERM MEMORY
  6 4   image1.nii.gz           # associated images (discovered by way of contrasts)
  7 4   image2.nii.gz
'''

# We need a dictionary to look up image lists by contrast ids
unique_contrasts = images.cognitive_contrast_cogatlas_id.unique().tolist()
image_lookup = dict()
for u in unique_contrasts:
   image_lookup[u] = images.image_id[images.cognitive_contrast_cogatlas_id==u].tolist()

output_triples_file = "%s/task_contrast_triples.tsv" % outfolder

# Create a data structure of tasks and contrasts for our analysis
concept_node_triples(image_dict=image_lookup,output_file=output_triples_file)
relationship_table = pandas.read_csv(output_triples_file,sep="\t")

# We want to give the concept categories as meta data so we produce category nodes
meta_data = get_concept_categories()




## STEP 2: VISUALIZATION WITH PYBRAINCOMPARE
from pybraincompare.ontology.tree import named_ontology_tree_from_tsv
from pybraincompare.template.visual import view

# First let's look at the tree structure
# output_json = "%s/task_contrast_tree.json" % outfolder
tree = named_ontology_tree_from_tsv(relationship_table,output_json=None,meta_data=meta_data)
html_snippet = make_ontology_tree_d3(tree)
view(html_snippet)



## STEP 3: REVERSE INFERENCE WITH PYBRAINCOMPARE
# The following steps should be run in a cluster environment
# this will show an example in a single batch script
from pybraincompare.ontology.inference import priors_groups_from_tree, save_priors_df, calculate_reverse_inference,
calculate_reverse_inference_threshes
from pybraincompare.compare.mrutils import get_images_df    

standard_mask = "/usr/share/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz"         # Brain mask
input_folder = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference/resampled" # Images folder
output_folder = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference/priors"
# Take a look at "image_pattern" and "node_pattern" inputs if not using NeuroVault and pybraincompare tree

###### 3.1 First generate priors groups
# If an output_folder is specified, will also save each to a pickle
priors_groups = priors_groups_from_tree(tree,standard_mask,input_folder)

###### 3.2 Calculate priors tables, and then reverse inference 
# P(node mental process|activation) = P(activation|mental process) * P(mental process)
# divided by
# P(activation|mental process) * P(mental process) + P(A|~mental process) * P(~mental process)
# P(activation|mental process): my voxelwise prior map

# If we want to make a custom range table (default is to use ranges defined by images)
range_table = make_range_table(mr,ranges=[[2.93,100]])

for group in priors_groups:
    range_table = group["range_table"]
    priors_df_files = save_priors_df(group["nid"],group["in"],group["out"],standard_mask,output_folder,range_table)
    priors_in = pandas.read_pickle(priors_df_files["in"])
    priors_out = pandas.read_pickle(priors_df_files["out"])
    in_count = len(group["in"])
    out_count = len(group["out"])
    scores = calculate_reverse_inference_threshes(priors_in,priors_out,in_count,out_count)
    # This is a reverse inference score, the p(cognitive process | activation in range [x1..x2])
    # p(cognitive process | an activation "level") based on an entire group of images with the tag

    # Now let's pretend we have a "query" image to calculate a score for, for now let's use the mean of our images  
    mrin = get_images_df(file_paths=group["in"],mask=standard_mask).mean()
    mrout = get_images_df(file_paths=group["out"],mask=standard_mask).mean()
    
    # This will use the priors tables as "lookups" to create a vector of probabilities (one per voxel)
    # matched to the appropriate probability [voxel,threshold] in the priors lookup tables. 
    # We can calculate a score using the "in" priors table (the images labeled with the concept)
    # and the "out" priors table (everything else) and could use this in some kind of classification framework
    ri_in = calculate_reverse_inference(mrin,range_table,priors_in,priors_out,in_count,out_count)
    ri_out = calculate_reverse_inference(mrout,range_table,priors_in,priors_out,in_count,out_count)

