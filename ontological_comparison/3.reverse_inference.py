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

from pybraincompare.ontology.inference import calculate_reverse_inference
from pybraincompare.ontology.tree import named_ontology_tree_from_tsv, 
from pybraincompare.template.visual import view

# First let's look at the tree structure
output_json = "%s/task_contrast_tree.json" % outfolder
tree = named_ontology_tree_from_tsv(relationship_table,output_json=None,meta=meta_data)
html_snippet = make_ontology_tree_d3(tree)
view(html_snippet)

## STEP 3: REVERSE INFERENCE WITH PYBRAINCOMPARE

# Here is top path for images defined above (resampled to MNI 2mm)
mrpath = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference/resampled"

# IN PROGRESS
