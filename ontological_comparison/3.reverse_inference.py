#!/usr/bin/python
### REVERSE INFERENCE

import pandas

# Read in images metadata
outfolder = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference"
images = pandas.read_csv("%s/contrast_defined_images.tsv" %outfolder,sep="\t")


## STEP 1: GENERATE TRIPLES DATA STRUCTURE
from cognitiveatlas.datastructure import concept_node_triples

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

data_structure = read_csv(output_triples_file,sep="\t")

## STEP 2: REVERSE INFERENCE WITH PYBRAINCOMPARE

# under development!

