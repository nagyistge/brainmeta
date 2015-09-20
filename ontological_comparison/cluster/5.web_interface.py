#!/usr/bin/python
### FINAL WEB REPORT

from cognitiveatlas.datastructure import concept_node_triples, get_concept_categories
from cognitiveatlas.api import get_concept
from vm import make_analysis_web_folder
from glob import glob
import pandas
import json
import os
import re

# For the VM: these paths will be environmental variables
base = "/share/PI/russpold/work/IMAGE_COMPARISON/ONTOLOGICAL_COMPARISON"
results = "%s/results" %base  # any kind of tsv/result file
data = "%s/data" %base        # mostly images
web = "%s/web" %base

### Step 1: Load meta data sources

# We will use image meta data
images = pandas.read_csv("%s/contrast_defined_images.tsv" %results,sep="\t")
collections = pandas.read_csv("%s/collections_with_dois.tsv" %results,sep="\t")
output_triples_file = "%s/task_contrast_triples.tsv" % results
relationship_table = pandas.read_csv(output_triples_file,sep="\t")

# We want to give the concept categories as meta data so we produce category nodes
categories = get_concept_categories()

# Get reverse inference scores from results
scores_df = pandas.read_csv("%s/reverse_inference_scores_compiled.tsv" %data,sep="\t",index_col="nid")
singles_ranges = pandas.read_csv("%s/reverse_inference_scores_ind_ranges.tsv" %data,sep="\t",index_col=0)
singles_binary = pandas.read_csv("%s/reverse_inference_scores_ind_binary.tsv" %data,sep="\t",index_col=0)

unique_nodes = relationship_table.id.unique().tolist()

# Output entire results table to html, in case we want it
scores_df.to_html("%s/reverse_inference_table.html" %web)

# We will store a data frame of meta data
# Lookup for meta_data is the id of the node!
meta_data = {}

for node in unique_nodes:
    meta_single = {}
    # This is an image node
    if re.search("node_",node):
        print "Found image node!"
        relationship_table_row = relationship_table[relationship_table.id==node]
        image_id = relationship_table_row.name.tolist()[0]
        # Reverse inference scores
        meta_single["category"] = ""
        meta_single["type"] = "nii"
        # NeuroVault metadata
        concepts = relationship_table.parent[relationship_table.name == image_id]
        concepts = [relationship_table.name[relationship_table.id==c].tolist()[0] for c in concepts]
        neurovault_row = images[images.image_id == int(image_id)]
        collection_row = collections[collections.collection_id == neurovault_row.collection.tolist()[0]]
        collection_meta = {"DOI":collection_row["DOI"].tolist()[0],
                           "authors":collection_row["authors"].tolist()[0],
                           "journal":collection_row["journal_name"].tolist()[0],
                           "url":collection_row["url"].tolist()[0],
                           "subjects":collection_row["number_of_subjects"].tolist()[0],
                           "smoothing_fwhm":str(collection_row["smoothing_fwhm"].tolist()[0]).encode("utf-8")}
        meta_single["url"] = neurovault_row["url"].tolist()[0]
        meta_single["thumbnail"] = neurovault_row["thumbnail"].tolist()[0]
        meta_single["images"] = neurovault_row["thumbnail"].tolist()
        meta_single["task"] = neurovault_row["cognitive_paradigm_cogatlas"].tolist()[0]
        meta_single["contrast"] = neurovault_row["cognitive_contrast_cogatlas"].tolist()[0]
        meta_single["download"] = neurovault_row["file"].tolist()[0]
        meta_single["concept"] = concepts
        if neurovault_row["description"].tolist()[0]:
            meta_single["description"] =  str(neurovault_row["description"].tolist()[0]).encode("utf-8")
        else:
            meta_single["description"] = ""
    else: # A node
        if node != "1":
            relationship_table_row = relationship_table[relationship_table.id==node]
            contrast_name = relationship_table_row.name.tolist()[0]
            concept = get_concept(id=node).json
            # Reverse inference scores - all images
            if node in singles_ranges.index: # a node with images below it
                meta_single["ri_range_scores"] = singles_ranges.loc[node,:].to_json()
                meta_single["ri_binary_scores"] = singles_binary.loc[node,:].to_json()
                image_ids =[int(x) for x in singles_binary.loc["trm_557b495cdde57",:].index.tolist()]
                meta_single["images"] = images["thumbnail"][images.image_id.isin(image_ids)].tolist()
            # Reverse inference scores - single
            if node in scores_df.index:
                node_scores = scores_df.loc[node]
                ri_scores = {}
                ri_scores["count_in"] = str(node_scores["count_in"])
                ri_scores["count_out"] = str(node_scores["count_out"])
                ri_scores["ri_score_binary"] = str(node_scores["scores_binary"])
                ri_scores["ri_binary_threshold"] = str(node_scores["threshold"])
                ri_scores["ri_range_score_in"] = str(node_scores["ri_in_ranges"])
                ri_scores["ri_range_score_out"] = str(node_scores["ri_out_ranges"])
                ri_scores["ri_binary_score_in"] = str(node_scores["ri_in_binary"])
                ri_scores["ri_binary_score_out"] = str(node_scores["ri_out_binary"])
                ri_scores["ri_range_bayes"] = str(node_scores["bayes_factor_ranges"])
                ri_scores["ri_binary_bayes"] = str(node_scores["bayes_factor_binary_%s" %node_scores["threshold"]])
                meta_single["scores"] = [ri_scores]
            # Cognitive Atlas meta data
            meta_single["url"] = "http://www.cognitiveatlas.org/term/id/%s" %node
            meta_single["type"] = "concept"
            meta_single["thumbnail"] = "http://www.cognitiveatlas.org/images/logo-front.png"
            meta_single["concept"] = [relationship_table.name[relationship_table.id==node].tolist()[0]]
            meta_single["task"] = ""
            meta_single["contrast"] = []
            meta_single["download"] = "http://www.cognitiveatlas.org/rdf/id/%s" %node
            meta_single["category"] = categories[node]["category"]
            if concept[0]["definition_text"]:
                meta_single["description"] = concept[0]["definition_text"].encode("utf-8")
            else:
                meta_single["description"] = ""
    meta_data[node] = meta_single

## STEP 2: VISUALIZATION WITH PYBRAINCOMPARE
from pybraincompare.ontology.tree import named_ontology_tree_from_tsv, make_ontology_tree_d3

# First let's look at the tree structure
# output_json = "%s/task_contrast_tree.json" % outfolder
tree = named_ontology_tree_from_tsv(relationship_table,output_json=None,meta_data=meta_data)
html_snippet = make_ontology_tree_d3(tree)
web_folder = "%s/tree" %web
make_analysis_web_folder(html_snippet,web_folder)

# To get a dump of just the tree (for use in more advanced custom web interface)
filey = open('%s/tree/reverseinference2.json' %web,'wb')
filey.write(json.dumps(tree, sort_keys=True,indent=4, separators=(',', ': ')))
filey.close()


## STEP 3: Export individual scores

### Images
single_scores_folder = "%s/data/individual_scores" %base  # any kind of tsv/result file
single_scores = glob("%s/*.pkl" %single_scores_folder)
scores_export_folder = "%s/indscores" %web
if not os.path.exists(scores_export_folder):
    os.mkdir(scores_export_folder)

for s in range(0,len(single_scores)):
    print "Parsing data for images %s of %s" %(s,len(single_scores))
    single_score_pkl = single_scores[s]
    ss = pickle.load(open(single_score_pkl,"rb"))
    meta_single = {}
    meta_single["scores"] = ss["single_scores"].to_dict(orient="records")
    meta_single["image_id"] = ss["image_id"]
    node = ss["image_id"]
    image_id = node
    # Again include meta data
    relationship_table_row = relationship_table[relationship_table.id==node]
    # Reverse inference scores
    meta_single["category"] = ""
    meta_single["type"] = "nii"
    concepts = relationship_table.parent[relationship_table.name == image_id]
    concepts = [relationship_table.name[relationship_table.id==c].tolist()[0] for c in concepts]
    neurovault_row = images[images.image_id == int(image_id)]
    collection_row = collections[collections.collection_id == neurovault_row.collection.tolist()[0]]
    collection_meta = {"DOI":collection_row["DOI"].tolist()[0],
                           "authors":collection_row["authors"].tolist()[0],
                           "journal":collection_row["journal_name"].tolist()[0],
                           "url":collection_row["url"].tolist()[0],
                           "subjects":collection_row["number_of_subjects"].tolist()[0],
                           "smoothing_fwhm":str(collection_row["smoothing_fwhm"].tolist()[0]).encode("utf-8")}
    meta_single["url"] = neurovault_row["url"].tolist()[0]
    meta_single["thumbnail"] = neurovault_row["thumbnail"].tolist()[0]
    meta_single["images"] = neurovault_row["thumbnail"].tolist()
    meta_single["task"] = neurovault_row["cognitive_paradigm_cogatlas"].tolist()[0]
    meta_single["contrast"] = neurovault_row["cognitive_contrast_cogatlas"].tolist()[0]
    meta_single["download"] = neurovault_row["file"].tolist()[0]
    meta_single["concept"] = concepts
    if neurovault_row["description"].tolist()[0]:
        description = str(neurovault_row["description"].tolist()[0]).encode("utf-8");
        if description != "nan":
            meta_single["description"] =  description
        else:
            meta_single["description"] = ""
    else:
        meta_single["description"] = ""
    output_file = "%s/ri_%s.json" %(scores_export_folder,meta_single["image_id"])
    filey = open(output_file,'wb')
    filey.write(json.dumps(meta_single, sort_keys=True,indent=4, separators=(',', ': ')))
    filey.close()
    
### Concepts
for node in unique_nodes:
    # This is a concept node
    if not re.search("node_",node):
        if node != "1":
            relationship_table_row = relationship_table[relationship_table.id==node]
            contrast_name = relationship_table_row.name.tolist()[0]
            concept = get_concept(id=node).json
            # Reverse inference scores? Otherwise, we don't care
            if node in scores_df.index:
                meta_single = {}
                node_scores = scores_df.loc[node]
                ri_scores = {}
                ri_scores["count_in"] = node_scores["count_in"]
                ri_scores["count_out"] = node_scores["count_out"]
                ri_scores["ri_score_binary"] = node_scores["scores_binary"]
                ri_scores["ri_binary_threshold"] = node_scores["threshold"]
                ri_scores["ri_range_score_in"] = node_scores["ri_in_ranges"]
                ri_scores["ri_range_score_out"] = node_scores["ri_out_ranges"]
                ri_scores["ri_binary_score_in"] = node_scores["ri_in_binary"]
                ri_scores["ri_binary_score_out"] = node_scores["ri_out_binary"]
                ri_scores["ri_range_bayes"] = node_scores["bayes_factor_ranges"]
                ri_scores["ri_binary_bayes"] = node_scores["bayes_factor_binary_%s" %node_scores["threshold"]]
                range_index_names = [x for x in node_scores.index if re.search("[[]",x)]
                ri_scores["ri_ranges"] =  node_scores[range_index_names].to_dict()
                meta_single["scores"] = ri_scores
                # Reverse inference scores - all images
                meta_single["ri_range_scores"] = singles_ranges.loc[node,:].to_dict()
                meta_single["ri_binary_scores"] = singles_binary.loc[node,:].to_dict()
                image_ids =[int(x) for x in singles_binary.loc["trm_557b495cdde57",:].index.tolist()]
                meta_single["images"] = images["thumbnail"][images.image_id.isin(image_ids)].tolist()
                # Cognitive Atlas meta data
                meta_single["url"] = "http://www.cognitiveatlas.org/term/id/%s" %node
                meta_single["type"] = "concept"
                meta_single["thumbnail"] = "http://www.cognitiveatlas.org/images/logo-front.png"
                meta_single["concept"] = [relationship_table.name[relationship_table.id==node].tolist()[0]]
                meta_single["task"] = ""
                meta_single["contrast"] = []
                meta_single["download"] = "http://www.cognitiveatlas.org/rdf/id/%s" %node
                meta_single["category"] = categories[node]["category"]
                if concept[0]["definition_text"]:
                    meta_single["description"] = concept[0]["definition_text"].encode("utf-8")
                else:
                    meta_single["description"] = ""
                output_file = "%s/ri_%s.json" %(scores_export_folder,node)
                filey = open(output_file,'wb')
                filey.write(json.dumps(meta_single, sort_keys=True,indent=4, separators=(',', ': ')))
                filey.close()


## STEP 4: Export image similarity

# Image similarity?
similarity = pandas.read_csv("%s/contrast_defined_images_pearsonpd_similarity.tsv" %results,sep="\t")

# Done!
