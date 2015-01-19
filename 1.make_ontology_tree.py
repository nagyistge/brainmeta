#!/usr/bin/env python2

from annotate import ontology
from template.visual import view

# This defines a "mock" ontology for the openfmri proof of concept set
relationship_table = "/home/vanessa/Documents/Work/BRAINMETA/mr/cogatlas_annotate_triples.tsv"

# Create a data structure for d3 visualization and data analysis (no output json, specify with output_json=)
data_structure = ontology.named_ontology_tree_from_tsv(relationship_table)

html_snippet = ontology.make_ontology_tree_d3(data_structure)

# View in browser!
view(html_snippet)


