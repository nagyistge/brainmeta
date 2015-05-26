#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Experiment 3 How to handle missing data - generate group maps on TACC

import sys
from clusterhcp.stats import run_randomise_local

groupA = sys.argv[1].split(",")
letter = sys.argv[2]
output_directory = sys.argv[3]
map_id = sys.argv[4
]

# Generate group maps with randomise, to output folder
run_randomise_local(groupA,output_directory,output_prefix="%s_group%s" %(map_id,letter))

