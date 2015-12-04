'''
4wordfish_sims.py

Your wordfish project home directory should be defined as WORDFISH_HOME

We have already generated vectors that project CNP questions onto different corpus (reddit,
neurosynth), and then we parsed CNP data and made a matrix of diagnoses. We also made
similarity matrices for the actual CNP question data (based on different disorder groups).
Now we want to generate similarity matrices for the CNP questions mapped to the wordfish 
vectors. THEN we can do RSA! :)

These analyses have not been integrated into the wordfish pipeline, but will be if
they are useful.

'''


# REPRESENTATIONAL SIMILARITY ANALYSIS ####################################################

from wordfish.utils import mkdir
from glob import glob
import pandas
import numpy
import pickle
import os
import re

base_dir = os.environ["WORDFISH_HOME"]

# Get analysis output directories
analysis_dir = mkdir("%s/analysis" %(base_dir))
scripts_dir = mkdir("%s/scripts" %(base_dir))
vectors_dir = mkdir("%s/vectors" %(analysis_dir))
vectors_files = glob("%s/cnp_vectors_wrt_*" %vectors_dir)

# We don't need paralell clusterizing for this!
for vectors_file in vectors_files:
    label = os.path.basename(vectors_file).split("_")[-1].replace(".tsv","")
    print "Parsing %s" %(label)
    data = pandas.read_csv(vectors_file,sep="\t",index_col=0)
    corr = data.T.corr()
    corr.to_csv("%s/cnp_vectors_corr_wrt_%s.tsv" %(vector_dir,label),sep="\t")
# boumboumboum!
