'''
0build_models.py

Your wordfish project home directory should be defined as WORDFISH_HOME

This script will build models based on reddit and neurosynth corpus, to
prepare for extraction of feature vectors for CNP questions based
on each corpus. General workflow is as follows:

     Generate word embeddings for 120K reddit corpus
     Map CNP questions into this space
     Do representational similarity analysis
     
Building of models could be run in parallel, if you get impatient.

'''

from wordfish.analysis import build_models, load_models, extract_similarity_matrix, export_vectors, featurize_to_corpus
from wordfish.models import build_svm
from wordfish.corpus import get_corpus, get_meta, subset_corpus
from wordfish.terms import get_terms
from wordfish.utils import mkdir
import pandas
import pickle
import os

base_dir = os.environ["WORDFISH_HOME"]

# Get analysis output directories
analysis_dir = mkdir("%s/analysis" %(base_dir))
model_dir = mkdir("%s/models" %(analysis_dir))
vector_dir = mkdir("%s/vectors" %(analysis_dir))

# Generate more specific corpus by way of file naming scheme
corpus = get_corpus(base_dir)

# Subset corpus to different boards
boards = subset_corpus(corpus["reddit"])

# Build missing models
tobuild = [b for b in boards.keys() if b not in models.keys()]
boards = dict((k, boards[k]) for k in tobuild)
newmodels = build_models(boards)

# Save all models
save_models(newmodels,base_dir)

# The next step is to do RSA analysis, and we will submit scripts to run in parallel.
# see rsa_prep.py, run_rsa_prep.py, rsa.py and run_rsa.py
