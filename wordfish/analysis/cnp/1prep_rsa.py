# REPRESENTATIONAL SIMILARITY ANALYSIS ####################################################

# The vectors that we just created describe the CNP questions projected onto a board, 
# meaning we have a vector for each question based on taking an average of words defined 
# by the model. We now can generate similarity matrices for each set of vectors (comparing 
# similarity of questionA to questionB in CNP) and then compare those matrices with
# representational similarity analysis. 

from wordfish.analysis import build_models, load_models, extract_similarity_matrix, export_vectors, featurize_to_corpus
from wordfish.models import build_svm
from wordfish.corpus import get_corpus, get_meta, subset_corpus
from wordfish.terms import get_terms
from wordfish.utils import mkdir
import pandas
import numpy
import pickle
import os
import re
import sys

model_id = sys.argv[1]
base_dir = sys.argv[2]

analysis_dir = mkdir("%s/analysis" %(base_dir))
vector_dir = mkdir("%s/vectors" %(analysis_dir))

# Load CNP questions
cnp = get_meta(base_dir)["cnp"]
models = load_models(base_dir)
model = models[model_id]

# Map CNP questions onto corpus vectors
vectors,labels = featurize_to_corpus(model,cnp)

# Put correct label back on question
label_list=[]
for label in labels.iterrows():
    idx=label[0]
    label_list.append(label[1][label[1]==1].index[0])
vectors.index = label_list

# One question is duplicated
duplicates = vectors.loc["EYSENCK2"].copy()
duplicates.index = ["EYSENCK2","DROP"]
vectors = vectors.drop("EYSENCK2")
vectors.loc["EYSENCK2"] = duplicates.loc["EYSENCK2"].copy()

# Get rid of the test questions
exp = re.compile("TEST_")
test_names = [x for x in vectors.index if exp.match(x)]
vectors = vectors.drop(test_names)

# Some questions have all zeros, not mappable to the space
vectors = vectors.loc[~(vectors==0).all(axis=1)]

# Save to file
vectors.to_csv("%s/cnp_vectors_wrt_%s.tsv" %(vector_dir,model_id),sep="\t")
