'''
run_analysis.py

USAGE:

python run_analysis.py

Your wordfish project home directory should be defined as WORDFISH_HOME

This script will extract feature vectors for CNP questions based
on a reddit corpus. General workflow is as follows:

     Generate word embeddings for 120K reddit corpus
     Map CNP questions into this space
     
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
reddit = corpus["reddit"]
cnp = get_meta(base_dir)["cnp"]

models = load_models(base_dir)
model = models["reddit"]

# Prepare a meta object for each sentences.txt

# Fit CNP questions to reddit model
vectors,labels = featurize_to_corpus(model,cnp)
# Put correct label back on question
label_list=[]
for label in labels.iterrows():
   idx=label[0]
   label_list.append(label[1][label[1]==1].index[0])
   
vectors.index = label_list
vectors.to_csv("%s/cnp_vectors.tsv" %(vector_dir),sep="\t")

# SUBSETS OF QUESTIONS ###########################################################
# We will want to make our own groupings of questions based on some context.
# Let's choose random terms and come up with sets!
random_words = ["anxiety","depression","bipolar","schizophrenia",
                "happy","impulsive","sad","decision","angry","woman",
                "man","young","old","helpless","dependent","medication",
                "troubled","eager","shy","quiet","lonely","choice","manic",
                "adhd","impaired","chronic","symptoms","anorexia","abuse"]

def make_match_df(model,vectors,words):
    # Let's try saving absolute value (meaning we want extremes) and not (only positive)
    matchdf = pandas.DataFrame(columns=vectors.index)
    for random_word in words:
        try:
            vectors_copy = vectors.copy()
            # Find most similar cnp vector
            vectors_copy.loc["QUERY"] = model[random_word]
            vectors_similar = vectors_copy.T.corr()
            vectors_similar.index = vectors_copy.index
            matchdf.loc[random_word,vectors.index] = vectors_similar.loc["QUERY",vectors.index]
        except:
            print "Skipping %s" %(random_word)
    return matchdf

matchdf = make_match_df(model,vectors,random_words)
matchdf.to_csv("%s/cnp_word_matches.tsv" %(vector_dir),sep="\t")
