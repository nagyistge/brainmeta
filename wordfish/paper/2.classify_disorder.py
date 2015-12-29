'''
2.classify_disorder.py

USAGE:

This script will use a specific word2vec model, corpus, and terminology to
build pairwise SVM classifiers to predict a disorder label (a reddit board)
based on the text of the board. Note that the same functionality has been
compressed into wordfish functions:

model = load_models(base_dir,"neurosynth")["neurosynth"]
meta = get_meta(base_dir)["neurosynth"]
vectors,labels = featurize_to_corpus(model,meta)
classifiers = build_svm(vectors=vectors,labels=labels,kernel="linear")

'''

from wordfish.analysis import build_models, save_models, export_models_tsv, load_models, extract_similarity_matrix, export_vectors, featurize_to_corpus
from wordfish.models import build_svm
from wordfish.corpus import get_corpus, get_meta, subset_corpus
from wordfish.terms import get_terms
from wordfish.utils import mkdir
import os
import pickle

base_dir = os.environ["WORDFISH_HOME"]

# CLASSIFICATION OF DISORDER with reddit text ############################
# Can we train a model to predict disorder based on text from reddit?
# Load meta data associated with corpus, this is where we have labels

meta = get_meta(base_dir)

# First we will generate a vector to describe each reddit post as
# an average of the words in it from our model. This should be ok
# to do as the word2vec model does not know anything about 
# the groupings.

reddit_corpus = get_corpus(base_dir)["reddit"]
#len(reddit_corpus)
#121862

# Let's generate a vector of labels
labels = [os.path.basename(x).split("_")[0] for x in reddit_corpus]
numbers = [os.path.basename(x).split("_")[1] for x in reddit_corpus]

# Use reddit to build a model
model = load_models(base_dir,"reddit")
analyzer = DeepTextAnalyzer(model)
    
vectors = pandas.DataFrame(columns=range(300))

for r in range(len(reddit_corpus)):
    post = reddit_corpus[r]
    label = "%s_%s" %(labels[r],numbers[r])
    # Build a model for everyone else
    if label not in vectors.index:
        try:
            print "Processing %s of %s" %(r,len(reddit_corpus))
            vectors.loc[label] = analyzer.text2mean_vector(post)
        except:
            pass
    if r%10000.0==0:
        # Save pickle of df foruse later
        pickle.dump(vectors,open("%s/analysis/models/classifier_reddit.pkl" %(base_dir),"wb"))

# Save pickle of df foruse later, this is a matrix of reddit posts with labels we can use to train model (classify_blog.py)
pickle.dump(vectors,open("%s/analysis/models/classifier_reddit.pkl" %(base_dir),"wb"))
#vectors = pickle.load(open("%s/analysis/models/classifier_reddit.pkl" %(base_dir),"rb"))

classifiers = build_svm(vectors=vectors,labels=labels,kernel="linear")
