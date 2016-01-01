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

from wordfish.analysis import train_word2vec_model, save_models, export_models_tsv, load_models, vocab_term_intersect, extract_similarity_matrix, extract_vectors, DeepTextAnalyzer
from wordfish.models import build_svm
from wordfish.corpus import get_corpus, get_meta, subset_corpus
from wordfish.terms import get_terms
from wordfish.utils import mkdir
import os
import numpy
import pandas
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

classifiers = build_svm(vectors=vectors,labels=labels,kernel="linear")


# Now we have extracted new reddit posts, let's see if we can predict the post
meta = get_meta(base_dir)["reddit"]
from sklearn.neighbors import NearestNeighbors

# Here are the reddit vectors (using different reddit data from months ago, same boards)
reddit_vectors = pickle.load(open("%s/analysis/models/classifier_reddit.pkl" %(base_dir),"rb"))
nbrs = NearestNeighbors(n_neighbors=100, algorithm='ball_tree').fit(reddit_vectors)
model = load_models(base_dir,"reddit")["reddit"]
test_labels = [l.split("/")[-1].replace("_meta.txt","") for l in meta]
train_labels = [l.split("_")[0] for l in reddit_vectors.index]

# We will save a data frame of predictions - counts for each test
predictions = pandas.DataFrame(columns=numpy.unique(train_labels).tolist())

analyzer = DeepTextAnalyzer(model)
for m in range(len(meta)):
    print "Processing %s of %s" %(m,len(meta))
    post = meta[m]
    post = json.load(open(post,"r"))
    board = post["labels"][0]
    label = meta[m].split("/")[-1]
    if label not in predictions.index:
        vector = analyzer.text2mean_vector(post["text"],read_file=False)
        if vector != None:
            distances, indices = nbrs.kneighbors(vector)
            rx = [train_labels[x] for x in indices[0]]
            # Diagnosis is the top label
            predict = pandas.DataFrame(rx)[0].value_counts()
            predictions.loc[label,predict.index] = predict

# How many did we get right?
correct_labels = [x.split("_")[0] for x in predictions.index]
predicted_labels = predictions.idxmax(axis=1).tolist()
accuracy = 0
for x in range(len(correct_labels)):
    if correct_labels[x]==predicted_labels[x]:
        accuracy +=1
accuracy/float(len(correct_labels))

# Let's now do KNN with a new vector. First map to the space
vectors = pickle.load(open("%s/analysis/models/classifier_reddit.pkl" %(base_dir),"rb"))
vectors[vectors.isnull()] = 0

model = load_models(base_dir,"reddit")["reddit"]
text = "I get super anxious about doing things at night. I tend to get increasingly tired after about 8, and I don't eat after 6. Generally, I try to be in the house by 7 p.m. But tonight was my friend's birthday and it was at a restaurant not far from me at 7 p.m. I actually walked there - also a success because I have an auto-immune disease that makes it painful to walk. And we ate dinner together and I didn't have to flee and I didn't get sick! I sound pretty dumb, eh. These are things that normal people do. But whatever. Suck it! I haven't been out to dinner past 7 p.m. in two years. I'm now exhausted. edit: too many words"
analyzer = DeepTextAnalyzer(model)
analyzer.text2mean_vector(text)
vector = analyzer.text2mean_vector(text,read_file=False)
# Do KNN

from sklearn.neighbors import NearestNeighbors
nbrs = NearestNeighbors(n_neighbors=15, algorithm='ball_tree').fit(vectors)
distances, indices = nbrs.kneighbors(vector)
rx = [labels[x] for x in indices[0]]


# Try mapping neurosynth abstracts to reddit, how many terms we get in top
model = load_models(base_dir,"neurosynth")["neurosynth"]
analyzer = DeepTextAnalyzer(model)
    
vectors = pandas.DataFrame(columns=range(300))
meta = get_meta(base_dir)["neurosynth"]
pmids = [x.split("/")[-1].replace("_meta.txt","") for x in meta]

for r in range(len(meta)):
    article = json.load(open(meta[r],"rb"))
    pmid = pmids[r]
    # Build a model for everyone else
    if pmid not in vectors.index:
        try:
            print "Processing %s of %s" %(r,len(meta))
            vectors.loc[pmid] = analyzer.text2mean_vector(article["text"],read_file=False)
        except:
            pass

# Make a labels lookup based on pmid
lookup = {}
allterms = []
for r in range(len(meta)):
    article = json.load(open(meta[r],"rb"))
    pmid = pmids[r]
    lookup[pmid] = article["labels"]    
    allterms = allterms + article["labels"]

# Get actual counts of all terms
term_counts = pandas.DataFrame(allterms)[0].value_counts()

# Save pickle of df foruse later, this is a matrix of pmids
from numpy.random import choice
vectors[vectors.isnull()==True] = 0
clf = {"vectors":vectors,"labels":lookup,"all_terms":allterms,"term_counts":term_counts}
pickle.dump(clf,open("%s/analysis/models/classifier_neurosynth.pkl" %(base_dir),"wb"))

# Calculate percentages
from scipy.stats import ttest_1samp

# Now use the vectors to find the most similar other vectors
for r in range(len(meta)):
    article = json.load(open(meta[r],"rb"))
    pmid = pmids[r]
    vector = vectors.loc[pmid]
    # Do KNN to find N most similar
    N = len(article["labels"])
    nbrs = NearestNeighbors(n_neighbors=N, algorithm='ball_tree').fit(vectors.loc[vectors.index!=pmid])
    distances, indices = nbrs.kneighbors(vector)
    pmid_matches = [pmids[x] for x in indices[0]]
    label_matches = []
    for p in pmid_matches:
        label_matches = label_matches + lookup[p]
    result_counts = pandas.DataFrame(label_matches)[0].value_counts()
    # Now we need to determine how the count compares to by chance?
    count_iters = pandas.DataFrame(columns=term_counts.index)
    pvalues = pandas.DataFrame(columns=term_counts.index)
    for iter in range(1000):
        print "Running permutation %s of 1000" %iter
        random_selection = choice([p for p in pmids if p != pmid],N)
        label_random_matches = [lookup[p] for p in random_selection]
        random_counts = pandas.DataFrame(label_random_matches)[0].value_counts()
        count_iters.loc[iter] = random_counts
    count_iters[count_iters.isnull()==True] = 0
    # Now for each of our actual terms, do a 1 sample T test to determine if count is different from mean
    for result in result_counts.index:
        actual_count = result_counts.loc[result]
        counts_permuted = count_iters[result]
        pval = ttest_1samp(counts_permuted,actual_count)[1]
        pvalues.loc[0,result] = pval    





