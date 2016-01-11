from gensim import corpora, models, similarities
from wordfish.corpus import get_corpus, get_meta
from wordfish.nlp import sentence2words
from wordfish.utils import mkdir
import pickle
import gensim
import pandas
import numpy
import json
import re
import os
base_dir = os.environ["WORDFISH_HOME"]
dict_dir = mkdir("%s/dict" %(base_dir))
models_dir = mkdir("%s/models" %(base_dir))
vectors_dir = mkdir("%s/analysis/vectors" %(base_dir))
rsa_dir = mkdir("%s/analysis/rsa" %(base_dir))

# First we need to build a dictionary
def build_dictionary(corpus_files):
    word_lists = []
    for c in range(len(corpus_files)):
        print "Parsing %s of %s" %(c,len(corpus_files))
        corpus_file = corpus_files[c]
        text = open(corpus_file,"rb").read().strip("\n")
        words = sentence2words(text)
        if len(words) != 0:
            word_lists.append(words)    
    dictionary = corpora.Dictionary(word_lists)     
    return dictionary

corpus_files = get_corpus(base_dir)["reddit"]
dictionary = build_dictionary(corpus_files) # This takes a while
dictionary.save('%s/reddit.dict' %dict_dir)
# I think this will be how to load it later
#dictionary = gensim.corpora.Dictionary.load_from_text('reddit.dict')
print(dictionary.token2id)

# Now we need to generate a bow vector for each document
def build_corpus(corpus_files,dictionary):
    corpus = []
    for c in range(len(corpus_files)):
        print "Parsing %s of %s" %(c,len(corpus_files))       
        corpus_file = corpus_files[c]
        text = open(corpus_file,"rb").read().strip("\n")
        words = sentence2words(text)
        corpus.append(dictionary.doc2bow(words))
    return corpus

corpus = build_corpus(corpus_files,dictionary)
corpora.MmCorpus.serialize('%s/reddit.mm' %models_dir, corpus)

# Latent semantic analysis
#lsi = gensim.models.lsimodel.LsiModel(corpus=corpus, id2word=dictionary, num_topics=400)
#lsi.print_topics(10)

# Latent Dirichlet Allocation
lda = gensim.models.ldamodel.LdaModel(corpus=corpus, id2word=dictionary, num_topics=100, update_every=1, chunksize=10000, passes=1)
lda10 = gensim.models.ldamodel.LdaModel(corpus=corpus, id2word=dictionary, num_topics=10, update_every=1, chunksize=10000, passes=1)

lda.print_topics(10)
pickle.dump(lda,open("%s/reddit_lda_model.pkl" %models_dir,"wb"))
pickle.dump(lda10,open("%s/reddit_lda_model_10.pkl" %models_dir,"wb"))

# Load CNP questions
cnp_questions = get_meta(base_dir)["cnp"]

# Here is how to transform new text to model
cnp_lda = dict()
cnp_lda_10 = dict()
for cnp_question in cnp_questions:
    print "Parsing %s" %cnp_question
    text = json.loads(open(cnp_question,"rb").read().strip("\n"))
    words = sentence2words(text["text"])
    doc_bow = dictionary.doc2bow(words)
    cnp_lda[text["labels"][0]] = lda[doc_bow]
    cnp_lda_10[text["labels"][0]] = lda10[doc_bow]
    
pickle.dump(cnp_lda,open("%s/cnp_lda.pkl" %models_dir,"wb"))
pickle.dump(cnp_lda_10,open("%s/cnp_lda_10.pkl" %models_dir,"wb"))

# Let's make a data frame of 
scores = pandas.DataFrame(columns=range(100))
scores10 = pandas.DataFrame(columns=range(10))
for cnp_question, topic_mapping in cnp_lda.iteritems():
    for single_mapping in topic_mapping:
        scores.loc[cnp_question,single_mapping[0]] = single_mapping[1]
scores.to_csv("%s/cnp_lda_scores.csv" %models_dir)

for cnp_question, topic_mapping in cnp_lda_10.iteritems():
    for single_mapping in topic_mapping:
        scores10.loc[cnp_question,single_mapping[0]] = single_mapping[1]
scores10.to_csv("%s/cnp_lda_scores_10.csv" %models_dir)

allscores=[]
for cnp_question, topic_mapping in cnp_lda_10.iteritems():
    for single_mapping in topic_mapping:
        allscores.append(single_mapping[1])

# Calculate similarity matrix of questions based on LDA
from scipy.stats import pearsonr
sims_lda = pandas.DataFrame(columns=cnp_lda_10.keys(),index=cnp_lda_10.keys())
for q1 in sims_lda.index.tolist():
    print "Parsing %s" %(q1)
    for q2 in sims_lda.index.tolist():
        if q1!=q2:
            vector1 = scores10.loc[q1]
            vector2 = scores10.loc[q2]
            score = pearsonr(vector1,vector2)[1]
            sims_lda.loc[q1,q2] = score
            sims_lda.loc[q2,q1] = score
        else:
            sims_lda.loc[q1,q2] = 1.0

sims_lda.to_csv("%s/cnp_sims_df_10.csv" %models_dir)

# Local machine work
scores = pandas.read_csv("lda/cnp_lda_scores_10.csv")
cnp_lda = pickle.load(open("lda/cnp_lda_10.pkl","rb"))
lda = pickle.load(open("lda/reddit_lda_model_10.pkl","rb"))
sims_lda = pandas.read_csv("lda/cnp_sims_df_10.csv",index_col=0)
cnp_sims = pickle.load(open("%s/rsa_matrix_lookup.pkl" %rsa_dir,"rb"))
sima_lda=sims_lda.fillna(0)

# Function for rsa
def rsa(mfile,sims_lda,question_filter=None):
    m1 = pandas.read_csv(mfile,sep="\t",index_col=0)
    m2 = sims_lda
    overlap = m1.index[m1.index.isin(m2.index)]
    m1 = m1.loc[overlap,overlap]
    m2 = m2.loc[overlap,overlap]
    if question_filter != None:
       ls = [x for x in m1.index.tolist() if re.search(question_filter,x)]
       m1 = m1.loc[ls,ls]
       m2 = m2.loc[ls,ls]
    vector1 = m1.mask(numpy.triu(numpy.ones(m1.shape)).astype(numpy.bool)).values.flatten()
    vector2 = m2.mask(numpy.triu(numpy.ones(m2.shape)).astype(numpy.bool)).values.flatten()
    # Find overlapping (non-nans)
    vector1_defined = numpy.where(~numpy.isnan(vector1))[0]
    vector2_defined = numpy.where(~numpy.isnan(vector2))[0]
    idx = numpy.intersect1d(vector1_defined,vector2_defined)
    return pearsonr(vector1[idx],vector2[idx])[0]

cnp_labels = [x for x in cnp_sims.keys() if not re.match("wordfish",x) and x not in ["question_scale","assessment"]]

results = pandas.DataFrame(columns=["RSA"])
for label in cnp_labels:
    print "Processing %s" %(label)
    mfile = cnp_sims[label] 
    res = rsa(mfile,sims_lda)
    results.loc[label,"RSA"] = res

# Now run for subsets of questionnaires
questions = numpy.unique([''.join([i for i in x if not i.isdigit()]) for x in scores.index.tolist()]).tolist()
questions = [x for x in questions if not re.search("ASRS",x) and not re.search("ASSESS",x)]
questions = questions + ["ASRS","ASSESS"]
questions.pop(questions.index("TEST_"))
for label in cnp_labels:
    for question in questions:
        print "Processing %s and %s" %(label,question)
        mfile = cnp_sims[label] 
        res = rsa(mfile,sims_lda,question)
        results.loc["%s_%s" %(label,question),"RSA"] = res
        
results.to_csv("%s/rsa_cnp_lda.tsv" %rsa_dir,sep="\t")
