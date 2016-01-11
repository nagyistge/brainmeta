from gensim import corpora, models, similarities
from wordfish.corpus import get_corpus, get_meta
from wordfish.nlp import sentence2words
from wordfish.utils import mkdir
import pickle
import gensim
import pandas
import json
import os
base_dir = os.environ["WORDFISH_HOME"]
dict_dir = mkdir("%s/dict" %(base_dir))
models_dir = mkdir("%s/models" %(base_dir))

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
lda.print_topics(10)

# Load CNP questions
cnp_questions = get_meta(base_dir)["cnp"]

# Here is how to transform new text to model
cnp_lda = dict()
for cnp_question in cnp_questions:
    print "Parsing %s" %cnp_question
    text = json.loads(open(cnp_question,"rb").read().strip("\n"))
    words = sentence2words(text["text"])
    doc_bow = dictionary.doc2bow(words)
    cnp_lda[text["labels"][0]] = lda[doc_bow]
    
pickle.dump(cnp_lda,open("%s/cnp_lda.pkl" %models_dir,"wb"))

# Let's make a data frame of 
scores = pandas.DataFrame(columns=range(100))
for cnp_question, topic_mapping in cnp_lda.iteritems():
    for single_mapping in topic_mapping:
        scores.loc[cnp_question,single_mapping[0]] = single_mapping[1]
scores.to_csv("%s/cnp_lda_scores.csv" %models_dir)
