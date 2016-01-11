from gensim import corpora, models, similarities
from wordfish.corpus import get_corpus
from wordfish.nlp import sentence2words
from wordfish.utils import mkdir
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
    for corpus_file in corpus_files:
        text = open(corpus_file,"rb").read().strip("\n")
        text = " ".join(sentence2words(text))
        corpus.append(dictionary.doc2bow(text))
    return corpus

corpus = build_corpus(corpus_files,dictionary)
corpora.MmCorpus.serialize('%s/reddit.mm' %models_dir, corpus)

# Latent semantic analysis
lsi = gensim.models.lsimodel.LsiModel(corpus=corpus, id2word=dictionary, num_topics=400)
