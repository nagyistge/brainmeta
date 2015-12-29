'''
1.build_model.py

USAGE:

This script will use corpus / terms defined from wordfish plugins, and build 
a dataframe (vectors) for each entity in the reddit corpus (rows) associated
with a particular label (reddit board). This data frame is then used as input
to build pairwise models to predict disorder labels from the wordfish vectors.

Your wordfish project home directory should be defined as WORDFISH_HOME

    - Merge terms and relationships into a common corpus
    - For all text extract features with deep learning (word2vec)
    - For each set of terms, parse over text and find occurrences

'''

# First train simple word2vec model with different corpus
from wordfish.analysis import train_word2vec_model, save_models, export_models_tsv, load_models, vocab_term_intersect, extract_similarity_matrix, extract_vectors, DeepTextAnalyzer
from wordfish.corpus import get_corpus, get_meta
from wordfish.terms import get_terms
from wordfish.utils import mkdir
import pandas
import pickle
import sys
import os

base_dir = os.environ["WORDFISH_HOME"]

# Setup analysis output directory
analysis_dir = mkdir("%s/analysis" %(base_dir))
model_dir = mkdir("%s/models" %(analysis_dir))
vector_dir = mkdir("%s/vectors" %(analysis_dir))

corpus = get_corpus(base_dir)

disorders = dict()
reddit = corpus["reddit"]
for red in reddit:
    topic = os.path.basename(red).split("_")[0]
    if topic in disorders:
        disorders[topic].append(red)
    else:
        disorders[topic] = [red]

corpus.update(disorders)


# Train corpus specific models
models = dict()
print "Training models..."
for corpus_id,sentences in corpus.iteritems():
    try:
        models[topic] = train_word2vec_model(sentences)
    except:
        print "Error building model for %s" %(topic)
        pass

# Export models to tsv, and save - will use later
save_models(models,base_dir)
export_models_tsv(models,base_dir)

# Here are all terminologies
terms = get_terms(base_dir,subset=True)
#terms.keys()
#['cattell', 'neurosynth', 'fsl', 'fma_nif', 'cognitiveatlas']

# We want to find overlapping terms between each terminology and word2vec model
for model_name,model in models.iteritems(): # word2vec
    intersects = vocab_term_intersect(terms,model) # terms are terminology
    #>>> intersects.keys() intersect of model with each terminology below
    #['cattell', 'neurosynth', 'fsl', 'fma_nif', 'cognitiveatlas']
    for tag,ints in intersects.iteritems():
        vs = numpy.unique([x[3] for x in ints]).tolist()
        # (755, 131, u'conversation', 'conversation')
        # And export the model using the vocabulary (vs) as a subset
        export_models_tsv({"%s_%s" %(tag,model_name):model},base_dir,vocabs=[vs])

# Export vectors - will use to compare contexts in later analysis
for model_name,model in models.iteritems():
    print "Processing %s" %(model_name)
    vecs = extract_vectors(model)
    vecs.to_csv("%s/%s.tsv" %(vector_dir,model_name),sep="\t")
