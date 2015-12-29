
from glob import glob
import pandas
relations_files = glob("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/wordfish/paper/models/*.tsv")
from scipy import stats
import numpy
import os

# First let's make a list of terms defined for each matrix
terms = dict()
allterms = []
for relation_file in relations_files:
    tmp = pandas.read_csv(relation_file,sep="\t",index_col=0)
    tmp.index = tmp.columns
    terms[relation_file] = tmp
    allterms = numpy.unique(allterms + tmp.columns.tolist()).tolist()

# For each term, find contexts for which it is defined
for t in range(len(allterms)):
    print "Parsing %s of %s" %(t,len(allterms))
    term = allterms[t] 
    contexts = [x for x in terms.keys() if term in terms[x].columns]

    # Get terms for those contexts
    shared_terms = []
    for context in contexts:
        shared_terms = shared_terms + terms[context].columns.tolist()
    shared_terms = numpy.unique(shared_terms).tolist()
    df = pandas.DataFrame(columns=shared_terms)
    for context in contexts:
        row = terms[context].loc[term]
        df.loc[context,row.index] = row.tolist()
    df = df.fillna(0)
    df.to_csv("%s_contexts.tsv" %(term),sep="\t")

labels =  [os.path.basename(x).replace(".tsv","") for x in relations_files]
enricheddf = pandas.DataFrame(index=labels,columns=allterms)

# For each term, calculate enriched contexts, and associated words
for t in range(len(allterms)):
    term = allterms[t] 
    contexts = [x for x in terms.keys() if term in terms[x].columns]

    # Get terms for those contexts
    shared_terms = []
    for context in contexts:
        shared_terms = shared_terms + terms[context].columns.tolist()
    shared_terms = numpy.unique(shared_terms).tolist()
    df = pandas.DataFrame(columns=shared_terms)
    for context in contexts:
        row = terms[context].loc[term]
        df.loc[context,row.index] = row.tolist()
    df = df.fillna(0)
  
    # TASK 1: OVERALL ENRICHMENT OF A WORD TO A CONTEXT
    # Overall enrichment for the term against a larger corpus 
    # is the normalized sum of similarities across all other terms,
    # A higher enrichment score means the term has a similar wordvector
    # to some subset of terms, which (since we are using word2vec)
    # means that it is commonly found / associated / near to the term
    enrichment = df.sum(axis=1) / df.sum().sum()
    enrichment_labels = [os.path.basename(x).replace(".tsv","") for x in df.index]
    enricheddf.loc[enrichment_labels,term] = enrichment.tolist()
  
    # For what words?
    # Difference, how about 3SD
    sdy = sd(as.matrix(tmp))
    meany = df.mean(axis=1)
    threesdup = mean+(3*df.std(axis=1))
    context_up = which(as.matrix(tmp) >= threesdup, arr.ind = TRUE)
  
  # We can then further look at the highest scores for each to
  # determine what those other enriched terms are...
  # These are terms and context that are most similar to query word
  enriched_terms = colnames(tmp)[context_up[,2]]
  enriched_context = rownames(tmp)[context_up[,1]]  

    # TASK 2: ENRICHMENT OF A TERM/TERM RELATIONSHIP FOR SOME CONTEXT
    word_means = df.sum(axis=0)  # mean similarity of all words to the given word
    word_stds = word_means.std()    


# Here we have different Neurosynth terms across contexts, and we would
# want to know which term/context pairs are significantly different from the rest.

# take differences between group matrices and look for significant differences - the differences are the features of importance to look for given some group (context)

Then we have to take some data from people with disorders, and figure out if it is enriched for those terms. It could also be the case that disorders overlap (and so we don't necessarily want the most significantly different relationships) but rather the ones that are more prevalent? something?

Each term 1 / term 2 is a distribution of values, and we can do two things. 1) determine which contexts have significantly same and differenat relatoinsips for the term. This could be a way to identify terms that are "definitive" of the context (differences) and contexts that are similar (same). We then would want to see the same similarities reflected in those terms in some real data... 2) "classify" someone into this space by projecting their data (turning it into a vector) based on each set of relationships, 
# Run classification across many
# models load_models(base)

When people describe themselves, do the terms map more onto one context than another? Aka, if someone gives me a list of terms, can I predict they have a disorder? Likelihood?

# with little text
1) For each term:
   Subset context matrix to other terms in description
   Rank based on similarity
   Highest similar terms 
1) Select terms
2) Select context matrices
3) 

# REQUIRES LOTS OF TEXT
# Take person entire blog, build word2vec
# For terms that are defined in their model, assess similarity to ecah context
# rank based on similarity


    


