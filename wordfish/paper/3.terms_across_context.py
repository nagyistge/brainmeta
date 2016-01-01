from scipy.stats import pearsonr
from glob import glob
import pandas
relations_files = glob("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/wordfish/paper/models/allwords/*.tsv")
from scipy import stats
import numpy
import os

# RSA ANALYSIS ################################################################################


# First let's make a list of terms defined for each matrix
terms = dict()
allterms = []
for relation_file in relations_files:
    tmp = pandas.read_csv(relation_file,sep="\t",index_col=0)
    tmp.index = tmp.columns
    terms[relation_file] = tmp
    allterms = numpy.unique(allterms + tmp.columns.tolist()).tolist()

# Set working directory to where we will save contexts
os.chdir("/home/vanessa/Documents/Work/WORDFISH/analysis/models/allwords")


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

# Make a lookup dataframe with stemmed versions
#lookup = dict()
#stemmed = []
#for relation_file in relations_files:
#    tmp = pandas.read_csv(relation_file,sep="\t",index_col=0)
#    tmp.index = do_stem(tmp.columns)
#    tmp.columns = tmp.index
#    lookup[relation_file] = tmp
#    allterms = numpy.unique(allterms + tmp.columns.tolist()).tolist()


from wordfish.nlp import do_stem, processText

# Try classification of new paragraph using context
input_vector = "I get super anxious about doing things at night. I tend to get increasingly tired after about 8, and I don't eat after 6. Generally, I try to be in the house by 7 p.m. But tonight was my friend's birthday and it was at a restaurant not far from me at 7 p.m. I actually walked there - also a success because I have an auto-immune disease that makes it painful to walk. And we ate dinner together and I didn't have to flee and I didn't get sick! I sound pretty dumb, eh. These are things that normal people do. But whatever. Suck it! I haven't been out to dinner past 7 p.m. in two years. I'm now exhausted. edit: too many words"

stemmed = processText(" ".join(allterms))
vector = processText(input_vector)
overlap_vector = [v for v in vector if v in stemmed]

# We want to know which context are these terms most enriched? Enriched means that the words are likely to be found together, so they have higher scores to the other terms.
enricheddf = pandas.DataFrame(index=relations_files,columns=overlap_vector) # preserve duplicated terms
for relation_file in relations_files:
    tmp = pandas.read_csv(relation_file,sep="\t",index_col=0)
    stemmed_columns = processText(" ".join(tmp.columns))
    tmp.index = stemmed_columns
    tmp.columns = tmp.index
    # Find the words that we have in the vocabulary
    overlap = [v for v in vector if v in stemmed_columns]
    unique_overlap = numpy.unique(overlap).tolist()
    rxdf = tmp.loc[[x for x in stemmed_columns if x in unique_overlap],[x for x in stemmed_columns if x in unique_overlap]]
    enrichment = rxdf.mean() # This is a mean pearson score
    enricheddf.loc[relation_file,unique_overlap] = enrichment.loc[unique_overlap]
    # Not being in the vocabulary means the word wasn't even included in the model ==> enrichment is 0
    enricheddf[enricheddf.isnull()==True] = 0
    # The total enrichment in a corpus is a mean score across all the words
    enricheddf.mean(axis=1)


# RSA ANALYSIS ################################################################################
# We can assess differences between the different contexts (eg disorders) by doing RSA analysis
# with the matrices

# Write a function to perform RSA
def rsa(mfile1,mfile2):
    m1 = pandas.read_csv(mfile1,sep="\t",index_col=0)
    m2 = pandas.read_csv(mfile2,sep="\t",index_col=0)
    m1.index = m1.columns
    m2.index = m2.columns
    # Reduce to common questions (only an issue if neurosynth), diff is one question
    overlap = m1.index[m1.index.isin(m2.index)]
    m1 = m1.loc[overlap,overlap]
    m2 = m2.loc[overlap,overlap]
    vector1 = m1.mask(numpy.triu(numpy.ones(m1.shape)).astype(numpy.bool)).values.flatten()
    vector2 = m2.mask(numpy.triu(numpy.ones(m2.shape)).astype(numpy.bool)).values.flatten()
    # Find overlapping (non-nans)
    vector1_defined = numpy.where(~numpy.isnan(vector1))[0]
    vector2_defined = numpy.where(~numpy.isnan(vector2))[0]
    idx = numpy.intersect1d(vector1_defined,vector2_defined)
    return pearsonr(vector1[idx],vector2[idx])[0]

results = pandas.DataFrame(columns=relations_files,index=relations_files)
for relation_file1 in relations_files:
    print "Processing %s" %(relation_file1)
    for relation_file2 in relations_files:
        res = rsa(relation_file1,relation_file2)
        results.loc[relation_file1,relation_file2] = res
        results.loc[relation_file2,relation_file1] = res

results.to_csv("rsa_disorders.tsv",sep="\t")

# COMPARE TO CNP RSA ##########################################################################
# We can assess differences between the different contexts (eg disorders) by doing RSA analysis
# with the matrices


# The above is running, then do below.
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
  
        

    results.to_csv("rsa_disorders.tsv",sep="\t")

    # STOPPED HERE - we would want to show that this matrix is similar to matrix derived from data
    # RSA matrix for vectors SIMILAR TO RSA Matrix derived for behavior groups CNP

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
    threesds = 3*word_means.std()  


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


    


