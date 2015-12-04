# This script will generate a similarity matrix for all CNP questions based on the actual
# data, for either all of the subjects/people, or a subset based on having a disorder label

from sklearn.preprocessing import scale
from scipy.stats import pearsonr
import pandas
import numpy
import pickle
import sys
import os
import re

r = sys.argv[1]
data_pkl = sys.argv[2]
disorder_pkl = sys.argv[3]
question_pkl = sys.argv[4]

sanity_check = False

data = pandas.read_pickle(data_pkl)
rx = pandas.read_pickle(disorder_pkl)
questions = pickle.load(open(question_pkl,"rb"))
behavior_dir = os.path.dirname(data_pkl)

# This means we don't filter down the data
if r != "all":
    # Find all subtypes of the disorder
    exp = re.compile("^%s" %r)
    disorders = [x for x in rx.columns if exp.match(x)]
    label =  rx[disorders].sum(axis=1)
    individuals = label.index[label==1].tolist()
    # Subset the data to those individuals
    data = data.loc[individuals,questions]
else:
    data = data.loc[:,questions]

# Keep track of how many we have
N = data.shape[0]

# Function to filter/parse a column
def parse_column(col):
    col = col[col.isnull()==False]
    col = col.astype(float)
    col = col[col.eq(-9999)==False]
    return col

# We need at least 2 people
if N > 1: 
    # Missing values of "." ... NO.
    data = data.replace(".",numpy.nan)
    subset = data.loc[:,~(data.isnull()).all()]
    if r == "all" and sanity_check == True:
        filey = open("%s/cnp_unique_values.txt" %(behavior_dir),"w")
        for c1 in subset.columns:
            col1 = subset[c1]
            col1 = parse_column(col1)
            filey.writelines("Unique values for %s are %s\n" %(c1,col1.unique().tolist()))
        filey.close()
    # Create a correlation matrix for subset that doesn't include missing or -9999 values
    corr = pandas.DataFrame(columns=subset.columns,index=subset.columns)
    count=1
    for c1 in subset.columns:
        print "Parsing %s, %s of %s" %(c1,count,subset.shape[1])
        col1 = subset[c1]
        col1 = parse_column(col1)
        for c2 in subset.columns:
            if c1!=c2: 
                col2 = subset[c2]
                col2 = parse_column(col2)
                if len(col2)>0 and len(col1)>0:
                    overlap = col1.index[col1.index.isin(col2.index)]
                    if len(overlap) > 1:
                        # Find the overlap
                        col1_overlap = col1[overlap]
                        col2_overlap = col2[overlap]
                        corr.loc[c1,c2] = pearsonr(col1_overlap,col2_overlap)[0]
                    else:
                        corr.loc[c1,c2] = 0
                else:
                    corr.loc[c1,c2] = 0    
            else:
                corr.loc[c1,c2] = 1    
        count+=1

# Save based on group number and N
corr.to_csv("%s/%s_%s_corr.tsv" %(behavior_dir,r,N),sep="\t")
