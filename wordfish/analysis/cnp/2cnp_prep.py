'''
2rsa.py

Your wordfish project home directory should be defined as WORDFISH_HOME

This script will parse the raw CNP data and generate a matrix of pairwise similarity scores
for questions, subsetting to different groups of individuals with disorder labels. We can
then compare these different matrices with RSA to the same disorder models derived from
the boards.


First we will generate vectors that describe the CNP questions projected onto a corpus (eg,
a reddit board), meaning we have a vector for each question based on taking an average of words defined 
by the model. We now can generate similarity matrices for each set of vectors (comparing 
similarity of questionA to questionB in CNP) and then compare those matrices with
representational similarity analysis. This script will submit jobs to:

     Map CNP questions into the space of corpus word embeddings
     Save vectors
     Generate similarity matrix for QuestionA vs. Question B
     
These analyses have not been integrated into the wordfish pipeline, but will be if
they are useful.

'''

from wordfish.utils import mkdir
import pandas
import numpy
import pickle
import os
import re
import sys

base_dir = os.environ["WORDFISH_HOME"]
analysis_dir = mkdir("%s/analysis" %(base_dir))
vector_dir = mkdir("%s/vectors" %(analysis_dir))
behavior_dir = mkdir("%s/behavior" %(analysis_dir))

# Load CNP data
data1 = pandas.read_csv("%s/HTAC_Qry_1.csv" %behavior_dir,low_memory=False)
data2 = pandas.read_csv("%s/HTAC_Qry_2.csv" %behavior_dir,low_memory=False)
data = pandas.DataFrame.merge(data1,data2,on="ptid")
data.index = data.ptid
data = data.loc[data.ptid.isnull()==False]
data.columns = [x.upper() for x in data.columns]

# These labels are just missing the ASRS
to_renames = ["ORGANIZE", "REMAPPOINTMENT", "AVIODSTART", "FIDGET", "OVERACTIVE" "FINALDETAIL"]
for to_rename in to_renames:
    data = data.rename(columns={to_rename:'ASRS_%s' %(to_rename)})
data = data.rename(columns={u'OVERACTIVE':'ASRS_OVERACTIVE'})
data = data.rename(columns={u'FINALDETAIL':'ASRS_FINALDETAIL'})

# Save filtered versions.
data.to_pickle("%s/CNP_raw.pkl" %behavior_dir)

def filter_columns(df,exp):
    exp = re.compile(exp)
    rx = [x for x in df.columns if exp.match(x)]    
    return df[rx]

# Generate table with diagnoses ################################################################
scid = filter_columns(data,"SCID")
# Only keep diagnoses columns (with DX)
scid = filter_columns(scid,"SCID_DX[0-9]")

# Get unique labels
unique_diagnoses = []
for col in scid.columns:
    unique_diagnoses = unique_diagnoses + scid[col].unique().tolist()

unique_diagnoses = numpy.unique(unique_diagnoses).tolist()
unique_diagnoses.pop(unique_diagnoses.index("nan"))
unique_diagnoses.sort()

# Save list 
filey = open("%s/rx_list.txt" %behavior_dir,"w")
filey.writelines("\n".join(unique_diagnoses))
filey.close()

# Now put into a matrix
rx = pandas.DataFrame(0,columns=unique_diagnoses,index=data.index)

# Go through people and find their troubles!
for row in rx.iterrows():
    ptid = row[0]
    diagnoses = scid.loc[ptid]
    diagnoses = diagnoses[diagnoses.isnull()==False].tolist()
    rx.loc[ptid,diagnoses] = 1


# In the next script, we will subset data (and create a correlation matrix of questions)
# based on each rx group (or everyone)
rx.to_csv("%s/cnp_rx_matrix.tsv" %behavior_dir,sep="\t")
rx.to_pickle("%s/cnp_rx_matrix.pkl" %behavior_dir)
