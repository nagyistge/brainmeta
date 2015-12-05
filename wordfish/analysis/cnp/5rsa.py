'''
4wordfish_sims.py

Your wordfish project home directory should be defined as WORDFISH_HOME

We have already generated vectors that project CNP questions onto different corpus (reddit,
neurosynth), and then we parsed CNP data and made a matrix of diagnoses. We then made similarity matrices for the actual CNP question data (based on different disorder groups).
along with similarity matrices for the CNP questions mapped to the wordfish 
vectors. Time for RSA!

These analyses have not been integrated into the wordfish pipeline, but will be if
they are useful.

'''


# REPRESENTATIONAL SIMILARITY ANALYSIS ####################################################

from wordfish.utils import mkdir
from glob import glob
import pandas
import numpy
import pickle
import os
import re

base_dir = os.environ["WORDFISH_HOME"]

# Get analysis output directories
analysis_dir = mkdir("%s/analysis" %(base_dir))
scripts_dir = mkdir("%s/scripts" %(base_dir))
vectors_dir = mkdir("%s/vectors" %(analysis_dir))
rsa_dir = mkdir("%s/rsa" %(analysis_dir))

# Let's collect our similarity matrices:
matrices = {"question_scale":"%s/cnp_scale_sim.tsv" %behavior_dir,
            "assessment":"%s/cnp_assessment_sim.tsv" %behavior_dir}

vectors_sims = glob("%s/cnp_vectors_corr_wrt_*" %vectors_dir)
for vector_sim in vectors_sims:
    label = os.path.basename(vector_sim).split("_")[-1].split(".")[0]
    matrices["wordfish_%s"%label] = vector_sim

# Now add the matrices to compare the cnp questions (format is disorder ID, N)
question_sims = glob("%s/*corr.tsv" %behavior_dir)
for question_sim in question_sims:
    label = os.path.basename(question_sim).split("_")[0]
    matrices[label] = question_sim

# We will want to look up what the number labels correspond to
lookup = {"291":"cnp_alcohol_induced_mood_or_anxiety_disorder",
          "292":"cnp_drug_related_mood_or_anxiety_disorder",
          "293":"cnp_mood_disorder_due_to_medical_condition",
          "295":"cnp_schizophrenia",
          "296":"cnp_major_depressive_disorder_or_bipolar_disorder",
          "300":"cnp_anxiety_ocd_body_disorder_or_phobia",
          "303":"cnp_alcohol_dependence",
          "304":"cnp_drug_dependence",
          "305":"cnp_drug_abuse",
          "307":"cnp_eating_or_pain_disorder",
          "309":"cnp_ptsd_or_adjustment_disorder",
          "311":"cnp_depressive_disorder_nos",
          "314":"cnp_adhd",
          "799":"cnp_no_diagnosis_on_axis_i"}
pickle.dump(lookup,open("%s/dsm_rx_lookup.pkl" %behavior_dir,"wb"))
sims = dict()
for label,file_name in matrices.iteritems():
    if label in lookup:
        sims[lookup[label]] = file_name
    else:
        sims[label] = file_name

pickle.dump(sims,open("%s/rsa_matrix_lookup.pkl" %rsa_dir,"wb"))
    

# Write a function to perform RSA
def rsa(mfile1,mfile2):
    m1 = pandas.read_csv(mfile1,sep="\t",index_col=0)
    m2 = pandas.read_csv(mfile2,sep="\t",index_col=0)
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

results = pandas.DataFrame(columns=sims.keys(),index=sims.keys())
for label1,mfile1 in sims.iteritems():
    print "Processing %s" %(label1)
    for label2,mfile2 in sims.iteritems():
        res = rsa(mfile1,mfile2)
        results.loc[label1,label2] = res
        results.loc[label2,label1] = res
        

results.to_csv("%s/rsa_all_settona.tsv" %rsa_dir,sep="\t")
