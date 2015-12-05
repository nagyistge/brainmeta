'''
3run_cnp_sims.py

Your wordfish project home directory should be defined as WORDFISH_HOME

We have already generated vectors that project CNP questions onto different corpus (reddit,
neurosynth), and then we parsed CNP data and made a matrix of diagnoses. Now we want to 
generate similarity matrices for the actual CNP question data, for everyone, and for
different subsets of individuals (with different diagnoses). This script will submit jobs to
do that.

These analyses have not been integrated into the wordfish pipeline, but will be if
they are useful.

'''


# REPRESENTATIONAL SIMILARITY ANALYSIS ####################################################

from wordfish.utils import mkdir
import pandas
import numpy
import pickle
import os
import re

base_dir = os.environ["WORDFISH_HOME"]

# Get analysis output directories
analysis_dir = mkdir("%s/analysis" %(base_dir))
behavior_dir = mkdir("%s/behavior" %(analysis_dir))
scripts_dir = mkdir("%s/scripts" %(base_dir))
vectors_dir = mkdir("%s/vectors" %(analysis_dir))

# We will filter down to behavioral questions defined in CNP that were mappable to either nsyn or reddit
vectors_reddit = pandas.read_csv("%s/cnp_vectors_wrt_reddit.tsv" %vectors_dir,sep="\t",index_col=0)
vectors_nsyn = pandas.read_csv("%s/cnp_vectors_wrt_neurosynth.tsv" %vectors_dir,sep="\t",index_col=0)
questions = numpy.unique(vectors_reddit.index.tolist() + vectors_nsyn.index.tolist()).tolist()
question_pkl = "%s/cnp_question_labels.pkl" %behavior_dir
pickle.dump(questions,open(question_pkl,"wb"))

# Read in the disorder rx data frame
disorder_pkl = "%s/cnp_rx_matrix.pkl" %behavior_dir
data_pkl = "%s/CNP_raw.pkl" %behavior_dir
rx = pandas.read_pickle(disorder_pkl) 

# Function to submit job
def submit_job(scripts_dir,data_pkl,disorder_pkl,question_pkl,r):
    filey = "%s/.job/wordfish_%s.job" %(scripts_dir,r)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=wordfish_%s\n" %(r))
    filey.writelines("#SBATCH --output=.out/%s.out\n" %(r))
    filey.writelines("#SBATCH --error=.out/%s.err\n" %(r))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python %s/3cnp_sims.py %s %s %s %s" %(scripts_dir,r,data_pkl,disorder_pkl,question_pkl))
    filey.close()
    os.system("sbatch -p russpold " + "%s/.job/wordfish_%s.job" %(scripts_dir,r)) 

# We will first not take subtypes into account, to maximize data, and do based on disorder groups
exp = re.compile("^[0-9]{3}")
groups = numpy.unique([y[0] for y in [exp.findall(x) for x in rx.columns if exp.match(x)]]).tolist()
# ['291', '292', '293', '295', '296', '300', '303', '304', '305', '307', '309', '311', '314', '799']
for g in groups:
    submit_job(scripts_dir,data_pkl,disorder_pkl,question_pkl,g)

# Now let's defined our own groups based on labels
groups = ["anxiety","depress","dependence","abuse","bipolar","schizo","attention","alcohol","cannabis"]
for g in groups:
    submit_job(scripts_dir,data_pkl,disorder_pkl,question_pkl,g)

# Run one script to do for all
submit_job(scripts_dir,data_pkl,disorder_pkl,question_pkl,"all")

# Let's also generate matrices to compare similarity of "meta stuffs" about the questions
meta = pandas.read_csv("%s/cogpheno_739.tsv" %behavior_dir,sep="\t")
meta = meta[meta.question_label.isin(questions)]
meta.index = meta.question_label

# We will look at the following:
#   Similar scales
#   assessment name

def make_meta_simmatrix(field_name):
    # Similar scales
    sim = pandas.DataFrame(columns=meta.index,index=meta.index)
    for c in range(meta.shape[0]):
        c1 = meta.index[c]
        print "Parsing %s, %s of %s" %(c1,c,meta.shape[0])
        for d in range(meta.shape[0]):
            c2 = meta.index[d]
            if c1 <= c2:
                if meta.loc[c1,field_name] == meta.loc[c2,field_name]:
                    sim.loc[c1,c2] = 1
                    sim.loc[c2,c1] = 1
                else:
                    sim.loc[c1,c2] = 0
                    sim.loc[c2,c1] = 0
    return sim

sim = make_meta_simmatrix("question_options")
sim.to_csv("%s/cnp_scale_sim.tsv" %behavior_dir,sep="\t")
sim = make_meta_simmatrix("assessment_name")
sim.to_csv("%s/cnp_assessment_sim.tsv" %behavior_dir,sep="\t")
