'''
1run_prep_rsa.py

Your wordfish project home directory should be defined as WORDFISH_HOME

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


# REPRESENTATIONAL SIMILARITY ANALYSIS ####################################################

from wordfish.analysis import load_models
from wordfish.utils import mkdir
import os

base_dir = os.environ["WORDFISH_HOME"]

# Get analysis output directories
analysis_dir = mkdir("%s/analysis" %(base_dir))
vector_dir = mkdir("%s/vectors" %(analysis_dir))
scripts_dir = mkdir("%s/scripts" %(base_dir))

models = load_models(base_dir)

for model_id, model in models.iteritems():
    filey = "%s/.job/wordfish_%s.job" %(scripts_dir,model_id)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=%s\n" %(model_id))
    filey.writelines("#SBATCH --output=.out/%s.out\n" %(model_id))
    filey.writelines("#SBATCH --error=.out/%s.err\n" %(model_id))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python %s/1prep_rsa.py %s %s" %(scripts_dir,model_id,base_dir))
    filey.close()
    os.system("sbatch -p russpold " + "%s/.job/wordfish_%s.job" %(scripts_dir,model_id)) 
