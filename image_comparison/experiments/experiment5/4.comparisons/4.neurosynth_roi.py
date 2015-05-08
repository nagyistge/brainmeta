#!/usr/bin/python

# metric: ROI based on neurosynth

import os
import sys
import pandas
import numpy
import random
import nibabel
from scipy.stats import pearsonr
from nilearn.masking import apply_mask
from similarity_metrics import weighted_correlation_coefficient

contrast = sys.argv[1]
data_outfile_prefix = sys.argv[2]
standard = sys.argv[3]
subject_path_pkl = sys.argv[4]
raw_data_path = sys.argv[5]
metric = sys.argv[6]

subjects = pandas.read_pickle(subject_path_pkl).index

# Here is the neurosynth topic map
topic_map_path = "/scratch/PI/russpold/data/NEUROSYNTH/topic_maps/maps/topic8_pFgA_z_FDR_0.05.nii.gz"
topic_map = nibabel.load(topic_map_path)
brain_mask = nibabel.load(standard)

# Brain mask the topic map to include same number of voxels
topic_map = apply_mask([topic_map],brain_mask,ensure_finite=False)[0]
indexer = numpy.where(topic_map!=0)[0]

# Read in the raw data, filter to subject subset
matrix = pandas.read_pickle(raw_data_path)
matrix = matrix.loc[subjects]

# Filter the matrix to the voxels defined in the topic map
matrix = matrix.loc[:,indexer]

# Create matrix of scores - for weighted and non-weighted
scores = pandas.DataFrame(columns=subjects,index=subjects)
scores_weighted = pandas.DataFrame(columns=subjects,index=subjects)

# Prepare weights for weighted correlation
weights = topic_map[indexer]
weights = (weights - numpy.min(weights)) / (numpy.max(weights) - numpy.min(weights))

for s in range(0,len(scores.index)):
    sub1 = scores.index[s]
    for sub2 in scores.column:
        if sub1==sub2:
            scores.loc[sub1,sub2] = 1
        elif numpy.isnan(scores.loc[sub1,sub2]):
            # Get complete case mask
            pdmask = numpy.zeros(len(matrix.columns))
            image1_data = matrix.loc[sub1].tolist()
            image2_data = matrix.loc[sub2].tolist()
            pdmask[(numpy.squeeze(image1_data != 0)) * (numpy.isnan(numpy.squeeze(image1_data)) == False)] += 1
            pdmask[(numpy.squeeze(image2_data != 0)) * (numpy.isnan(numpy.squeeze(image2_data)) == False)] += 1
            pdmask[pdmask != 2] = 0
            pdmask[pdmask == 2] = 1
            image1_data = matrix.loc[sub1,numpy.where(pdmask==1)]
            image2_data = matrix.loc[sub2,numpy.where(pdmask==1)]
            # 1) Standard pearson correlation
            scores.loc[sub1,sub2] = pearsonr(image1_data,image2_data)[0]
            scores.loc[sub2,sub1] = pearsonr(image1_data,image2_data)[0]
            # 2) Weighted correlation
            score = weighted_correlation_coefficient(image1_data,image2_data,weights)
            scores_weighted.loc[sub1,sub2] = score
            scores_weighted.loc[sub2,sub1] = score
           
    # Save the similarity matrix to file, delete last version
    scores.to_csv("%s_%s_neurosynth_motor%s.tsv" %(data_outfile_prefix,metric,s),sep="\t")
    scores_weighted.to_csv("%s_%s_neurosynth_motor_weighted%s.tsv" %(data_outfile_prefix,metric,s),sep="\t")
    if os.path.exists("%s_%s_neurosynth_motor%s.tsv" %(data_outfile_prefix,metric,(s-1))):
        os.remove("%s_%s_neurosynth_motor%s.tsv" %(data_outfile_prefix,metric,(s-1)))
    if os.path.exists("%s_%s_neurosynth_motor_weighted%s.tsv" %(data_outfile_prefix,metric,(s-1))):
        os.remove("%s_%s_neurosynth_motor_weighted%s.tsv" %(data_outfile_prefix,metric,(s-1)))
