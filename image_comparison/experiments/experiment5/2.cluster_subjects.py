#!/usr/bin/python

# Here we want to, for each motor contrast, find a group of subjects that are most closely related. We will use kmeans clustering (and feature selection, if necessary)

from similarity_metrics import weighted_correlation_coefficient
from scipy.stats import pearsonr
from nilearn.masking import apply_mask
import numpy
import pandas
import nibabel
import sys
import os

contrast = sys.argv[1]
data_outfile_prefix = sys.argv[2]
standard = sys.argv[3]
subject_path_pkl = sys.argv[4]
do_data_prep = bool(sys.argv[5])

# This is a matrix of all cope image paths, cols are contrasts, rows subjects
motor_ss_maps = pandas.read_pickle(subject_path_pkl)

# We will create a brain mask
brain_mask = nibabel.load(standard)
number_voxels = len(brain_mask.get_data()[brain_mask.get_data()==1])
paths = motor_ss_maps[contrast]

# Data preparation
if do_data_prep == True:
    matrix = pandas.DataFrame(index=motor_ss_maps.index,columns=range(0,number_voxels))
    for p in range(0,len(paths)):
        print "Processing %s of %s" %(p,len(paths))
        path = paths[p]
        subid = paths.index[p]
        mr = nibabel.load(path)
        masked = apply_mask([mr],brain_mask,ensure_finite=False)[0]
        matrix.loc[subid] = masked

    # Save to data directory
    matrix.to_pickle("%s_brainmask.pkl" %(data_outfile_prefix))
else:
    matrix = pandas.read_pickle("%s_brainmask.pkl" %(data_outfile_prefix))    

# Now calculate similarity between all subjects
sim_matrix = pandas.DataFrame(index=motor_ss_maps.index,columns=motor_ss_maps.index) # using all voxels
sim_matrix_weighted_mean = pandas.DataFrame(index=motor_ss_maps.index,columns=motor_ss_maps.index) # weighted voxels, mean of 2 images
sim_matrix_weighted_min = pandas.DataFrame(index=motor_ss_maps.index,columns=motor_ss_maps.index) # weighted voxels, min of two images

for s in range(0,len(motor_ss_maps.index)):
    sub1 = motor_ss_maps.index[s]
    print "Processing %s of %s" %(s,len(motor_ss_maps.index)) 
    for sub2 in motor_ss_maps.index:
        if sub1==sub2:
            sim_matrix.loc[sub1,sub2] = 1
            sim_matrix_weighted_mean.loc[sub1,sub2] = 1
            sim_matrix_weighted_min.loc[sub1,sub2] = 1
        elif numpy.isnan(sim_matrix.loc[sub1,sub2]):
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
            # 1) Standard pearson correlation for all voxels
            sim_matrix.loc[sub1,sub2] = pearsonr(image1_data,image2_data)[0]
            sim_matrix.loc[sub2,sub1] = pearsonr(image1_data,image2_data)[0]
            # 2) Weighted correlation - mean
            weights = (image1_data+image2_data) / 2
            # Normalize between 0 and 1
            weights = (weights - numpy.min(weights)) / (numpy.max(weights) - numpy.min(weights))
            score = weighted_correlation_coefficient(image1_data,image2_data,weights)
            sim_matrix_weighted_mean.loc[sub1,sub2] = score
            sim_matrix_weighted_mean.loc[sub2,sub1] = score
            # 3) Weighted correlation - min            
            pair = pandas.DataFrame()                       
            pair[sub1] = image1_data
            pair[sub2] = image2_data   
            weights = pair.min(axis=1)
            # Normalize between 0 and 1
            weights = (weights - numpy.min(weights)) / (numpy.max(weights) - numpy.min(weights))
            score = weighted_correlation_coefficient(image1_data,image2_data,weights)
            sim_matrix_weighted_min.loc[sub1,sub2] = score
            sim_matrix_weighted_min.loc[sub2,sub1] = score
    # Save the similarity matrices to file, delete last version
    sim_matrix.to_csv("%s_pearsonr%s.tsv" %(data_outfile_prefix,s),sep="\t")
    sim_matrix_weighted_min.to_csv("%s_weightmin%s.tsv" %(data_outfile_prefix,s),sep="\t")
    sim_matrix_weighted_mean.to_csv("%s_weightmean%s.tsv" %(data_outfile_prefix,s),sep="\t")
    if os.path.exists("%s_pearsonr%s.tsv" %(data_outfile_prefix,(s-1))):
        os.remove("%s_pearsonr%s.tsv" %(data_outfile_prefix,(s-1)))
    if os.path.exists("%s_weightmin%s.tsv" %(data_outfile_prefix,(s-1))):
        os.remove("%s_weightmin%s.tsv" %(data_outfile_prefix,(s-1)))
    if os.path.exists("%s_weightmean%s.tsv" %(data_outfile_prefix,(s-1))):
        os.remove("%s_weightmean%s.tsv" %(data_outfile_prefix,(s-1)))
