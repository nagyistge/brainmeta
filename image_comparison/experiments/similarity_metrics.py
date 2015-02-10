#/usr/bin/python2

'''Similarity Metrics

For experiment1: all image input are nibabel, Z maps, already registered to brain mask

@author vsoch 
@data 2/2015
'''

import numpy as np
import nibabel as nib
import image_transformations as IT
from nilearn.masking import apply_mask
from scipy.spatial.distance import cdist
from nipy.algorithms.registration.histogram_registration import HistogramRegistration

# Mean Absolute Differences
# Activation Scores

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist

# METRIC EXTRACTION
def run_all(image1,image2,image1_label,image2_label,brain_mask):
  '''Extract all similarity metrics for image1 and image2. image2 should be the "gold standard", if applicable (the image we register to). Returns a dictionary of metric scores.'''

  strategies = dict()
  metrics = dict()

  # Strategy 1: intersection of nonzero,non-nan voxels (pairwise deletion)
  strategies["pairwise_deletion_mask"] = IT.get_pairwise_deletion_mask(image1,image2,brain_mask)

  # Strategy 2: everything in whole brain mask (use mask)
  strategies["inside_brain_mask"] = brain_mask
  
  # For each of the strategies, calculate similarity metrics
  for strategy,strategy_mask in strategies.iteritems():
    data = apply_mask([image1,image2],strategy_mask)
    label1 = "%s_%s" %(image1_label,strategy) 
    label2 = "%s_%s" %(image2_label,strategy) 
    
    # Comparison metrics [continuous]
    metrics["covariance_%s_vs_%s" %(label1,label2)] = covariance(data[0],data[1])
    metrics["correlation_coefficient_%s_vs_%s" %(label1,label2)] = correlation_coefficient(image1,image2,strategy_mask)
    metrics["correlation_ratio_%s_vs_%s" %(label1,label2)] = correlation_ratio(image1,image2,strategy_mask)
    metrics["correlation_ratio_norm_%s_vs_%s" %(label1,label2)] = correlation_ratio_norm(image1,image2,strategy_mask)
    metrics["mutual_information_norm_%s_vs_%s" %(label1,label2)] = mutual_information_norm(image1,image2,strategy_mask)
    metrics["mutual_information_%s_vs_%s" %(label1,label2)] = mutual_information_norm(image1,image2,strategy_mask)
    metrics["supervised_ll_ratio_%s_vs_%s" %(label1,label2)] = supervised_ll_ratio(image1,image2,strategy_mask)
    metrics["cosine_%s_vs_%s" %(label1,label2)] = cosine_metric(data[0],data[1])
    metrics["activation_differences_%s_vs_%s" %(label1,label2)] = activation_differences(data[0],data[1])
    distances = ["euclidean","minkowski","cityblock","seuclidean","sqeuclidean",
                 "kulsinki","chebyshev","canberra","braycurtis","mahalanobis",
                 "wminkowski"]  
    for dist in distances:
        metrics["%s_%s_vs_%s" %(dist,label1,label2)] = pdist(data[0],data[1],dist)

    # Individual metrics
    metrics["variance_%s" %(label1)] = variance(data[0])   
    metrics["variance_%s" %(label2)] = variance(data[1])   
    metrics["std_%s" %(label1)] = variance(data[0])   
    metrics["std_%s" %(label2)] = variance(data[1])   

    # Comparison metrics [boolean]
    data = data.astype(bool).astype(int)
    distances = ["hamming","yule","matching","dice","kulsinski","rogerstanimoto",
                 "russellrao","sokalmichener"]
    for dist in distances:
        metrics["%s_%s_vs_%s" %(dist,label1,label2)] = pdist(data[0],data[1],dist)

  return metrics  

# SIMILARITY METRICS

# Covariance
def covariance(image1,image2):  
  '''Covariance is basically summing up the product of the differences of each value to the sample mean, and then dividing by N-1. In the case of a population, we would divide by N.'''
  return np.cov(image1,image2)[0][1]

# Variance
def variance(image1):
  '''Variance is average deviation from the mean.'''
  return np.var(image1)

# Standard Deviation
def standard_deviation(image1):
  '''Standard deviation is (a human interpretable) version of average deviation from the mean (eg, back in same scale as data)'''
  return np.std(image1)

# Correlation Coefficient
def correlation_coefficient(image1,image2):
  '''Correlation coefficient is ratio between covariance and product of standard deviations'''
  #return covariance(image1,image2) / (standard_deviation(image1) * standard_deviation(image2))
  mi = HistogramRegistration(image2, image1, similarity='cc')  
  T = mi.optimize("affine")
  return mi.explore(T)[0][0]

# Correlation Coefficient Without Centering (cosine)
def cosine_metric(image1,image2):
  X = [image1,image2]
  return pdist(X,"cosine")

# Correlation Ratio
def correlation_ratio(image1,image2,mask):
  mi = HistogramRegistration(image2, image1, similarity='cr',from_mask=mask,to_mask=mask)  
  T = mi.optimize("affine")
  return mi.explore(T)[0][0]

def correlation_ratio_norm(image1,image2,mask):
  '''L1 based correlation ratio'''
  mi = HistogramRegistration(image2, image1, similarity='crl1',from_mask=mask,to_mask=mask)  
  T = mi.optimize("affine")
  return mi.explore(T)[0][0]

# Mutual Information
def mutual_information(image1,image2,mask):
  '''mutual information'''
  mi = HistogramRegistration(image2, image1, similarity='mi',from_mask=mask,to_mask=mask)  
  T = mi.optimize("affine")
  return mi.explore(T)[0][0]

# Normalized Mutual Information
def mutual_information_norm(image1,image2,mask):
  '''normalized mutual information'''
  mi = HistogramRegistration(image2, image1, similarity='nmi',from_mask=mask,to_mask=mask)  
  T = mi.optimize("affine")
  return mi.explore(T)[0][0]
  
def supervised_ll_ratio(image1,image2,mask):
  '''supervised log likihood ratio'''
  mi = HistogramRegistration(image2, image1, similarity='slr',from_mask=mask,to_mask=mask)  
  T = mi.optimize("affine")
  return mi.explore(T)[0][0]

# Activation Scores
def activation_differences(image1,image2):
  '''
  http://vsoch.com/wiki/doku.php?id=summer_2011#section_1_-_overview_of_method (see end proposal)
  image2 is the original
  '''
  template_in_mask = image2.astype(bool).astype(int)
  template_out_mask =  np.abs(template_mask - 1)
  image_mask = image1.astype(bool).astype(int)  
  activation_in_roi = np.sum(np.abs(image1)*template_in_mask) 
  activation_out_roi = np.sum(np.abs(image1)*template_out_mask)
  return (activation_in_roi/float(np.sum(template_in_mask))) - (activation_out_roi/float(np.sum(template_out_mask)))
