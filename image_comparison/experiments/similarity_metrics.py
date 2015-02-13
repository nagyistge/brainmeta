#/usr/bin/python2

'''Similarity Metrics

For experiment1: all image input are nibabel, Z maps, already registered to brain mask

@author vsoch 
@data 2/2015
'''

import os
import sys,ctypes
_old_rtld = sys.getdlopenflags()
sys.setdlopenflags(_old_rtld|ctypes.RTLD_GLOBAL)
import numpy as np
import nibabel as nib
import sklearn.metrics as skm
import image_transformations as IT
from nilearn.masking import apply_mask
from scipy.spatial.distance import pdist
from nipype.interfaces.nipy.utils import Similarity
#--end other packages that need MKL
sys.setdlopenflags(_old_rtld)

# Mean Absolute Differences
# Activation Scores

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist

# METRIC EXTRACTION
def run_all(image1,image2,label1,label2,brain_mask,tmpdir):
  '''Extract all similarity metrics for image1 and image2. image2 should be the "gold standard", if applicable (the image we register to). Returns a dictionary of metric scores.'''

  # We will use a pairwise deletion mask
  pairwise_deletion_mask = IT.get_pairwise_deletion_mask(image1,image2,brain_mask)
  # If there is overlap:
  if len(np.unique(pairwise_deletion_mask.get_data())) != 1:
    data = apply_mask([image1,image2],pairwise_deletion_mask)
    pairwise_metrics = run_pairwise(data=data,image1=image1,image2=image2,brain_mask=pairwise_deletion_mask,
                                  label1=label1,label2=label2,tmpdir=tmpdir)    
    single_metrics = dict()
    single_metrics[label1] = run_single(data[0])
    single_metrics[label2] = run_single(data[1])
  
  else: 
    pairwise_metrics = get_none_metrics_pairwise()
    single_metrics = dict()
  return pairwise_metrics,single_metrics


def get_column_labels():
  return ["covariance","correlation_coefficient","correlation_ratio","correlation_ratio_norm",
          "mutual_information_norm","mutual_information","cosine",
          "activation_differences","euclidean","minkowski","cityblock","seuclidean",
          "sqeuclidean","chebyshev","canberra","braycurtis","skl_linear_kernel",
          "skl_l1","skl_l2","skl_sigmoid_kernel","skl_polynomial_kernel","skl_rbf_kernel",
          "hamming","yule","matching","dice","kulsinski","rogerstanimoto",
          "russellrao","sokalmichener"]

'''In the case that there is no overlap in the pairwise deletion mask, we return None for all'''
def get_none_metrics_pairwise():
  metric_names = get_column_labels()
  metrics = dict()
  for m in metric_names:
    metrics[m] = None
  return metrics

def run_single(image1):

  metrics = dict()
  metrics["standard_deviation"] = standard_deviation(image1)
  metrics["variance"] = variance(image1)

def run_pairwise(data,image1,image2,brain_mask,label1,label2,tmpdir):

  metrics = dict()

  # We need to generate temporary images to use Similarity module
  image1_tmp = IT.make_tmp_nii(image1,tmp_file_prefix="%s/%s_vs_%s1" %(tmpdir,label1,label2))
  image2_tmp = IT.make_tmp_nii(image2,tmp_file_prefix="%s/%s_vs_%s2" %(tmpdir,label1,label2))
  mask_tmp = IT.make_tmp_nii(brain_mask,tmp_file_prefix="%s/%s_%s_mask" %(tmpdir,label1,label2))
  
  # Comparison metrics [continuous]
  metrics["covariance"]  = covariance(data[0],data[1])
  metrics["correlation_coefficient"] = correlation_coefficient(image1_tmp,image2_tmp,mask_tmp)
  metrics["correlation_ratio"] = correlation_ratio(image1_tmp,image2_tmp,mask_tmp)
  metrics["correlation_ratio_norm"] = correlation_ratio_norm(image1_tmp,image2_tmp,mask_tmp)
  metrics["mutual_information_norm"] = mutual_information_norm(image1_tmp,image2_tmp,mask_tmp)
  metrics["mutual_information"] = mutual_information_norm(image1_tmp,image2_tmp,mask_tmp)
  metrics["cosine"] = cosine_metric(data[0],data[1])
  metrics["activation_differences"] = activation_differences(data[0],data[1])
  distances = ["euclidean","minkowski","cityblock","seuclidean","sqeuclidean",
               "chebyshev","canberra","braycurtis"]  
  for dist in distances:
    print "Calculating %s" %(dist)
    metrics["%s" %(dist)] = pdist(data,dist)[0]

  metrics["skl_linear_kernel"] = skm.pairwise.linear_kernel(data[0],data[1])[0][0]
  metrics["skl_l1"] = skm.pairwise_distances(data[0],data[1],metric="l1")[0][0]
  metrics["skl_l2"] = skm.pairwise_distances(data[0],data[1],metric="l2")[0][0]
  metrics["skl_sigmoid_kernel"] = skm.pairwise_kernels(data[0],data[1],metric="sigmoid")[0][0]
  metrics["skl_polynomial_kernel"] = skm.pairwise_kernels(data[0],data[1],metric="polynomial")[0][0]
  metrics["skl_rbf_kernel"] = skm.pairwise_kernels(data[0],data[1],metric="rbf")[0][0]

  # Comparison metrics [boolean]
  data = data.astype(bool).astype(int)
  distances = ["hamming","yule","matching","dice","kulsinski","rogerstanimoto",
              "russellrao","sokalmichener"]
  for dist in distances:
      metrics["%s" %(dist)] = pdist(data,dist)[0]

  if os.path.exists(image1_tmp): os.remove(image1_tmp)
  if os.path.exists(image2_tmp): os.remove(image2_tmp)
  if os.path.exists(mask_tmp): os.remove(mask_tmp)
  return metrics  


# SIMILARITY METRICS

# Covariance
def covariance(image1,image2):  
  '''Covariance is basically summing up the product of the differences of each value to the sample mean, and then dividing by N-1. In the case of a population, we would divide by N.'''
  print "Calculating covariance for %s and %s" %(image1,image2)
  return np.cov(image1,image2)[0][1]


# Variance
def variance(image1):
  '''Variance is average deviation from the mean.'''
  print "Calculating variance for %s" %(image1)
  return np.var(image1)


# Standard Deviation
def standard_deviation(image1):
  '''Standard deviation is (a human interpretable) version of average deviation from the mean (eg, back in same scale as data)'''
  print "Calculating standard deviation for %s" %(image1)
  return np.std(image1)


# Correlation Coefficient
def correlation_coefficient(image1_file,image2_file,mask_file):
  '''Correlation coefficient is ratio between covariance and product of standard deviations'''
  print "Calculating correlation coefficient for %s and %s" %(image1_file,image2_file)
  #return covariance(image1,image2) / (standard_deviation(image1) * standard_deviation(image2))
  try:
    similarity = Similarity()
    similarity.inputs.volume1 = image1_file
    similarity.inputs.volume2 = image2_file
    similarity.inputs.mask1 = mask_file
    similarity.inputs.mask2 = mask_file
    similarity.inputs.metric = 'cc'
    res = similarity.run() # doctest: +SKIP  
    return res.outputs.similarity
  except: return None

# Correlation Coefficient Without Centering (cosine)
def cosine_metric(image1,image2):
  print "Calculating cosine for %s and %s" %(image1,image2)
  X = [image1,image2]
  return pdist(X,"cosine")[0]


# Correlation Ratio
def correlation_ratio(image1_file,image2_file,mask_file):
  print "Calculating correlation ratio for %s and %s" %(image1_file,image2_file)
  try:
    similarity = Similarity()
    similarity.inputs.volume1 = image1_file
    similarity.inputs.volume2 = image2_file
    similarity.inputs.mask1 = mask_file
    similarity.inputs.mask2 = mask_file
    similarity.inputs.metric = 'cr'
    res = similarity.run() # doctest: +SKIP  
    return res.outputs.similarity
  except: return None

def correlation_ratio_norm(image1_file,image2_file,mask_file):
  '''L1 based correlation ratio'''
  print "Calculating correlation ratio norm for %s and %s" %(image1_file,image2_file)
  try:
    similarity = Similarity()
    similarity.inputs.volume1 = image1_file
    similarity.inputs.volume2 = image2_file
    similarity.inputs.mask1 = mask_file
    similarity.inputs.mask2 = mask_file
    similarity.inputs.metric = 'crl1'
    res = similarity.run() # doctest: +SKIP 
    return res.outputs.similarity
  except: return None

# Mutual Information
def mutual_information(image1_file,image2_file,mask_file):
  '''mutual information'''
  print "Calculating mutual information for %s and %s" %(image1_file,image2_file)
  try:
    similarity = Similarity()
    similarity.inputs.volume1 = image1_file
    similarity.inputs.volume2 = image2_file
    similarity.inputs.mask1 = mask_file
    similarity.inputs.mask2 = mask_file
    similarity.inputs.metric = 'mi'
    res = similarity.run() # doctest: +SKIP  
    return res.outputs.similarity
  except: return None


# Normalized Mutual Information
def mutual_information_norm(image1_file,image2_file,mask_file):
  '''normalized mutual information'''
  print "Calculating normalized mutual information for %s and %s" %(image1_file,image2_file)
  try:
    similarity = Similarity()
    similarity.inputs.volume1 = image1_file
    similarity.inputs.volume2 = image2_file
    similarity.inputs.mask1 = mask_file
    similarity.inputs.mask2 = mask_file
    similarity.inputs.metric = 'nmi'
    res = similarity.run() # doctest: +SKIP  
    return res.outputs.similarity
  except: return None

# Activation Scores
def activation_differences(image1,image2):
  '''
  http://vsoch.com/wiki/doku.php?id=summer_2011#section_1_-_overview_of_method (see end proposal)
  image2 is the original
  '''
  print "Calculating activation scores for %s and %s" %(image1,image2)
  try:
    template_in_mask = image2.astype(bool).astype(int)
    template_out_mask =  np.abs(template_in_mask - 1)
    image_mask = image1.astype(bool).astype(int)  
    activation_in_roi = np.sum(np.abs(image1)*template_in_mask) 
    activation_out_roi = np.sum(np.abs(image1)*template_out_mask)
    return (activation_in_roi/float(np.sum(template_in_mask))) - (activation_out_roi/float(np.sum(template_out_mask)))
  except: return None
