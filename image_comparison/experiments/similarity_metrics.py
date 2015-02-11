#/usr/bin/python2

'''Similarity Metrics

For experiment1: all image input are nibabel, Z maps, already registered to brain mask

@author vsoch 
@data 2/2015
'''

import os
import numpy as np
import nibabel as nib
import image_transformations as IT
from nilearn.masking import apply_mask
from scipy.spatial.distance import cdist
from nipype.interfaces.nipy.utils import Similarity

# Mean Absolute Differences
# Activation Scores

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist

# METRIC EXTRACTION
def run_all(image1,image2,label1,label2,brain_mask,tmpdir):
  '''Extract all similarity metrics for image1 and image2. image2 should be the "gold standard", if applicable (the image we register to). Returns a dictionary of metric scores.'''

  # We will use a pairwise deletion mask
  pairwise_deletion_mask = IT.get_pairwise_deletion_mask(image1,image2,brain_mask)
  data = apply_mask([image1,image2],pairwise_deletion_mask)
  pairwise_metrics = run_pairwise(data=data,image1=image1,image2=image2,brain_mask=pairwise_deletion_mask,
                                  label1=label1,label2=label2,tmpdir=tmpdir)    
  single_metrics = dict()
  single_metrics[label1] = run_single(data[0],label1)
  single_metrics[label2] = run_single(data[1],label2)
  return pairwise_metrics,single_metrics


def get_column_labels():
  return ["covariance","correlation_coefficient","correlation_ratio","correlation_ratio_norm",
          "mutual_information_norm","mutual_information","supervised_ll_ratio","cosine",
          "activation_differences","euclidean","minkowski","cityblock","seuclidean",
          "sqeuclidean","kulsinki","chebyshev","canberra","braycurtis","mahalanobis",
          "wminkowski","hamming","yule","matching","dice","kulsinski","rogerstanimoto",
          "russellrao","sokalmichener"]


def run_single(image1):

  metrics = dict()
  metrics["standard_deviation"] = standard_deviation(image1)
  metrics["variance"] = variance(image1)

def run_pairwise(data,image1,image2,brain_mask,label1,label2,tmpdir):

  metrics = dict()

  # We need to generate temporary images to use Similarity module
  image1_tmp = IT.make_tmp_nii(image1,"%s/%s.nii" %(tmpdir,label1))
  image2_tmp = IT.make_tmp_nii(image2,"%s/%s.nii" %(tmpdir,label2))
  mask_tmp = IT.make_tmp_nii(brain_mask,"%s/%s_%s_mask" %(tmpdir,label1,label2))
  
  # Comparison metrics [continuous]
  metrics["covariance"]  = covariance(data[0],data[1])
  metrics["correlation_coefficient"] = correlation_coefficient(image1_tmp,image2_tmp,mask_tmp)
  metrics["correlation_ratio"] = correlation_ratio(image1_tmp,image2_tmp,mask_tmp)
  metrics["correlation_ratio_norm"] = correlation_ratio_norm(image1_tmp,image2_tmp,mask_tmp)
  metrics["mutual_information_norm"] = mutual_information_norm(image1_tmp,image2_tmp,mask_tmp)
  metrics["mutual_information"] = mutual_information_norm(image1_tmp,image2_tmp,mask_tmp)
  metrics["supervised_ll_ratio"] = supervised_ll_ratio(image1_tmp,image2_tmp,mask_tmp)
  metrics["cosine"] = cosine_metric(data[0],data[1])
  metrics["activation_differences"] = activation_differences(data[0],data[1])
  distances = ["euclidean","minkowski","cityblock","seuclidean","sqeuclidean",
               "kulsinki","chebyshev","canberra","braycurtis","mahalanobis",
               "wminkowski"]  
  for dist in distances:
    metrics["%s" %(dist)] = pdist(data[0],data[1],dist)

  # Comparison metrics [boolean]
  data = data.astype(bool).astype(int)
  distances = ["hamming","yule","matching","dice","kulsinski","rogerstanimoto",
              "russellrao","sokalmichener"]
  for dist in distances:
      metrics["%s" %(dist)] = pdist(data[0],data[1],dist)

  os.remove(image1_tmp)
  os.remove(image2_tmp)
  os.remove(mask_tmp)
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
  similarity = Similarity()
  similarity.inputs.volume1 = image1_file
  similarity.inputs.volume2 = image2_file
  similarity.inputs.mask1 = mask_file
  similarity.inputs.mask2 = mask_file
  similarity.inputs.metric = 'cc'
  res = similarity.run() # doctest: +SKIP  
  return res.outputs.similarity


# Correlation Coefficient Without Centering (cosine)
def cosine_metric(image1,image2):
  print "Calculating cosine for %s and %s" %(image1,image2)
  X = [image1,image2]
  return pdist(X,"cosine")


# Correlation Ratio
def correlation_ratio(image1_file,image2_file,mask_file):
  print "Calculating correlation ratio for %s and %s" %(image1_file,image2_file)
  similarity = Similarity()
  similarity.inputs.volume1 = image1_file
  similarity.inputs.volume2 = image2_file
  similarity.inputs.mask1 = mask_file
  similarity.inputs.mask2 = mask_file
  similarity.inputs.metric = 'cr'
  res = similarity.run() # doctest: +SKIP  
  return res.outputs.similarity


def correlation_ratio_norm(image1_file,image2_file,mask_file):
  '''L1 based correlation ratio'''
  print "Calculating correlation ratio norm for %s and %s" %(image1_file,image2_file)
  similarity = Similarity()
  similarity.inputs.volume1 = image1_file
  similarity.inputs.volume2 = image2_file
  similarity.inputs.mask1 = mask_file
  similarity.inputs.mask2 = mask_file
  similarity.inputs.metric = 'crl1'
  res = similarity.run() # doctest: +SKIP  
  return res.outputs.similarity


# Mutual Information
def mutual_information(image1_file,image2_file,mask_file):
  '''mutual information'''
  print "Calculating mutual information for %s and %s" %(image1_file,image2_file)
  similarity = Similarity()
  similarity.inputs.volume1 = image1_file
  similarity.inputs.volume2 = image2_file
  similarity.inputs.mask1 = mask_file
  similarity.inputs.mask2 = mask_file
  similarity.inputs.metric = 'mi'
  res = similarity.run() # doctest: +SKIP  
  return res.outputs.similarity



# Normalized Mutual Information
def mutual_information_norm(image1_file,image2_file,mask_file):
  '''normalized mutual information'''
  print "Calculating normalized mutual information for %s and %s" %(image1_file,image2_file)
  similarity = Similarity()
  similarity.inputs.volume1 = image1_file
  similarity.inputs.volume2 = image2_file
  similarity.inputs.mask1 = mask_file
  similarity.inputs.mask2 = mask_file
  similarity.inputs.metric = 'nmi'
  res = similarity.run() # doctest: +SKIP  
  return res.outputs.similarity

  
def supervised_ll_ratio(image1_file,image2_file,mask_file):
  '''supervised log likihood ratio'''
  print "Calculating supervised ll ratio for %s and %s" %(image1_file,image2_file)
  similarity = Similarity()
  similarity.inputs.volume1 = image1_file
  similarity.inputs.volume2 = image2_file
  similarity.inputs.mask1 = mask_file
  similarity.inputs.mask2 = mask_file
  similarity.inputs.metric = 'slr'
  res = similarity.run() # doctest: +SKIP  
  return res.outputs.similarity


# Activation Scores
def activation_differences(image1,image2):
  '''
  http://vsoch.com/wiki/doku.php?id=summer_2011#section_1_-_overview_of_method (see end proposal)
  image2 is the original
  '''
  print "Calculating activation scores for %s and %s" %(image1,image2)
  template_in_mask = image2.astype(bool).astype(int)
  template_out_mask =  np.abs(template_mask - 1)
  image_mask = image1.astype(bool).astype(int)  
  activation_in_roi = np.sum(np.abs(image1)*template_in_mask) 
  activation_out_roi = np.sum(np.abs(image1)*template_out_mask)
  return (activation_in_roi/float(np.sum(template_in_mask))) - (activation_out_roi/float(np.sum(template_out_mask)))
