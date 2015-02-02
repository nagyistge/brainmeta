#/usr/bin/python2

'''Similarity Metrics

@author vsoch 
@data 2/2015
'''

import numpy as np
from nipy.algorithms.registration import histogram_registration

def covariance(from_img,to_img):  
  '''Covariance is basically summing up the product of the differences of each value to the sample mean, and then dividing by N-1. In the case of a population, we would divide by N.'''
  return np.cov(from_img,to_img)

def variance(from_img):
  '''Variance is average deviation from the mean.'''
  return np.var(from_img)

def standard_deviation(from_img):
  '''Standard deviation is (a human interpretable) version of average deviation from the mean (eg, back in same scale as data)'''
  return np.std(from_img)

def correlation_coefficient(from_img,to_img):
  '''Correlation coefficient is ratio between covariance and product of standard deviations'''
  # return covariance(from_img,to_img) / (standard_deviation(from_img) * standard_deviation(to_img))
  # Note that this "implementation" is different from numpy, the formula above
  # returns something different.
  return histogram_registration(from_img, to_img, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='cc', interp='pv')

def correlation_ratio(from_img,to_img):
  return histogram_registration(from_img, to_img, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='cr', interp='pv')

def correlation_ratio_norm(from_img,to_img):
  '''L1 based correlation ratio'''
  return histogram_registration(from_img, to_img, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='crl1', interp='pv')

def mutual_information(from_img,to_img):
  '''mutual information'''
  return histogram_registration(from_img, to_img, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='mi', interp='pv')

def mutual_information_norm(from_img,to_img):
  '''normalized mutual information'''
  return histogram_registration(from_img, to_img, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='nmi', interp='pv')

def supervised_ll_ratio(from_img,to_img):
  '''supervised log likihood ratio'''
  return histogram_registration(from_img, to_img, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='nmi', interp='pv')
