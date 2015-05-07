#!/usr/bin/python

import numpy as np

def covariance(image1_vector,image2_vector):  
  '''Covariance is basically summing up the product of the differences of each value to the sample mean, and then dividing by N-1. In the case of a population, we would divide by N.'''
  return np.cov(image1_vector,image2_vector)[0][1]


# Variance
def variance(image1_vector):
  '''Variance is average deviation from the mean.'''
  return np.var(image1_vector)


# Standard Deviation
def standard_deviation(image1_vector):
  '''Standard deviation is (a human interpretable) version of average deviation from the mean (eg, back in same scale as data)'''
  return np.std(image1_vector)

# Weighted mean
def weighted_mean(image1_vector,weights):
    return np.dot(image1_vector,weights)/numpy.sum(weights)
   
# Weighted covariance
def weighted_covariance(image1_vector,image2_vector,weights):
     weighted_mean1 = weighted_mean(image1_vector,weights)
     weighted_mean2 = weighted_mean(image2_vector,weights)
     numerator = ((image1_vector-weighted_mean1)*(image2_vector-weighted_mean2)*weights).sum()
     denominator = numpy.sum(weights)
     return numerator/denominator

# Correlation Coefficient
def correlation_coefficient(image1_vector,image2_vector):
    '''Correlation coefficient is ratio between covariance and product of standard deviations'''
    return covariance(image1_vector,image2_vector) / (standard_deviation(image1_vector) * standard_deviation(image2_vector))


# Weighted Correlation Coefficient
def weighted_correlation_coefficient(image1_vector,image2_vector,weights):
    '''Correlation coefficient is ratio between covariance and product of standard deviations'''
    weighted_cov = weighted_covariance(image1_vector,image2_vector,weights)
    image1_cov = weighted_covariance(image1_vector,image1_vector,weights)
    image2_cov = weighted_covariance(image2_vector,image2_vector,weights)   
    return (weighted_cov / (numpy.sqrt((image1_cov*image2_cov))))
