#/usr/bin/python2

'''Similarity Metrics
These are hand coded, standard metrics, mostly for my own learning.
Yes, there are functions to do these things automatically (and more efficiently) 
I don't learn anything that way.

@author vsoch 
@data 2/2015
'''

import numpy as np


def covariance(X,Y,dof=1):  
  '''Covariance is basically summing up the product of the differences of each value to the sample mean, and then dividing by N-1. In the case of a population, we would divide by N.'''
  X = X.flatten(); Y=Y.flatten()
  if len(X) == len(Y):
    #cov = 0
    #for i in range(0,len(X)):
    #  cov += ( ( (X[i] - np.mean(X)) * (Y[i] - np.mean(Y) )) )
    #return cov / (len(X) - 1)
    return np.dot(X-np.mean(X),Y-np.mean(Y)) / (len(X)-dof)

def variance(X,dof=1):
  '''Variance is average deviation from the mean.'''
  #var = 0
  #for i in range(0,len(X)):
  #  var += np.square(X[i] - np.mean(X))
  #return var / len(X)
  #Matrix version of the above
  #X = X.flatten()
  return np.sum(np.square(X-np.mean(X))) / len(X-dof)

def standard_deviation(X):
  '''Standard deviation is (a human interpretable) version of average deviation from the mean (eg, back in same scale as data)'''
  return np.sqrt(variance(X))

def correlation_coefficient(X,Y):
  '''Correlation coefficient is ratio between covariance and product of standard deviations'''
  # return covariance(X,Y) / (standard_deviation(X) * standard_deviation(Y))
  # Note that this "implementation" is different from numpy, the formula above
  # returns something different.
  return np.corrcoef(X,Y)[0][1]

def correlation_ratio(X,Y):
  return covariance(X,Y) / np.sqrt(variance(X),variance(Y))

