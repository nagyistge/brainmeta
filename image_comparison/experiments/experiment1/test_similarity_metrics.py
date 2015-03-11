#!/usr/bin/env python2

# This script is to test similarity metrics. It will be in the format:
# Pybraincompare
# Some standard library

# Rest assured these are NOT used in the package, this is for learning purposes only!
import numpy as np
from similarity_metrics import covariance, variance, standard_deviation

X = np.array([1,2,3])
Y = np.array([4,5,6])

# Covariance: how much two independent random variables move together
covariance(X,Y)
np.cov(X,Y)[0][1]

# Variance: average deviation for a random variable from the mean
variance(X)
np.var(X)

# Standard deviation is also a measure of average variation from the mean, but back in the space of the data (so its more human interpretable)
standard_deviation(X)
np.std(X)

# Correlation coefficient is the ratio between covariance and product of standard deviations.
np.corrcoef(X,Y)[0][1]
correlation_coefficient(X,Y)


