#/usr/bin/python2

'''Similarity Metrics

@author vsoch 
@data 2/2015
'''

import numpy as np
import nibabel as nib
from scipy.spatial.distance import pdist
from nipy.algorithms.registration import histogram_registration

# Mean Absolute Differences
# Activation Scores

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist

# IMAGE MANIPULATION
def load_image_data(image_file):
  mr = nib.load(image_file)
  return mr.get_data()

def flatten_image(image):
  return 

# METRIC EXTRACTION
def run_all(image1,image2,brainmask):
  '''Extract all similarity metrics for image1 and image2. image2 should be the "gold standard", if applicable (the image we register to). Returns a dictionary of metric scores.'''
  data1 = load_image_data(image1)
  data2 = load_image_data(image2)

  # Strategy 1: intersection of nonzero,non-nan voxels (pairwise deletion)

  # Strategy 2: everything in whole brain mask (brain mask)


# (intersection corresponds to an "AND conjunction")
bin_p_values = (log_p_values != 0)
mask_vt_filename = haxby_dataset.mask_vt[0]
vt = nibabel.load(mask_vt_filename).get_data().astype(bool)
bin_p_values_and_vt = np.logical_and(bin_p_values, vt)

plot_roi(nibabel.Nifti1Image(bin_p_values_and_vt.astype(np.int),
         fmri_img.get_affine()),
         mean_img, title='Intersection with ventral temporal mask',
         cut_coords=cut_coords)

  # Get intersection of nonzero voxels
  similarity_scores = dict()
  similarity_scores["covariance"] = covariance(image1,image2)
  similarity_scores["variance"] = variance(image1,image2)
   

# SIMILARITY METRICS

# Covariance
def covariance(image1,image2):  
  '''Covariance is basically summing up the product of the differences of each value to the sample mean, and then dividing by N-1. In the case of a population, we would divide by N.'''
  return np.cov(image1,image2)

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
  # return covariance(image1,image2) / (standard_deviation(image1) * standard_deviation(image2))
  # Note that this "implementation" is different from numpy, the formula above
  # returns something different.
  return histogram_registration(image1, image2, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='cc', interp='pv')

# Correlation Coefficient Without Centering (cosine ratio)
def cosine_ratio(from_image,image2):
  X = [from_image.flatten()]
  return pdist(

# Correlation Ratio
def correlation_ratio(image1,image2):
  return histogram_registration(image1, image2, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='cr', interp='pv')

def correlation_ratio_norm(image1,image2):
  '''L1 based correlation ratio'''
  return histogram_registration(image1, image2, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='crl1', interp='pv')

# Mutual Information
def mutual_information(image1,image2):
  '''mutual information'''
  return histogram_registration(image1, image2, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='mi', interp='pv')

# Normalized Mutual Information
def mutual_information_norm(image1,image2):
  '''normalized mutual information'''
  return histogram_registration(image1, image2, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='nmi', interp='pv')

def supervised_ll_ratio(image1,image2):
  '''supervised log likihood ratio'''
  return histogram_registration(image1, image2, from_bins=256, to_bins=None, from_mask=None, to_mask=None, similarity='nmi', interp='pv')

