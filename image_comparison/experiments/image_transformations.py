#/usr/bin/python2

'''Image Transformations

The image should already be a Z score, nibabel object

@author vsoch 
@data 2/2015
'''

import sys,ctypes
_old_rtld = sys.getdlopenflags()
sys.setdlopenflags(_old_rtld|ctypes.RTLD_GLOBAL)
import numpy as np
import nibabel as nib
from nilearn.masking import apply_mask
from scipy.spatial.distance import pdist
from nipy.algorithms.registration import histogram_registration
#--end other packages that need MKL
sys.setdlopenflags(_old_rtld)

# Thresholding and Segmentation (return entire images) ------------------------------------------------

# Positive and negative thresholding
def threshold_abs(image1,thresholds=[0.0,0.5,1.0,1.5,1.65,1.7,1.75,1.8,1.85,1.9,1.96,2,2.58,3,3.5,4.0]):
  '''threshold the image at a range of Z score thresholds, including high positive and negative.'''
  thresholded = dict()
  data = image1.get_data()
  for thresh in thresholds:
    tmp = np.zeros(image1.shape)  
    tmp[np.abs(data) >= thresh] = data[np.abs(data) >= thresh]  
    new_image = nib.Nifti1Image(tmp,affine = image1.get_affine(),header=image1.get_header())
    thresholded[thresh] = new_image
  return thresholded


# Segment to group of anatomical ROIs
def anatomical_rois(image1,atlas):
  '''segment image for all regions defined by some atlas'''
  print "WRITE ME"


# Masking ---------------------------------------------------------------------------------------------

def get_pairwise_deletion_mask(image1,image2,mask):
  '''return pandas data frame with only intersection of brain masked, non zero voxels'''
  if image1.shape == image2.shape:
    pdmask = np.zeros(image1.shape)
    pdmask[(np.squeeze(image1.get_data() != 0)) * (np.isnan(np.squeeze(image1.get_data())) == False)] += 1
    pdmask[(np.squeeze(image2.get_data() != 0)) * (np.isnan(np.squeeze(image2.get_data())) == False)] += 1
    pdmask[pdmask != 2] = 0
    pdmask[pdmask == 2] = 1
    pdmask = np.logical_and(pdmask, mask.get_data()).astype(int)
    pdmask_img = nib.Nifti1Image(pdmask,affine=mask.get_affine(),header=mask.get_header())
    return pdmask_img    

# File Operations
def make_tmp_nii(image1,tmp_file_prefix):
  tmp_file = "%s.nii" %(tmp_file_prefix.replace(".","pt"))
  image1_tmp = nib.Nifti1Image(image1.get_data(),affine=image1.get_affine(),header=image1.get_header())
  nib.save(image1_tmp,tmp_file)
  return tmp_file
