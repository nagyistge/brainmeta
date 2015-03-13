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
#--end other packages that need MKL
sys.setdlopenflags(_old_rtld)

# Convert to Z Scores (return entire images) ------------------------------------------------
# Note: not comfortable doing this, needs discussion
def to_Z(image1):
  data = image1.get_data()
  mask = np.zeros(data.shape)
  mask[data!=0] = 1
  themean = np.mean(data[mask==1])
  thesd = np.std(data[mask==1])
  Z = np.zeros(data.shape)
  Z[mask==1] = (data[mask==1] - themean) / thesd
  Z_nii = nib.nifti1.Nifti1Image(Z,affine=image1.get_affine(),header=image1.get_header())
  return Z_nii

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

def threshold_pos(image1,thresholds=[0.0,0.5,1.0,1.5,1.65,1.7,1.75,1.8,1.85,1.9,1.96,2,2.58,3,3.5,4.0]):
  '''threshold the image at a range of Z score thresholds, including only high positive values.'''
  thresholded = dict()
  data = image1.get_data()
  for thresh in thresholds:
    tmp = np.zeros(image1.shape)  
    tmp[data >= thresh] = data[data >= thresh]  
    new_image = nib.Nifti1Image(tmp,affine = image1.get_affine(),header=image1.get_header())
    thresholded[thresh] = new_image
  return thresholded


# Segment to only include some region of interest
def get_masked_images(images,roi):
  '''segment image for all regions defined by some roi
     roi should be the data matrix, not a nibabel image
  '''
  if not isinstance(images,list): images = [images]
  masked = []
  for image in images:
    tmp = np.zeros(image.shape)
    tmp[roi==1] = image.get_data()[roi==1]
    new_img = nib.Nifti1Image(tmp,header=image.get_header(),affine=image.get_affine())
    masked.append(new_img)
  return masked


# Masking ---------------------------------------------------------------------------------------------

# Mask includes intersection of nonzero, non-nan voxels
def get_pairwise_deletion_mask(image1,image2,mask):
  '''return pandas data frame with only intersection of brain masked, non zero voxels'''
  if image1.shape == image2.shape:
    image1_data = image1.get_data()
    image2_data = image2.get_data()
    pdmask = np.zeros(image1.shape)
    pdmask[(np.squeeze(image1_data != 0)) * (np.isnan(np.squeeze(image1_data)) == False)] += 1
    pdmask[(np.squeeze(image2_data != 0)) * (np.isnan(np.squeeze(image2_data)) == False)] += 1
    pdmask[pdmask != 2] = 0
    pdmask[pdmask == 2] = 1
    pdmask = np.logical_and(pdmask, mask.get_data()).astype(int)
    pdmask_img = nib.Nifti1Image(pdmask,affine=mask.get_affine(),header=mask.get_header())
    return pdmask_img    


# Mask includes union of nonzero voxels in both images
def get_pairwise_inclusion_mask(image1,image2,mask,absolute_value=False):
  '''return pandas data frame with union of brain masked, non zero voxels'''
  if image1.shape == image2.shape:
    image1_data = image1.get_data()
    image2_data = image2.get_data()
    pimask = np.zeros(image1.shape)
    pimask[(np.squeeze(image1_data != 0)) * (np.isnan(np.squeeze(image1_data)) == False)] += 1
    pimask[(np.squeeze(image2_data != 0)) * (np.isnan(np.squeeze(image2_data)) == False)] += 1
    pimask[pimask != 0] = 1
    pimask = np.logical_and(pimask, mask.get_data()).astype(int)
    pimask_img = nib.Nifti1Image(pimask,affine=mask.get_affine(),header=mask.get_header())
    return pimask_img    


# File Operations
def make_tmp_nii(image1,tmp_file_prefix):
  tmp_file = "%s.nii" %(tmp_file_prefix.replace(".","pt"))
  image1_tmp = nib.Nifti1Image(image1.get_data(),affine=image1.get_affine(),header=image1.get_header())
  nib.save(image1_tmp,tmp_file)
  return tmp_file
