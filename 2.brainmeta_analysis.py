#!/usr/bin/env python2

# DATA STRUCTURES -------------------------------------------------------

import json
import math
import nibabel
import pickle
import numpy as np
import pandas as pd
from annotate import ontology
from nilearn.masking import apply_mask

# This defines a "mock" ontology for the openfmri proof of concept set
relationship_table = "/home/vanessa/Documents/Work/BRAINMETA/mr/cogatlas_annotate_triples.tsv"

# Create a data structure for data analysis
data_structure = ontology.named_ontology_tree_from_tsv(relationship_table)

# Here is top path for images defined above (resampled to MNI 2mm)
mrpath = "/home/vanessa/Documents/Work/BRAINMETA/mr/resampled"

# Our data structure goes into a json object!
myjson = json.loads(data_structure)


# HIGH LEVEL IDEA ------------------------------------------------------

# We want to generate a reverse inference score for EVERY SINGLE NODE in this map.  Each node is a grouping for the maps defined below it! Let's first write a bunch of functions for working with our data structures, calculating priors, and calculating reverse inference scores for some group of brain maps. I think Russ will be happier with voxelwise analysis, (vs regional) so I'll start with that.

# Steps ----------------------------
# calculate priors for each node, meaning p(activation in threshold [1..n] | region (voxel))
# We must do this for images within (under) node [priors_in], and for all others [priors_out]
# Calculate reverse inference scores using these priors!
# Add scores to our visualization

# Researcher Workflow --------------
# I have a new result map!
# [I don't give a label]: Find best matching (most similar image) node in map
# [I do give a label]: Direct me to reverse inference score for that node

# Database Workflow ----------------
# priors_in and priors_out must be updated, and reverse_inference values recalculated at some frequency when new data is added

# Need to think about --------------
# What if some maps are thresholded? (because they ARE)
# Need to check for any redundancy in groups (if we can assign an image to more than one place in ntology, and it's in both groups, how to choose?)
# What is most efficient data structure storage / calculation strategy for these?


# FUNCTIONS ------------------------------------------------------------


# This function will "grab" a node by the name ( a subset of the tree )
def locate_by_name(myjson,name):
    if myjson.get('name',None) == name:
      return myjson
    for child in myjson.get('children',[]):
      result = locate_by_name(child,name)
      if result is not None:
        return result
    return None

# This function will get all the names of the nodes
def get_node_names(myjson,names=[]):
  if myjson.get('name',None) == None:
    return names
  else: names.append(myjson.get('name'))
  for child in myjson.get('children',[]):
    names = get_node_names(child,names)
    if not names: return None
  return names

# This function will return all nifti images for some json tree
def get_nifti_names(myjson,nii_files=[]):
  import re
  expression = re.compile(".nii.gz")
  if myjson.get('name',None) == None:
    return nii_files
  else: 
    if expression.search(myjson.get('name')):
      nii_files.append(myjson.get('name'))
  for child in myjson.get('children',[]):
    nii_files = get_nifti_names(child,nii_files)
    if not nii_files: return None
  return nii_files

# Function to calculate priors from a regionally-based df for ranges of values
# ranges_df should have columns ["start","stop"]
# ranges_df.index (row names) should correspond to region name!
def calculate_regional_priors_in_ranges(region_df,ranges_df):
  # A table of priors, columns --> thresholds, rows --> regions)
  priors = pd.DataFrame(columns=ranges_df.index)  # defined for images in nodeg
  for row in ranges_df.iterrows(): 
    # Nan means that value does not pass
    bool_df = region_df[(region_df >= row[1].start) & (region_df <=row[1].stop)]    
    # [Numerator] Count the non-NaN values in each column, add 1 for laplace smoothing
    # [Denominator] sum(numerator with no laplace) + V (words in vocabulary --> regions/voxels)
    # [Overall] probability that the image voxel (region) is in the range given entire image set!
    numerator = (bool_df.shape[0] - bool_df.isnull().sum())
    numerator_laplace_smoothed = numerator + 1
    denominator = np.sum(numerator) + bool_df.shape[0]
    priors[row[0]] = numerator_laplace_smoothed / denominator
  return priors


# Function to save brain images for each level of priors
def priors_to_brain_image(priors_df,output_name,mask_image):
  import numpy as np
  mask_image = nibabel.load(mask_image)
  roi = mask_image.get_data()
  # This will hold a 4D object for the data
  new_dim = mask_image.get_shape()
  new_dim = [new_dim[0],new_dim[1],new_dim[2],len(priors_df.columns)]
  data = np.zeros(shape=new_dim)
  for t in range(0,len(priors_df.columns)):
    thresh = priors_df.columns[t]
    new_image = np.zeros(shape=mask_image.get_shape())
    new_image[np.where(roi!=0)] = priors_df[thresh]
    data[:,:,:,t] = new_image
  final_img = nibabel.Nifti1Image(data, affine=mask_image.get_affine())
  final_img.set_filename("%s" %(output_name))
  nibabel.save(final_img,output_name)


# Function to return reverse inference value for each threshold in priors matrix
def calculate_reverse_inference_threshes(p_in,p_out,num_in,num_out):
  import numpy as np
  total = in_count + out_count # total number of nifti images
  p_process_in = float(in_count) / total   # percentage of niftis in
  p_process_out = float(out_count) / total # percentage out
  # If we multiply, we will get 0, so we take sum of logs
  p_in_log = np.log(p_in)
  p_out_log = np.log(p_out)
  numerators = p_in_log.sum(axis=0) * p_process_in
  denominators = (p_in_log.sum(axis=0) * p_process_in) + (p_out_log.sum(axis=0) * p_process_out)
  return (numerators / denominators)

# Function to return reverse inference value based on particular thresholds of a brain stat map (or average of the node set) - this will return one value!
def calculate_reverse_inference(stat_map,ranges_df,p_in,p_out,num_in,num_out):
  import numpy as np, pandas as pd
  total = in_count + out_count # total number of nifti images
  p_process_in = float(in_count) / total   # percentage of niftis in
  p_process_out = float(out_count) / total # percentage out
  # If we multiply, we will get 0, so we take sum of logs
  # For each ACTUAL voxel value, assign to its threshold (all thresholds should be represented)
  stat_map_levels = stat_map.copy()
  for row in ranges_df.iterrows():
    stat_map_levels.loc[(stat_map >= row[1].start) & (stat_map <=row[1].stop)] = row[0]
  # BUG WILL BE CHANGED: For now we will put values slightly above threshold in upper/lower - this needs
  idx_upper =  [x for x in range(0,len(stat_map_levels)) if (stat_map_levels[x] not in range_table.index) and stat_map_levels[x] > 3.0]
  idx_lower =  [x for x in range(0,len(stat_map_levels)) if (stat_map_levels[x] not in range_table.index) and stat_map_levels[x] < -8.0]
  stat_map_levels.loc[idx_upper] = "[2.5,3.0]"
  stat_map_levels.loc[idx_lower] = "[-8.0,-7.5]"
  # Now use actual threshold labels to choose the appropriate probability values for each voxel
  # I need advice how to make this faster, and this should be a separate function
  p_in_vector = []; p_out_vector = [];
  for v in range(0,len(stat_map_levels)):
    level = stat_map_levels[v] 
    p_in_vector.append(p_in.loc[v,level])
    p_out_vector.append(p_out.loc[v,level]) 
  p_in_log = np.log(p_in_vector)
  p_out_log = np.log(p_out_vector)
  numerator = p_in_log.sum(axis=0) * p_process_in
  denominator = (p_in_log.sum(axis=0) * p_process_in) + (p_out_log.sum(axis=0) * p_process_out)
  return (numerator / denominator)

# DATA PREPARATION ------------------------------------------------------

# Get the names of the nodes, our "groups" of images
names = get_node_names(myjson,[])

# Brain mask
brain_mask = "/usr/share/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz"

# First create table of voxels for all images, the top node, BRAINMETA
node = locate_by_name(myjson, names[0])
nifti_files = get_nifti_names(node,[]) # This is all the images
full_paths = ["%s/%s" %(mrpath,x.split(":")[1]) for x in nifti_files]
mr = apply_mask(full_paths, brain_mask, dtype='f', smoothing_fwhm=None, ensure_finite=True)  
mrtable = pd.DataFrame(mr)
mrtable.index = nifti_files	
# Need to double check this with Russ, filling Na == 0
# It makes sense because a NaN in an image is likely some division by 0?
mrtable.fillna(0) 
# Now calculate the intervals 
mins = [np.min(x) for x in mr]
maxs = [np.max(x) for x in mr] 
min = math.floor(np.min(mins)); max = math.ceil(np.min(maxs))
steps = ((abs(min)+max)*2)+1
breaks = np.linspace(min,max,num=steps,retstep=True)
breaks = breaks[0]

# Make a table of the ranges, as numeric
# Row name corresponds to the threshold name
range_table = pd.DataFrame(columns=["start","stop"])
for s in range(0,len(breaks)-1):
  start = breaks[s]
  stop = breaks[s+1]
  name = "[%s,%s]" %(start,stop)
  range_table.loc[name] = [start,stop]


# PRIORS CALCULATION ----------------------------------------------------
 
# Calculate priors for image sets at each node (**node names must be unique) 
# QUESTION FOR RUSS AND NOLAN: this is images in nodes vs EVERYTHING else, not just those
# on the same level. We may just want those on the same level

# We will save a dictionary of priors, columns --> thresholds, rows --> voxels)
priors_in = dict()  # priors for images included in node
priors_out = dict()  # priors for images NOT included

# Now iterate through nodes (names), and calculate prior table for each
for n in range(0,len(names)):
  name = names[n]
  print "Processing %s, %s of %s" %(name,n+1,len(names))
  # Filter table to those images
  in_node = locate_by_name(myjson, name)
  in_nifti_subset = get_nifti_names(in_node,[])
  out_nifti_subset = [i for i in nifti_files if i not in in_nifti_subset]  
  in_subset = mrtable.loc[in_nifti_subset] 
  out_subset = mrtable.loc[out_nifti_subset]
  # Add to dict of priors tables: columns --> thresholds, rows --> voxels)
  if in_subset.shape[0] != 0:
    priors_in[name] = calculate_regional_priors_in_ranges(in_subset,range_table)
  else: priors_in[name] = []  # This should never happen
  if out_subset.shape[0] != 0:
    priors_out[name] = calculate_regional_priors_in_ranges(out_subset,range_table)
  else: priors_out[name] = [] # This will happen for head node   


# Boom! Now we have, for each node, a table of prior probabilities, p[activation in range | region (voxel)] and we can calculate reverse inference!

# NOTE: saving all priors in one pickle is not feasible!

# SAVING DATA ------------------------------------------------------------------------------------
# Generate priors image for each of "in" images
priors_folder = "/home/vanessa/Documents/Work/BRAINMETA/mr/priors"
priors_keys = priors_in.keys()
pickle_folder = "/home/vanessa/Documents/Work/BRAINMETA/data"

for name in priors_keys:
  inny = priors_in[name]
  out = priors_out[name]
  print "Processing %s" %(name)
  in_map = "%s/priors_in_%s.nii.gz" %(priors_folder,name.replace(".nii.gz",""))
  out_map = "%s/priors_out_%s.nii.gz" %(priors_folder,name.replace(".nii.gz",""))
  pickle.dump(inny,open("%s/priors_in_%s.pkl" %(pickle_folder,name),"wb"))
  pickle.dump(out,open("%s/priors_out_%s.pkl" %(pickle_folder,name),"wb"))
  priors_to_brain_image(inny,in_map,brain_mask)
  priors_to_brain_image(out,out_map,brain_mask)
  

# Reverse Inference Calculation ------------------------------------------------------------------

# P(node mental process|activation) = P(activation|mental process) * P(mental process)

# divided by

# P(activation|mental process) * P(mental process) + P(A|~mental process) * P(~mental process)

# P(activation|mental process): my voxelwise prior map

# Given how I set this up, I'm going to start by calculating the probability of mental process given activation in a THRESHOLD. This means that I will have a reverse inference score for each level of activation! Or I can combine the priors maps differently to do something else.

# Now iterate through nodes (names), and calculate reverse inferences for each
reverse_inferences_thresh = dict()
for n in range(0,len(names)):
  name = names[n]
  print "Processing %s, %s of %s" %(name,n+1,len(names))
  # Filter table to those images
  in_node = locate_by_name(myjson, name)
  in_nifti_subset = get_nifti_names(in_node,[])
  out_nifti_subset = [i for i in nifti_files if i not in in_nifti_subset]  
  in_count = len(in_nifti_subset)
  out_count = len(out_nifti_subset)
  reverse_inferences_thresh[name] = calculate_reverse_inference_threshes(priors_in[name],priors_out[name],in_count,out_count)

# Save to pickle
pickle.dump(reverse_inferences_thresh,open("/home/vanessa/Documents/Work/BRAINMETA/data/brainmeta_reverse_inf_threshes.pkl","wb"))

# Load from pickle!
#reverse_inferences = pickle.load(open("/home/vanessa/Documents/Work/BRAINMETA/data/brainmeta_reverse_inf_threshes.pkl","rb"))

inference_df = pd.DataFrame(columns=range_table.index)
for name,inferences in reverse_inferences_thresh.iteritems():
  inference_df.loc[name] = inferences
inference_df.to_csv("/home/vanessa/Documents/Work/BRAINMETA/data/brainmeta_reverse_inf.tsv",sep="\t")

# Ok, so what we have now is a reverse inference score, the p(cognitive process | activation in range [x1..x2]) for all voxels in an image.  The way I'm thinking about this is the p(cognitive process | an activation "level") This is a good start, and the next step is bringing together these distributions - eg, if you think about an ACTUAL brain map result, we have many levels for each region (voxel) and so the same calculation should be done, but sampling from the appropriate prior for the level that we find in the image.  Let me try to write a new function.

# We could (and should) use a different brain map, but for now I'll just take an average over
# the maps in the node set. Let's just do this for the main node sets (eg, eliminate any that are just a single image!)

import re
expression = re.compile(".nii.gz")
names_subset = [x for x in names if not expression.search(x)]
names_subset.pop(0)

reverse_inferences = dict()
for n in range(17,len(names_subset)):
  name = names_subset[n]
  print "Processing %s, %s of %s" %(name,n+1,len(names_subset))
  # Filter table to those images
  in_node = locate_by_name(myjson, name)
  in_nifti_subset = get_nifti_names(in_node,[])
  out_nifti_subset = [i for i in nifti_files if i not in in_nifti_subset]  
  in_count = len(in_nifti_subset)
  out_count = len(out_nifti_subset)
  average_stat_map = mrtable.loc[in_nifti_subset].mean()
  reverse_inferences[name] = calculate_reverse_inference(average_stat_map,range_table,priors_in[name],priors_out[name],in_count,out_count)
