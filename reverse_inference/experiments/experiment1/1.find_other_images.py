# Find other images with a DOI that we can guess task / contrast for

from pyneurovault import api
import pandas as pd
import numpy as np
import pylab as plt
import nibabel as nb
from nilearn.plotting.img_plotting import plot_anat
from pyneurovault.nsynth import get_neurosynth_terms
from pyneurovault.analysis import get_frequency_map

# Use a joblib memory, to avoid depending on an Internet connection
from joblib import Memory
mem = Memory(cachedir='/tmp/neurovault_analysis/cache')

# Will extract all collections and images in one query to work from
nv = api.NeuroVault()

# Get combined data frame
combined_df = nv.get_images_with_collections_df()
filtered_df = combined_df[-combined_df.DOI.isnull()]

# Remove HCP and openfmri
filtered_df = filtered_df[-filtered_df["description_collection"].str.contains("OpenfMRI")]
# Remove HCP and openfmri
filtered_df = filtered_df[-filtered_df["description_collection"].str.contains("HCP")]

descriptions = []
descriptions_images = []
authors = []
# In description, remove newlines
for row in filtered_df.iterrows():
  tmpd = row[1].description_collection.replace("\n","").replace("\r","")
  tmpa = row[1].authors.replace("\n","").replace("\r","")
  tmpi = row[1].description_image.replace("\n","").replace("\r","")
  if tmpd != None:
      descriptions.append(tmpd)
  else:
      descriptions.append("None")
  if tmpi != None:
      descriptions_images.append(tmpi)
  else:    
      descriptions_images.append("None")
  if tmpa != None:
      authors.append(tmpa)
  else:
      authors.append("None")
  
filtered_df.description_collection = descriptions
filtered_df.authors = authors
filtered_df.description_image = descriptions_images
filtered_df = filtered_df.drop("images",1)
filtered_df.to_csv("/home/vanessa/Desktop/neurovault_doi.tsv",sep="\t",encoding="utf-8")

