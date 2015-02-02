#!/usr/bin/env python2

# Here is downloading openFmri data for the pilot / proof of concept

from pyneurovault import api

# Will extract all collections and images in one query to work from
nv = api.NeuroVault()

# Get unique cognitive atlas contrasts and counts
contrasts = nv.get_contrasts()

# Download images, collections, or both
nv.export_images_tsv("/home/vanessa/Documents/Work/BRAINMETA/doc/images.tsv")
nv.export_collections_tsv("/home/vanessa/Documents/Work/BRAINMETA/doc/collections.tsv")

# Download images from collections 102 and 106
outfolder = "/home/vanessa/Documents/Work/BRAINMETA/mr"
standard = "/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
download_table = nv.download_and_resample(outfolder,standard, collection_ids=[42,98,39])
download_table.to_csv("/home/vanessa/Documents/Work/BRAINMETA/doc/download_neurovault.tsv")
