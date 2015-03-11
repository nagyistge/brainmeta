#/usr/bin/python Experiment 1:
# Individual images compared to transformations of themselves.

#1. Define "gold standard" as a comparison of image A with itself (using brain mask I assume?).
#2. Define a range of thresholds to apply to image (from unthresholded --> coordinate)
#3. Define a set of similarity metrics to use
#4. For each of pairwise comparison, pairwise inclusion, and brain mask, and for each similarity metric, test how scores change at the different levels of threshold.

from pyneurovault import api
import nibabel as nib
import numpy as np
from glob import glob
import similarity_metrics
import pandas

# IMAGE DOWNLOAD -----------

# Directory for analysis
outdir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON"

# Find openfmri images
nv = api.NeuroVault()
df = nv.get_images_with_collections_df()
result = nv.search(df=df,column_name="description_collection",search_string="OpenfMRI")

# Download and resample
collection_ids = result["collection_id"].unique().tolist()
outfolder = "%s/mr" %(outdir)
standard = "/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
download_table = nv.download_and_resample(outfolder,standard, collection_ids=collection_ids)
download_table.to_pickle("%s/openfmri_images.pkl" %(outfolder))

# Prepare meta data (labels) for images - we want study, image ID
meta = pandas.DataFrame()
meta["ID"] =  download_table.image_id
meta["COLLECTION_ID"] = download_table["collection_id"]
meta["MAP_TYPE"] = download_table["map_type"]
meta["DOI"] = download_table["DOI"]
meta["ORIGINAL_FILE"] = download_table["file"]
meta["DESCRIPTION"] = download_table["description_image"]
meta["FILE"] = ["/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/mr/resampled/000%s.nii.gz" %(x) for x in download_table.image_id]
meta.to_pickle("%s/openfmri_labels.pkl" %(outfolder))
meta.to_csv("%s/openfmri_labels.tsv" %(outfolder),sep="\t")
