#!/usr/bin/env python2

# Use neurovault, cognitive atlas APIs, and pybraincompare to calculate similarity (complete case analysis with pearson score) for images in NeuroVault tagged with a cognitive atlas contrast (N=93)

import nibabel
import pandas
from pyneurovault import api
from cognitiveatlas.api import get_task, get_concept

## STEP 1: DOWNLOAD OF NV IMAGES ######################################################

# Set up work folders for data
outfolder = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference"

# Get all collections
collections = api.get_collections()

# Filter images to those that have a DOI
collections = collections[collections.DOI.isnull()==False]
collections.to_csv("%s/collections_with_dois.tsv" %outfolder,encoding="utf-8",sep="\t")

# Get image meta data for collections (N=1023)
images = api.get_images(collection_pks=collections.collection_id.tolist())

# Filter images to those with contrasts defined (N=93)
images = images[images.cognitive_contrast_cogatlas_id.isnull()==False]
images.to_csv("%s/contrast_defined_images.tsv" %outfolder,encoding="utf-8",sep="\t")

standard = "/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
api.download_images(dest_dir = outfolder,images_df=images,target=standard)
standard = nibabel.load(standard)

## STEP 2: IMAGE SIMILARITY
######################################################

from pybraincompare.compare.mrutils import make_binary_deletion_mask
from pybraincompare.compare.maths import calculate_correlation
""" Usage:
calculate_correlation(images,mask=None,atlas=None,summary=False,corr_type="pearson"):
make_binary_deletion_mask(images)
"""

# Function to pad ID with appropriate number of zeros
def pad_zeros(the_id,total_length=6):
    return "%s%s" %((total_length - len(str(the_id))) * "0",the_id)

# Calculate image similarity with pearson correlation
image_ids = images.image_id.tolist()
simmatrix = pandas.DataFrame(columns=image_ids,index=image_ids)
for id1 in image_ids:
    print "Processing %s..." %id1
    mr1_id = pad_zeros(id1)
    mr1_path = "%s/resampled/%s.nii.gz" %(outfolder,mr1_id)
    mr1 = nibabel.load(mr1_path)
    for id2 in image_ids:
        mr2_id = pad_zeros(id2)
        mr2_path = "%s/resampled/%s.nii.gz" %(outfolder,mr2_id)
        mr2 = nibabel.load(mr2_path)
        # Make a pairwise deletion / complete case analysis mask
        pdmask = make_binary_deletion_mask([mr1,mr2])
        pdmask = nibabel.Nifti1Image(pdmask,affine=standard.get_affine())
        score = calculate_correlation([mr1,mr2],mask=pdmask)
        simmatrix.loc[id1,id2] = score
        simmatrix.loc[id2,id1] = score

simmatrix.to_csv("%s/contrast_defined_images_pearsonpd_similarity.tsv" %outfolder,sep="\t")
