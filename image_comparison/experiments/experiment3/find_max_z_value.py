#!/usr/bin/env python

# This script will read in the complete set of maps, convert each to a z score map, and determine the max and min zscores across all maps. We will also save a directory of images
# for each map, including unthresholded and thresholded

import pickle
import pandas
import image_transformations as IT
from nilearn.plotting import plot_stat_map
from pybraincompare.template.futils import make_dir
import nibabel as nib

basedir = "/scratch/users/vsochat/DATA/BRAINMETA/experiment3"
input_file = "%s/doc/hcp_groupmaps_filter.tsv" %(basedir)
imgdir = "%s/img"
inputs = pandas.read_csv(input_file,sep="\t")
filenames = inputs.files.tolist()

# Groups file tells us the dof for each group
groups_pkl = "%s/doc/hcp_10groups460_alltasks.pkl" %(basedir)
groups = pickle.load(open(groups_pkl,"rb"))

minval = 0
maxval = 0

# Calculate max and min z score values
for f in range(0,len(filenames)):
    print "Processing %s of %s" %(f,len(filenames))
    filename = filenames[f]
    # Load the image, convert to Z, and calculate the thresholded version
    mr = nib.load(filename)
    # Convert to Z score map
    group1 = inputs.groups[inputs.files==filename].tolist()[0]
    dof = len(groups[group1]) - 2
    z = IT.t_to_z(mr,dof)
    if z.get_data().min() < minval:
        minval = z.get_data().min()
    if z.get_data().max() > maxval:
        maxval = z.get_data().max()

# Actual min
# >>> minval
# -12.269853691998309
# 11.183601479336019

thresholds = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0,12.0,13.0]

# Create images of maps
# Now we want to do the same thing, but threshold the images, and save an image at each threshold 
for f in inputs.iterrows():
    print "Processing %s of %s" %(f,len(filenames))
    uid = f[1].uid
    mr = nib.load(f[1].files)
    subject_directory = "%s/img/%s" %(basedir,uid)
    make_dir(subject_directory)
    mrthresh_posneg = IT.threshold_abs(mr,thresholds=thresholds)
    mrthresh_pos = IT.threshold_pos(mr,thresholds=thresholds)
    for thresh in thresholds:
        plotposneg = plot_stat_map(mrthresh_posneg[thresh],
                                title="%s_thresh%s_posneg" %(uid,thresh),
                                display_mode="z",
                                cut_coords=5,vmax=maxval)
        plotposneg.savefig("%s/%s_thresh%s_posneg.png" %(subject_directory,uid,thresh))
        plotpos = plot_stat_map(mrthresh_pos[thresh],
                                title="%s_thresh%s_pos" %(uid,thresh),
                                display_mode="z",
                                cut_coords=5,vmax=maxval)
        plotpos.savefig("%s/%s_thresh%s_pos.png" %(subject_directory,uid,thresh))
        plotpos.close()
        plotposneg.close()

