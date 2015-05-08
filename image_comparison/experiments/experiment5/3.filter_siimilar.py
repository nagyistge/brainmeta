#!/usr/bin/python

subject_path_pkl = "/scratch/users/vsochat/DATA/BRAINMETA/experiment5/doc/HCP_465_motor_cope_input.pkl"
output_path_pkl = "/scratch/users/vsochat/DATA/BRAINMETA/experiment5/doc/HCP_465_motor_cope_input_similar.pkl"

# This is a matrix of all cope image paths, cols are contrasts, rows subjects
maps = pandas.read_pickle(subject_path_pkl)

# Here is the list of filtered subjects
filtered_list = list(pandas.read_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment5/doc/HCP_motor_similar.txt")[0])

# Do the filter
maps = maps.loc[filtered_list]

# Save to file
maps.to_pickle(output_path_pkl)

