#/usr/bin/python

import pandas
import os

# We will be using all subjects from HCP that have data across tasks
ids = pandas.read_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_465_with_all_tasks.tsv",sep="\t")
ids = [str(id) for id in ids.id.tolist()]

contrasts = pandas.read_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment3/doc/hcp_contrasts_id.tsv",sep="\t")

# Filter down to just motor
contrasts = contrasts[contrasts.task=="tfMRI_MOTOR"]

# For each contrast, get subject data paths
group_list_directory = "/scratch/PI/russpold/work/HCP/group_maps/copes"

motor_ss_maps = pandas.DataFrame(columns=contrasts.contrasts,index=ids)

for con in contrasts.iterrows():
    copefile = "%s/%s_%s_copes.txt" %(group_list_directory,con[1].task,con[1].contrasts)
    paths = pandas.read_csv(copefile,header=None)    
    path_list = paths[0].tolist()
    # Get subject IDS
    subids = [p.split("/")[7] for p in path_list]
    paths.index = subids
    data_paths = paths.loc[ids]
    motor_ss_maps.loc[data_paths.index,con[1].contrasts] = data_paths[0]
    
motor_ss_maps.to_csv("/scratch/users/vsochat/DATA/BRAINMETA/experiment5/doc/HCP_465_motor_cope_input.tsv",sep="\t")
motor_ss_maps.to_pickle("/scratch/users/vsochat/DATA/BRAINMETA/experiment5/doc/HCP_465_motor_cope_input.pkl")
