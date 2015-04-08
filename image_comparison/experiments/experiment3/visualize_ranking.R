# This script will generate a web report to visualize which images we "got wrong" 
# as a sanity check to the ranking algorithms(s)

library(rjson)
library(plyr)

# First, read in all input files
datadir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3"
setwd("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/image_comparison/experiments/experiment3")
source("helper_functions.R")
setwd(datadir)

input_data_files = list.files(datadir,pattern="*n_pd_v2.tsv|*n_bm_v2.tsv|*n_pi_v2.tsv")
pearsons_pd = remove_columns(read.csv(input_data_files[2],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
pearsons_pi = remove_columns(read.csv(input_data_files[3],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
pearsons_bm = remove_columns(read.csv(input_data_files[1],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
spearman_pd = remove_columns(read.csv(input_data_files[5],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
spearman_pi = remove_columns(read.csv(input_data_files[6],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
spearman_bm = remove_columns(read.csv(input_data_files[4],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
thresholds = sort(unique(pearsons_pd$thresh))
results = list(pds=spearman_pd,pis=spearman_pi,bms=spearman_bm,pdp=pearsons_pd,pip=pearsons_pi,bmp=pearsons_bm)

# Setup paths for outdirectory
outdir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/img/experiment3/thresholded"
setwd(outdir)

# Index page will have links to:
# Complete Case Analysis (link to first subject)
# Pairwise Inclusion (link to first subject)
# Brain Mask (link to first subject)

# Sub-results folders:
outdir_pd = paste(outdir,"/pd",sep="")
outdir_pi = paste(outdir,"/pi",sep="")
dir.create(outdir_pd)
dir.create(outdir_pi)

# Single results page, one per subject in each folder
# QUERY IMAGE AT TOP (unthresholded)

## HEATMAP MATRICES: green for correct, red for incorrect, mouseover shows other brain map

# Heatmap matrix positive and negative:
# thresh0.0
# .
# .
# .
# thresh4.0

# Heatmap matrix positive:
# thresh0.0
# .
# .
# .
# thresh4.0

# Arrows at bottom go to next map

# CHUNK INDIFFERENT RANKING ALGORITHM

# Single subject page templates
template = "templates/template.html"

# Here are the "last links" that we start with
prev = "../index.html"
  
# For each image, for each threshold, we assess distance from the "gold standard" ordering based on task/contrast
image_ids = unique(results[["pdp"]]$UID)
labels = c("intersect.pearson","union.pearson","intersect.spearman","union.spearman")
thresh_labels = as.character(thresholds)
thresh_labels[1] = "0.0"
thresh_labels[3] = "1.0"
thresh_labels[12] = "2.0"
thresh_labels[16] = "4.0"
names(thresh_labels) = thresholds
strategies = c("pairwise deletion","pairwise inclusion","pairwise deletion","pairwise inclusion")

for (i in 1:length(image_ids)){

  cat("Processing",i,"of",length(image_ids),"\n")
    
  # Calculate gold standard, and actual rankings
  image_id = image_ids[i]
  other_ids = image_ids[-which(image_ids==image_id)]
  gs = make_gold_standard_ranking(image_id,other_ids)
  # We only care about CONTRAST for now
  gs = gs[,c(1,2)]
  # For each of pearson, spearman [union, intersect, brainmask] get scores for image_id
  dfposneg = get_single_result(image_id,results,direction="posneg")   
  dfpos = get_single_result(image_id,results,direction="pos")   
  for (thresh in thresholds){
    for (l in 1:length(labels)){

      label = labels[l]
      strategy = strategies[l]
      metric = strsplit(label,"[.]")[[1]][2]
      thresh_label = thresh_labels[as.character(thresh)]
      
      # Subject specific output pages
      if (strategy=="pairwise deletion"){
        ss_file = paste(outdir_pd,"/pd",i,"_",thresh_label,"_",metric,".html",sep="")
        ss = paste("pd",i,"_",thresh_label,"_",metric,".html",sep="")        

        # Next pages
        if (i==length(image_ids)){
          pd_next = paste("pd1_",thresh_label,"_",metric,".html",sep="")
        } else {
          pd_next = paste("pd",i+1,"_",thresh_label,"_",metric,".html",sep="")
        }  
        
        } else {
          ss_file = paste(outdir_pi,"/pi",i,"_",thresh_label,"_",metric,".html",sep="")        
          ss = paste("pi",i,"_",thresh_label,"_",metric,".html",sep="")      
      
        # Next pages
        if (i==length(image_ids)){
          pd_next = paste("pi1_",thresh_label,"_",metric,".html",sep="")
        } else {
          pd_next = paste("pi",i+1,"_",thresh_label,"_",metric,".html",sep="")          
        }  
      
      }  
    
      # Get ordering based on actual scores
      sorted_posneg = filter_single_result(dfposneg,thresh,label,other_ids,image_id)
      sorted_pos = filter_single_result(dfpos,thresh,label,other_ids,image_id)
      ranking_posneg = get_ranking(gs,sorted_posneg)
      ranking_pos = get_ranking(gs,sorted_pos)
      
      # Add predicted order based on similarity scores
      sorted_pos = sorted_pos[rownames(ranking_pos$CONTRAST)]
      sorted_pos = seq(1,length(sorted_pos))
      sorted_posneg = sorted_posneg[rownames(ranking_posneg$CONTRAST)]
      sorted_posneg = seq(1,length(sorted_posneg))
      
      # For each, format image names
      img_pos = paste(rownames(ranking_pos$CONTRAST),"_thresh",as.integer(thresh),".0_pos.png",sep="")
      img_posneg = paste(rownames(ranking_posneg$CONTRAST),"_thresh",as.integer(thresh),".0_posneg.png",sep="")
      ranking_posneg$CONTRAST$img = img_posneg
      ranking_pos$CONTRAST$img = img_pos
      ranking_posneg$CONTRAST$ranking = sorted_posneg
      ranking_pos$CONTRAST$ranking = sorted_pos
      
      # Write each to json
      json_posneg = format_json(ranking_posneg$CONTRAST)
      json_pos = format_json(ranking_pos$CONTRAST)
      
      # FILL IN TEMPLATE
      tmp = readLines(template)
      tmp = gsub("repNEXTrep",pd_next,tmp,perl=TRUE)
      tmp = gsub("repPREVIOUSrep",prev,tmp,perl=TRUE)
      tmp = gsub("repSTRATEGYrep",paste(strategy,", ",metric,": ",thresh_label,sep=""),tmp,perl=TRUE)
      tmp = gsub("repQUERY_IMAGErep",paste(image_id,"/",image_id,"_thresh0.0_posneg.png",sep=""),tmp,perl=TRUE)
      if (metric == "spearman"){
        tmp = gsub("repSTRATEGY_COLORrep","#33CC99",tmp,perl=TRUE)
      } else {  
        tmp = gsub("repSTRATEGY_COLORrep","#FF3366",tmp,perl=TRUE)
      }
      tmp = gsub("repPOSrep",json_pos,tmp,perl=TRUE)
      tmp = gsub("repPOSNEGrep",json_posneg,tmp,perl=TRUE)
      
      filey = file(ss_file)
      writeLines(tmp, filey)
      close(filey)
      
      # Set new previous page
      if (i==length(image_ids)){
        prev = "../index.html"
      }
    }
  }
}
