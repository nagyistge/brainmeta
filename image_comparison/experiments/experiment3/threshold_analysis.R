library(plyr)

# First, read in all input files
datadir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3"
input_files = list.files(datadir,pattern="*.tsv")
input_data = list.files(datadir,pattern="*spearman*|*pearson*") # data files have extra column, parse differently
other_files = input_files[-which(input_files %in% input_data)]
setwd("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/image_comparison/experiments/experiment3")
source("helper_functions.R")
setwd(datadir)

# We first will filter data to remove redundant contrasts (negative and inverse subtractions)
for (input_file in other_files){
  cat("Filtering",input_file,"\n")
  output_file = paste(datadir,"/",gsub(".tsv","_filter.tsv",input_file),sep="")
  filter_data(input_file,output_file,extra_cols=3)  
}

# Data files have an extra column, df_mr
for (input_file in input_data){
  cat("Filtering",input_file,"\n")
  output_file = paste(datadir,"/",gsub(".tsv","_filter.tsv",input_file),sep="")
  filter_data(input_file,output_file,extra_cols=4)  
}

# Read in all filtered data
input_size_files = list.files(datadir,pattern="*sizes.tsv")
input_data_files = list.files(datadir,pattern="*n_cca.tsv|*n_svi.tsv")
input_nanlog_files = list.files(datadir,pattern="*g_cca.tsv|*g_svi.tsv")
pearsons_pd = remove_columns(read.csv(input_data_files[1],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
pearsons_pi = remove_columns(read.csv(input_data_files[2],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
spearman_pd = remove_columns(read.csv(input_data_files[3],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))
spearman_pi = remove_columns(read.csv(input_data_files[4],stringsAsFactors=FALSE,sep="\t"),c("dof_mr1","mr_df"))

thresholds = sort(unique(pearsons_pd$thresh))

# Part I: Handling Missing Values to Optimize Similarity Search Ranking
# ? Visualize: How do distributions of scores change? ###################################################################

# NOTE: R can get buggy and not correctly return plots from a function, so this was run manually.
savedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/img/experiment3/distributions"
results = list(pds=spearman_pd,pis=spearman_pi,pdp=pearsons_pd,pip=pearsons_pi)
for (thresh in thresholds){
  outfile = paste(savedir,"/density_posneg",thresh,".png",sep="")
  #plot_result(results,thresh,"posneg",outfile)  
}

for (thresh in thresholds){
  outfile = paste(savedir,"/density_pos",thresh,".png",sep="")
  #plot_result(results,thresh,"pos",outfile)
}

# ? Visualize: How do mean scores change? ##############################################################################

# Now do for each of pearson, posneg
IP = flatten_data(pearsons_pd,direction="posneg",label="intersect.pearson")
UP = flatten_data(pearsons_pi,direction="posneg",label="union.pearson")
IS = flatten_data(spearman_pd,direction="posneg",label="intersect.spearman")
US = flatten_data(spearman_pi,direction="posneg",label="union.spearman")

plot.new()
# Put them all together and plot
ALL = as.data.frame(rbind(IP,UP,IS,US),stringsAsFactors=FALSE)
ALL$score = as.numeric(ALL$score)
ALL$thresh = as.numeric(ALL$thresh)
ALL$strategy = as.character(ALL$strategy)
ALL$strategy[ALL$strategy=="intersect.pearson"] = "cca pearson"
ALL$strategy[ALL$strategy=="intersect.spearman"] = "cca spearman"
ALL$strategy[ALL$strategy=="union.pearson"] = "svi pearson"
ALL$strategy[ALL$strategy=="union.spearman"] = "svi spearman"

save(ALL,file=paste(datadir,"/all_scores_posneg.Rda",sep=""))
#load(file=paste(datadir,"/all_scores_posneg.Rda",sep=""))
ALLSUM = ddply(ALL, c("strategy","thresh"), summarise, mscore = mean(score), up=get_ci(score,"upper"), down=get_ci(score,"lower"))
ggplot(ALLSUM, aes(x=thresh, group=strategy,y=mscore,ymin=down,ymax=up,fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  facet_wrap(~strategy) +
  ylab("Mean Score") + title("Scores with Different Strategies for Handling Missing Data")
ggsave(paste(savedir,"/meanscores_with_ci_posneg.png",sep=""))

ggplot(ALL[which(ALL$thresh<=8),],aes(x=score, fill=strategy)) +
  geom_density(alpha=0.25) + 
  ylab("Density") + 
  xlab("Score") +
  xlim(-1,1) +
  facet_wrap(~thresh)
ggsave(paste(savedir,"/scores_densities_posneg_smaller.png",sep=""))

# We want to visualize them separately - too crowded in the plot above
pdp_plot = ALL[ALL$strategy == "cca pearson",]
pds_plot = ALL[ALL$strategy == "cca spearman",]
pip_plot = ALL[ALL$strategy == "svi pearson",]
pis_plot = ALL[ALL$strategy == "svi spearman",]

# We will need ggplot colors
pdp_gg = plot_distribution(pdp_plot,"complete case analysis, pearson",ymax=2.8)
pip_gg = plot_distribution(pip_plot,"single-value imputation, pearson",ymax=2.8)
pds_gg = plot_distribution(pds_plot,"complete case analysis, spearman",ymax=2.8)
pis_gg = plot_distribution(pis_plot,"single-value imputation, spearman",ymax=2.8)
g = arrangeGrob(pdp_gg, pds_gg, pip_gg, pis_gg, ncol=1)
ggsave(file=paste(datadir,"scores_distributions_posneg.png",sep=""),g)    


# Now do for each of pearson, pos
IP = flatten_data(pearsons_pd,direction="pos",label="intersect.pearson")
UP = flatten_data(pearsons_pi,direction="pos",label="union.pearson")
IS = flatten_data(spearman_pd,direction="pos",label="intersect.spearman")
US = flatten_data(spearman_pi,direction="pos",label="union.spearman")

plot.new()
# Put them all together and plot
ALL = as.data.frame(rbind(IP,UP,IS,US),stringsAsFactors=FALSE)
ALL$score = as.numeric(ALL$score)
ALL$thresh = as.numeric(ALL$thresh)
ALL$strategy = as.character(ALL$strategy)
ALL$strategy[ALL$strategy=="intersect.pearson"] = "cca pearson"
ALL$strategy[ALL$strategy=="intersect.spearman"] = "cca spearman"
ALL$strategy[ALL$strategy=="union.pearson"] = "svi pearson"
ALL$strategy[ALL$strategy=="union.spearman"] = "svi spearman"

save(ALL,file=paste(datadir,"/all_scores_pos.Rda",sep=""))
#load(file=paste(datadir,"/all_scores_pos.Rda",sep=""))
ALLSUM = ddply(ALL, c("strategy","thresh"), summarise, mscore = mean(score), up=get_ci(score,"upper"), down=get_ci(score,"lower"))
ggplot(ALLSUM, aes(x=thresh, group=strategy,y=mscore,ymin=down,ymax=up,fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +") +
#  facet_wrap(~strategy) +
  ylab("Mean Score") + title("Scores with Different Strategies for Handling Missing Data")
ggsave(paste(savedir,"/meanscores_with_ci_pos_singleplot.png",sep=""))

ggplot(ALL,aes(x=score, fill=strategy)) + 
  geom_density(alpha=0.25) + 
  ylab("Density") + 
  xlab("Score") +
  ylim(0,2.8) +
  xlim(-1,1) +
  facet_wrap(~thresh)
ggsave(paste(savedir,"/scores_densities_pos.png",sep=""))


# ? Visualize: How do sizes change? ##############################################################################

pd_sizes = read.csv(input_size_files[1],sep="\t")
pi_sizes = read.csv(input_size_files[2],sep="\t")

# Now we will look at the distributions of values themselves!
ncomps = c()
direction = "pos"
thresholds = sort(unique(pd_sizes$thresh))
for (thresh in thresholds){
    cat("  THRESH",thresh,"\n")
    tmp = c()
    # Get mask size distributions
    pd = unlist(get_masksize_dist(pd_sizes,thresh,direction))
    pi = unlist(get_masksize_dist(pi_sizes,thresh,direction))
    tmp[["cca"]] = pd
    tmp[["svi"]] = pi
    ncomps[[as.character(thresh)]] = tmp
}
save(ncomps,file=paste(datadir,"/dist_masksizes",direction,".Rda",sep=""))

# Plot the number of comparisons - first put into a big data frame
df = c()
ts = names(ncomps)
for (t in ts){
    tsubset = ncomps[[t]]
    tmp1 = cbind(tsubset$cca,rep(t,length(tsubset$cca)),rep(direction,length(tsubset$cca)),rep("complete case analysis",length(tsubset$cca)))
    tmp2 = cbind(tsubset$svi,rep(t,length(tsubset$svi)),rep(direction,length(tsubset$svi)),rep("single value imputation",length(tsubset$svi)))
    tmp = rbind(tmp1,tmp2)
    df = rbind(df,tmp)
}  
colnames(df) = c("mask.size","thresh","direction","strategy")
df = as.data.frame(df)
df$mask.size = as.numeric(as.character(df$mask.size))
df$thresh = as.numeric(as.character(df$thresh))
df$strategy = as.character(df$strategy)
save(df,file=paste(datadir,"/dist_masksizes_flat_",direction,"_df.Rda",sep=""))
df = df[,-3] # remove posneg
tmp = ddply(df, c("thresh","strategy"), summarise, mask.size.mean=mean(mask.size),mask.size.min = min(mask.size), mask.size.max=max(mask.size))

# Here we are looking at the mean pearsons (+/- one standard deviation)
ggplot(tmp, aes(x=thresh, y=mask.size.mean, ymin=mask.size.min, ymax=mask.size.max, fill=strategy,group=strategy,color=strategy)) + 
  geom_line(size=1) + 
  title("Mask Sizes at Different Thresholds +") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +/-") +
  ylab("Mask Size (voxels)") +
  scale_x_continuous(breaks = round(seq(0, 13, by = 1.0),1))
ggsave(paste(savedir,"/mask_sizes_voxels_",direction,".png",sep=""))

# Now plot percentage of voxels
total_size = max(df$mask.size)
df2 = df
df2$mask.size = df2$mask.size / total_size
tmp = ddply(df2, c("thresh","strategy"), summarise, mask.size.mean=mean(mask.size),mask.size.min = min(mask.size), mask.size.max=max(mask.size))
ggplot(tmp, aes(x=thresh, y=mask.size.mean, ymin=mask.size.min, ymax=mask.size.max, fill=strategy,group=strategy,color=strategy)) + 
  geom_line(size=1) + 
  title("Mask Sizes at Different Thresholds +") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +/-") +
  ylab("size of comparison set as % of brain mask") + 
  scale_x_continuous(breaks = round(seq(0, 13, by = 1.0),1))
ggsave(paste(savedir,"/mask_percsizes_atthresh_",direction,".png",sep=""))

# Visualize - the number of nans (and why) ######################################################################################################

pd_nan = read.csv(input_nanlog_files[1],sep="\t")
pi_nan = read.csv(input_nanlog_files[2],sep="\t")

thresholds = sort(unique(pd_nan$thresh))
direction = "posneg"
# For each threshold, and strategy, count the number of successful vs nans
nanlog = data.frame(nan_mrthresh_empty=1,success=1,nan_fewer_3_values=1,nan_no_overlap=1,thresh=1,strategy=1,direction=1,stringsAsFactors=FALSE)
count=1
for (thresh in thresholds) {
  pdsub = pd_nan[pd_nan$thresh==thresh,-which(colnames(pd_nan) %in% c("X","thresh","dof_mr1"))]    
  pdsub = pdsub[pdsub$direction==direction,-which(colnames(pdsub) == "direction")]
  pdsub = table(unlist(pdsub))
  pisub = pi_nan[pi_nan$thresh==thresh,-which(colnames(pi_nan) %in% c("X","thresh","dof_mr1"))]    
  pisub = pisub[pisub$direction==direction,-which(colnames(pisub) == "direction")]
  pisub = table(unlist(pisub))
  nanlog[count,c(which(colnames(nanlog)%in% names(pdsub)),5,6,7)] = c(pdsub,thresh,"cca",direction)
  nanlog[count+1,c(which(colnames(nanlog)%in% names(pisub)),5,6,7)] = c(pisub,thresh,"svi",direction)
  count=count+2
}
direction = "pos"
# For each threshold, and strategy, count the number of successful vs nans
for (thresh in thresholds) {
  pdsub = pd_nan[pd_nan$thresh==thresh,-which(colnames(pd_nan) %in% c("X","thresh","dof_mr1"))]    
  pdsub = pdsub[pdsub$direction==direction,-which(colnames(pdsub) == "direction")]
  pdsub = table(unlist(pdsub))
  pisub = pi_nan[pi_nan$thresh==thresh,-which(colnames(pi_nan) %in% c("X","thresh","dof_mr1"))]    
  pisub = pisub[pisub$direction==direction,-which(colnames(pisub) == "direction")]
  pisub = table(unlist(pisub))
  nanlog[count,c(which(colnames(nanlog)%in% names(pdsub)),5,6,7)] = c(pdsub,thresh,"cca",direction)
  nanlog[count+1,c(which(colnames(nanlog)%in% names(pisub)),5,6,7)] = c(pisub,thresh,"svi",direction)
  count=count+2
}
nanlog[is.na(nanlog)] = 0
save(nanlog,file=paste(datadir,"/nanlog_cca-svi.Rda",sep=""))
write.table(nanlog,file=paste(datadir,"/nanlog_cca-svi.tsv",sep=""),sep="\t",row.names=FALSE)


# ? Quantify: How often are rankings significantly different? ##############################################################################
# We will use the Chunk Indifferent Ranking Algorithm, which compares a gold standard ranking against the actual
# ranking, and assigns accuracy for each group ("chunk") within the ranking (a chunk corresponds to a group of images of the same task, contrast, etc, 
# A correct prediction is given full credit (an accuracy of 1/N, where N is the number of members of the chunk), and an incorrect
# prediction is given a weighted accuracy, 1/N * the distance it falls between the best and worst case position. 
# [A] For group == the first chunk: 
#  The best case position is 1 position after the last member of the gold standard chunk, and the worst is the end of the list.
# [B] For group == the last chunk: 
#  The best case position is the index of the first member of the chunk, and the worst is the first index of the list (1).
# [C] For group == middle chunk: 
#  Members with actual positions up the list use algorithm A, above, and members with actual positions down in the list use algorithm B


# Specify direction
direction = "posneg"

# For each image, for each threshold, we assess distance from the "gold standard" ordering based on task/contrast
image_ids = unique(results[["pdp"]]$UID)
labels = c("intersect.pearson","union.pearson","intersect.spearman","union.spearman")
thresholds = thresholds[1:11]
acc_df = c()
for (i in 1:length(image_ids)){
  cat("Processing",i,"of",length(image_ids),"\n","Threshold:","\n")
  image_id = image_ids[i]
  other_ids = image_ids[-which(image_ids==image_id)]
  gs = make_gold_standard_ranking(image_id,other_ids)
  # For each of pearson, spearman [union, intersect, brainmask] get scores for image_id
  df = get_single_result(image_id,results,direction=direction)   
  for (thresh in thresholds){
    cat(thresh,",",sep="")
    for (label in labels){
      # Get ordering based on actual scores
      sorted = filter_single_result(df,thresh,label,other_ids,image_id)
      acc = calculate_accuracy(gs,sorted,9)
      
      # For higher thresholds, for smaller groups we may not have an entire group, 
      # so the accuracy will not be calculatable, and we do not append these
      # this is only true for threshold of 4 for "pos" data
      acc_df_single = c()
      if (!is.null(acc$CONTRAST[[1]]) && !is.null(acc$CONTRAST[[2]])){
        contrast = rbind(c(image_id,thresh,label,acc$CONTRAST[[1]],1,"contrast"),
                         c(image_id,thresh,label,acc$CONTRAST[[2]],2,"contrast"))
        acc_df_single = rbind(acc_df_single,contrast)
      }
      if (!is.null(acc$TASK[[1]]) && !is.null(acc$TASK[[2]])){
        task = rbind(c(image_id,thresh,label,acc$TASK[[1]],1,"task"),
                     c(image_id,thresh,label,acc$TASK[[2]],2,"task"))
        acc_df_single = rbind(acc_df_single,task)
      }

      if (!is.null(acc$GROUP[[1]]) && !is.null(acc$GROUP[[2]])){
        group = rbind(c(image_id,thresh,label,acc$GROUP[[1]],1,"group"),
                    c(image_id,thresh,label,acc$GROUP[[2]],2,"group"))
        acc_df_single = rbind(acc_df_single,group)
      }

      if (!is.null(acc$TASK_CONTRAST[[1]]) && !is.null(acc$TASK_CONTRAST[[2]]) &&
          !is.null(acc$TASK_CONTRAST[[3]]) && !is.null(acc$TASK_CONTRAST[[4]])){
          task_contrast = rbind(c(image_id,thresh,label,acc$TASK_CONTRAST[[1]],1,"task_contrast"),
                            c(image_id,thresh,label,acc$TASK_CONTRAST[[2]],2,"task_contrast"),
                            c(image_id,thresh,label,acc$TASK_CONTRAST[[3]],3,"task_contrast"))      
          acc_df_single = rbind(acc_df_single,task_contrast)
      }
      acc_df = rbind(acc_df,acc_df_single)
    }
  }
}    
colnames(acc_df) = c("image_id","thresh","strategy","accuracy","group","standard")
save(acc_df,file=paste(datadir,"/accuracy_df_",direction,"9.Rda",sep=""))
load(paste(datadir,"/accuracy_df_",direction,".Rda",sep=""))
write.table(acc_df,file=paste(datadir,"/accuracy_df_",direction,"9.tsv",sep=""),sep="\t")
acc_df = as.data.frame(acc_df,stringsAsFactors=FALSE)

# ? Visualize: How often are rankings significantly different? ##############################################################################
# We don't care about group 2 for task, contrast or group
subset = acc_df
contrast = subset[subset$standard=="contrast",]
contrast = contrast[which(contrast$group==1),]
task = subset[subset$standard=="task",]
task = task[which(task$group==1),]
task_contrast = subset[subset$standard=="task_contrast",]
# For task_contrast, calculate weighted accuracy as 2/3,1/3
task_contrast_grp1 = task_contrast[which(task_contrast$group==1),]
task_contrast_grp2 = task_contrast[which(task_contrast$group==2),]
tc_accuracy = (2/3)*as.numeric(task_contrast_grp1$accuracy) + (1/3)*as.numeric(task_contrast_grp2$accuracy)
task_contrast = task_contrast_grp1
task_contrast$accuracy = tc_accuracy
group = subset[subset$standard=="group",]
group = group[which(group$group==1),]
# Not including group - accuracy is essentially 0
subset = rbind(contrast,task,task_contrast)
subset$accuracy = as.numeric(subset$accuracy)
subset$thresh = as.numeric(subset$thresh)

subset$strategy[subset$strategy=="intersect.pearson"] = "cca pearson"
subset$strategy[subset$strategy=="intersect.spearman"] = "cca spearman"
subset$strategy[subset$strategy=="union.pearson"] = "svi pearson"
subset$strategy[subset$strategy=="union.spearman"] = "svi spearman"

subset_summary = ddply(subset, c("strategy","thresh","standard"), summarise, mscore = mean(accuracy), up=get_ci(accuracy,"upper"), down=get_ci(accuracy,"lower"))
ggplot(subset_summary, aes(x=thresh,y=mscore,ymax=up,ymin=down, fill=strategy,colour=strategy,standard=standard)) + 
  geom_line(size=1) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +") +
  facet_wrap(~standard) +
  ylab("Accuracy")
ggsave(paste(savedir,"/chunk_ranking_accuracy",direction,"_withCI_9.png",sep=""))
write.table(subset_summary,file=paste(datadir,"/chunkrank_summary_",direction,".tsv",sep=""),sep="\t",row.names=FALSE)

# Let's just look at complete case analysis
load(paste(datadir,"/accuracy_df_posneg.Rda",sep=""))
acc_df = as.data.frame(acc_df,stringsAsFactors=FALSE)
subset = as.data.frame(acc_df)
contrast = subset[subset$standard=="contrast",]
contrast = contrast[which(contrast$group==1),]
task = subset[subset$standard=="task",]
task = task[which(task$group==1),]
task_contrast = subset[subset$standard=="task_contrast",]
# For task_contrast, calculate weighted accuracy as 2/3,1/3
task_contrast_grp1 = task_contrast[which(task_contrast$group==1),]
task_contrast_grp2 = task_contrast[which(task_contrast$group==2),]
tc_accuracy = (2/3)*as.numeric(task_contrast_grp1$accuracy) + (1/3)*as.numeric(task_contrast_grp2$accuracy)
task_contrast = task_contrast_grp1
task_contrast$accuracy = tc_accuracy
subset = rbind(contrast,task,task_contrast)
subset$accuracy = as.numeric(subset$accuracy)
subset$thresh = as.numeric(subset$thresh)

load(paste(datadir,"/accuracy_df_pos.Rda",sep=""))
subset1 = as.data.frame(acc_df,stringsAsFactors=FALSE)
contrast = subset1[subset1$standard=="contrast",]
contrast = contrast[which(contrast$group==1),]
task = subset1[subset1$standard=="task",]
task = task[which(task$group==1),]
task_contrast = subset1[subset1$standard=="task_contrast",]
# For task_contrast, calculate weighted accuracy as 2/3,1/3
task_contrast_grp1 = task_contrast[which(task_contrast$group==1),]
task_contrast_grp2 = task_contrast[which(task_contrast$group==2),]
tc_accuracy = (2/3)*as.numeric(task_contrast_grp1$accuracy) + (1/3)*as.numeric(task_contrast_grp2$accuracy)
task_contrast = task_contrast_grp1
task_contrast$accuracy = tc_accuracy
subset1 = rbind(contrast,task,task_contrast)
subset1$accuracy = as.numeric(subset1$accuracy)
subset1$thresh = as.numeric(subset1$thresh)

posneg = subset[subset$strategy=="intersect.pearson",]
posneg$strategy = "cca pearson +/-"
pos = subset1[subset1$strategy=="intersect.pearson",]
pos$strategy = "cca pearson, +"
data = rbind(posneg,pos)
data$thresh = as.numeric(data$thresh)
data$accuracy = as.numeric(data$accuracy)

subset_summary = ddply(data, c("strategy","thresh","standard"), summarise, mscore = mean(accuracy), up=get_ci(accuracy,"upper"), down=get_ci(accuracy,"lower"))
ggplot(subset_summary, aes(x=thresh,y=mscore,ymax=up,ymin=down, fill=strategy,colour=strategy,standard=standard)) + 
  geom_line(size=1) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Comparing +/- and +") +
  facet_wrap(~standard,nrow=1) +
  ylab("Accuracy")
ggsave(paste(savedir,"/chunk_ranking_accuracy_direction_comparison_withCI.png",sep=""))


group$thresh = as.numeric(group$thresh)
group$accuracy = as.numeric(group$accuracy)
group_summary = ddply(group, c("strategy","thresh","standard"), summarise, mscore = mean(accuracy), up=get_ci(accuracy,"upper"), down=get_ci(accuracy,"lower"))
ggplot(group_summary, aes(x=thresh,y=mscore,ymax=up,ymin=down, fill=strategy,colour=strategy,standard=standard)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  facet_wrap(~standard) +
  ylab("Accuracy")
ggsave(paste(savedir,"/chunk_ranking_group_accuracy",direction,"_withCI.png",sep=""))


# Part II: Image Classification
unique_groups = c(paste("GRP0",seq(1,9),sep=""),"GRP10")

# For each pairwise group, A and B
#  Subset data to groups A and B, remove GRP# from labels
#   For each handling missing data strategy
#    For each threshold 0:4
#      Define gold standard labels as labels from A
#      For every map in B assign it a label (TASK_CON) corresponding to the most similar in A
#      Calculate (and save) accuracy

# We will keep track of which we have compared, so there are no repeats (eg, group1 vs group2 == group2 vs group1)
repeat_check = c()
all_acc = c()
direction = "pos"

# For each pairwise group, A and B
for (group1 in unique_groups){
  for (group2 in unique_groups){
    if (group1!=group2){
      #  Subset data to groups A and B, remove group from the labels
      # We will create a label with sorted names to eliminate repeats
      label = paste(sort(c(group1,group2)),collapse="vs")
      if (!(label %in% repeat_check)){
        acc_scores = get_group_result(group1,group2,results,direction=direction)
        acc_scores$label = label
        all_acc = rbind(all_acc,acc_scores)      
        repeat_check = c(repeat_check,label)
      }
    }
  }  
}

# Finally, we can calculate a mean accuracy for each threshold, strategy
all_acc$strategy[all_acc$strategy=="intersect.pearson"] = "cca pearson"
all_acc$strategy[all_acc$strategy=="intersect.spearman"] = "cca spearman"
all_acc$strategy[all_acc$strategy=="union.pearson"] = "svi pearson"
all_acc$strategy[all_acc$strategy=="union.spearman"] = "svi spearman"

acc_summary = ddply(all_acc, c("strategy","thresh"), summarise, 
                    accuracy_mean = mean(accuracy), accuracy_upper=get_ci(accuracy,"upper"), accuracy_down=get_ci(accuracy,"lower"),
                    sensitivity_mean = mean(sensitivity), sensitivity_upper = get_ci(sensitivity,"upper"),sensitivty_lower = get_ci(sensivity,"lower"),
                    specificity_mean = mean(specificity), specificity_upper = get_ci(specificity,"upper"),specificity_down = get_ci(specificity,"lower"),                   
                    roc_mean = mean(roc), roc_upper = get_ci(roc,"upper"),roc_lower = get_ci(roc,"lower"))

ggplot(acc_summary, aes(x=thresh,y=mscore,ymax=up,ymin=down, fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  ylab("Accuracy") +
  ylim(0,1) +
  scale_x_continuous(breaks = round(seq(0, 13, by = 1.0),1))
ggsave(paste(savedir,"/ml_accuracy",direction,"_withCI.png",sep=""))
write.table(acc_summary,file=paste(datadir,"/ml_accuracy_",direction,".tsv",sep=""),sep="\t")

acc_summary = read.table(file=paste(datadir,"/ml_accuracy_",direction,".tsv",sep=""),sep="\t")
# Woot!

# Finally, if we match NeuroVault percent voxels in brainmap to our result, what accuracy can we expect to get?
load(file=paste(datadir,"/dist_masksizes_flat_",direction,"_df.Rda",sep=""))
mask_sizes = df 

# Read in the threshold report for all neurovault images!
brainmask_size = 228235
subset = df[df$strategy=="complete case analysis",]
subset$value = as.numeric(as.character(subset$mask.size/brainmask_size))
subset = subset[,-c(1,3)]
# Get the mean for each size
ss = unique(subset$thresh)
df=c()
for (s in ss){
  sub = subset[which(subset$thresh==s),]
  df=rbind(df,cbind(s,mean(sub$value)))
}
colnames(df) = c("thresh","perc_nonzero_voxels")

# Here is where we can load accuracy values for posneg, or pos
load(paste(datadir,"/accuracy_df_posneg.Rda",sep="")) #(acc_df)
acc_df = as.data.frame(acc_df,stringsAsFactors=FALSE)
acc_df = acc_df[which(acc_df$standard=="contrast"),]

# Here we will save lists of accuracies based on the real thresholding, we only care about contrast
pd_q=c()
pi_q=c()

thresh_report = read.csv("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/nv_all_774_threshold_report.tsv",sep="\t")
for (it in 1:1000){
  cat("iteration",it,"\n")
  for (t in 1:nrow(thresh_report)){
    percent_nonzero = thresh_report$percent_nonzero[t]
    # Find the most similar sized value (smallest squared distance)
    tmp = abs(df[,2] - percent_nonzero)
    idx = which(tmp==min(tmp))[1]
    thresh = df[idx,1]
    # Now we sample an accuracy at that threshold
    acc_pd = acc_df[which(acc_df$strategy=="intersect.pearson"),]
    if (thresh %in% acc_pd$thresh) {
      pd_q=c(pd_q,sample(acc_pd$accuracy[which(acc_pd$thresh==thresh)],1))  
    }    
    acc_pi = acc_df[which(acc_df$strategy=="intersect.spearman"),]
    if (thresh %in% acc_pi$thresh) {
      pi_q=c(pi_q,sample(acc_pi$accuracy[which(acc_pi$thresh==thresh)],1))  
    }
  }
}

# Calculate mean accuracy
mean(as.numeric(pi_q))
mean(as.numeric(pd_q))