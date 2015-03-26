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
input_size_files = list.files(datadir,pattern="*sizes_filter.tsv")
input_data_files = list.files(datadir,pattern="*n_pd_filter.tsv|*n_bm_filter.tsv|*n_pi_filter.tsv")
input_nanlog_files = list.files(datadir,pattern="*g_pd_filter.tsv|*g_bm_filter.tsv|*g_pi_filter.tsv")
pearsons_pd = read.table(input_data_files[2],stringsAsFactors=FALSE)
pearsons_pi = read.table(input_data_files[3],stringsAsFactors=FALSE)
pearsons_bm = read.table(input_data_files[1],stringsAsFactors=FALSE)
spearman_pd = read.table(input_data_files[5],stringsAsFactors=FALSE)
spearman_pi = read.table(input_data_files[6],stringsAsFactors=FALSE)
spearman_bm = read.table(input_data_files[4],stringsAsFactors=FALSE)

thresholds = unique(pearsons_pd$thresh)

# Part I: Handling Missing Values to Optimize Similarity Search Ranking
# ? Visualize: How do distributions of scores change? ###################################################################

# NOTE: R can get buggy and not correctly return plots from a function, so this was run manually.
savedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/img/experiment3/distributions"
results = list(pds=spearman_pd,pis=spearman_pi,bms=spearman_bm,pdp=pearsons_pd,pip=pearsons_pi,bmp=pearsons_bm)
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
MP = flatten_data(pearsons_bm,direction="posneg",label="brain.mask.pearson")
IS = flatten_data(spearman_pd,direction="posneg",label="intersect.spearman")
US = flatten_data(spearman_pi,direction="posneg",label="union.spearman")
MS = flatten_data(spearman_bm,direction="posneg",label="brain.mask.spearman")

plot.new()
# Put them all together and plot
ALL = as.data.frame(rbind(IP,UP,MP,IS,US,MS),stringsAsFactors=FALSE)
ALL$score = as.numeric(ALL$score)
#save(ALL,file=paste(datadir,"/all_scores_posneg.Rda",sep=""))
load(file=paste(datadir,"/all_scores_posneg.Rda",sep=""))
ALLSUM = ddply(ALL, c("strategy","thresh"), summarise, mscore = mean(score), up=get_ci(score,"upper"), down=get_ci(score,"lower"))
ggplot(ALLSUM, aes(x=thresh, group=strategy,y=mscore,ymin=down,ymax=up,fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  facet_wrap(~strategy) +
  ylab("Mean Score") + title("Scores with Different Strategies for Handling Missing Data")
ggsave(paste(savedir,"/meanscores_with_ci_posneg.png",sep=""))

ggplot(ALL,aes(x=score, fill=strategy)) +
  geom_density(alpha=0.25) + 
  ylab("Density") + 
  xlab("Score") +
  ylim(0,2.5) +
  xlim(-1,1) +
  facet_wrap(~thresh)
ggsave(paste(savedir,"/scores_densities_posneg.png",sep=""))

# Now do for each of pearson, pos
IP = flatten_data(pearsons_pd,direction="pos",label="intersect.pearson")
UP = flatten_data(pearsons_pi,direction="pos",label="union.pearson")
MP = flatten_data(pearsons_bm,direction="pos",label="brain.mask.pearson")
IS = flatten_data(spearman_pd,direction="pos",label="intersect.spearman")
US = flatten_data(spearman_pi,direction="pos",label="union.spearman")
MS = flatten_data(spearman_bm,direction="pos",label="brain.mask.spearman")

plot.new()
# Put them all together and plot
ALL = as.data.frame(rbind(IP,UP,MP,IS,US,MS),stringsAsFactors=FALSE)
ALL$score = as.numeric(ALL$score)
#save(ALL,file=paste(datadir,"/all_scores_pos.Rda",sep=""))
load(file=paste(datadir,"/all_scores_pos.Rda",sep=""))
ALLSUM = ddply(ALL, c("strategy","thresh"), summarise, mscore = mean(score), up=get_ci(score,"upper"), down=get_ci(score,"lower"))
ggplot(ALLSUM, aes(x=thresh, group=strategy,y=mscore,ymin=down,ymax=up,fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  facet_wrap(~strategy) +
  ylab("Mean Score") + title("Scores with Different Strategies for Handling Missing Data")
ggsave(paste(savedir,"/meanscores_with_ci_pos.png",sep=""))

ggplot(ALL,aes(x=score, fill=strategy)) + 
  geom_density(alpha=0.25) + 
  ylab("Density") + 
  xlab("Score") +
  ylim(0,2.5) +
  xlim(-1,1) +
  facet_wrap(~thresh)
ggsave(paste(savedir,"/scores_densities_pos.png",sep=""))



# We will read in size data later
pd_sizes = read.csv(input_size_files[2],sep="\t")
pi_sizes = read.csv(input_size_files[3],sep="\t")
bm_sizes = read.csv(input_size_files[1],sep="\t")


# ? Quantify: How often are distributions significantly different? ##############################################################################
# we want to know "how often" it looks different, or we get weird results.

# We will be doing wilcox tests to assess for differences in means for each strategy for handling missing data, for each threshold, vs. the gold standard
# We will save all results into a data frame, and correct for multiple comparisons within thresholds and strategies
# Then we can calculate the percentage of the time, for each strategy and threshold, that we have sig. different result
# We could then represent an actual database based on the content and estimate the % of time / threshold we get different results from gs
direction="pos"
image_ids = unique(results[["pdp"]]$UID)
plot_single = FALSE
wilcox_tests = c()

for (i in 1:length(image_ids)){
  cat("Processing",i,"of",length(image_ids),"\n")
  image_id = image_ids[i]
  df = get_single_result(image_id,results,direction=direction) 
  
  if (plot_single==TRUE){
    # Plot the distribution of values for each strategy and threshold
    ggplot(df,aes(y=value,x=strategy, fill=strategy)) +
      facet_wrap(~thresh) +
      geom_violin(alpha=0.25) +
      title("Strategies for Handling Missing Data Influence on Scores for Image",image_id) +
      ylab("Density") +
      xlab("Score")
    
    ggplot(df,aes(x=value, fill=strategy)) +
      facet_wrap(~thresh) +
      geom_density(alpha=0.25) +
      ylab("Density") +
      xlab("Pearsons R") 
  
    ggsave(paste(savedir,"/single_map_posneg.png",sep=""))
  }
  # We will do wilcox test to determine when we have significantly different means
  gs = df[df$thresh==0.0,]
  # We will compare threshold of 0.0 to itself as sanity check
  for (thresh in thresholds){
    wt_pdp = wilcox_test(df,thresh,gs,"intersect.pearson")
    wt_pip = wilcox_test(df,thresh,gs,"union.pearson")
    wt_bmp = wilcox_test(df,thresh,gs,"brain.mask.pearson")
    wt_pds = wilcox_test(df,thresh,gs,"intersect.spearman")
    wt_pis = wilcox_test(df,thresh,gs,"union.spearman")
    wt_bms = wilcox_test(df,thresh,gs,"brain.mask.spearman")    
    tmp = cbind(wt_pdp$p.value,wt_pip$p.value,wt_bmp$p.value,
                wt_pds$p.value,wt_pis$p.value,wt_bms$p.value,
                wt_pdp$conf.int[1],wt_pdp$conf.int[2],wt_pip$conf.int[1],
                wt_pip$conf.int[2],wt_bmp$conf.int[1],wt_bmp$conf.int[2],
                wt_pds$conf.int[1],wt_pds$conf.int[2],wt_pis$conf.int[1],
                wt_pis$conf.int[2],wt_bms$conf.int[1],wt_bms$conf.int[2],
                thresh,image_id)  
    wilcox_tests = rbind(wilcox_tests,tmp)
  }
}
colnames(wilcox_tests) = c("pdp.pval","pip.pval","bmp.pval","pds.pval",
                           "pis.pval","bms.pval","pdp.ci.lower","pdp.ci.upper","pip.ci.lower","pip.ci.upper","bmp.ci.lower","pdp.ci.upper",
                           "pds.ci.lower","pds.ci.upper","pis.ci.lower","pis.ci.upper","bms.ci.lower","pds.ci.upper",
                           "thresh","image.id")
test = as.data.frame(wilcox_tests,stringsAsFactors = FALSE)
for (col in 1:(ncol(test)-1)){
  test[,col] = as.numeric(test[,col])
}
save(test,file=paste(datadir,"/wilcox_tests_pos_nocorrection.Rda",sep=""))
testadjust = test # Here we will save adjusted (FDR corrected) q values
save(wilcox_tests,file=paste(datadir,"/wilcox_tests_distributions_pos_raw.Rda",sep=""))

# Now we will keep track of percentage of significantly different for each
per_sigdiff = c()

# Now correct FDR 0.05 for each masking strategy
# we are correcting within thresholds
for (thresh in unique(test$thresh)){
  pdp_adj = p.adjust(test$pdp.pval[test$thresh==thresh],method="fdr")
  pip_adj = p.adjust(test$pip.pval[test$thresh==thresh],method="fdr")
  bmp_adj = p.adjust(test$bmp.pval[test$thresh==thresh],method="fdr")
  pds_adj = p.adjust(test$pds.pval[test$thresh==thresh],method="fdr")
  pis_adj = p.adjust(test$pis.pval[test$thresh==thresh],method="fdr")
  bms_adj = p.adjust(test$bms.pval[test$thresh==thresh],method="fdr")
  testadjust$pdp.pval[test$thresh==thresh] = pdp_adj
  testadjust$pip.pval[test$thresh==thresh] = pip_adj
  testadjust$bmp.pval[test$thresh==thresh] = bmp_adj
  testadjust$pds.pval[test$thresh==thresh] = pds_adj
  testadjust$pis.pval[test$thresh==thresh] = pis_adj
  testadjust$bms.pval[test$thresh==thresh] = bms_adj
  pdp_per = length(which(pdp_adj <= 0.05)) / length(pdp_adj)
  pip_per = length(which(pip_adj <= 0.05)) / length(pip_adj)
  bmp_per = length(which(bmp_adj <= 0.05)) / length(bmp_adj)
  pds_per = length(which(pds_adj <= 0.05)) / length(pds_adj)
  pis_per = length(which(pis_adj <= 0.05)) / length(pis_adj)
  bms_per = length(which(bms_adj <= 0.05)) / length(bms_adj)
  per_sigdiff = rbind(per_sigdiff,cbind(pdp_per,pip_per,bmp_per,pds_per,pis_per,bms_per,thresh))
}

save(testadjust,file=paste(datadir,"/wilcox_tests_pos_fdr.Rda",sep=""))
colnames(per_sigdiff) = c("intersect.pearson","union.pearson","brain.mask.pearson","intersect.spearman","union.spearman","brain.mask.spearman","thresh")
save(per_sigdiff,file=paste(datadir,"/wilcox_tests_per_sigdiff_pos_fdr05.Rda",sep=""))
write.table(per_sigdiff,file=paste(datadir,"/wilcox_tests_per_sigdiff_pos_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)
tmp = melt(as.data.frame(per_sigdiff),id.vars=c("thresh"))
colnames(tmp) = c("threshold","masking.strategy","value")
ggplot(tmp,aes(x=threshold,y=value,group=masking.strategy,colour=masking.strategy)) + 
  geom_line(size=1) + 
  ylab("% different from threshold of 0") +
  xlab("Threshold +/-") +
  theme(text = element_text(size=20),legend.position="none") +
  facet_wrap(~masking.strategy)
ggsave(paste(savedir,"/percent_means_sigdiff_pos_fdr05.png",sep=""))

# For which thresholds are masking strategies pearson score distributions sig. dif?
# Let's get an overall pvalue for each threshold
wilcox_thresh = c()
for (thresh in thresholds){
  cat("Threshold:",thresh,"\n")
  df = get_direction_result(results,"pos")
  gs = df[df$thresh==0.0,]
  wt_pdp = wilcox_test(df,thresh,gs,"intersect.pearson")
  wt_pip = wilcox_test(df,thresh,gs,"union.pearson")
  wt_bmp = wilcox_test(df,thresh,gs,"brain.mask.pearson")
  wt_pds = wilcox_test(df,thresh,gs,"intersect.spearman")
  wt_pis = wilcox_test(df,thresh,gs,"union.spearman")
  wt_bms = wilcox_test(df,thresh,gs,"brain.mask.spearman")    
  tmp = rbind(
      cbind("intersect.pearson",wt_pdp$p.value,wt_pdp$conf.int[1],wt_pdp$conf.int[2],thresh),
      cbind("union.pearson",wt_pip$p.value,wt_pip$conf.int[1],wt_pip$conf.int[2],thresh),
      cbind("brain.mask.pearson",wt_bmp$p.value,wt_bmp$conf.int[1],wt_bmp$conf.int[2],thresh),
      cbind("intersect.spearman",wt_pds$p.value,wt_pds$conf.int[1],wt_pds$conf.int[2],thresh),
      cbind("union.spearman",wt_pis$p.value,wt_pis$conf.int[1],wt_pis$conf.int[2],thresh),
      cbind("brain.mask.spearman",wt_bms$p.value,wt_bms$conf.int[1],wt_bms$conf.int[2],thresh)      
    )
  wilcox_thresh = rbind(wilcox_thresh,tmp)
}
colnames(wilcox_thresh) = c("strategy","p.value","ci.low","ci.high","thresh")
save(wilcox_thresh,file="wilcox_thresh_pos_all_uncorrected.Rda")
fdr = p.adjust(as.numeric(wilcox_thresh[,2]),method="fdr")
wilcox_thresh = cbind(wilcox_thresh,fdr)
tmp = as.data.frame(wilcox_thresh,stringsAsFactors=FALSE)
tmp$p.value = as.numeric(tmp$p.value)
tmp$ci.low = as.numeric(tmp$ci.low)
tmp$ci.high = as.numeric(tmp$ci.high)
tmp$fdr = as.numeric(tmp$fdr)
write.table(tmp,file=paste(datadir,"/wilcox_tests_distributions_corrected_pos.tsv",sep=""),sep="\t",row.names=FALSE)
save(tmp,file=paste(datadir,"/wilcox_tests_distributions_corrected_pos.Rda",sep=""))
tmp = melt(tmp,id.vars=c("strategy","thresh","ci.low","ci.high","fdr"))
write.table(tmp,file=paste(datadir,"/wilcox_tests_distributions_corrected_pos_flat.tsv",sep=""),sep="\t",row.names=FALSE)

# Significantly different means
tmp[which(tmp$fdr<0.01),]


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
direction = "pos"

# For each image, for each threshold, we assess distance from the "gold standard" ordering based on task/contrast
image_ids = unique(results[["pdp"]]$UID)
labels = c("intersect.pearson","union.pearson","brain.mask.pearson","intersect.spearman","union.spearman","brain.mask.spearman")

acc_df = c()
for (i in 1:length(image_ids)){
  #cat("Processing",i,"of",length(image_ids),"\n","Threshold:","\n")
  image_id = image_ids[i]
  other_ids = image_ids[-which(image_ids==image_id)]
  gs = make_gold_standard_ranking(image_id,other_ids)
  # For each of pearson, spearman [union, intersect, brainmask] get scores for image_id
  df = get_single_result(image_id,results,direction=direction)   
  for (thresh in thresholds){
    #cat(thresh,",",sep="")
    for (label in labels){
      # Get ordering based on actual scores
      sorted = filter_single_result(df,thresh,label,other_ids,image_id)
      acc = calculate_accuracy(gs,sorted)
      
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
save(acc_df,file=paste(datadir,"/accuracy_df_",direction,".Rda",sep=""))
write.table(acc_df,file=paste(datadir,"/accuracy_df_",direction,".tsv",sep=""),sep="\t")
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

subset_summary = ddply(subset, c("strategy","thresh","standard"), summarise, mscore = mean(accuracy), up=get_ci(accuracy,"upper"), down=get_ci(accuracy,"lower"))
ggplot(subset_summary, aes(x=thresh,y=mscore,ymax=up,ymin=down, fill=strategy,colour=strategy,standard=standard)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  facet_wrap(~standard) +
  ylab("Accuracy")
ggsave(paste(savedir,"/chunk_ranking_accuracy",direction,"_withCI.png",sep=""))

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
acc_summary = ddply(all_acc, c("strategy","thresh"), summarise, mscore = mean(acc), up=get_ci(acc,"upper"), down=get_ci(acc,"lower"))
ggplot(acc_summary, aes(x=thresh,y=mscore,ymax=up,ymin=down, fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  ylab("Accuracy")
ggsave(paste(savedir,"/ml_accuracy",direction,"_withCI.png",sep=""))

# Woot!