library(plyr)

# First, read in all input files
datadir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3_unrelated"
setwd(datadir)

# We will take an average for each sample, for each threshold, strategy, score
pearsons_pd = parse_input(read.csv("pearsons_cca.tsv",sep="\t",row.names=1))
pearsons_pi = parse_input(read.csv("pearsons_svi.tsv",sep="\t",row.names=1))
spearmans_pd = parse_input(read.csv("spearmans_cca.tsv",sep="\t",row.names=1)) 
spearmans_pi = parse_input(read.csv("spearmans_svi.tsv",sep="\t",row.names=1))

setwd("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/image_comparison/experiments/experiment3")
source("helper_functions.R")
setwd(datadir)

thresholds = sort(unique(pearsons_pi$thresh))

# Part I: Handling Missing Values to Optimize Similarity Search Ranking
# ? Visualize: How do distributions of scores change? ###################################################################

# NOTE: R can get buggy and not correctly return plots from a function, so this was run manually.
savedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/img/experiment3_unrelated"
results = list(pds=spearmans_pd,pis=spearmans_pi,pdp=pearsons_pd,pip=pearsons_pi)
for (thresh in thresholds){
  outfile = paste(savedir,"/density_posneg",thresh,".png",sep="")
  #plot_result(results,thresh,"posneg",outfile)  
}

for (thresh in thresholds){
  outfile = paste(savedir,"/density_pos",thresh,".png",sep="")
  #plot_result(results,thresh,"pos",outfile)
}

# ? Visualize: How do mean scores change? ##############################################################################

plot.new()
# Put them all together and plot
pearsons_pd$strategy = "cca.pearson"
pearsons_pi$strategy = "svi.pearson"
spearmans_pd$strategy = "cca.spearman"
spearmans_pi$strategy = "svi.spearman"
ALL = as.data.frame(rbind(pearsons_pd,pearsons_pi,spearmans_pi,spearmans_pd),stringsAsFactors=FALSE)
ALL$score = as.numeric(ALL$value)
ALL$thresh = as.numeric(ALL$thresh)
ALL$strategy = as.character(ALL$strategy)
save(ALL,file=paste(datadir,"/all_scores.Rda",sep=""))
#load(file=paste(datadir,"/all_scores.Rda",sep=""))

direction = "pos"
ALLSUM = ALL[ALL$direction==direction,]
ALLSUM = ALLSUM[-which(is.na(ALLSUM$score)),]
ALLSUM = ddply(ALLSUM, c("strategy","thresh"), summarise, mscore = mean(score), up=get_ci(score,"upper"), down=get_ci(score,"lower"))
ggplot(ALLSUM, aes(x=thresh, group=strategy,y=mscore,ymin=down,ymax=up,fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +") +
  ylab("Mean Score") + title("Scores with Different Strategies for Handling Missing Data")
ggsave(paste(savedir,"/meanscores_with_ci_pos.png",sep=""))


# ? Visualize: How do sizes change? ##############################################################################


pd_sizes = parse_input(read.csv("sizes_cca.tsv",sep="\t",row.names=1))
pi_sizes = parse_input(read.csv("sizes_svi.tsv",sep="\t",row.names=1))
pd_sizes$value = as.numeric(pd_sizes$value)
pd_sizes$strategy = "cca"
pi_sizes$strategy = "svi"
pi_sizes$thresh = as.numeric(pi_sizes$thresh)
pd_sizes$thresh = as.numeric(pd_sizes$thresh)
sizes = rbind(pd_sizes,pi_sizes)
direction = "pos"

tmp = ddply(sizes[sizes$direction==direction,], c("thresh","strategy"), summarise, mask.size.mean=mean(value),mask.size.min = min(value), mask.size.max=max(value))


# Here we are looking at the mean pearsons (+/- one standard deviation)
ggplot(tmp, aes(x=thresh, y=mask.size.mean, ymin=mask.size.min, ymax=mask.size.max, fill=strategy,group=strategy,color=strategy)) + 
  geom_line(size=1) + 
  title("Mask Sizes at Different Thresholds +") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +") +
  ylab("Mask Size (voxels)") +
  scale_x_continuous(breaks = round(seq(0, 13, by = 1.0),1))
ggsave(paste(savedir,"/mask_sizes_voxels_",direction,".png",sep=""))

# Now plot percentage of voxels
total_size = max(sizes$value)
df2 = sizes
df2 = df2[df2$direction==direction,]
df2$value = df2$value / total_size
df2 = ddply(df2, c("thresh","strategy"), summarise, mask.size.mean=mean(value),mask.size.min = min(value), mask.size.max=max(value))
ggplot(df2, aes(x=thresh, y=mask.size.mean, ymin=mask.size.min, ymax=mask.size.max, fill=strategy,group=strategy,color=strategy)) + 
  geom_line(size=1) + 
  title("Mask Sizes at Different Thresholds +") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +") +
  ylab("size of comparison set as % of brain mask") + 
  scale_x_continuous(breaks = round(seq(0, 13, by = 1.0),1))
ggsave(paste(savedir,"/mask_percsizes_atthresh_",direction,".png",sep=""))


# Part II: Image Classification
iters = 0:500

# For each of 500 permutations:
#  Subset data to unrelated groups A and B
#     For each strategy for handling missingness
#       For each unthresholded map, Ai in A
#         Apply each threshold in Z=+/-0:13, and Z=+0:13 to all of B
#           Calculate similarity for each of B to Ai                
#           Assign correct classification if contrast Ai most similar to equivalent contrast in B
# Calculate (and save) accuracy

direction = "posneg"

setwd("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3_unrelated/permutations")

# We will save a data frame of accuracies
directions = c("pos","posneg")

allres = c()
for (thresh in thresholds){
  for (direction in directions){
    cat("THRESHOLD:",thresh,"DIRECTION",direction,"\n")
    res = c()
    for (iter in iters){
      # Read in each of the input files
      cat(iter,"\n")
      pi = parse_single_input(read.csv(paste(iter,"_pearson_svi.tsv",sep=""),sep="\t",row.names=1))  
      pd = parse_single_input(read.csv(paste(iter,"_pearson_cca.tsv",sep=""),sep="\t",row.names=1))  
      si = parse_single_input(read.csv(paste(iter,"_spearman_svi.tsv",sep=""),sep="\t",row.names=1))  
      sd = parse_single_input(read.csv(paste(iter,"_spearman_cca.tsv",sep=""),sep="\t",row.names=1))  
      queryimages = unique(pi$queryimage) 
      thresholds = unique(pi$thresh)  
      piacc = calculate_accuracy(pi,queryimages,thresh,direction,"svi.pearson")
      pdacc = calculate_accuracy(pd,queryimages,thresh,direction,"cca.pearson")
      siacc = calculate_accuracy(si,queryimages,thresh,direction,"svi.spearman")
      sdacc = calculate_accuracy(sd,queryimages,thresh,direction,"cca.spearman")
      allacc = rbind(piacc,pdacc,siacc,sdacc)
      res = rbind(res,allacc)
    }
    # Calculate the mean across thresholds and directions
    res = ddply(res,c("direction","thresh","strategy"), summarise, accuracy_mean=mean(acc),up=get_ci(acc,"upper"),down=get_ci(acc,"lower"))
    allres = rbind(allres,res)
  }
}

# STOPPED HERE - this needs to finish running!

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