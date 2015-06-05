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
#save(ALL,file=paste(datadir,"/all_scores.Rda",sep=""))
load(file=paste(datadir,"/all_scores.Rda",sep=""))

# Plot distributions separately
ALL = ALL[ALL$thresh<=7.0,]
ALL = ALL[ALL$strategy=="cca.pearson",]
ggplot(ALL[ALL$direction=="posneg",],aes(x=value, fill=strategy, thresh=thresh)) + 
  geom_density(alpha=0.25) + 
  facet_wrap(~thresh,ncol=2) +
  ylab("Density") + 
  xlab("Threshold +/-") + 
  xlim(-1,1) +
  theme(legend.position="none")
ggsave(paste(savedir,"/density_ccapearson_posneg.png",sep=""))

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
iters = 0:499

# For each of 500 permutations:
#  Subset data to unrelated groups A and B
#     For each strategy for handling missingness
#       For each unthresholded map, Ai in A
#         Apply each threshold in Z=+/-0:13, and Z=+0:13 to all of B
#           Calculate similarity for each of B to Ai                
#           Assign correct classification if contrast Ai most similar to equivalent contrast in B
# Calculate (and save) accuracy


setwd("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3_unrelated/permutations")

# We will save a data frame of accuracies
directions = c("pos","posneg")
thresholds=1:13
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
    save(allres,file="accuracy_df_finished.Rda")
  }
}

direction="posneg"
ggplot(all[all$direction==direction,], aes(x=thresh,y=accuracy_mean,ymax=up,ymin=down, fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  ylab("Accuracy") +
  ylim(0,1) +
  scale_x_continuous(breaks = round(seq(0, 13, by = 1.0),1))
ggsave(paste(savedir,"/ml_accuracy",direction,"_withCI.png",sep=""))
write.table(allres,file=paste(datadir,"/ml_accuracy.tsv",sep=""),sep="\t")

# Woot!

# Finally, we want to create a confusion matrix - and get a sense for which contrasts are misclassified
iters=0:499
setwd("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3_unrelated/permutations")

# We will save a data frame of accuracies
directions = c("pos","posneg")
thresholds=0:13
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
      piacc = get_classifications(pi,queryimages,thresh,direction,"svi.pearson")
      pdacc = get_classifications(pd,queryimages,thresh,direction,"cca.pearson")
      siacc = get_classifications(si,queryimages,thresh,direction,"svi.spearman")
      sdacc = get_classifications(sd,queryimages,thresh,direction,"cca.spearman")
      allacc = rbind(piacc,pdacc,siacc,sdacc)
      res = rbind(res,allacc)
    }
    # Calculate the mean across thresholds and directions
    res$direction = as.character(res$direction)
    res = ddply(res,c("direction","thresh","strategy"), summarise, accuracy_mean=mean(acc),up=get_ci(acc,"upper"),down=get_ci(acc,"lower"))
    allres = rbind(allres,res)
    save(allres,file="accuracy_df_finished.Rda")
  }
}

direction="posneg"
ggplot(all[all$direction==direction,], aes(x=thresh,y=accuracy_mean,ymax=up,ymin=down, fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  ylab("Accuracy") +
  ylim(0,1) +
  scale_x_continuous(breaks = round(seq(0, 13, by = 1.0),1))
ggsave(paste(savedir,"/ml_accuracy",direction,"_withCI.png",sep=""))
write.table(allres,file=paste(datadir,"/ml_accuracy.tsv",sep=""),sep="\t")
