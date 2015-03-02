# First, read in all input files

input_files =list.files("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores_all",pattern="144*")
datadir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores_all"
setwd(datadir)
gs_file = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/atlas_scores_all/gs_comparisons.tsv"

# IMPORTANT: the pos_only variable corresponds to "absolute_value" in the running scripts, and "True" means
# that we include both negative and positive values. "False" would mean including only positive values.
# For most / all analyses, we want pos_only == "True"
pearsons_gs = read.csv(gs_file,sep="\t",row.names=1)
pearsons_gs = as.matrix(pearsons_gs)
diag(pearsons_gs) = NA
pearsons_pd = read.csv(input_files[3],sep="\t")
pearsons_pi = read.csv(input_files[4],sep="\t")
pearsons_bm = read.csv(input_files[2],sep="\t")
pd_sizes = read.csv(input_files[5],sep="\t")
pi_sizes = read.csv(input_files[6],sep="\t")
bm_sizes = read.csv(input_files[1],sep="\t")
#ts_means = read.csv(input_files[8],sep="\t")
#ts_sds = read.csv(input_files[9],sep="\t")

thresholds = unique(pearsons_pd$thresh)
setwd("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/image_comparison/experiments")
source("experiment2_functions.R")
results = list(gs=pearsons_gs,pd=pearsons_pd[,-1],pi=pearsons_pi[,-1],bm=pearsons_bm[,-1])

# PART I: Assessing significant differences in DISTRIBUTIONS / VARIANCE / MEANS of scores

# QUALITATIVE ASSESSMENT (look at change in mean scores for all images over thresholds)
# First we will look at density plots of the pearson scores
savedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/img/all"
for (thresh in thresholds){
  plot_result(results,thresh,"True")
  ggsave(paste(savedir,"/pearson_density_posonly",thresh,".png",sep=""))
}

# Next we will look at boxplots of pearson scores
for (thresh in thresholds){
  plot_result(results,thresh,plot_type="boxplot")
  ggsave(paste(savedir,"/pearson_boxplot_",thresh,".png",sep=""))
}

# Now we want to look at change in means over time with a timeseries

# - START: did not use this portion -----------------------------------------
# We need to format data as follows:
# pearsonr pearsonr.low pearsonr.up thresh masking
maskings = colnames(ts_means)[1:4]
TS = c()
for (mask in maskings){
  mean_subset = ts_means[,c(mask,"thresh")]
  sd_subset = ts_sds[,c(mask,"thresh")]
  high_subset = sd_subset[,1] + mean_subset[,1]
  masking = rep(mask,nrow(mean_subset))
  low_subset = mean_subset[,1] - sd_subset[,1]
  tmp = cbind(mean_subset,high_subset,low_subset,masking)
  colnames(tmp)[1] = "pearsonr"
  TS = rbind(TS,tmp)
}

# Here we are looking at the mean pearsons (+/- one standard deviation)
ggplot(TS, aes(x=thresh, y=pearsonr, ymin=low_subset, ymax=high_subset, colour=masking,fill=masking)) + 
      geom_line(size=1) + 
      geom_ribbon(alpha=0.1) +
      xlab("Threshold +/-") +
      ylab("Pearsons R")
# - END: did not use this portion -------------------------------------------

# This isn't really what we want - we want to see real confidence intervals!
# We need ALL score data in format of:
# pearsonr thresh masking

# First we need to generate the same thing for the gold standard at every threshold
GS = c()
for (thresh in thresholds){
  vector = as.vector(as.matrix(pearsons_gs))
  vector = vector[-which(is.na(vector))]
  GS = rbind(GS,cbind(vector,rep(thresh,length(vector))))
}
GS = cbind(GS,rep("base.standard",nrow(GS)))
colnames(GS) = c("pearsonr","thresh","masking")

# Now do for each of pearson
PD = flatten_data(pearsons_pd,pos_only="True",label="pairwise.deletion")
PI = flatten_data(pearsons_pi,pos_only="True",label="pairwise.inclusion")
BM = flatten_data(pearsons_bm,pos_only="True",label="brain.mask")

library(plyr)
plot.new()
# Put them all together and plot
ALL = as.data.frame(rbind(GS,PD,PI,BM),stringsAsFactors=FALSE)
ALL$pearsonr = as.numeric(ALL$pearsonr)
save(ALL,file=paste(datadir,"/all_pearsons_posonly.Rda",sep=""))
ALLSUM = ddply(ALL, c("masking","thresh"), summarise, mpearsonr = mean(pearsonr), up=get_ci(pearsonr,"upper"), down=get_ci(pearsonr,"lower"))
ggplot(ALLSUM, aes(x=thresh, group=masking,y=mpearsonr,ymin=down,ymax=up,fill=masking,colour=masking)) + 
  geom_line(size=1.5) + 
  geom_ribbon(alpha=0.15,linetype=0) +
  xlab("Threshold +/-") +
  ylab("Pearsons R") + title("Pearson Scores with Different Masking Strategies")
ggsave(paste(savedir,"/pearson_timeseries_posonly.png",sep=""))
save(ALLSUM,file=paste(datadir,"/allsum_pearsons_posonly.Rda",sep=""))
ggplot(ALL,aes(x=pearsonr, fill=masking)) + geom_density(alpha=0.25) + ylab("Density") + xlab("Pearsons R") + facet_wrap(~thresh)
ggsave(paste(savedir,"/pearson_densities_posonly.png",sep=""))

# Save a version that only goes up to 3.02 threshold
ggplot(ALL[which(ALL$thresh<3.03),],aes(x=pearsonr, fill=masking)) + geom_density(alpha=0.25) + ylab("Density") + xlab("Pearsons R") + facet_wrap(~thresh)
ggsave(paste(savedir,"/pearson_densities_lt3.png",sep=""))

# QUANTITATIVE ASSESSMENT (look at change in distribution of scores for each image)
# we want to know "how often" it looks different, or we get weird results.

# Part A: use thresholding as a proxy for coverage - when I threshold an image, number of voxels decreases

# We will be doing wilcox tests to assess for differences in means for each masking strategy, for each threshold, vs. the gold standard
# We will save all results into a data frame, and correct for multiple comparisons within thresholds and masking strategies
# Then we can calculate the percentage of the time, for each strategy and threshold, that we have sig. different result
# We could then represent an actual database based on the content and estimate the % of time / threshold we get different results from gs
pos_only="False"
results = list(gs=pearsons_gs,pd=pearsons_pd,pi=pearsons_pi,bm=pearsons_bm)
image_ids = rownames(results[["gs"]])
plot_single = FALSE

wilcox_tests = c()

for (image_id in image_ids){
  df = get_single_result(image_id,results,pos_only=pos_only) 
  
  if (plot_single==TRUE){
    # Plot the distribution of values for each strategy and threshold
    ggplot(df,aes(y=value,x=strategy, fill=strategy)) +
      facet_wrap(~thresh) +
      geom_violin(alpha=0.25) +
      title("Masking Strategy Influence on Pearson R for Image",image_id) +
      ylab("Density") +
      xlab("Pearsons R")
    
    ggplot(df,aes(x=value, fill=strategy)) +
      facet_wrap(~thresh) +
      geom_density(alpha=0.25) +
      title("Masking Strategy Influence on Pearson R for Image",image_id) +
      ylab("Density") +
      xlab("Pearsons R") 
  
    ggsave(paste(savedir,"/single_map_posneg.png",sep=""))
  }
  # We will do wilcox test to determine when we have significantly different means
  df_gs = df[df$strategy=="GS",]
  df_gs = df_gs$value[df_gs$thresh==0.00]
  for (thresh in thresholds){
    df_pd = df[df$strategy=="PD",]
    df_pd = df_pd$value[df_pd$thresh==thresh]  
    wt_pd = wilcox.test(df_gs, y = df_pd,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
    df_pi = df[df$strategy=="PI",]
    df_pi = df_pi$value[df_pi$thresh==thresh]  
    wt_pi = wilcox.test(df_gs, y = df_pi,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
    df_bm = df[df$strategy=="BM",]
    df_bm = df_bm$value[df_bm$thresh==thresh]  
    wt_bm = wilcox.test(df_gs, y = df_bm,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
    tmp = cbind(wt_pd$p.value,wt_pi$p.value,wt_bm$p.value,wt_pd$conf.int[1],wt_pd$conf.int[2],wt_pi$conf.int[1],wt_pi$conf.int[2],wt_bm$conf.int[1],wt_bm$conf.int[2],thresh,image_id)  
    wilcox_tests = rbind(wilcox_tests,tmp)
  }
}
colnames(wilcox_tests) = c("pd.pval","pi.pval","bm.pval","pd.ci.lower","pd.ci.upper","pi.ci.lower","pi.ci.upper","bm.ci.lower","pd.ci.upper","thresh","image.id")
test = as.data.frame(wilcox_tests,stringsAsFactors = FALSE)
for (col in 1:ncol(test)){
  test[,col] = as.numeric(test[,col])
}
save(test,file=paste(datadir,"/wilcox_tests_all_nocorrection.Rda",sep=""))
testadjust = test # Here we will save adjusted (FDR corrected) q values

# For which thresholds are masking strategies pearson score distributions sig. dif?
# Let's get an overall pvalue for each threshold
wilcox_thresh = c()
for (thresh in thresholds){
  cat("Threshold:",thresh,"\n")
  df_gs = df[df$strategy=="GS",]
  df_gs = df_gs$value[df_gs$thresh==thresh]
  df_pd = df[df$strategy=="PD",]
  df_pd = df_pd$value[df_pd$thresh==thresh]
  wt_pd = wilcox.test(df_gs, y = df_pd,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
  df_pi = df[df$strategy=="PI",]
  df_pi = df_pi$value[df_pi$thresh==thresh]
  wt_pi = wilcox.test(df_gs, y = df_pi,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
  df_bm = df[df$strategy=="BM",]
  df_bm = df_bm$value[df_bm$thresh==thresh]
  wt_bm = wilcox.test(df_gs, y = df_bm,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
  tmp = rbind(
      cbind("pd",wt_pd$p.value,wt_pd$conf.int[1],wt_pd$conf.int[2],thresh),
      cbind("pi",wt_pi$p.value,wt_pi$conf.int[1],wt_pi$conf.int[2],thresh),
      cbind("bm",wt_bm$p.value,wt_bm$conf.int[1],wt_bm$conf.int[2],thresh)
    )
  wilcox_thresh = rbind(wilcox_thresh,tmp)
}
colnames(wilcox_thresh) = c("strategy","p.value","ci.low","ci.high","thresh")
fdr = p.adjust(as.numeric(wilcox_thresh[,2]),method="fdr")
wilcox_thresh = cbind(wilcox_thresh,fdr)
write.table(wilcox_thresh,file=paste(datadir,"/wilcox_tests_corrected_all.tsv",sep=""),sep="\t",row.names=FALSE)
tmp = as.data.frame(wilcox_thresh,stringsAsFactors=FALSE)
tmp$p.value = as.numeric(tmp$p.value)
tmp$ci.low = as.numeric(tmp$ci.low)
tmp$ci.high = as.numeric(tmp$ci.high)
tmp = tmp[,-6]
tmp = melt(tmp,id.vars=c("strategy","thresh","ci.low","ci.high"))
write.table(tmp,file=paste(datadir,"/wilcox_tests_uncorrected_withci_all.tsv",sep=""),sep="\t",row.names=FALSE)

# Significantly different means
tmp[which(tmp$value<0.05),]

# Now we will keep track of percentage of significantly different for each
per_sigdiff = c()

# Now correct FDR 0.001 for each masking strategy
# we are not correcting within thresholds
for (thresh in unique(test$thresh)){
  pd_adj = p.adjust(test$pd.pval[test$thresh==thresh],method="fdr")
  pi_adj = p.adjust(test$pi.pval[test$thresh==thresh],method="fdr")
  bm_adj = p.adjust(test$bm.pval[test$thresh==thresh],method="fdr")
  testadjust$pd.pval[test$thresh==thresh] = pd_adj
  testadjust$pi.pval[test$thresh==thresh] = pi_adj
  testadjust$bm.pval[test$thresh==thresh] = bm_adj
  pd_per = length(which(pd_adj <= 0.05)) / length(pd_adj)
  pi_per = length(which(pi_adj <= 0.05)) / length(pi_adj)
  bm_per = length(which(bm_adj <= 0.05)) / length(bm_adj)
  per_sigdiff = rbind(per_sigdiff,cbind(pd_per,pi_per,bm_per,thresh))
}

save(testadjust,file=paste(datadir,"/wilcox_tests_all_fdr.Rda",sep=""))
colnames(per_sigdiff) = c("pd.percent_diff","pi.percent_diff","bm.percent_diff","thresh")
per_sigdiff = per_sigdiff[ order(per_sigdiff[,4]), ]
save(per_sigdiff,file=paste(datadir,"/wilcox_tests_per_sigdiff_fdr05.Rda",sep=""))
write.table(per_sigdiff,file=paste(datadir,"/wilcox_tests_per_sigdiff_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)
tmp = melt(as.data.frame(per_sigdiff),id.vars=c("thresh"))
colnames(tmp) = c("threshold","masking.strategy","value")
ggplot(tmp,aes(x=threshold,y=value,group=masking.strategy,colour=masking.strategy)) + geom_line(size=1) + ylab("% different from base standard")
ggsave(paste(savedir,"/percent_pearsonmeans_sigdiff_fdr05.png",sep=""))

# Part B: use region parcellation as a proxy for coverage - when I threshold an image, number of voxels decreases

# TODO: need to do this analysis!

# PART II: Assessing significant differences in ORDERING / RANKING of similar images
# We will save our results in a data frame with format
# imageid strategy rho rho-pvalue tau tau-pvalue
all_df = c()
abs_values = c("True","False")

for (abs_v in abs_values){
  df = c()
  for (thresh in thresholds){  
    # For each *gsr* and each masking strategy *pd*,*pi*,and *bm*
    for (g in 1:nrow(pearsons_gs)){
      rowname = rownames(pearsons_gs)[g]
      gsr = seq(1,ncol(pearsons_gs)-1)
      # Names of gsr correspond with the image id associated with the ranking
      names(gsr) = gsub("X","",(names(sort(abs(pearsons_gs[g,]),decreasing=TRUE))))
      # For each of pdr,pir,bmr, we get order based on the names in gsr
      pdr = get_similar(pearsons_pd,rowname,thresh,abs_v)
      pdr = gsr[pdr]      
      pir = get_similar(pearsons_pi,rowname,thresh,abs_v)
      pir = gsr[pir]
      bmr = get_similar(pearsons_bm,rowname,thresh,abs_v)
      bmr = gsr[bmr]
      # PAIRWISE DELETION
      s = get_sets(gsr,pdr)
      rho = cor.test(s$set1,s$set2, method = c("spearman"), conf.level = 0.95)
      tau = cor.test(s$set1,s$set2, method = c("kendall"), conf.level = 0.95)
      row = cbind(rowname,"PD",rho$p.value,rho$estimate,tau$p.value,tau$estimate,thresh)
      df = rbind(df,row)
      # PAIRWISE INCLUSION
      s = get_sets(gsr,pir)
      rho = cor.test(s$set1,s$set2, method = c("spearman"), conf.level = 0.95)
      tau = cor.test(s$set1,s$set2, method = c("kendall"), conf.level = 0.95)
      row = cbind(rowname,"PI",rho$p.value,rho$estimate,tau$p.value,tau$estimate,thresh)
      df = rbind(df,row)
      # BRAIN MASK
      s = get_sets(gsr,bmr)
      rho = cor.test(s$set1,s$set2, method = c("spearman"), conf.level = 0.95)
      tau = cor.test(s$set1,s$set2, method = c("kendall"), conf.level = 0.95)
      row = cbind(rowname,"BM",rho$p.value,rho$estimate,tau$p.value,tau$estimate,thresh)
      df = rbind(df,row)
    }
  }
  colnames(df) = c("imageid","strategy","rho_pvalue","rho","tau_pvalue","tau","thresh")
  rownames(df) = seq(1,nrow(df))
  df = as.data.frame(df,stringsAsFactors=FALSE)
  df$rho_pvalue = as.numeric(df$rho_pvalue)
  df$tau_pvalue = as.numeric(df$tau_pvalue)
  all_df[[abs_v]] = df
}
save(all_df,file=paste(datadir,"/all_taurho_scores_uncorrected.Rda",sep=""))

# This function returns FDR corrected qvalues (c) for each threshold
counts = c()
for (abs_v in abs_values){
  df = all_df[[abs_v]]
  for (thresh in thresholds){
    subset = df[which(df$thresh==thresh),]
    c = plot_pval(subset,thresh,savedir,0.05)
    c = cbind(c,rep(abs_v,nrow(c)),rep(thresh,nrow(c)))
    counts = rbind(counts,c)
  }
}
colnames(counts)[5:6] = c("abs_value","thresh")
save(counts,file=paste(datadir,"/counts_taurho_fdr0.05.Rda",sep=""))
write.table(counts,file=paste(datadir,"/counts_taurho_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)

# Which one are rho vs tau?
rho = counts[grep("RHO",counts$STRATEGY),]
tau = counts[grep("TAU",counts$STRATEGY),]
rho = rho[order(-rho$perc_diff),]
tau = tau[order(-tau$perc_diff),]

# Save to table
write.table(rho,file=paste(datadir,"/rho_table_all_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)
write.table(tau,file=paste(datadir,"/tau_table_all_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)

# Show the tables
rho
tau

# Finally plot the percentages
subset = tau[tau$abs_value=="False",]
subset$STRATEGY[subset$STRATEGY == "BM_TAU"] = "BM"
subset$STRATEGY[subset$STRATEGY == "PD_TAU"] = "PD"
subset$STRATEGY[subset$STRATEGY == "PI_TAU"] = "PI"
ggplot(subset, aes(x=thresh,y=perc_diff,color=STRATEGY)) +
   geom_line(size=2,alpha=0.25,stat="identity") +
   ylab("% significantly different rankings") + xlab("threshold") + facet_wrap(~STRATEGY)
ggsave(paste(savedir,"/tau_sigdiff_posonly_fdr05.png",sep=""))

subset = rho[rho$pos_only=="True",]
ggplot(subset, aes(x=thresh,y=perc_diff,color=STRATEGY)) + geom_line(size=2,alpha=0.25,stat="identity") + title("Percentage Significantly Different Images, RHO") + ylab("percent") + xlab("") + facet_wrap(~ STRATEGY) + scale_x_continuous(breaks=unique(subset$thresh),labels=unique(subset$thresh))
ggsave(paste(savedir,"/rho_sigdiff_posonly.png",sep=""))


# Finally plot the percentages
subset = tau[tau$pos_only=="False",]
ggplot(subset, aes(x=thresh,y=perc_diff,color=STRATEGY)) + geom_line(size=2,alpha=0.25,stat="identity") + title("Percentage Significantly Different Images, RHO") + ylab("percent") + xlab("") + facet_wrap(~ STRATEGY) + scale_x_continuous(breaks=unique(subset$thresh),labels=unique(subset$thresh))
ggsave(paste(savedir,"/tau_sigdiff_all.png",sep=""))

subset = tau[tau$pos_only=="True",]
ggplot(subset, aes(x=thresh,y=perc_diff,color=STRATEGY)) + geom_line(size=2,alpha=0.25,stat="identity") + title("Percentage Significantly Different Images, RHO") + ylab("percent") + xlab("") + facet_wrap(~ STRATEGY) + scale_x_continuous(breaks=unique(subset$thresh),labels=unique(subset$thresh))
ggsave(paste(savedir,"/tau_sigdiff_posonly.png",sep=""))

# Finally, we want to see how the mask sizes change
# First we will look at the means
ncomps_all = c()
for (abs_v in abs_values){
  ncomps_means = c()
  for (thresh in thresholds){
    pdmean = mean(get_number_comparisons(pd_sizes,thresh,abs_v))
    pimean = mean(get_number_comparisons(pi_sizes,thresh,abs_v))
    bmmean = mean(get_number_comparisons(bm_sizes,thresh,abs_v))
    tmp = cbind(pdmean,pimean,bmmean)
    ncomps_means = rbind(ncomps_means,tmp)
  }
  rownames(ncomps_means) = thresholds
  colnames(ncomps_means) = c("pairwise.deletion","pairwise.inclusion","brain.mask")
  ncomps_all[[abs_v]] = ncomps_means
}
save(ncomps_all,file=paste(datadir,"/mean_comparison_counts.Rda",sep=""))

# Now we will look at the distributions of values themselves!
ncomps_all = c()
for (abs_v in abs_values){
  ncomps = c()
  cat("ABSOLUTE VALUE",abs_v,"\n")
  for (thresh in thresholds){
    cat("  THRESH",thresh,"\n")
    tmp = c()
    pd = unlist(get_masksize_comparisons(pd_sizes,thresh,abs_v))
    pi = unlist(get_masksize_comparisons(pi_sizes,thresh,abs_v))
    bm = unlist(get_masksize_comparisons(bm_sizes,thresh,abs_v))
    tmp[["pairwise.deletion"]] = pd
    tmp[["pairwise.inclusion"]] = pi
    tmp[["brain.mask"]] = bm
    ncomps[[as.character(thresh)]] = tmp
  }
  ncomps_all[[abs_v]] = ncomps
}
save(ncomps_all,file=paste(datadir,"/dist_masksizes_all.Rda",sep=""))

# Plot the number of comparisons - first put into a big data frame
df = c()
for (abs_v in abs_values){
  subset = ncomps_all[[abs_v]] 
  ts = names(subset)
  for (t in ts){
    tsubset = subset[[t]]
    tmp1 = cbind(tsubset$pairwise.deletion,rep(t,length(tsubset$pairwise.deletion)),rep(abs_v,length(tsubset$pairwise.deletion)),rep("PD",length(tsubset$pairwise.deletion)))
    tmp2 = cbind(tsubset$pairwise.inclusion,rep(t,length(tsubset$pairwise.inclusion)),rep(abs_v,length(tsubset$pairwise.inclusion)),rep("PI",length(tsubset$pairwise.inclusion)))
    tmp3 = cbind(tsubset$brain.mask,rep(t,length(tsubset$brain.mask)),rep(abs_v,length(tsubset$brain.mask)),rep("BM",length(tsubset$brain.mask)))
    tmp = rbind(tmp1,tmp2,tmp3)
    df = rbind(df,tmp)
  }  
}
colnames(df) = c("mask.size","thresh","pos_only","strategy")
df = as.data.frame(df)
df$mask.size = as.numeric(as.character(df$mask.size))
save(df,file=paste(datadir,"/dist_masksizes_flat_df.Rda",sep=""))

subset = df[df$pos_only=="True",-which(colnames(df)=="pos_only")]
tmp = ddply(subset, c("thresh","strategy"), summarise, mask.size.mean=mean(mask.size),mask.size.min = min(mask.size), mask.size.max=max(mask.size))

# Here we are looking at the mean pearsons (+/- one standard deviation)
ggplot(tmp, aes(x=thresh, y=mask.size.mean, ymin=mask.size.min, ymax=mask.size.max,colour=strategy,fill=strategy,group=strategy)) + 
  geom_line(size=1) + 
  title("Mask Sizes at Different +/- Thresholds") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +/-") +
  ylab("Mask Size (voxels)")
ggsave(paste(savedir,"/mask_sizes_atthresh_posneg.png",sep=""))

subset = df[df$pos_only=="True",-which(colnames(df)=="pos_only")]
tmp = ddply(subset, c("thresh","strategy"), summarise, mask.size.mean=mean(mask.size),mask.size.min = min(mask.size), mask.size.max=max(mask.size))

# Here we are looking at the mean pearsons (+/- one standard deviation)
ggplot(tmp, aes(x=thresh, y=mask.size.mean, ymin=mask.size.min, ymax=mask.size.max,colour=strategy,fill=strategy,group=strategy)) + 
  geom_line(size=1) + 
  title("Mask Sizes at Different + Thresholds") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +") +
  ylab("Mask Size")
ggsave(paste(savedir,"/mask_sizes_atthresh_posonly.png",sep=""))
