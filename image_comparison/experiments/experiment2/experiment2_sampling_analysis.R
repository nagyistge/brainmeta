# First, read in all input files

input_files =list.files("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/random_sampling",pattern="spearman*")
datadir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/random_sampling"
setwd(datadir)

library(plyr)
library(reshape2)
library(ggplot2)
library(gridExtra)

pearsons = read.csv(input_files[1],sep="\t")
sizes = read.csv(input_files[3],sep="\t")
gs = read.csv(input_files[2],sep="\t")
samples = sort(unique(pearsons$percent_sample))

setwd("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/image_comparison/experiments")
source("experiment2_functions.R")
gs_flat = melt(gs,id.vars=c("X"))
# Remove values of 1 - when comparing an image to itself
gs_flat = gs_flat[-which(gs_flat$value==1),]
gs_flat = cbind(gs_flat,rep("GS",nrow(gs_flat)),rep(NA,nrow(gs_flat)))
colnames(gs_flat)[c(1,4,5)] = c("X","strategy","percent_sample")
sizes_flat = melt(sizes,id.vars=c("X","strategy","percent_sample"))
# Remove values of NA - when comparing an image to itself
sizes_flat = sizes_flat[-which(is.na(sizes_flat$value)),]
pearsons_flat = melt(pearsons,id.vars=c("X","strategy","percent_sample"))
# Remove values of NA - when comparing an image to itself
pearsons_flat = pearsons_flat[-which(is.na(pearsons_flat$value)),]

# PART I: How does coverage influence pearson scores for each masking strategy?
# We are basically sampling N ROIs at a time, a proxy for increasing the amount of coverage, to generate a parcellation with which to test each masking strategy. This resulted in a distribution of pearson r scores for each of the Ui, 
# for each parcellation P, for each masking strategy S.
savedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/img/spearman_sampling"
plot.new()
for (n in samples){
  # First subset the data to only include the nroi of interest
  subset = pearsons_flat[which(pearsons_flat$percent_sample==n),]
  tmp = gs_flat[,colnames(subset)]
  both = rbind(subset,tmp)
  plot1 = ggplot(both,aes(x=value, fill=strategy)) + 
    geom_density(alpha=0.25) + 
    #ylim(0,3) +
    ylab("Density") + 
    xlab(paste(n,"percent random sampling"))
  # We also need the size data in the next plot
  size_subset = sizes_flat[which(sizes_flat$percent_sample==n),]
  plot2 = ggplot(size_subset,aes(y=value,x=strategy, fill=strategy)) + 
    geom_boxplot(alpha=0.25) + 
    ylab("number voxels") + 
    ylim(0,max(sizes_flat$value)) +
    xlab("total size of sample") +
    guides(fill=FALSE)
  g = arrangeGrob(plot2, plot1, ncol=2)
  ggsave(file=paste(savedir,"/density_plot/spearman_sampling_dens_",n,".png",sep=""),g)
}

# Try making one big plot (this did not work)
#ALL = melt(pearsons_flat[,-which(colnames(pearsons_flat) %in% c("X","variable"))],id.vars=c("n_roi","strategy"),meanr=mean(value))
#ggplot(ALL, aes(group=strategy,y=value,fill=strategy,colour=strategy)) + 
#  geom_density(alpha=0.25) +
#  facet_wrap(~n_roi) +
#  ylab("Pearson R (mean)") +
#  xlab("Number of Regions (Coverage)")
#ggsave(paste(savedir,"/pearson_timeseries_posonly.png",sep=""))

# Now we can assess when the means are significantly different
wilcox_tests = c()
count=1
for (n in nroi){
  cat(count,"of",length(nroi),"\n")
  subset = pearsons_flat[which(pearsons_flat$n_roi==n),]
  # PAIRWISE DELETION
  h_gs = gs_flat$value
  h_pd = pearsons_flat$value[pearsons_flat$strategy=="PD"]
  h_pd = sample(h_pd,length(h_gs))
  wt_pd = wilcox.test(h_gs, y = h_pd,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
  # PAIRWISE INCLUSION
  h_gs = gs_flat$value
  h_pi = pearsons_flat$value[pearsons_flat$strategy=="PI"]
  h_pi = sample(h_pi,length(h_gs))
  wt_pi = wilcox.test(h_gs, y = h_pi,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
  # BRAIN MASK
  h_gs = gs_flat$value
  h_bm = pearsons_flat$value[pearsons_flat$strategy=="BM"]
  h_bm = sample(h_bm,length(h_gs))
  wt_bm = wilcox.test(h_gs, y = h_bm,alternative="two.sided",conf.int=TRUE,conf.level=0.95)  
  tmp = cbind(wt_pd$p.value,wt_pi$p.value,wt_bm$p.value,wt_pd$conf.int[1],wt_pd$conf.int[2],wt_pi$conf.int[1],wt_pi$conf.int[2],wt_bm$conf.int[1],wt_bm$conf.int[2],n)  
  wilcox_tests = rbind(wilcox_tests,tmp)
  count = count+1
}

colnames(wilcox_tests) = c("pd.pval","pi.pval","bm.pval","pd.ci.lower","pd.ci.upper","pi.ci.lower","pi.ci.upper","bm.ci.lower","bm.ci.upper","n_roi")
test = as.data.frame(wilcox_tests,stringsAsFactors = FALSE)
for (col in 1:ncol(test)){
  test[,col] = as.numeric(test[,col])
}
save(test,file=paste(datadir,"/wilcox_regional_tests_all_nocorrection.Rda",sep=""))
testadjust = test # Here we will save adjusted (FDR corrected) q values
testadjust$pd.pval = p.adjust(testadjust$pd.pval,method="fdr")
testadjust$pi.pval = p.adjust(testadjust$pi.pval,method="fdr")
testadjust$bm.pval = p.adjust(testadjust$bm.pval,method="fdr")
save(testadjust,file=paste(datadir,"/wilcox_regional_tests_all_fdr.Rda",sep=""))
write.table(testadjust,file=paste(datadir,"/wilcox_regional_tests_all_fdr.tsv",sep=""),sep="\t",row.names=FALSE)

# Put them all together and plot
PD = cbind("pd",testadjust$pd.pval,testadjust$pd.ci.lower,testadjust$pd.ci.upper,test$n_roi)
PI = cbind("pi",testadjust$pi.pval,testadjust$pi.ci.lower,testadjust$pi.ci.upper,test$n_roi)
BM = cbind("bm",testadjust$bm.pval,testadjust$bm.ci.lower,testadjust$bm.ci.upper,test$n_roi)

ALL = as.data.frame(rbind(PD,PI,BM),stringsAsFactors=FALSE)
colnames(ALL) = c("strategy","qval","ci.lower","ci.upper","n_roi")
ALL$qval = as.numeric(ALL$qval)
ggplot(ALL, aes(x=n_roi, group=strategy,y=qval,fill=strategy,colour=strategy)) + 
  geom_line(size=1.5) 
ggsave(paste(savedir,"/wilcox_tests_all_sig_fdr.png",sep=""))

# Now we want to do the same, but calculate a frequency based on individual tests
wilcox_single = c()
image_ids = unique(pearsons_flat$X)
count=1
for (image_id in image_ids){
  cat(count,"of",length(nroi),"\n")
  subset = pearsons_flat[which(pearsons_flat$X==image_id),]
  subset_gs = gs_flat$value[which(gs_flat$X==image_id)]
  for (n in nroi){
    subset_n = subset[which(subset$n_roi==n),]
    # PAIRWISE DELETION
    h_pd = subset_n$value[subset_n$strategy=="PD"]
    wt_pd = wilcox.test(subset_gs, y = h_pd,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
    # PAIRWISE INCLUSION
    h_pi = subset_n$value[subset_n$strategy=="PI"]
    wt_pi = wilcox.test(subset_gs, y = h_pi,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
    # BRAIN MASK
    h_bm = subset_n$value[subset_n$strategy=="BM"]
    wt_bm = wilcox.test(subset_gs, y = h_bm,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
    tmp = rbind(
      cbind("pd",wt_pd$p.value,wt_pd$conf.int[1],wt_pd$conf.int[2],image_id,n),
      cbind("pi",wt_pi$p.value,wt_pi$conf.int[1],wt_pi$conf.int[2],image_id,n),
      cbind("bm",wt_bm$p.value,wt_bm$conf.int[1],wt_bm$conf.int[2],image_id,n)
    )
    wilcox_single = rbind(wilcox_single,tmp)
  }
  count = count+1
}
colnames(wilcox_single) = c("strategy","p.value","ci.low","ci.high","image_id","n_roi")
test = as.data.frame(wilcox_single,stringsAsFactors = FALSE)
for (col in 2:ncol(test)){
  test[,col] = as.numeric(test[,col])
}
save(test,file=paste(datadir,"/wilcox_ind_regional_tests_all_nocorrection.Rda",sep=""))
testadjust = test # Here we will save adjusted (FDR corrected) q values
testadjust$p.value = p.adjust(testadjust$p.value,method="fdr")
save(testadjust,file=paste(datadir,"/wilcox_ind_regional_tests_all_fdr.Rda",sep=""))
write.table(testadjust,file=paste(datadir,"/wilcox_ind_regional_tests_all_fdr.tsv",sep=""),sep="\t",row.names=FALSE)
testadjust$ci.low = as.numeric(testadjust$ci.low)
testadjust$ci.high = as.numeric(testadjust$ci.high)

# Now we will keep track of percentage of significantly different for each
per_sigdiff = c()

# Now correct FDR 0.001 for each masking strategy
# we are not correcting within thresholds
for (n in unique(testadjust$n_roi)){
  pd = testadjust[which(testadjust$strategy=="pd"),]
  pi = testadjust[which(testadjust$strategy=="pi"),]
  bm = testadjust[which(testadjust$strategy=="bm"),]
  pd_per = length(which(pd$p.value[which(pd$n_roi==n)] <= 0.05)) / length(pd$p.value[which(pd$n_roi==n)])
  pi_per = length(which(pi$p.value[which(pi$n_roi==n)] <= 0.05)) / length(pi$p.value[which(pi$n_roi==n)])
  bm_per = length(which(bm$p.value[which(bm$n_roi==n)] <= 0.05)) / length(bm$p.value[which(bm$n_roi==n)])  
  tmp = rbind(cbind("pd",pd_per,n),
    cbind("pi",pi_per,n),
    cbind("bm",bm_per,n))
  per_sigdiff = rbind(per_sigdiff,tmp)
}
colnames(per_sigdiff) = c("strategy","percent_diff","n_roi")
for (col in 2:ncol(per_sigdiff)){
  per_sigdiff[,col] = as.numeric(per_sigdiff[,col])
}
per_sigdiff = per_sigdiff[ order(per_sigdiff[,3]), ]
per_sigdiff = as.data.frame(per_sigdiff)
per_sigdiff$percent_diff = as.numeric(as.character(per_sigdiff$percent_diff))
per_sigdiff$n_roi = as.numeric(as.character(per_sigdiff$n_roi))
save(per_sigdiff,file=paste(datadir,"/wilcox_regional_ind_per_sigdiff_fdr05.Rda",sep=""))
write.table(per_sigdiff,file=paste(datadir,"/wilcox_regional_ind_per_sigdiff_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)
ggplot(per_sigdiff,aes(x=n_roi,y=percent_diff,group=strategy,colour=strategy)) + 
  geom_line(size=1) + 
  ylab("% different from base standard (fdr 0.05)") +
  xlab("number of rois (coverage)")
ggsave(paste(savedir,"/percent_regional_sigdiff_fdr05.png",sep=""))

# Part B: use region parcellation as a proxy for coverage - when I threshold an image, number of voxels decreases
# Finally we need to address the ordering of the results!
image_ids = gs$X
df=c()
count=1
for (image_id in image_ids){
  cat(count,"of",length(image_ids),"\n")
  gsr = seq(1,ncol(gs)-2)
  # Names of gsr correspond with the image id associated with the ranking
  sorted_gs = sort(gs[which(gs$X==image_id),-c(1)],decreasing=TRUE)
  sorted_gs = sorted_gs[-which(names(sorted_gs)==paste("X",image_id,sep=""))]
  names(gsr) = gsub("X","",names(sorted_gs))  
  for (n in samples){
    # For each of pdr,pir,bmr, we get order based on the names in gsr
    # PAIRWISE DELETION
    pdr = pearsons[which(pearsons$strategy=="PD"),-which(colnames(pearsons)=="strategy")]
    pdr = pdr[which(pdr$X==image_id),-which(colnames(pdr)=="X")]
    pdr = pdr[which(pdr$percent_sample==n),-which(colnames(pdr)=="percent_sample")]
    pdr = gsub("X","",names(sort(pdr,decreasing=TRUE)))
    pdr = gsr[pdr]      
    s = get_sets(gsr,pdr)
    tau = cor.test(s$set1,s$set2, method = c("kendall"), conf.level = 0.95)
    row = cbind(image_id,"PD",tau$p.value,tau$estimate,n)
    df = rbind(df,row)
    # PAIRWISE INCLUSION
    pir = pearsons[which(pearsons$strategy=="PI"),-which(colnames(pearsons)=="strategy")]
    pir = pir[which(pir$X==image_id),-which(colnames(pir)=="X")]
    pir = pir[which(pir$percent_sample==n),-which(colnames(pir)=="percent_sample")]
    pir = gsub("X","",names(sort(pir,decreasing=TRUE)))
    pir = gsr[pir]      
    s = get_sets(gsr,pir)
    tau = cor.test(s$set1,s$set2, method = c("kendall"), conf.level = 0.95)
    row = cbind(image_id,"PI",tau$p.value,tau$estimate,n)
    df = rbind(df,row)
    # BRAIN MASK
    bmr = pearsons[which(pearsons$strategy=="BM"),-which(colnames(pearsons)=="strategy")]
    bmr = bmr[which(bmr$X==image_id),-which(colnames(bmr)=="X")]
    bmr = bmr[which(bmr$percent_sample==n),-which(colnames(bmr)=="percent_sample")]
    bmr = gsub("X","",names(sort(bmr,decreasing=TRUE)))
    bmr = gsr[bmr]      
    s = get_sets(gsr,bmr)
    tau = cor.test(s$set1,s$set2, method = c("kendall"), conf.level = 0.95)
    row = cbind(image_id,"BM",tau$p.value,tau$estimate,n)
    df = rbind(df,row)
  }
  count=count+1
}
    
colnames(df) = c("imageid","strategy","tau_pvalue","tau","percent_sample")
rownames(df) = seq(1,nrow(df))
tmp = as.data.frame(df,stringsAsFactors=FALSE)
tmp$tau_pvalue = as.numeric(tmp$tau_pvalue)
tmp$percent_sample = as.numeric(tmp$percent_sample)
save(tmp,file=paste(datadir,"/single_tau_sampleperc_scores_uncorrected.Rda",sep=""))

# This function returns FDR corrected qvalues (c) for each threshold
counts = c()
for (n in samples){
    subset = tmp[which(tmp$percent_sample==n),]
    c = plot_pval_region(subset,n,savedir,0.05)
    c = cbind(c,rep(n,nrow(c)))
    counts = rbind(counts,c)
}

colnames(counts)[5] = c("n_roi")
save(counts,file=paste(datadir,"/counts_tau_fdr0.05.Rda",sep=""))
write.table(counts,file=paste(datadir,"/counts_tau_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)

# Which one are rho vs tau?
tau = counts[grep("TAU",counts$STRATEGY),]
tau = tau[order(-tau$perc_diff),]

# Save to table
write.table(tau,file=paste(datadir,"/tau_table_all_fdr05.tsv",sep=""),sep="\t",row.names=FALSE)

# Show the tables
tau

# Finally plot the percentages
tau$STRATEGY[tau$STRATEGY == "BM_TAU"] = "BM"
tau$STRATEGY[tau$STRATEGY == "PD_TAU"] = "PD"
tau$STRATEGY[tau$STRATEGY == "PI_TAU"] = "PI"
ggplot(tau, aes(x=n_roi,y=perc_diff,color=STRATEGY)) +
   geom_line(size=2,alpha=0.25,stat="identity") +
   ylab("% significantly different rankings") + xlab("number of regions (coverage)") + facet_wrap(~STRATEGY)
ggsave(paste(savedir,"/tau_sigdiff_fdr05.png",sep=""))


# Finally, we want to see how the mask sizes change, so we can ascribe a size with each roi
# and strategy
head(sizes_flat)
sizes_flat$value = as.numeric(as.character(sizes_flat$value))
sizes_flat$n_roi = as.numeric(as.character(sizes_flat$n_roi))
ggplot(sizes_flat, aes(x=n_roi,y=value,color=strategy)) +
  geom_line(size=2,alpha=0.25,stat="identity") +
  ylab("total size of parcellation (voxels)") + 
  xlab("number of regions (coverage)") + 
  facet_wrap(~strategy)
ggsave(paste(savedir,"/tau_sigdiff_masksizes_fdr05.png",sep=""))

# Let's express the number of voxels as a percentage of the brain mask, so we can apply globally to different sizes / sampled images
brainmask_size = unique(sizes_flat$value[sizes_flat$strategy=="BM"])
sizes_flat$value = as.numeric(as.character(sizes_flat$value/brainmask_size))
ggplot(sizes_flat, aes(x=n_roi,y=value,color=strategy)) +
  geom_line(size=2,alpha=0.25,stat="identity") +
  ylab("total size of parcellation (% of brain mask)") + 
  xlab("number of regions (coverage)") + 
  facet_wrap(~strategy)
ggsave(paste(savedir,"/tau_sigdiff_maskpersizes_fdr05.png",sep=""))

# Read in the threshold report for all neurovault images!
brainmask_size = unique(sizes_flat$value[sizes_flat$strategy=="BM"])
subset = sizes_flat[sizes_flat$strategy=="PD",]
subset$value = as.numeric(as.character(subset$value/brainmask_size))
subset = subset[,-c(1,2,4)]
# Get the mean for each size
ss = unique(subset$n_roi)
df=c()
for (s in ss){
  sub = subset[which(subset$n_roi==s),]
  df=rbind(df,cbind(s,mean(sub$value)))
}

# Here is where we can load pvalues associated with each threshold (n_roi) in "testadjust" variable
load(paste(datadir,"/wilcox_ind_regional_tests_all_fdr.Rda",sep=""))
adjusted_bm = testadjust[which(testadjust$strategy=="bm"),] 
adjusted_pd = testadjust[which(testadjust$strategy=="pd"),]
adjusted_pi = testadjust[which(testadjust$strategy=="pi"),]

# Here we will save data frames of single fdr for each based on the real data
bm_q=c()
pd_q=c()
pi_q=c()

thresh_report = read.csv("/home/vanessa/Desktop/Z/threshold_report.tsv",sep="\t")
for (it in 1:1000){
  for (t in 1:nrow(thresh_report)){
    percent_nonzero = thresh_report$percent_nonzero[t]
    # Find the most similar sized value (smallest squared distance)
    tmp = abs(df[,2] - percent_nonzero)
    idx = which(tmp==min(tmp))[1]
    nroi = df[idx,1]
    # Now we sample a pvalue from the distribution of 144 images with masking type and n_roi
    bm_q=c(bm_q,sample(adjusted_bm$p.value[which(adjusted_bm$n_roi==nroi)],1))  
    pd_q=c(pd_q,sample(adjusted_pd$p.value[which(adjusted_pd$n_roi==nroi)],1))  
    pi_q=c(pi_q,sample(adjusted_pi$p.value[which(adjusted_pi$n_roi==nroi)],1))  
  }
}

# Calculate percentage of times:
bm_perc_diff = length(which(bm_q<=0.05)) / length(bm_q)
pd_perc_diff = length(which(pd_q<=0.05)) / length(pd_q)
pi_perc_diff = length(which(pi_q<=0.05)) / length(pi_q)

# What if we remove the bad ones?
good_maps = thresh_report[thresh_report$percent_nonzero >= 0.5,]
# Here we will save data frames of single fdr for each based on the real data
bm_q=c()
pd_q=c()
pi_q=c()

for (it in 1:1000){
  for (t in 1:nrow(good_maps)){
    percent_nonzero = good_maps$percent_nonzero[t]
    # Find the most similar sized value (smallest squared distance)
    tmp = abs(df[,2] - percent_nonzero)
    idx = which(tmp==min(tmp))[1]
    nroi = df[idx,1]
    # Now we sample a pvalue from the distribution of 144 images with masking type and n_roi
    bm_q=c(bm_q,sample(adjusted_bm$p.value[which(adjusted_bm$n_roi==nroi)],1))  
    pd_q=c(pd_q,sample(adjusted_pd$p.value[which(adjusted_pd$n_roi==nroi)],1))  
    pi_q=c(pi_q,sample(adjusted_pi$p.value[which(adjusted_pi$n_roi==nroi)],1))  
  }
}
length(which(bm_q<=0.05)) / length(bm_q)
length(which(pd_q<=0.05)) / length(pd_q)
length(which(pi_q<=0.05)) / length(pi_q)

# Now let's look at thresholding for ALL neurovault images!
thresh_report = read.csv("/home/vanessa/Desktop/NV/774_threshold_report.tsv",sep="\t")
mean(thresh_report$percent_nonzero)
sd(thresh_report$percent_nonzero)

bm_q=c()
pd_q=c()
pi_q=c()

for (it in 1:1000){
  for (t in 1:nrow(good_maps)){
    percent_nonzero = good_maps$percent_nonzero[t]
    # Find the most similar sized value (smallest squared distance)
    tmp = abs(df[,2] - percent_nonzero)
    idx = which(tmp==min(tmp))[1]
    nroi = df[idx,1]
    # Now we sample a pvalue from the distribution of 144 images with masking type and n_roi
    bm_q=c(bm_q,sample(adjusted_bm$p.value[which(adjusted_bm$n_roi==nroi)],1))  
    pd_q=c(pd_q,sample(adjusted_pd$p.value[which(adjusted_pd$n_roi==nroi)],1))  
    pi_q=c(pi_q,sample(adjusted_pi$p.value[which(adjusted_pi$n_roi==nroi)],1))  
  }
}
length(which(bm_q<=0.05)) / length(bm_q)
length(which(pd_q<=0.05)) / length(pd_q)
length(which(pi_q<=0.05)) / length(pi_q)


