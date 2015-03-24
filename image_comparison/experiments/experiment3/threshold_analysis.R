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
  image_id = image_ids[1]
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
per_sigdiff = per_sigdiff[ order(per_sigdiff[,4]), ]
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
  df = get_direction_result(results,"posneg")
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
fdr = p.adjust(as.numeric(wilcox_thresh[,2]),method="fdr")
wilcox_thresh = cbind(wilcox_thresh,fdr)
tmp = as.data.frame(wilcox_thresh,stringsAsFactors=FALSE)
tmp$p.value = as.numeric(tmp$p.value)
tmp$ci.low = as.numeric(tmp$ci.low)
tmp$ci.high = as.numeric(tmp$ci.high)
tmp$fdr = as.numeric(tmp$fdr)
write.table(tmp,file=paste(datadir,"/wilcox_tests_distributions_corrected_posneg.tsv",sep=""),sep="\t",row.names=FALSE)
save(tmp,file=paste(datadir,"/wilcox_tests_distributions_corrected_posneg.Rda",sep=""))
tmp = melt(tmp,id.vars=c("strategy","thresh","ci.low","ci.high","fdr"))

# Significantly different means
tmp[which(tmp$fdr<0.01),]


# ? Quantify: How often are rankings significantly different? ##############################################################################

# THIS PART IS NOT DONE YET 3/23/2015
# PART II: Assessing significant differences in ORDERING / RANKING of similar images
# We will save our results in a data frame for each of the contrast, task, group, and task_contrast gold standards
# with format imageid strategy rho rho-pvalue tau tau-pvalue
contrast_df = c()
task_df = c()
group_df = c()
task_group_df = c()

directions = c("posneg","pos")

for (abs_v in abs_values){
  df = c()
  for (thresh in thresholds){  
    # For each *gsr* and each masking strategy *pd*,*pi*,and *bm*
    for (g in 1:nrow(pearsons_gs)){
      
      gs = make_gold_standard_ranking(image_id,image_ids)
      
      #TODO here:
      # which to use, tau or spearman?
      # For each gs ordering, sort data based on that
      # then calculate tau or rho
      # Save in different dataframes
      
      rowname = rownames(pearsons_gs)[g]
      gsr = seq(1,ncol(pearsons_gs))
      # Names of gsr correspond with the image id associated with the ranking
      names(gsr) = c(rowname,gsub("X","",(names(sort(abs(pearsons_gs[g,]),decreasing=TRUE)))))
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
   ylab("% significantly different rankings") + 
   xlab("threshold +/-") + 
  facet_wrap(~STRATEGY) +
  theme(text = element_text(size=20))
ggsave(paste(savedir,"/tau_sigdiff_posonly_fdr05.png",sep=""))

# Finally, we want to see how the mask sizes change
# First we will look at the means
# DID NOT USE THIS IN PAPER ----------------
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
save(ncomps_all,file=paste(datadir,"/mean_sizes_comparison_counts.Rda",sep=""))
# END DID NOT USE THIS IN PAPER ----------------

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
  title("Mask Sizes at Different Thresholds +") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +") +
  ylab("Mask Size (voxels)") 
ggsave(paste(savedir,"/mask_sizes_atthresh_posonly.png",sep=""))

# Now plot percentage of voxels
total_size = max(df$mask.size)
subset$mask.size = subset$mask.size / total_size
tmp = ddply(subset, c("thresh","strategy"), summarise, mask.size.mean=mean(mask.size),mask.size.min = min(mask.size), mask.size.max=max(mask.size))
ggplot(tmp, aes(x=thresh, y=mask.size.mean, ymin=mask.size.min, ymax=mask.size.max,colour=strategy,fill=strategy,group=strategy)) + 
  geom_line(size=1) + 
  title("Mask Sizes at Different Thresholds +/-") +
  geom_ribbon(alpha=0.25,linetype=0) +
  xlab("Threshold +") +
  ylab("% Non-missing") 
ggsave(paste(savedir,"/mask_percsizes_atthresh_posneg.png",sep=""))
