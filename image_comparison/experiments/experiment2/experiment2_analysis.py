# Redoing this analysis in python because R is terrible.

from scipy import stats
import pylab as plt
import numpy as np
import seaborn as sns

# Generate pearson correlated data with approximately cor(X, Y) = r
def generate_vectors(length = 200, r = 0.5):
    data = np.random.multivariate_normal([0, 0], [[1, r], [r, 1]], size=length)
    X, Y = data[:,0], data[:,1]

    # That's it! Now let's take a look at the actual correlation:
    import scipy.stats as stats
    return (X,Y)

    print 'r=', stats.pearsonr(X, Y)[0]
In [2]:
def test_correlations(vector_length = 200, r = 0.5, mask="first_n", mask_len = 20):
    #find correlation between full vectors
    X, Y = generate_vectors(vector_length, r)
    full_sample_r = stats.pearsonr(X, Y)[0]
    
    
    if mask == "small_values":
        #mask out small values (between -1 and 1 sigma)
        mask = (X > 2) | (X < -2)
    elif mask == "first_n":
        #mask out first n elements
        mask = np.zeros(X.shape[0])
        mask[-mask_len:] = 1
        mask = (mask == 1)
        assert mask.sum() == mask_len
    
    #find correlation only withing the mask
    masked_sample_r = stats.pearsonr(X[mask], Y[mask])[0]
    masked_value_error = full_sample_r - masked_sample_r
    masked_perc_error = masked_value_error/full_sample_r*100
    
    #zero voxels of one map and recalculate the correlation
    zeroed_X = X.copy()
    zeroed_X[mask == False] = 0
    #zeroed_Y = Y.copy()
    #zeroed_Y[mask == False] = 0
    zeroed_sample_r = stats.pearsonr(zeroed_X, Y)[0]
    zeroed_value_error = full_sample_r - zeroed_sample_r
    zeroed_perc_error = zeroed_value_error/full_sample_r*100
    return masked_perc_error, zeroed_perc_error, full_sample_r, masked_sample_r, zeroed_sample_r
In [3]:
test_correlations()
Out[3]:
(2.9589286656536831,
 59.253174080905922,
 0.46383203792100114,
 0.45010757879047097,
 0.1889968330486568)
Removing NAs vs zeroeing NAs when NAs happen at random
In [13]:
masked_errors = []
zeroed_errors = []
full_rs = []
masked_rs = []
zeroed_rs = []
for _ in range(200):
    masked_error, zeroed_error, full_r, masked_r, zeroed_r = test_correlations(mask="first_n")
    masked_errors.append(masked_error)
    zeroed_errors.append(zeroed_error)
    full_rs.append(full_r)
    masked_rs.append(masked_r)
    zeroed_rs.append(zeroed_r)
In [14]:
sns.violinplot([full_rs, masked_rs, zeroed_rs], names=["full correlation", "masked NAs", "zeroed NAs"], color="coral");
Results are not surprising. When we choose a random supsample we get noisier estimate of the correlation hence "masked NA" distribution is wider. When we impute zeros at random the correlation will decreaese (akin to systematic error).

Removing NAs vs zeroeing NAs when NAs are happen to be values between -2sigma and 2sigma
In [15]:
masked_errors = []
zeroed_errors = []
full_rs = []
masked_rs = []
zeroed_rs = []
for _ in range(200):
    masked_error, zeroed_error, full_r, masked_r, zeroed_r = test_correlations(mask="small_values")
    masked_errors.append(masked_error)
    zeroed_errors.append(zeroed_error)
    full_rs.append(full_r)
    masked_rs.append(masked_r)
    zeroed_rs.append(zeroed_r)
In [16]:
sns.violinplot([full_rs, masked_rs, zeroed_rs], names=["full correlation", "masked NAs", "zeroed NAs"], color="coral");
However if the NAs happen to be small values the using pairwise deletion (masking out NAs) inflates the correlation. Zeroing NAs stil underestimates the correlation.

In [ ]:
 
Back to top
This web site does not host notebooks, it only renders notebooks available on other websites.
Delivered by Fastly, Rendered by Rackspace
nbviewer GitHub repository.
nbviewer version: 852295b
IPython version: 2.4.1-maint
Rendered 5 minutes ago

# 1. Prepare data matrices
source("../experiments/experiment2_functions.R")
suppressPackageStartupMessages=TRUE
indir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores"  
setwd(indir)

thresholds = c(1.96,2.58)
abs_values = c("pos","pos_neg")

results = list()
for (thresh in thresholds){
  for (abs_v in abs_values){
    indir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/analysis/masking_scores"  
    res = read_inputs(indir,thresh,abs_v)
    for (r in names(res)){
      res[[r]] = format_df(res[[r]])
    }
    results[[paste(abs_v,thresh,sep="_")]] = res
  }
}
#save(results,file=paste(indir,"/masking_matrices_all.Rda",sep=""))

# 2. Make plots
savedir = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/img/all"

for (thresh in thresholds){
  for (abs_v in abs_values){
    plot_result(results[[paste(abs_v,thresh,sep="_")]],thresh)
    ggsave(paste(savedir,"/pearson_density_",thresh,"_",abs_v,".png",sep=""))
    plot_result(results[[paste(abs_v,thresh,sep="_")]],thresh,"violin")
    ggsave(paste(savedir,"/pearson_violin_",thresh,"_",abs_v,".png",sep=""))
    plot_result(results[[paste(abs_v,thresh,sep="_")]],thresh,"boxplot")
    ggsave(paste(savedir,"/pearson_boxplot_",thresh"_",abs_v,,".png",sep=""))
  }
}

# We will save our results in a data frame with format
# imageid strategy rho rho-pvalue tau tau-pvalue
all_df = c()

for (thresh in thresholds){
  for (abs_v in abs_values){
    df = c()
    res = results[[paste(abs_v,thresh,sep="_")]]
    # For each *gsr* and each masking strategy *pd*,*pi*,and *bm*
    for (g in 1:nrow(res$pearson_gs)){
      rowname = rownames(res$pearson_gs)[g]
      gsr = as.numeric(names(sort(abs(res$pearson_gs[g,]),decreasing=TRUE)))
      pdr = as.numeric(names(sort(abs(res$pearson_pd[g,]),decreasing=TRUE)))
      pir = as.numeric(names(sort(abs(res$pearson_pi[g,]),decreasing=TRUE)))
      bmr = as.numeric(names(sort(abs(res$pearson_bm[g,]),decreasing=TRUE)))
      # PAIRWISE DELETION
      rho = cor.test(gsr, pdr, method = c("spearman"), conf.level = 0.95)
      tau = cor.test(gsr, pdr, method = c("kendall"), conf.level = 0.95)
      row = cbind(rowname,"PD",rho$p.value,rho$estimate,tau$p.value,tau$estimate)
      df = rbind(df,row)
      # PAIRWISE INCLUSION
      rho = cor.test(gsr, pir, method = c("spearman"), conf.level = 0.95)
      tau = cor.test(gsr, pir, method = c("kendall"), conf.level = 0.95)
      row = cbind(rowname,"PI",rho$p.value,rho$estimate,tau$p.value,tau$estimate)
      df = rbind(df,row)
      # BRAIN MASK
      rho = cor.test(gsr, bmr, method = c("spearman"), conf.level = 0.95)
      tau = cor.test(gsr, bmr, method = c("kendall"), conf.level = 0.95)
      row = cbind(rowname,"BM",rho$p.value,rho$estimate,tau$p.value,tau$estimate)
      df = rbind(df,row)
    }
    colnames(df) = c("imageid","strategy","rho_pvalue","rho","tau_pvalue","tau")
    rownames(df) = seq(1,nrow(df))
    df = as.data.frame(df,stringsAsFactors=FALSE)
    all_df[[paste(abs_v,thresh,sep="_")]] = df
  }
}
# Do the plotting
counts = list()
for (thresh in thresholds){
  for (abs_v in abs_values){
    df = all_df[[paste(abs_v,thresh,sep="_")]]
    c = plot_pval(df,thresh,savedir)
    c = cbind(c,paste(c$STRATEGY,rep(abs_v,nrow(c)),rep(thresh,nrow(c)),sep="_"))
    counts = rbind(counts,c)
  }
}
colnames(counts)[4] = "perc_diff"
colnames(counts)[5] = "threshold"

# Which one are rho vs tau?
rho = counts[grep("RHO",counts$STRATEGY),]
tau = counts[grep("TAU",counts$STRATEGY),]
rho = rho[order(-rho$perc_diff),]
tau = tau[order(-tau$perc_diff),]

# Save to table
write.table(counts,file=paste(indir,"/result_table_all.tsv",sep=""),sep="\t",row.names=FALSE)

# Show the tables
rho
tau

# Finally plot the percentages
rho$threshold = gsub("_RHO","",rho$threshold)
ggplot(rho, aes(x=threshold,y=perc_diff,fill=STRATEGY)) + geom_bar(alpha=0.25,stat="identity") + title("Percentage Significantly Different Images, RHO") + ylab("percent") + xlab("masking strategy") + ylim(0,20)
ggsave(paste(savedir,"/rho_sigdiff.png",sep=""))

# Finally plot the percentages
tau$threshold = gsub("_TAU","",tau$threshold)
ggplot(tau, aes(x=threshold,y=perc_diff,fill=STRATEGY)) + geom_bar(alpha=0.25,stat="identity") + title("Percentage Significantly Different Images, RHO") + ylab("percent") + xlab("masking strategy") + ylim(0,20)
ggsave(paste(savedir,"/tau_sigdiff.png",sep=""))
