# The R Markdown was buggy!

# 1. Prepare data matrices
source("../experiments/experiment2_functions.R")
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
    ggsave(paste(savedir,"/pearson_density_",thresh,".png",sep=""))
    plot_result(results[[paste(abs_v,thresh,sep="_")]],thresh,"violin")
    ggsave(paste(savedir,"/pearson_violin_",thresh,".png",sep=""))
    plot_result(results[[paste(abs_v,thresh,sep="_")]],thresh,"boxplot")
    ggsave(paste(savedir,"/pearson_boxplot_",thresh,".png",sep=""))
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
