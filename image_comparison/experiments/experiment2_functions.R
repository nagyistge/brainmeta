library(reshape2)
library(ggplot2)

# Write a function to read in the files for a particular threshold.
read_inputs = function(indir,thresh) {
  
  # Similarity Scores
  pearson_pd = read.csv(paste(indir,"/144_masking_pd_",thresh,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pearson_pi = read.csv(paste(indir,"/144_masking_pi_",thresh,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pearson_bm = read.csv(paste(indir,"/144_masking_bm_",thresh,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pearson_gs = read.csv(paste(indir,"/144_masking_gs_",thresh,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  
  # Mask Size Differences
  pd_vs_pi = read.csv(paste(indir,"/144_pd_vs_pi_sizediff_",thresh,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pd_vs_bm = read.csv(paste(indir,"/144_pd_vs_bm_sizediff_",thresh,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  bm_vs_pi = read.csv(paste(indir,"/144_pi_vs_bm_sizediff_",thresh,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)  
  
  result = list(bm_vs_pi=bm_vs_pi,pd_vs_pi=pd_vs_pi,pd_vs_bm=pd_vs_bm,pearson_gs=pearson_gs,pearson_bm=pearson_bm,pearson_pi=pearson_pi,pearson_pd=pearson_pd)
  return(result)
}

# Write a function to format the data frames
format_df = function(df){
  colnames(df)[1] = "ID"
  colnames(df) = gsub("X","",colnames(df))
  rownames(df) = df[,1]
  df = df[,-1]
  return(df)
}

plot_result = function(res,thresh,type="density") {
  
  gs_vector = as.vector(as.matrix(res$pearson_gs))
  pd_vector = as.vector(as.matrix(res$pearson_pd))
  pi_vector = as.vector(as.matrix(res$pearson_pi))
  bm_vector = as.vector(as.matrix(res$pearson_bm))
  
  pearsons_df = data.frame(GS=gs_vector,PD=pd_vector,PI=pi_vector,BM=bm_vector)
  pearsons_flat =  melt(pearsons_df)
  
  # Mask Sizes
  pd_vs_bm_vector = as.vector(as.matrix(res$pd_vs_bm))
  pd_vs_pi_vector = as.vector(as.matrix(res$pd_vs_pi))
  bm_vs_pi_vector = as.vector(as.matrix(res$bm_vs_pi))
  
  mask_df = data.frame(pd_vs_bm_vector,pd_vs_pi_vector,bm_vs_pi_vector)
  mask_flat = melt(mask_df)
  
  if (type == "boxplot"){
    return(ggplot(pearsons_flat, aes(variable,value,fill=variable)) + geom_boxplot(alpha=0.25) + ylab("Pearsons R") + title(paste("Pearson Correlations for Different Masking Strategies, Threshold",thresh,sep="")))  
  } else if (type == "violin"){
    return(ggplot(pearsons_flat, aes(variable,value,fill=variable)) + geom_violin(alpha=0.25) + ylab("Pearsons R") + title(paste("Pearson Correlations for Different Masking Strategies, Threshold",thresh,sep=""))) 
  } else {
    return(ggplot(pearsons_flat,aes(x=value, fill=variable)) + geom_density(alpha=0.25) + title("Comparing Density Plots of Pearson Correlations") + ylab("Density") + xlab("Pearsons R"))
  }
}

plot_pval = function(df,thresh,savedir){
  
  # Rho
  df[which(df$strategy=="BM"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="BM"),3]),method="fdr"))
  df[which(df$strategy=="PI"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PI"),3]),method="fdr"))
  df[which(df$strategy=="PD"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PD"),3]),method="fdr"))
  
  rho = data.frame(variable=df$strategy,value=as.numeric(df$rho_pvalue),stringsAsFactors=FALSE)
  ggplot(rho[-which(rho$value>.25),],aes(x=value, fill=variable)) + geom_density(alpha=0.25) + title(paste("Significantly different orderings based on masking, threshold",thresh,sep="")) + ylab("Density") + xlab("RHO: fdr q-value")
  ggsave(paste(savedir,"/rho_",thresh,".png",sep=""))
  
  # Tau
  df[which(df$strategy=="BM"),5] = p.adjust(as.numeric(as.character(df[which(df$strategy=="BM"),5]),method="fdr"))
  df[which(df$strategy=="PI"),5] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PI"),5]),method="fdr"))
  df[which(df$strategy=="PD"),5] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PD"),5]),method="fdr"))
  
  tau = data.frame(variable=df$strategy,value=as.numeric(df$tau_pvalue),stringsAsFactors=FALSE)
  # Get rid of anything over -.25 too much clutter
  ggplot(tau[-which(tau$value>.25),],aes(x=value, fill=variable)) + geom_density(alpha=0.25) + title(paste("Significantly different orderings based on masking, threshold",thresh,sep="")) + ylab("Density") + xlab("TAU: fdr q-value")
  ggsave(paste(savedir,"/tau_",thresh,".png",sep=""))
  
  counts = rbind(table(rho$value[which(rho$variable=="PD")] <= 0.05),
                 table(tau$value[which(rho$variable=="PD")] <= 0.05),
                 table(rho$value[which(rho$variable=="PI")] <= 0.05),
                 table(tau$value[which(rho$variable=="PI")] <= 0.05),
                 table(rho$value[which(rho$variable=="BM")] <= 0.05),
                 table(tau$value[which(rho$variable=="BM")] <= 0.05))
  counts = cbind(c("PD_RHO","PD_TAU","PI_RHO","PI_TAU","BM_RHO","BM_TAU"),counts)
  colnames(counts) = c("STRATEGY","NOT_SIG","SIG_DIFF")
  counts = as.data.frame(counts,stringsAsFactors=FALSE)
  counts$SIG_DIFF = as.numeric(counts$SIG_DIFF)
  counts$NOT_SIG = as.numeric(counts$NOT_SIG)
  counts = cbind(counts,counts$SIG_DIFF / (counts$SIG_DIFF + counts$NOT_SIG))
  colnames(counts)[3] = "perc_diff"
  return(counts)
}

# Write a function to format the data frames
format_df = function(df){
  colnames(df)[1] = "ID"
  colnames(df) = gsub("X","",colnames(df))
  rownames(df) = df[,1]
  df = df[,-1]
  return(df)
}
