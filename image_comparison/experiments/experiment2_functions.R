library(reshape2)
library(pheatmap)
library(ggplot2)

# Write a function to read in the files for a particular threshold.
read_inputs = function(indir,thresh,abs_value) {
  
  # Similarity Scores
  pearson_pd = read.csv(paste(indir,"/144_masking_pd_",thresh,"_",abs_value,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pearson_pi = read.csv(paste(indir,"/144_masking_pi_",thresh,"_",abs_value,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pearson_bm = read.csv(paste(indir,"/144_masking_bm_",thresh,"_",abs_value,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pearson_gs = read.csv(paste(indir,"/144_masking_gs_",thresh,"_",abs_value,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  
  # Mask Size Differences
  pd_vs_pi = read.csv(paste(indir,"/144_pd_vs_pi_sizediff_",thresh,"_",abs_value,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  pd_vs_bm = read.csv(paste(indir,"/144_pd_vs_bm_sizediff_",thresh,"_",abs_value,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)
  bm_vs_pi = read.csv(paste(indir,"/144_pi_vs_bm_sizediff_",thresh,"_",abs_value,".tsv",sep=""),sep="\t",stringsAsFactors=FALSE)  
  
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

# A function to return a flat data frame for a single image against all others, for all thresholds
get_single_result = function(image_id,results,pos_only="False") {
  gs_base = results[["gs"]][which(rownames(results[["gs"]])==image_id),]
  # We need to define the gs for each threshold
  gs = c()
  thresholds = unique(results[["pd"]]$thresh)
  for (thresh in thresholds){
    tmp = cbind(gs_base,thresh)
    names(tmp)[ncol(gs_base)+1] = c("thresh")    
    gs = rbind(gs,tmp)
  }
  gs = cbind(gs,rep("GS",nrow(gs))); colnames(gs)[ncol(gs)] = "strategy"
  pd = results[["pd"]][which(results[["pd"]][,"X"]==image_id),-which(colnames(results[["pd"]])=="X")]
  pi = results[["pi"]][which(results[["pi"]][,"X"]==image_id),-which(colnames(results[["pi"]])=="X")]
  bm = results[["bm"]][which(results[["bm"]][,"X"]==image_id),-which(colnames(results[["bm"]])=="X")]  
  # Get rid of pos_only column
  pd = pd[which(pd$pos_only==pos_only),-which(colnames(pd)=="pos_only")]
  pi = pi[which(pi$pos_only==pos_only),-which(colnames(pi)=="pos_only")]
  bm = bm[which(bm$pos_only==pos_only),-which(colnames(bm)=="pos_only")]
  # The column corresponding to the image id will be all NA- eliminate
  gs = gs[-which(is.na(gs[,1])),]
  pd = pd[,-which(colnames(pd)==paste("X",image_id,sep=""))]
  pi = pi[,-which(colnames(pi)==paste("X",image_id,sep=""))]
  bm = bm[,-which(colnames(bm)==paste("X",image_id,sep=""))]
  # Add on strategy
  pd = cbind(pd,rep("PD",nrow(pd))); colnames(pd)[ncol(pd)] = "strategy"
  pi = cbind(pi,rep("PI",nrow(pi))); colnames(pi)[ncol(pi)] = "strategy"
  bm = cbind(bm,rep("BM",nrow(bm))); colnames(bm)[ncol(bm)] = "strategy"
  all = rbind(pd,pi,bm)
  all_flat = melt(all,id.vars=c("thresh","strategy"),factorsAsStrings=TRUE)
  gs = as.data.frame(gs,stringsAsFactors=FALSE)
  gs$gs_base = as.numeric(gs$gs_base)
  gs_flat = melt(gs,id.vars=c("thresh","strategy"),factorsAsStrings=TRUE)
  gs_flat$thresh = as.numeric(gs_flat$thresh)
  all_flat$strategy = as.character(all_flat$strategy)
  gs_flat$strategy = as.character(gs_flat$strategy)
  all_flat$thresh = as.numeric(all_flat$thresh)
  all = rbind(gs_flat,all_flat)
  return(all)
}


plot_result = function(res,thresh,direction="False",plot_type="density") {
  
  pd = res$pd[res$pd$pos_only==direction,-which(colnames(res$pd)=="pos_only")]
  pi = res$pi[res$pd$pos_only==direction,-which(colnames(res$pi)=="pos_only")]
  bm = res$bm[res$pd$pos_only==direction,-which(colnames(res$bm)=="pos_only")]
  gs = as.vector(as.matrix(res$gs))
  
  # Filter down to threshold of interest
  pd = as.vector(as.matrix(pd[pd$thresh==thresh,-which(colnames(pd)=="thresh")]))
  pi = as.vector(as.matrix(pi[pi$thresh==thresh,-which(colnames(pi)=="thresh")]))
  bm = as.vector(as.matrix(bm[bm$thresh==thresh,-which(colnames(bm)=="thresh")]))
  
  # Remove nans (comparisons to self)
  pd = pd[-which(is.na(pd))]
  gs = gs[-which(is.na(gs))]
  pi = pi[-which(is.na(pi))]
  bm = bm[-which(is.na(bm))]
  
  pd = flatten(pd,"PD")
  gs = flatten(gs,"GS")
  pi = flatten(pi,"PI")
  bm = flatten(bm,"BM")

  flat = as.data.frame(rbind(pd,gs,pi,bm),stringsAsFactors=FALSE)
  flat$value = as.numeric(flat$value)
 
  if (plot_type=="density"){
    return(ggplot(flat,aes(x=value, fill=variable)) + geom_density(alpha=0.25) + ylab("Density") + xlab(paste("Pearsons R, Threshold",thresh)) + ylim(0,9))
  } else {
    return(ggplot(flat,aes(x=variable, y=value,fill=variable)) + geom_boxplot(alpha=0.25) + ylab("Pearsons R") + xlab(paste("Threshold",thresh)) + ylim(-1,1))    
  }
}

# Function to flatten so we can have different numbers of elements
flatten = function(values,label){
  df = cbind(variable=rep(label,length(values)),value=values)
  return(df)
}

flatten_data = function(df,pos_only="False",label){
  subset = df[which(df$pos_only==pos_only),]
  ts = unique(subset$thresh)
  ALL = c()
  for (t in ts){
    sub = subset[which(subset$thresh==t),-which(colnames(subset)%in%c("pos_only","thresh","X"))]
    vector = as.vector(as.matrix(sub))
    vector = vector[-which(is.na(vector))]
    ALL = rbind(ALL,cbind(vector,rep(t,length(vector))))    
  }
  ALL = cbind(ALL,rep(label,nrow(ALL)))
  colnames(ALL) = c("pearsonr","thresh","masking")  
  return(ALL)
}

# Function to calculate confidence intervals
get_ci = function(dat,direction="upper"){
  error = qnorm(0.975)*sd(dat)/sqrt(length(dat))
  if (direction=="upper"){
    return(mean(dat)+error)
  } else {
    return(mean(dat)-error)    
  }
}

# Function to return sorted ordering for a particular threshold and direction
get_similar = function(df,id,thresh,pos_only="False"){
  subset = df[which(df$X==id),]
  subset = subset[which(subset$pos_only==pos_only),-which(colnames(subset)%in%c("X","pos_only"))]  
  subset = subset[which(subset$thresh==thresh),-which(colnames(subset)%in%c("thresh"))]  
  subset = subset[-which(is.na(subset))]
  labels = gsub("X","",names(subset))
  subset = as.numeric(subset)
  names(subset) = labels
  ordered = names(sort(abs(subset),decreasing=TRUE))
  return(ordered)
}

# Function to return intersection of two sets, in the case a high
# threshold makes comparison between two images impossible, we need
# to remove that comparison from the gold standard list
get_sets = function(set1,set2){
  tonix = setdiff(names(set1),names(set2)) 
  if (length(tonix)!=0){
    set1=as.numeric(gsr[-which(names(gsr)%in%tonix)])
    set2=as.numeric(pdr[-which(names(pdr)%in%tonix)])    
  }
  return(list(set1=set1,set2=set2))
}

# Function to plot the change for a single ordering at every threshold
get_similar_single_data = function(df,gs,id,pos_only="False"){
  tmp = as.numeric(gs) # gets rid of names
  thresholds = unique(df$thresh)
  for (thresh in thresholds){
    subset = get_similar(df,id,thresh,pos_only)
    tmp = rbind(tmp,gsr[subset])
  }
  rownames(tmp) = c(0.0,thresholds)
  return(tmp)  
}  

# Plot the image!
plot_single_similar = function(df,gs,id,pos_only="False"){
  data = get_similar_single_data(df,gs,id,pos_only=pos_only)
  #TODO: for now we save static image, but later need to export data to some file... json?
  pheatmap(data)  
}

# This function will return the number of comparisons done (eg, not NA values)
get_number_comparisons = function(df,thresh,pos_only="False"){
  image_ids = unique(df$X)
  nc = c()
  for (id in image_ids){
    subset = df[which(df$X==id),]
    subset = subset[which(subset$thresh==thresh),-which(colnames(subset)=="thresh")]
    subset = subset[which(subset$pos_only==pos_only),-which(colnames(subset)=="pos_only")]  
    nc = c(nc,length(which(!is.na(subset)))-1) # subtract 1 for ID
  }
  return(nc)
}

# This function will return the mask sizes
get_masksize_comparisons = function(df,thresh,pos_only="False"){
  image_ids = unique(df$X)
  nc = c()
  for (id in image_ids){
    subset = df[which(df$X==id),]
    subset = subset[which(subset$thresh==thresh),-which(colnames(subset)=="thresh")]
    subset = subset[which(subset$pos_only==pos_only),-which(colnames(subset)=="pos_only")]  
    subset = subset[!is.na(subset)]
    subset = subset[-1]
    if (!length(subset)==0){
      nc[[id]] = subset      
    }
  }
  return(nc)
}


plot_pval = function(df,thresh,savedir){
  
  # Rho
  df[which(df$strategy=="BM"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="BM"),3]),method="fdr"))
  df[which(df$strategy=="PI"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PI"),3]),method="fdr"))
  df[which(df$strategy=="PD"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PD"),3]),method="fdr"))
  
  rho = data.frame(variable=df$strategy,value=as.numeric(df$rho_pvalue),stringsAsFactors=FALSE)
  rho$value = as.numeric(rho$value)
  #ggplot(rho,aes(x=variable, y=value,fill=variable)) + geom_boxplot(alpha=0.25) + title(paste("Significantly different orderings based on masking, threshold",thresh,sep="")) + ylab("Q-Value") + xlab("RHO: fdr q-value")
  #ggsave(paste(savedir,"/rho_",thresh,".png",sep=""))
  
  # Tau
  df[which(df$strategy=="BM"),5] = p.adjust(as.numeric(as.character(df[which(df$strategy=="BM"),5]),method="fdr"))
  df[which(df$strategy=="PI"),5] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PI"),5]),method="fdr"))
  df[which(df$strategy=="PD"),5] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PD"),5]),method="fdr"))
  
  tau = data.frame(variable=df$strategy,value=as.numeric(df$tau_pvalue),stringsAsFactors=FALSE)
  tau$value = as.numeric(tau$value)
  #ggplot(tau,aes(x=variable, y=value, fill=variable)) + geom_boxplot(alpha=0.25) + title(paste("Significantly different orderings based on masking, threshold",thresh,sep="")) + ylab("Density") + xlab("TAU: fdr q-value")
  #ggsave(paste(savedir,"/tau_",thresh,".png",sep=""))
  
  generate_count_table = function(dat,label,fdr_thresh){
    tabley = table(dat$value[which(dat$variable==label)] <= fdr_thresh)
    if (!("TRUE" %in% names(tabley))) {
      tabley["TRUE"] = 0
    }
    if (!("FALSE" %in% names(tabley))) {
      tabley["FALSE"] = 0
    }    
    return(tabley)
  }
  
  counts = rbind(generate_count_table(rho,"PD",0.001),
                 generate_count_table(tau,"PD",0.001),
                 generate_count_table(rho,"PI",0.001),
                 generate_count_table(tau,"PI",0.001),
                 generate_count_table(rho,"BM",0.001),
                 generate_count_table(tau,"BM",0.001))
                 
  counts = cbind(c("PD_RHO","PD_TAU","PI_RHO","PI_TAU","BM_RHO","BM_TAU"),counts)
  counts = as.data.frame(counts,stringsAsFactors=FALSE)
  colnames(counts) = c("STRATEGY","SIG_DIFF","NOT_SIG")
  counts[,"perc_diff"] = as.numeric(counts$SIG_DIFF) / (as.numeric(counts$SIG_DIFF) + as.numeric(counts$NOT_SIG))
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
