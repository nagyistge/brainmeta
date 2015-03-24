library(plyr)
library(reshape2)
library(pheatmap)
library(gridExtra)
library(ggplot2)

# FILTER FUNCTIONS ######################################################

filter_columns = function(tmp,extra_cols){
  column_names = gsub("GRP[0-9]+_","",colnames(tmp)[1:(ncol(tmp)-extra_cols)])
  keepers = get_keepers()
  # First we are filtering columns
  filtered = tmp[,which(column_names %in% keepers)]
  filtered$direction = tmp$direction
  filtered$thresh = tmp$thresh
  filtered$UID = tmp$X
  return(filtered)
}

# Function to filter data down to non-redundant contrasts
filter_data = function(input_file,output_file,extra_cols){
  tmp = read.csv(input_file,sep="\t")
  filtered = filter_columns(tmp,extra_cols)
  # Now we filter rows
  row_names = gsub("GRP[0-9]+_","",tmp$X)
  filtered = filtered[which(row_names %in% keepers),]
  write.table(filtered,file=output_file,row.names=TRUE,col.names=TRUE)
}


# Get list to keep
get_keepers = function() {
  return(c("TASK06_CON01","TASK06_CON02","TASK06_CON06","TASK01_CON07","TASK01_CON08","TASK01_CON09",
         "TASK05_CON13","TASK05_CON14","TASK05_CON15","TASK03_CON19","TASK03_CON20","TASK03_CON22",
         "TASK02_CON25","TASK02_CON26","TASK02_CON27","TASK07_CON31","TASK07_CON32","TASK07_CON33",
         "TASK07_CON34","TASK07_CON35","TASK07_CON36","TASK07_CON37","TASK07_CON38","TASK07_CON39",
         "TASK07_CON40","TASK07_CON41","TASK07_CON45","TASK07_CON46","TASK07_CON47","TASK07_CON48",
         "TASK07_CON49","TASK07_CON50","TASK07_CON51","TASK07_CON52","TASK04_CON61","TASK04_CON62",
         "TASK04_CON63","TASK04_CON64","TASK04_CON65","TASK04_CON66","TASK04_CON67","TASK04_CON68",
         "TASK04_CON69","TASK04_CON70","TASK04_CON71","TASK04_CON72","TASK04_CON73"))
}

# A function to return a flat data frame for a single image against all others, for all thresholds
get_single_result = function(image_id,results,direction="posneg") {
  
  pdp = results[["pdp"]][which(results[["pdp"]][,"UID"]==image_id),-which(colnames(results[["pdp"]])=="UID")]
  pip = results[["pip"]][which(results[["pip"]][,"UID"]==image_id),-which(colnames(results[["pip"]])=="UID")]
  bmp = results[["bmp"]][which(results[["bmp"]][,"UID"]==image_id),-which(colnames(results[["bmp"]])=="UID")]  
  pds = results[["pds"]][which(results[["pds"]][,"UID"]==image_id),-which(colnames(results[["pds"]])=="UID")]
  pis = results[["pis"]][which(results[["pis"]][,"UID"]==image_id),-which(colnames(results[["pis"]])=="UID")]
  bms = results[["bms"]][which(results[["bms"]][,"UID"]==image_id),-which(colnames(results[["bms"]])=="UID")]  
  # Get rid of direction
  pdp = pdp[which(pdp$direction==direction),-which(colnames(pdp)=="direction")]
  pip = pip[which(pip$direction==direction),-which(colnames(pip)=="direction")]
  bmp = bmp[which(bmp$direction==direction),-which(colnames(bmp)=="direction")]
  pds = pds[which(pds$direction==direction),-which(colnames(pds)=="direction")]
  pis = pis[which(pis$direction==direction),-which(colnames(pis)=="direction")]
  bms = bms[which(bms$direction==direction),-which(colnames(bms)=="direction")]
  # The column corresponding to the image id will be all 1 eliminate
  pdp = pdp[,-which(colnames(pdp)==image_id)]
  pip = pip[,-which(colnames(pip)==image_id)]
  bmp = bmp[,-which(colnames(bmp)==image_id)]
  pds = pds[,-which(colnames(pds)==image_id)]
  pis = pis[,-which(colnames(pis)==image_id)]
  bms = bms[,-which(colnames(bms)==image_id)]
  # Add on strategy
  pdp = cbind(pdp,rep("intersect.pearson",nrow(pdp))); colnames(pdp)[ncol(pdp)] = "strategy"
  pip = cbind(pip,rep("union.pearson",nrow(pip))); colnames(pip)[ncol(pip)] = "strategy"
  bmp = cbind(bmp,rep("brain.mask.pearson",nrow(bmp))); colnames(bmp)[ncol(bmp)] = "strategy"
  pds = cbind(pds,rep("intersect.spearman",nrow(pds))); colnames(pds)[ncol(pds)] = "strategy"
  pis = cbind(pis,rep("union.spearman",nrow(pis))); colnames(pis)[ncol(pis)] = "strategy"
  bms = cbind(bms,rep("brain.mask.spearman",nrow(bms))); colnames(bms)[ncol(bms)] = "strategy"
  
  all = rbind(pdp,pip,bmp,pds,pis,bms)
  all_flat = melt(all,id.vars=c("thresh","strategy"),factorsAsStrings=FALSE)
  all_flat$strategy = as.character(all_flat$strategy)
  all_flat$thresh = as.numeric(all_flat$thresh)
  return(all_flat)
}


# A function to return a flat data frame for a single threshold, direction, across all people
get_direction_result = function(results,direction="posneg") {
  
  # Get rid of direction
  pdp = results$pdp[which(results$pdp$direction==direction),-which(colnames(results$pdp)%in%c("direction","UID"))]
  pip = results$pip[which(results$pip$direction==direction),-which(colnames(results$pip)%in%c("direction","UID"))]
  bmp = results$bmp[which(results$bmp$direction==direction),-which(colnames(results$bmp)%in%c("direction","UID"))]
  pds = results$pds[which(results$pds$direction==direction),-which(colnames(results$pds)%in%c("direction","UID"))]
  pis = results$pis[which(results$pis$direction==direction),-which(colnames(results$pis)%in%c("direction","UID"))]
  bms = results$bms[which(results$bms$direction==direction),-which(colnames(results$bms)%in%c("direction","UID"))]
  # The column corresponding to the image id will be all 1 eliminate
  # Add on strategy
  pdp = cbind(pdp,rep("intersect.pearson",nrow(pdp))); colnames(pdp)[ncol(pdp)] = "strategy"
  pip = cbind(pip,rep("union.pearson",nrow(pip))); colnames(pip)[ncol(pip)] = "strategy"
  bmp = cbind(bmp,rep("brain.mask.pearson",nrow(bmp))); colnames(bmp)[ncol(bmp)] = "strategy"
  pds = cbind(pds,rep("intersect.spearman",nrow(pds))); colnames(pds)[ncol(pds)] = "strategy"
  pis = cbind(pis,rep("union.spearman",nrow(pis))); colnames(pis)[ncol(pis)] = "strategy"
  bms = cbind(bms,rep("brain.mask.spearman",nrow(bms))); colnames(bms)[ncol(bms)] = "strategy"
  
  all = rbind(pdp,pip,bmp,pds,pis,bms)
  all_flat = melt(all,id.vars=c("thresh","strategy"),factorsAsStrings=FALSE)
  all_flat$strategy = as.character(all_flat$strategy)
  all_flat$thresh = as.numeric(all_flat$thresh)
  all_flat$value = as.numeric(as.character(all_flat$value))
  # We need to remove comparisons of images to themselves at no threshold
  # We are keeping comparisons of images to themselves at different thresholds, of course
  if (any(all_flat$value ==1)){
      all_flat = all_flat[-which(all_flat$value==1),]
  }
  return(all_flat)
}


# PLOTTING FUNCTIONS ####################################################
plot_result = function(res,thresh,direction="posneg",outfile) {

  pdp = res$pdp[res$pdp$direction==direction,-which(colnames(res$pdp)=="direction")]
  pip = res$pip[res$pip$direction==direction,-which(colnames(res$pip)=="direction")]
  bmp = res$bmp[res$bmp$direction==direction,-which(colnames(res$bmp)=="direction")]
  pds = res$pds[res$pds$direction==direction,-which(colnames(res$pds)=="direction")]
  pis = res$pis[res$pis$direction==direction,-which(colnames(res$pis)=="direction")]
  bms = res$bms[res$bms$direction==direction,-which(colnames(res$bms)=="direction")]
  
  # Filter down to threshold of interest
  pdp = as.vector(as.matrix(pdp[pdp$thresh==thresh,-which(colnames(pdp)%in%c("thresh","UID"))]))
  pip = as.vector(as.matrix(pip[pip$thresh==thresh,-which(colnames(pip)%in%c("thresh","UID"))]))
  bmp = as.vector(as.matrix(bmp[bmp$thresh==thresh,-which(colnames(bmp)%in%c("thresh","UID"))]))
  pds = as.vector(as.matrix(pds[pds$thresh==thresh,-which(colnames(pds)%in%c("thresh","UID"))]))
  pis = as.vector(as.matrix(pis[pis$thresh==thresh,-which(colnames(pis)%in%c("thresh","UID"))]))
  bms = as.vector(as.matrix(bms[bms$thresh==thresh,-which(colnames(bms)%in%c("thresh","UID"))]))
    
  pdp = flatten(pdp,"intersect.pearson")
  pip = flatten(pip,"union.pearson")
  bmp = flatten(bmp,"brain.mask.pearson")
  pds = flatten(pds,"intersect.spearman")
  pis = flatten(pis,"union.spearman")
  bms = flatten(bms,"brain.mask.spearman")
  
  flat = as.data.frame(rbind(pdp,pip,bmp,pds,pis,bms),stringsAsFactors=FALSE)
  flat$value = as.numeric(flat$value)

  # Remove values of 1 - comparing an image to itself
  if (any(flat$value==1)) {
    flat = flat[-which(flat$value==1.0),]
  }
  # NAs mean the comparison was not possible - separate the data to plot
  number_na = length(which(is.na(flat$value)))
  if (number_na>0){
    nans = flat[which(is.na(flat$value)),]
    flat = flat[-which(is.na(flat$value)),]
    nans$value = 1
  } else {
    values = c(0,0,0,0,0,0)
    variables = unique(flat$variable)
    # Make something tiny no one will see for empty plot
    nans = data.frame(value = values, variable=variables)
  }
  
  densityplot = ggplot(flat,aes(x=value, fill=variable)) + 
    geom_density(alpha=0.25) + 
    ylab("Density") + 
    xlab(paste("Threshold",thresh)) + 
    ylim(0,4) +
    facet_wrap(~variable) +
    theme(legend.position="none")
  
  nanplot = ggplot(nans,aes(value, fill=variable)) + 
    geom_bar(alpha=0.25,binwidth=10) + 
    xlab("NaN Count") + 
    ylab("") +
    facet_wrap(~variable,nrow=1) +
    ylim(0,5000) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    guides(fill=guide_legend(title=NULL)) +
    scale_x_discrete(breaks=NULL)
  
  g = arrangeGrob(densityplot, nanplot, ncol=2)
  ggsave(file=outfile,g)    
}



# DATA TRANSFORMATION FUNCTIONS ####################################################

# Function to flatten so we can have different numbers of elements
flatten = function(values,label){
  df = cbind(variable=rep(label,length(values)),value=values)
  return(df)
}

flatten_data = function(df,direction="posneg",label){
  subset = df[which(df$direction==direction),]
  ts = unique(subset$thresh)
  ALL = c()
  for (t in ts){
    sub = subset[which(subset$thresh==t),-which(colnames(subset)%in%c("direction","thresh","UID"))]
    vector = as.vector(as.matrix(sub))
    if (length(which(is.na(vector))!=0)){
      vector = vector[-which(is.na(vector))]
    }
    vector = vector[which(vector!=1)]
    ALL = rbind(ALL,cbind(vector,rep(t,length(vector))))    
  }
  ALL = cbind(ALL,rep(label,nrow(ALL)))
  colnames(ALL) = c("score","thresh","strategy")  
  return(ALL)
}

# GOLD STANDARD FUNCTIONS ####################################################

make_gold_standard_labels = function(image_ids){
  # We will make gold standard lists for:
  CONTRAST = gsub("GRP[0-9]+_TASK[0-9]+_","",image_ids)    
  TASK = gsub("GRP[0-9]+_|_CON[0-9]+","",image_ids)        
  GROUP =  gsub("_TASK[0-9]+|_CON[0-9]+","",image_ids)       
  TASK_CONTRAST =  gsub("GRP[0-9]+_","",image_ids)   
  # For each of the above, we will try to predict labels from scores!
  # Which threshold does the best?
  gs = data.frame(UID=image_ids,CONTRAST=CONTRAST,TASK=TASK,GROUP=GROUP,TASK_CONTRAST=TASK_CONTRAST,stringsAsFactors=FALSE)
  return(gs)  
}

# Return gold standard ranked lists for an image (image_id) based on other image_ids
# eg 1 1 1 1 2 2 2 2 3 3 3 3 3
make_gold_standard_ranking = function(image_id,image_ids){
  # We will make gold standard lists for:
  contrast_gs = array(dim=length(image_ids))   
  task_gs = array(dim=length(image_ids))       
  group_gs = array(dim=length(image_ids))      
  task_contrast = array(dim=length(image_ids))
  
  tmp = strsplit(image_id,"_")[[1]]
  group = tmp[1]
  task = tmp[2]
  con = tmp[3]

  contrast_gs[grep(con,image_ids)] = 1
  contrast_gs[is.na(contrast_gs)] = 2
  task_gs[grep(task,image_ids)] = 1
  task_gs[is.na(task_gs)] = 2
  group_gs[grep(group,image_ids)] = 1
  group_gs[is.na(group_gs)] = 2
  task_contrast[grep(group,image_ids)] = 3
  task_contrast[grep(task,image_ids)] = 2
  task_contrast[grep(con,image_ids)] = 1
  task_contrast[is.na(task_contrast)] = 4
  
  gs = data.frame(UID=image_ids,CONTRAST=contrast_gs,TASK=task_gs,GROUP=group_gs,TASK_CONTRAST=task_contrast,stringsAsFactors=FALSE)
  return(gs)  
}


# STATISTICS FUNCTIONS ####################################################

# Function to calculate confidence intervals
get_ci = function(dat,direction="upper"){
  error = qnorm(0.975)*sd(dat)/sqrt(length(dat))
  if (direction=="upper"){
    return(mean(dat)+error)
  } else {
    return(mean(dat)-error)    
  }
}

# Function to do wilcox test to compare distributions at thresh 0.0
# to others
wilcox_test = function(df,thresh,gs,label){
  gs_vec = as.numeric(gs$value[gs$strategy==label])
  df_vec = df[df$strategy==label,]
  df_vec = as.numeric(df_vec$value[df_vec$thresh==thresh])  
  wt_vec = wilcox.test(gs_vec, y = df_vec,alternative="two.sided",conf.int=TRUE,conf.level=0.95)
  return(wt_vec)
}










# Write a function to format the data frames
format_df = function(df){
  colnames(df)[1] = "ID"
  colnames(df) = gsub("X","",colnames(df))
  rownames(df) = df[,1]
  df = df[,-1]
  return(df)
}



# Function to return sorted ordering for a particular threshold and direction
get_similar = function(df,id,thresh,pos_only="False"){
  subset = df[which(df$X==id),]
  subset = subset[which(subset$pos_only==pos_only),-which(colnames(subset)%in%c("X","pos_only"))]  
  subset = subset[which(subset$thresh==thresh),-which(colnames(subset)%in%c("thresh"))]  
  if (length(which(is.na(subset)))!=0) {
    subset = subset[-which(is.na(subset))]
  }
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


plot_pval = function(df,thresh,savedir,fdr){
  
  # A value exactly == 0 actually means no significant difference
  df$rho_pvalue[df$rho_pvalue==0] = 1
  df$tau_pvalue[df$tau_pvalue==0] = 1
  
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
    return(tabley[c("TRUE","FALSE")])
  }
  
  counts = rbind(generate_count_table(rho,"PD",fdr),
                 generate_count_table(tau,"PD",fdr),
                 generate_count_table(rho,"PI",fdr),
                 generate_count_table(tau,"PI",fdr),
                 generate_count_table(rho,"BM",fdr),
                 generate_count_table(tau,"BM",fdr))
                 
  counts = cbind(c("PD_RHO","PD_TAU","PI_RHO","PI_TAU","BM_RHO","BM_TAU"),counts)
  counts = as.data.frame(counts,stringsAsFactors=FALSE)
  colnames(counts) = c("STRATEGY","SIG_DIFF","NOT_SIG")
  counts[,"perc_diff"] = as.numeric(counts$SIG_DIFF) / (as.numeric(counts$SIG_DIFF) + as.numeric(counts$NOT_SIG))
  return(counts)
}

plot_pval_region = function(df,n_roi,savedir,fdr){
  
  # A value exactly == 0 actually means no significant difference
  df$tau_pvalue[df$tau_pvalue==0] = 1
  
  # Tau
  df[which(df$strategy=="BM"),] = p.adjust(as.numeric(as.character(df[which(df$strategy=="BM"),3]),method="fdr"))
  df[which(df$strategy=="PI"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PI"),3]),method="fdr"))
  df[which(df$strategy=="PD"),3] = p.adjust(as.numeric(as.character(df[which(df$strategy=="PD"),3]),method="fdr"))
  
  tau = data.frame(variable=df$strategy,value=as.numeric(df$tau_pvalue),stringsAsFactors=FALSE)
  tau$value = as.numeric(tau$value)
  
  generate_count_table = function(dat,label,fdr_thresh){
    tabley = table(dat$value[which(dat$variable==label)] <= fdr_thresh)
    if (!("TRUE" %in% names(tabley))) {
      tabley["TRUE"] = 0
    }
    if (!("FALSE" %in% names(tabley))) {
      tabley["FALSE"] = 0
    }    
    return(tabley[c("TRUE","FALSE")])
  }
  
  counts = rbind(generate_count_table(tau,"PD",fdr),
                 generate_count_table(tau,"PI",fdr),
                 generate_count_table(tau,"BM",fdr))
  
  counts = cbind(c("PD_TAU","PI_TAU","BM_TAU"),counts)
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
