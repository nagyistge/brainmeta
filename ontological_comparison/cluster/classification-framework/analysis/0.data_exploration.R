library(plyr)
library(dplyr)
library(ggplot2)

setwd("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/ontological_comparison/cluster/classification-framework/analysis")

# Reading in the result data

ri_score = read.csv("data/reverse_inference_scores.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1) # reverse inference scores
counts_in = read.csv("data/reverse_inference_counts_in.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1) # bayes for query images, ranges in
counts_out = read.csv("data/reverse_inference_counts_out.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1) # bayes for query images, ranges in

# Let's look at the overall posterior scores
hist(as.matrix(ri_score),main="Reverse Inference Scores",col="orange",xlab="posterior probability")

# Read in all groups
groups = read.csv("data/groups/all_groups.tsv",sep="\t",stringsAsFactors=FALSE)
image_ids = c()
for (image in groups$image){
  image = strsplit(image,"/")[[1]]
  image = as.numeric(strsplit(image[length(image)],"[.]")[[1]][1])
  image_ids = c(image_ids,image)
}
groups$image_ids = image_ids

# Make a lookup table for the node name
nodes = unique(groups$group)
node_lookup = c()
for (node in nodes){
  node_name = unique(groups$name[groups$group==node])
  node_lookup = c(node_lookup,node_name)
}
length(node_lookup) == length(nodes)
names(node_lookup) = nodes

# For each group, calculate an accuracy across thresholds
df = c()

# Calculate accuracy for each node group
# Count evidence for (meaning bayes_in > bayes_out or against (bayes_out > bayes_in)) each concept
# for each of ranges and bin data
for (threshold in seq(0,1,by=0.05)){
  cat("Parsing",threshold,"\n")  
  accuracies = c()
  for (node in nodes){
    # Find in group
    group = groups[groups$group==node,]
    in_group = group$image_ids[which(group$direction=="in")]
    out_group = group$image_ids[which(group$direction=="out")]
    # Get reverse inference scores
    if (node %in% colnames(ri_score)){
      scores = ri_score[,node]
      names(scores) = rownames(ri_score)
      # This case happens when N=1 for the node in question, since we removed the image from the group. The score should be 1.
      scores[is.na(scores)] = 1
      # Image index will have 1 for belonging to class, 0 otherwise
      real = array(0,dim=length(unique(image_ids)))
      predicted = array(0,dim=length(unique(image_ids)))
      names(real) = unique(image_ids)
      names(predicted) = unique(image_ids)
      predicted[names(which(scores>=threshold))] = 1
      real[as.character(in_group)] = 1
      # Calculate metrics
      # c("TP","FP","TN","FN","accuracy","in_count","out_count")
      TP = sum(real*predicted)  
      TN = length(intersect(names(which(real==0)),names(which(predicted==0))))
      FP = length(intersect(names(which(real==0)),names(which(predicted==1))))
      FN = length(intersect(names(which(real==1)),names(which(predicted==0))))
      if (TP+FP==0){
          sens = 0
      } else {
          sens = TP / (TP + FN)
      }
      if (TN+FN==0){
         spec = 0        
      } else {
         spec = TN / (TN + FP)
      }
      accuracy = (TP + TN)/ (TP + TN + FP + FN)
      accuracies = rbind(accuracies,c(node,TP,FP,TN,FN,sens,spec,accuracy,length(in_group),length(out_group),threshold))
    }
  
  }
  df = rbind(df,accuracies)
}

# Now look at accuracies for each threshold!
rownames(df) = seq(1,nrow(df))
colnames(df) = c("nid","TP","FP","TN","FN","sensitivity","specificity","accuracy","in_count","out_count","threshold")
df = as.data.frame(df,stringsAsFactors=FALSE)
save(df,file="accuracies_df_all.rda")

# Plot a basic ROC for each class
pdf("roc_all.pdf")
nodes = unique(df$nid)
pdf("roc_gr30.pdf")
for (node in nodes){
  subset = df[df$nid==node,]
  N = unique(subset$in_count)
  if (as.numeric(N)>30){
    title = paste("ROC Curve ",as.character(node_lookup[node])," N=(",N,")",sep="")
    plot(1-as.numeric(subset$specificity),as.numeric(subset$sensitivity),
         xlab="1-specificity",ylab="sensitivity",main=title,
         xlim=c(0,1),ylim=c(0,1),type="n")
    lines(1-as.numeric(subset$specificity),as.numeric(subset$sensitivity),col="blue",lwd=2,xlim=c(0,1),ylim=c(0,1))
    lines(seq(0,1,0.05),seq(0,1,0.05),col="red",lwd=2,xlim=c(0,1),ylim=c(0,1))
  }
}
dev.off()



# Now we want to assess a multilabel confusion matrix for each threshold.

# Here is a function to assess a multilabel confusion matrix
# From Sanmi Koyejo
# computes the confusion between labels Zt and predictions Ze.
# Assumes that Zt is coded as 0/1 0r -1/+1
# Assumes that the threshold has already been applied to Ze, so sign(Ze) corresponds to a decision
# Includes optional normalization wrt the rows

multilabel_confusion = function(Zt, Ze, normalize=FALSE) {
  L = dim(Zt)[2]
  M = array(0,dim=c(L, L))
  
  for (ix in seq(1,L)) {
    t = Zt[,ix]>0
    for (jx in seq(1,L)){
      p = Ze[,jx]>0
      M[ix, jx] = length(which((p & t)==TRUE))
    }
  }
  
  # To normalize, we divide by number of images tagged with the concept
  if (normalize==TRUE){
    # The images are in rows, concepts in columns
    # so getting a sum for each column --> total number of images tagged
    msums = as.numeric(colSums(Zt))
    for (ix in nrow(M)){
      # Dividing the number we got right by the possible gives an accuracy
      M[ix,] = M[ix,]/msums
    }
  }
  return(M)
}

unique_images = unique(image_ids)

# Our inputs are Zt, the labels, and Ze, the predictions.
# Each is an N by M matrix of N images and M contrasts
# A 1 at index Ze[N,M] indicates image N is predicted to be concept M
# A 1 at index Zt[N,M] means that this is actually the case 

# First let's build our "actual label" matrix, Zt
Zt = array(0,dim=c(length(unique_images),length(nodes)))
rownames(Zt) = unique_images
colnames(Zt) = nodes
# 1 means labeled == YES, -1 means NO
for (node in nodes){
  # Find in group
  group = groups[groups$group==node,]
  in_group = group$image_ids[which(group$direction=="in")]
  Zt[which(rownames(Zt)%in% in_group),node] = 1
}

# Now let's build a matrix to compare, for each threshold
# We will save a list of score matrices.
pdf("multilabel_confusions.pdf")
for (threshold in seq(0,1,by=0.05)){
  cat("Parsing",threshold,"\n")  
  Ze = array(0,dim=c(length(unique_images),length(nodes)))
  rownames(Ze) = unique_images
  colnames(Ze) = nodes
  for (node in nodes){
    # Find in group
    group = groups[groups$group==node,]
    in_group = group$image_ids[which(group$direction=="in")]
    out_group = group$image_ids[which(group$direction=="out")]
    # Get reverse inference scores
    if (node %in% colnames(ri_score)){
      scores = ri_score[,node]
      names(scores) = rownames(ri_score)
      # This case happens when N=1 for the node in question, since we removed the image from the group. The score should be 1.
      scores[is.na(scores)] = 1
      # Image index will have 1 for belonging to class, 0 otherwise
      real = array(0,dim=length(unique(image_ids)))
      predicted = array(0,dim=length(unique(image_ids)))
      names(real) = unique(image_ids)
      names(predicted) = unique(image_ids)
      correct_predictions = names(which(scores>=threshold))
      Ze[which(rownames(Ze) %in% correct_predictions),node] = 1 
    }
  }
  # Calculate multilabel confusion score
  mat = multilabel_confusion(Zt,Ze,FALSE)
  rownames(mat) = node_lookup[nodes]
  colnames(mat) = node_lookup[nodes]
  pheatmap(mat,cluster_rows=FALSE,cluster_cols=FALSE,main=paste("Multi-class confusion for threshold",threshold),fontsize=4)
}
dev.off()


  # Look at bayes for range and bin given "in" group
  for (image in in_group){
    count_range = count_for_against(image,node,bayes_in_ranges,bayes_out_ranges,count_range,"gt")    
    count_bin = count_for_against(image,node,bayes_in_bin,bayes_out_bin,count_bin,"gt")    
  }
  for (image in out_group){
    count_range = count_for_against(image,node,bayes_in_ranges,bayes_out_ranges,count_range,"lt")    
    count_bin = count_for_against(image,node,bayes_in_bin,bayes_out_bin,count_bin,"lt")    
  }
}


count_bin = matrix(0,nrow=ncol(bayes_in_bin),ncol=2)
rownames(count_bin) = colnames(bayes_in_bin)
colnames(count_bin) = c("for","against")
count_range = matrix(0,nrow=ncol(bayes_in_ranges),ncol=2)
rownames(count_range) = colnames(bayes_in_ranges)
colnames(count_range) = c("for","against")


# Write a function to generate evidence "for" or "against" a concept
count_for_against = function(image,node,bayesin,bayesout,count_df,direction){
  image_id = strsplit(image,"/")[[1]]
  image_id = as.character(as.numeric(sub(".nii.gz","",image_id[length(image_id)])))
  
  score_in = bayesin[image_id,node]
  score_out = bayesout[image_id,node]
  
  if (!is.na(score_in) && (!is.na(score_out))) {
    if (direction=="gt"){
      if (score_in > score_out){
        count_df[node,"for"] =  count_df[node,"for"] + 1
      } else {
        count_df[node,"against"] =  count_df[node,"against"] + 1
      }
    } else {
     
      if (score_in < score_out){
        count_df[node,"for"] =  count_df[node,"for"] + 1
      } else {
        count_df[node,"against"] =  count_df[node,"against"] + 1
      }
      
    }
  }
      return(count_df)  
}

# Count evidence for (meaning bayes_in > bayes_out or against (bayes_out > bayes_in)) each concept
# for each of ranges and bin data
for (node in nodes){
   cat("Parsing",node,"\n")  
   # Find in group
   group = groups[groups$group==node,]
   in_group = group$image[which(group$direction=="in")]
   out_group = group$image[which(group$direction=="out")]
   # Look at bayes for range and bin given "in" group
   for (image in in_group){
      count_range = count_for_against(image,node,bayes_in_ranges,bayes_out_ranges,count_range,"gt")    
      count_bin = count_for_against(image,node,bayes_in_bin,bayes_out_bin,count_bin,"gt")    
  }
   for (image in out_group){
     count_range = count_for_against(image,node,bayes_in_ranges,bayes_out_ranges,count_range,"lt")    
     count_bin = count_for_against(image,node,bayes_in_bin,bayes_out_bin,count_bin,"lt")    
    }
}

### STEP 1: VISUALIZATION #########################################################################
# Note - this does not completely coincide with order of google site

cr = melt(count_bin)
colnames(cr) = c("node","direction","value")

# Evidence for and against the concepts, not normalized
pdf("img/evidence_for_concepts.pdf")
for (node in nodes){
  subset = cr[cr$node==node,] 
  p = ggplot(subset,aes(x=direction,y=value,fill=direction)) + 
  geom_histogram(alpha=0.25,stat="identity",binwidth=1) +
  labs(title = paste("Evidence for/against",node_lookup[node]))
  print(p)
}
dev.off()

# Idea 1: If the "in" group images provide evidence for the concept, on the level of the node
# Evidence for and against the concepts, but now we only want to consider the images that are tagged AT the node:
count_bin_in = matrix(0,nrow=ncol(bayes_in_bin),ncol=2)
rownames(count_bin_in) = colnames(bayes_in_bin)
colnames(count_bin_in) = c("for","against")

for (node in nodes){
  cat("Parsing",node,"\n")  
  # Find in group
  group = groups[groups$group==node,]
  in_group = group$image[which(group$direction=="in")]
  # Look at bayes for range and bin given "in" group
  for (image in in_group){
    count_bin_in = count_for_against(image,node,bayes_in_bin,bayes_out_bin,count_bin_in,"gt")    
  }
}

crin = melt(count_bin_in)
colnames(crin) = c("node","direction","value")

pdf("img/evidence_for_concepts_ins.pdf")
for (node in nodes){
  subset = crin[crin$node==node,] 
  if (sum(subset$value)>0) {
  p = ggplot(subset,aes(x=direction,y=value,fill=direction)) + 
    geom_histogram(alpha=0.25,stat="identity",binwidth=1) +
    scale_y_continuous(limits = c(0, 93)) +
    labs(title = paste("Evidence for/against",node_lookup[node]))
    print(p)
  }
}
dev.off()

# Basic visualization of reverse inference scores
riranges$image_id = rownames(ri_ranges)
ribinary$image_id = rownames(ri_binary)
riranges = as.matrix(ri_ranges)
ribinary = as.matrix(ri_binary)
par(mfrow=c(1,2))
hist(riranges,col="purple",main="Reverse Inference Scores, Range Method")
hist(ribinary,col="blue",main="Reverse Inference Scores, Binary Method")

# They are essentially equivalent, but we will look at them in detail anyway
# Reverse inference scores by task
rir = melt(riranges)
colnames(rir) = c("image","task","value")
rib = melt(ribinary)
colnames(rib) = c("image","task","value")
rir$task = as.character(rir$task)
rir$image = as.character(rir$image)
rir$value = as.numeric(rir$value)
rib$task = as.character(rib$task)
rib$image = as.character(rib$image)
rib$value = as.numeric(rib$value)

# Let's look at the distributions individually for each concept!
# we will write to pdf
# We will also save a data frame with number in, number out, and average RI scores
# (not not used)
df=c("countin","countout","meanin","meanout")

par(mfrow=c(1,1))
pdf(file="img/node-ri-scores.pdf")
for (node in nodes){
  subset1 = rir$value[rir$task==node]
  subset2 = rib$value[rib$task==node]
  minvalue = min(min(subset1),min(subset2))
  maxvalue = max(max(subset1),max(subset2))
  #hist(subset,col="purple",main="RI Scores",xlim=c(minvalue,maxvalue),xlab="ranges method")
  #hist(subset,col="blue",main=as.character(node_lookup[node]),xlim=c(minvalue,maxvalue),xlab="binary method")
  # Now count the in vs out group
  group = groups[groups$group == node,]
  rir$direction = group$direction[as.character(group$image_ids) %in% rir$image]
  rib$direction = group$direction[as.character(group$image_ids) %in% rib$image]
  ingroup = length(group$image[group$direction=="in"])
  outgroup = length(group$image[group$direction=="out"])
  meanin = mean(rib$value[rib$direction=="in"])
  meanout = mean(rib$value[rib$direction=="out"])
  df = rbind(df,c(ingroup,outgroup,meanin,meanout))
  p = ggplot(rib,aes(x=value,fill=direction)) + 
  geom_histogram(alpha=0.25,stat="bin",binwidth=0.1) +
  facet_wrap(~direction) +
  labs(title = node_lookup[node])
  print(p)
}
dev.off()

# Function to calculate confidence intervals
get_ci = function(dat,direction="upper"){
  error = qnorm(0.975)*sd(dat)/sqrt(length(dat))
  if (direction=="upper"){
    return(mean(dat)+error)
  } else {
    return(mean(dat)-error)    
  }
}

# Let's look at mean RI scores for each concept node
ribsum = ddply(rib, c("task"), summarise, mscore = mean(value))
ribsum$name = as.character(node_lookup[ribsum$task])

# Sort by meanscore
tmp = ribsum[with(ribsum, order(-mscore)), ]
rownames(tmp) = seq(1,nrow(tmp))
tmp$sort = as.numeric(rownames(tmp))
ggplot(tmp, aes(x=sort,y=mscore,task=task,colour=mscore)) + 
  geom_bar(stat="identity") + 
  xlab("concept") +
  ylab(paste("mean reverse inference score")) +
  scale_x_discrete(limits=tmp$sort,labels=tmp$name) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

countin = c()
countout = c()
for (t in tmp$task){
  subset= groups[groups$group==t,]
  countin = c(countin,nrow(subset[subset$direction=="in",]))
  countout = c(countout,nrow(subset[subset$direction=="out",]))
}
tmp$countin = countin
tmp$countout = countout

# Look at overall mean reverse inference scores
ggplot(tmp, aes(x=sort,y=mscore,task=task,countin=countin,countout=countout,colour=mscore)) + 
  geom_bar(stat="identity") + 
  xlab("concept") +
  ylab(paste("mean reverse inference score")) +
  scale_x_discrete(limits=tmp$sort,labels=tmp$name) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Look at count in and out of set
ggplot(tmp, aes(x=sort,task=task,countin=countin,countout=countout)) + 
  geom_point(aes(y = countin,colour='pink')) +
  geom_point(aes(y = countout)) +
  scale_x_discrete(limits=tmp$sort,labels=tmp$name) +
  xlab("concept") +
  ylab(paste("count in (pink) and out")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")


# Idea 2: Can reverse inference scores be used to predict image labels
# First make a matrix of images by labels, with 1 if the image is in the "in group"
unique_ids = unique(image_ids)
labels = array(0,dim=c(length(unique_ids),length(nodes)))
rownames(labels) = unique_ids
colnames(labels) = nodes
for (node in nodes){
  subset = groups[groups$group==node,]
  ingroup = subset$image_ids[subset$direction=="in"]
  labels[which(rownames(labels)%in%ingroup),node] = 1
}

# Let's first be stupid and see if we can cluster concepts based on image scores
library(pheatmap)

ribname = ri_binary
colnames(ribname) = as.character(node_lookup[colnames(ribname)])
disty = dist(ribname)
dmat = as.matrix(disty)

# Correlation of concepts by image scores
corr = cor(ribname)
pheatmap(corr,fontsize_row=8)

# Correlation of images by image scores
corr = cor(t(ribname))
pheatmap(corr,fontsize_row=8)

# Get contrast names for the images
images = read.csv("/home/vanessa/Documents/Work/BRAINMETA/reverse_inference/contrast_defined_images.tsv",sep="\t")
contrast_names = as.character(images$cognitive_contrast_cogatlas[as.character(images$image_id)%in%rownames(ribname)])
contrast_names = strtrim(contrast_names, 25)
rownames(corr) = contrast_names
colnames(corr) = contrast_names
pheatmap(corr,fontsize_row=8)

# NEXT: What we would want to do is look at the change in bayes score as we add concepts KNOWN to be in the set.
# want to get feedback on this first before doing it
