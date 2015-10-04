library(plyr)
library(dplyr)
library(ggplot2)

setwd("/home/vanessa/Documents/Dropbox/Code/Python/brainmeta/ontological_comparison/cluster/classification-framework/analysis")

# Reading in the result data

ri_ranges = read.csv("data/reverse_inference_scores_ranges.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1) # reverse inference ranges scores  
ri_binary = read.csv("data/reverse_inference_scores_binary.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1) # binary scores
ri_priors_in = read.csv("data/reverse_inference_priors_in.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1) # priors for in and out of node sets
ri_priors_out = read.csv("data/reverse_inference_priors_out.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1)
bayes_in_ranges = read.csv("data/reverse_inference_bayes_in_ranges",sep="\t",stringsAsFactors=FALSE,row.names=1) # bayes for query images, ranges in
bayes_out_ranges = read.csv("data/reverse_inference_bayes_out_ranges.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1) # ranges out
bayes_in_bin = read.csv("data/reverse_inference_bayes_in_binary.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1)      # bayes binary bin
bayes_out_bin = read.csv("data/reverse_inference_bayes_out_binary.tsv",sep="\t",stringsAsFactors=FALSE,row.names=1)    # bayes binary out

# Read in all groups
groups = read.csv("data/groups/all_groups.tsv",sep="\t",stringsAsFactors=FALSE)
image_ids = c()
for (image in groups$image){
  image = strsplit(image,"/")[[1]]
  image = as.numeric(strsplit(image[length(image)],"[.]")[[1]][1])
  image_ids = c(image_ids,image)
}
groups$image_ids = image_ids

count_bin = matrix(0,nrow=ncol(bayes_in_bin),ncol=2)
rownames(count_bin) = colnames(bayes_in_bin)
colnames(count_bin) = c("for","against")
count_range = matrix(0,nrow=ncol(bayes_in_ranges),ncol=2)
rownames(count_range) = colnames(bayes_in_ranges)
colnames(count_range) = c("for","against")

# Make a lookup table for the node name
nodes = unique(groups$group)
node_lookup = c()
for (node in nodes){
  node_name = unique(groups$name[groups$group==node])
  node_lookup = c(node_lookup,node_name)
}
length(node_lookup) == length(nodes)
names(node_lookup) = nodes

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
# We will also safe a data frame with number in, number out, and average RI scores

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



# NEXT: What we would want to do is look at the change in bayes score as we add concepts KNOWN to be in the set.
