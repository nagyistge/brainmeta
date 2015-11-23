# Step 0: Data Visualization ################################################################

questions = read.csv("cnp_739.tsv",sep="\t")
rownames(questions) = questions$question_label

# Remove test that I exported
questions = questions[-which(questions$assessment_name=="Test"),]
assessments = unique(questions$assessment_name)

# Let's do MDS for all questions, color by assessment
d = dist(vectors) # euclidean distances between the rows
fit = cmdscale(d,eig=TRUE, k=2) # k is the number of dim
ques_text = as.character(questions$question_text[questions$question_label%in%rownames(vectors)])
ques_assessments = as.character(questions$assessment_name[questions$question_label%in%rownames(vectors)])

# Make some colors!
colors = topo.colors(length(unique(ques_assessments)))
lookup = list()
for (c in seq(1,length(unique(ques_assessments)))){
  lookup[[unique(ques_assessments)[c]]] = colors[c]
}
color_vector = c()
for (ques_assessment in ques_assessments) {
  color_vector = c(color_vector,lookup[[ques_assessment]])
}

# Try two dimensions
pdf("CNP_MDS_ALL.pdf",width=30,height=30)
plot(fit$points[,1],fit$points[,2],pch=19,col=color_vector,main="All CNP Questions",xlab="Dimension 1", ylab="Dimension 2") # view results
text(fit$points[,1],fit$points[,2],ques_text,cex=.5) # view results
plot.new()
legend(0,1,pch=19,legend=names(lookup),col=colors,cex=2)
plot(fit$points[,1],fit$points[,2],pch=19,col=color_vector,main="All CNP Questions",xlab="Dimension 1", ylab="Dimension 2") # view results
dev.off()

# Try three dimensions
library(scatterplot3d)
fit = cmdscale(d,eig=TRUE, k=3) # k is the number of dim
pdf("CNP_MDS_ALL_3D.pdf",width=30,height=30)
scatterplot3d(fit$points[,1],fit$points[,2],fit$points[,3],pch=19,color=color_vector,main="All CNP Questions",xlab="Dimension 1", ylab="Dimension 2") # view results
plot(fit$points[,1],fit$points[,2],pch=19,col=color_vector,main="Dimension 1 vs 2",xlab="Dimension 1", ylab="Dimension 2") # view results
plot(fit$points[,1],fit$points[,3],pch=19,col=color_vector,main="Dimension 1 vs 3",xlab="Dimension 1", ylab="Dimension 3") # view results
plot(fit$points[,2],fit$points[,3],pch=19,col=color_vector,main="Dimension 2 vs 3",xlab="Dimension 2", ylab="Dimension 3") # view results
plot.new()
legend(0,1,pch=19,legend=names(lookup),col=colors,cex=2)
dev.off()

# For each assessment, save MDS coordinates and text to visualize
pdf("CNP_MDS_ASSESSMENTS.pdf",width=10,height=10)
for (assessment in assessments){
  # Get all questions for assessment
  ques = as.character(questions$question_label[questions$assessment_name==assessment])
  ques_text = as.character(questions$question_text[questions$assessment_name==assessment])
  # Find vectors
  ques_vectors = vectors[which(rownames(vectors)%in%ques),]
  color = lookup[[assessment]]
  d = dist(ques_vectors) # euclidean distances between the rows
  fit = cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  plot(fit$points[,1],fit$points[,2],pch=19,col=color,main=assessment,xlab="Dimension 1",ylab="Dimension 2") # view results
  text(fit$points[,1],fit$points[,2],ques_text,cex=.8) # view results
}
dev.off()

# Step 1: Data Exploration #################################################################

# Labels
# First we need to parse labels for people. Unfortunately this is what we have to work with.
rx = read.csv("data/CNP_raw.csv",sep=",",stringsAsFactors=FALSE)
scid = rx[,grep("scid",colnames(rx))]
scid = scid[,3:ncol(scid)]

# Get unique labels
unique_diagnoses = c()
for (col in colnames(scid)){
  unique_diagnoses=c(unique_diagnoses,unique(scid[[col]]))
}
unique_diagnoses = unique(unique_diagnoses[-which(unique_diagnoses%in%c("",1,2,3,4,NA,7,-9998,-9999,"RULE OUT PTSD"))])
cat(unique_diagnoses,sep="\n")

diagnoses = array(0,dim=c(nrow(data),length(unique_diagnoses)))
rownames(diagnoses) = rownames(data)
colnames(diagnoses) = unique_diagnoses

# Go through people and find their troubles!
for (row in 1:nrow(data)){
  rx = as.character(scid[row,])
  ptid = rownames(data)[row]
  rx_cols = which(colnames(diagnoses)%in%rx[rx %in% colnames(diagnoses)])
  diagnoses[which(rownames(diagnoses)==ptid),rx_cols] = 1
}

# Save diagnoses matrix
write.table(diagnoses,file="data/cnp_rx_table_.tsv",sep="\t")
save(diagnoses,file="data/cnp_rx_table.Rda")

# Summarize
summy = sort(colSums(diagnoses),decreasing=TRUE)
cat(paste(names(summy),summy,sep=": "),sep="\n")


# Generate vector for each assessment
# Lets now make a very big assumption that the different sets of questions were selected to b
# meaningful for something, and so our first level of operations should be over assessment question sets
# Let's try using the question answers as a weight.

# We need to keep track of how many questions each person has had included
question_count = array(0,dim=nrow(data))
names(question_count) = rownames(data)

# Let's also keep track of the sum of "weights"
question_weights = array(0,dim=nrow(data))
names(question_weights) = rownames(data)

# Functions to repeat rows and columns
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

cnp_phenotypes = list()

# First let's say that a value of 0 is equivalent to saying "NA" or we don't know
for (assessment in assessments){
  cat("Parsing",assessment,"\n")
  ques = as.character(questions$question_label[questions$assessment_name==assessment])
  ques_text = as.character(questions$question_text[questions$assessment_name==assessment])
  ques_vectors = vectors[which(rownames(vectors)%in%ques),]
  # Now get individual data
  answers = data[,which(colnames(data)%in%ques)]
  phenotypes = array(0,dim=c(nrow(answers),ncol(vectors)))
  rownames(phenotypes) = rownames(data)
  # For each question, scale the vector by the answer
  for (question in ques){
    if(question %in% colnames(answers)){
      question_data = answers[,which(colnames(answers)%in%question)]
      names(question_data) = rownames(phenotypes)
      # This is to assume a value of "0" is equivalent to not knowing
      question_data = question_data[!is.na(question_data)]
      question_data = question_data[which(question_data!=-9999)]
      if (length(question_data)>0){ # some are all NA
        new_vectors = rep.row(as.numeric(ques_vectors[question,]),length(question_data))
        weightby = rep.col(as.numeric(question_data),ncol(phenotypes))
        new_vectors = new_vectors * weightby
        phenotypes[names(question_data),] = phenotypes[names(question_data),] + new_vectors  
        # Update counts of questions, and weights
        nonzero_peeps = names(question_data[question_data!=0])
        nonzero_weights = as.numeric(question_data[question_data!=0])
        question_count[nonzero_peeps] = question_count[nonzero_peeps] + 1 
        question_weights[nonzero_peeps] = question_weights[nonzero_peeps] + nonzero_weights 
      }
    }
  }
  result = list("phenotypes"=phenotypes,"question_count"=question_count,"question_weights"=question_weights)
  cnp_phenotypes[[assessment]] = result
}

# Save our progress
save(cnp_phenotypes,file="data/cnp_phenotypes.Rda")

# Let's make diagnosis groups just based on the "high level" disorder
diagnoses = as.data.frame(diagnoses)
disorder_groups = c("Bipolar","Depress","Schizo","Attention-Deficit","Abuse","Dependence","Alcohol","Cannabis")
for (disorder_group in disorder_groups){
  bipolar = rowSums(diagnoses[,grep(disorder_group,colnames(diagnoses))])
  bipolar[bipolar!=0] = 1
  diagnoses[,disorder_group] = bipolar
}

# Let's go through disorders based on N
sorted_rx = names(sort(colSums(diagnoses),decreasing=TRUE))

# For each assessment, let's try normalizing by weight, and then clustering 
# and applying disorder groups.
skip_these = c("Adult ADHD Clinical Diagnostic Scale")
for (assessment in names(cnp_phenotypes)){
  if (!(assessment %in% skip_these)){
    pdf(gsub(" ","_",paste("data/",assessment,"_phenos.pdf",sep="")))
    phenos = cnp_phenotypes[[assessment]]$phenotypes
    divideby = cnp_phenotypes[[assessment]]$question_weights
    phenos_norm = phenos / divideby  
    d = dist(phenos_norm) # euclidean distances between the rows
    fit = cmdscale(d,eig=TRUE, k=2) # k is the number of dim
    # Make a plot for sorted_rxrder label
    for (diagnosis in sorted_rx){
      color_vector = array("#000000",dim=nrow(phenos))
      names(color_vector) = rownames(phenos)
      peeps = rownames(diagnoses)[which(diagnoses[,which(colnames(diagnoses)==diagnosis)]==1)]
      color_vector[peeps] = "#FF3300"
      title = paste(assessment,"\n",diagnosis,"\nN=",length(peeps),sep="")
      plot(fit$points[,1],fit$points[,2],pch=19,col=color_vector,main=title,xlab="Dimension 1",ylab="Dimension 2") # view results
    }
 dev.off()
  }
}

# Holy crap. Signal!

# Experiment 2: Coming up with custom question groupings
# I would next want to derive a set of questions based on selecting by context only. We
# can let the researcher find vectors that are most similar to "anxiety" for example. The idea 
# is that we want to be able to get a set of questions that reflects a specific trait.

word_matches = read.csv("data/cnp_word_matches.tsv",sep="\t",header=TRUE,row.names=1)

# Set random threshold of .5
threshold = 0.5

# Derive new phenotypes
for (word in rownames(word_matches)){
  cat("Parsing",word,"\n")
  # Let's try taking both absolute value, and not
  word_vector = word_matches[word,which(!is.na(word_matches[word,]))]
  topN = sort(word_vector)
  topN = topN[which(topN>=threshold)]
  topNabs = sort(abs(word_vector))
  topNabs = topNabs[which(topNabs>=threshold)]
  ques_vectors = vectors[which(rownames(vectors)%in%names(topN)),]
  # Now get individual data
  answers = data[,which(colnames(data)%in%names(topN))]
  if (length(topN)>1){
      phenotypes = array(0,dim=c(max(nrow(answers),1),ncol(vectors)))
      rownames(phenotypes) = rownames(data)
      # For each question, scale the vector by the answer
      for (question in names(topN)){
        if(question %in% colnames(answers)){
          question_data = answers[,which(colnames(answers)%in%question)]
          names(question_data) = rownames(phenotypes)
          # This is to assume a value of "0" is equivalent to not knowing
          question_data = question_data[!is.na(question_data)]
          question_data = question_data[which(question_data!=-9999)]
          if (length(question_data)>0){ # some are all NA
            new_vectors = rep.row(as.numeric(ques_vectors[question,]),length(question_data))
            weightby = rep.col(as.numeric(question_data),ncol(phenotypes))
            new_vectors = new_vectors * weightby
            phenotypes[names(question_data),] = phenotypes[names(question_data),] + new_vectors  
            # Update counts of questions, and weights
            nonzero_peeps = names(question_data[question_data!=0])
            nonzero_weights = as.numeric(question_data[question_data!=0])
            question_count[nonzero_peeps] = question_count[nonzero_peeps] + 1 
            question_weights[nonzero_peeps] = question_weights[nonzero_peeps] + nonzero_weights 
          }
        }
      }
      result = list("phenotypes"=phenotypes,"question_count"=question_count,"question_weights"=question_weights)
      cnp_phenotypes[[word]] = result
    }
}

for (word in rownames(word_matches)){
    if (word %in% names(cnp_phenotypes)){
      pdf(gsub(" ","_",paste("data/",word,"_phenos.pdf",sep="")))
      phenos = cnp_phenotypes[[word]]$phenotypes
      divideby = cnp_phenotypes[[assessment]]$question_weights
      phenos_norm = phenos / divideby  
      d = dist(phenos_norm) # euclidean distances between the rows
      fit = cmdscale(d,eig=TRUE, k=2) # k is the number of dim
      # Make a plot for sorted_rxrder label
      for (diagnosis in sorted_rx){
        color_vector = array("#000000",dim=nrow(phenos))
        names(color_vector) = rownames(phenos)
        peeps = rownames(diagnoses)[which(diagnoses[,which(colnames(diagnoses)==diagnosis)]==1)]
        color_vector[peeps] = "#FF3300"
        title = paste(assessment,"\n",diagnosis,"\nN=",length(peeps),sep="")
        plot(fit$points[,1],fit$points[,2],pch=19,col=color_vector,main=title,xlab="Dimension 1",ylab="Dimension 2") # view results
      }
      dev.off()
    }
}
