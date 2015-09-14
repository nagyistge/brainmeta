options(java.parameters = "-Xmx4g") # This ensures we don't run out of memory
library(CogatSimilar) # https://github.com/CognitiveAtlas/cogat-similaR

args <- commandArgs(TRUE)
i = as.numeric(args[1])
output_file = args[2]
image_file = args[3]

# Read in table with images
images = read.csv(image_file,sep="\t",head=TRUE,stringsAsFactors=FALSE)

# A vector to hold similarities
similarities = array(dim=nrow(images))

cat("Processing",i,"of",nrow(images),"\n")
mr1 = images[i,]
CAID1 = mr1$cognitive_contrast_cogatlas_id
for (j in 1:nrow(images)){
  mr2 = images[j,]
  CAID2 = mr2$cognitive_contrast_cogatlas_id
  score = CogatSimilar(CAID1,CAID2)
  similarities[j] = score
}
names(similarities) = images$image_id

result = list()
result$scores = similarities
result$row = i

# Export to file
save(result,filename=output_file)