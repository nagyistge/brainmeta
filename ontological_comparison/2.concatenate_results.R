# Concatenate scores into one matrix

output_folder = "/scratch/users/vsochat/DATA/BRAINMETA/ontological_comparison/wang_scores"
input_files = list.files(output_folder,pattern="*.Rda")

# We will put our results in a data frame
similarities = matrix(nrow=nrow(images),ncol=nrow(images))
rownames(similarities) = images$image_id
colnames(similarities) = images$image_id

for (file in input_files){
  load(file)
  similarities[result$row1,result$row2] = result$score
}

output_file = cat(output_folder,"/contrast_defined_images_wang.tsv",sep="")
write.csv(similarities,filename=output_file,sep="\t")