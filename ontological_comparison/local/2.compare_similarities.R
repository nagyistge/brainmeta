library(pheatmap)

# Read in spatial and wang (graph-based) similarities
spatial_file = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference/contrast_defined_images_pearsonpd_similarity.tsv"
spatial = read.csv(spatial_file,sep="\t",head=TRUE,stringsAsFactors=FALSE,row.names=1)
graph_file = "/home/vanessa/Documents/Work/BRAINMETA/reverse_inference/contrast_defined_images_wang.tsv"
graph = read.csv(graph_file,sep=",",head=TRUE,stringsAsFactors=FALSE,row.names=1)

contrasts = c()
for (con in images$cognitive_contrast_cogatlas){
  contrasts=c(contrasts,strtrim(con,30))
}

colnames(spatial) = contrasts
colnames(graph) = contrasts

# Just plot as is
pheatmap(graph,cluster_rows=FALSE,cluster_cols=FALSE)
pheatmap(spatial,cluster_rows=FALSE,cluster_cols=FALSE)

# With clustering
pheatmap(graph)
pheatmap(spatial)
