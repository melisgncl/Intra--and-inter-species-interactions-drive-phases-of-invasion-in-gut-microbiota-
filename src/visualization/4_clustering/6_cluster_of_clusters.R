## CODE BLOCK
## Louis Gauthier &Melis Gencel
###Clustering of clone dynamics within cohort####
source("/Users/melis/Dropbox/mouse_gut/github/src/visualization/0_config/0_config.R")

library(dendextend)
library(ggdendro)
library(gplots)
library(viridis)
library(readr)
library(dplyr)
library(phylogram)

args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")
rootpath = args[1]

my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 100)

# reorder columns based on cluster number
reorderCols <- function(df,pattern){
  return(df[,order(sapply(strsplit(names(df),pattern),function(x) as.numeric(unlist(x)))[2,])])
}

# pad cluster dfs with rows of NA
padNA <- function(df,ign){
  df[nrow(df)+1,] <- NA
  return(df)
}

# allows to cbind with empty dataframe
# credited to https://stackoverflow.com/questions/7962267/
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# bind all samples of a cohort with padding
collateSamples <- function(samples.loess){
  maxrange = max(sapply(X = samples.loess,FUN = nrow)) + 1
  mat = data.frame()
  for(df in samples.loess){
    pad = maxrange - nrow(df)
    df = Reduce(padNA,1:pad,init=df,accumulate=FALSE)
    mat = cbind.fill(mat,df)
  }
  return(mat)
}

# plot (1) a dendrogram of all clusters and (2) a distance matrix with said dendrogram
plotDendrograms <- function(mat,group.name,clust.method="average"){
  dist = as.matrix(1 - cor(mat, use = "pairwise.complete.obs", method = "pearson"))
  clust <- hclust(as.dist(dist),method="average")
  #shuffle(as.dendrogram(clust)) -> dend
  as.dendrogram(clust) -> dend
  write.dendrogram(dend, file = paste(outdir, "/clusters_of_clusters/",group.name,".tree",sep=""), append = FALSE, edges = TRUE)
  write_csv(as.data.frame(dist),file = paste(outdir, "/clusters_of_clusters/",group.name,".matrix",sep=""),col_names = TRUE)
  
  
  cluster.colors=c("im1.C1"="#3cb44b","im1.C2"= "#4363d8","im1.C3"="#e6194B","im1.C4"="#e8ca00","im1.C5"="#911eb4",
                   "im1.C6"="#f58231","im1.C7"="#22766dff","im1.C8"="#42d4f4","im1.C9"="#f032e6","im1.C10"="#9A6324",
                   "im2.C1"="#3cb44b","im2.C2"= "#4363d8","im2.C3"="#e6194B","im2.C4"="#e8ca00","im2.C5"="#911eb4",
                   "im2.C6"="#f58231","im2.C7"="#22766dff","im2.C8"="#42d4f4","im2.C9"="#f032e6","im2.C10"="#9A6324",
                   "im3.C1"="#3cb44b","im3.C2"= "#4363d8","im3.C3"="#e6194B","im3.C4"="#e8ca00","im3.C5"="#911eb4",
                   "im3.C6"="#f58231","im3.C7"="#22766dff","im3.C8"="#42d4f4","im3.C9"="#f032e6","im3.C10"="#9A6324",
                   "im4.C1"="#3cb44b","im4.C2"= "#4363d8","im4.C3"="#e6194B","im4.C4"="#e8ca00","im4.C5"="#911eb4",
                   "im4.C6"="#f58231","im4.C7"="#22766dff","im4.C8"="#42d4f4","im4.C9"="#f032e6","im4.C10"="#9A6324",
                   "rm1.C1"="#3cb44b","rm1.C2"= "#4363d8","rm1.C3"="#e6194B","rm1.C4"="#e8ca00","rm1.C5"="#911eb4",
                   "rm1.C6"="#f58231","rm1.C7"="#22766dff","rm1.C8"="#42d4f4","rm1.C9"="#f032e6","rm1.C10"="#9A6324",
                   "rm2.C1"="#3cb44b","rm2.C2"= "#4363d8","rm2.C3"="#e6194B","rm2.C4"="#e8ca00","rm2.C5"="#911eb4",
                   "rm2.C6"="#f58231","rm2.C7"="#22766dff","rm2.C8"="#42d4f4","rm2.C9"="#f032e6","rm2.C10"="#9A6324",
                   "rm3.C1"="#3cb44b","rm3.C2"= "#4363d8","rm3.C3"="#e6194B","rm3.C4"="#e8ca00","rm3.C5"="#911eb4",
                   "rm3.C6"="#f58231","rm3.C7"="#22766dff","rm3.C8"="#42d4f4","rm3.C9"="#f032e6","rm3.C10"="#9A6324",
                   "rm4.C1"="#3cb44b","rm4.C2"= "#4363d8","rm4.C3"="#e6194B","rm4.C4"="#e8ca00","rm4.C5"="#911eb4",
                   "rm4.C6"="#f58231","rm4.C7"="#22766dff","rm4.C8"="#42d4f4","rm4.C9"="#f032e6","rm4.C10"="#9A6324",
                   "gf1.C1"="#3cb44b","gf1.C2"= "#4363d8","gf1.C3"="#e6194B","gf1.C4"="#e8ca00","gf1.C5"="#911eb4",
                   "gf1.C6"="#f58231","gf1.C7"="#22766dff","gf1.C8"="#42d4f4","gf1.C9"="#f032e6","gf1.C10"="#9A6324",
                   "gf2.C1"="#3cb44b","gf2.C2"= "#4363d8","gf2.C3"="#e6194B","gf2.C4"="#e8ca00","gf2.C5"="#911eb4",
                   "gf2.C6"="#f58231","gf2.C7"="#22766dff","gf2.C8"="#42d4f4","gf2.C9"="#f032e6","gf2.C10"="#9A6324",
                   "gf3.C1"="#3cb44b","gf3.C2"= "#4363d8","gf3.C3"="#e6194B","gf3.C4"="#e8ca00","gf3.C5"="#911eb4",
                   "gf3.C6"="#f58231","gf3.C7"="#22766dff","gf3.C8"="#42d4f4","gf3.C9"="#f032e6","gf3.C10"="#9A6324",
                   "gf4.C1"="#3cb44b","gf4.C2"= "#4363d8","gf4.C3"="#e6194B","gf4.C4"="#e8ca00","gf4.C5"="#911eb4",
                   "gf4.C6"="#f58231","gf4.C7"="#22766dff","gf4.C8"="#42d4f4","gf4.C9"="#f032e6","gf4.C10"="#9A6324")
  
  labels_colors(dend) <- cluster.colors[colnames(dist)][order.dendrogram(dend)]
  col_labels<-cluster.colors[colnames(dist)]
  png(paste(outdir, "/clusters_of_clusters/clusters_dendromatrix_", group.name,".png",sep=""),res = 300,width = 5.5,height = 5, units="in")
  par(mar = c(2,2,2,2))
  heatmap.2(dist,Rowv = dend,Colv = dend,col=rev(my_palette),density.info = "none",trace = "none",
            key.xlab="(1 - r)",cexRow = 0.5,cexCol = 0.5,colRow = col_labels, colCol = col_labels)
  dev.off()
  return(list(dist,dend))
}

########################################
rm1_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/rm1/rm1_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))
rm2_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/rm2/rm2_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))
rm3_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/rm3/rm3_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))
rm4_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/rm4/rm4_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))

rm1_clustered_loess.log10 = reorderCols(rm1_clustered_loess.log10,"C")
rm2_clustered_loess.log10 = reorderCols(rm2_clustered_loess.log10,"C")
rm3_clustered_loess.log10 = reorderCols(rm3_clustered_loess.log10,"C")
rm4_clustered_loess.log10 = reorderCols(rm4_clustered_loess.log10,"C")

colnames(rm1_clustered_loess.log10)=paste("rm1.",colnames(rm1_clustered_loess.log10),sep="")
colnames(rm2_clustered_loess.log10)=paste("rm2.",colnames(rm2_clustered_loess.log10),sep="")
colnames(rm3_clustered_loess.log10)=paste("rm3.",colnames(rm3_clustered_loess.log10),sep="")
colnames(rm4_clustered_loess.log10)=paste("rm4.",colnames(rm4_clustered_loess.log10),sep="")


########################################
gf1_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/gf1/gf1_clustered_loess_log10.csv"), 
    col_types = cols(time = col_skip()))
gf2_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/gf2/gf2_clustered_loess_log10.csv"), 
    col_types = cols(time = col_skip()))
gf3_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/gf3/gf3_clustered_loess_log10.csv"), 
    col_types = cols(time = col_skip()))
gf4_clustered_loess.log10 <- read_csv(paste0(outdir, "/clustering/gf4/gf4_clustered_loess_log10.csv"), 
    col_types = cols(time = col_skip()))

gf1_clustered_loess.log10 = reorderCols(gf1_clustered_loess.log10,"C")
gf2_clustered_loess.log10 = reorderCols(gf2_clustered_loess.log10,"C")
gf3_clustered_loess.log10 = reorderCols(gf3_clustered_loess.log10,"C")
gf4_clustered_loess.log10 = reorderCols(gf4_clustered_loess.log10,"C")

colnames(gf1_clustered_loess.log10)=paste("gf1.",colnames(gf1_clustered_loess.log10),sep="")
colnames(gf2_clustered_loess.log10)=paste("gf2.",colnames(gf2_clustered_loess.log10),sep="")
colnames(gf3_clustered_loess.log10)=paste("gf3.",colnames(gf3_clustered_loess.log10),sep="")
colnames(gf4_clustered_loess.log10)=paste("gf4.",colnames(gf4_clustered_loess.log10),sep="")

#############PLOT#######################


samples.loess = list(rm1_clustered_loess.log10,rm2_clustered_loess.log10,rm3_clustered_loess.log10,rm4_clustered_loess.log10)
rm_new.log10.mat = collateSamples(samples.loess)

plotDendrograms(rm_new.log10.mat,"rm_log10")
########################################
samples.loess = list(gf1_clustered_loess.log10,gf2_clustered_loess.log10,gf3_clustered_loess.log10,gf4_clustered_loess.log10)
gf.log10.mat = collateSamples(samples.loess)

plotDendrograms(gf.log10.mat,"gf_log10")



samples.loess = list(rm1_clustered_loess.log10,rm2_clustered_loess.log10,rm3_clustered_loess.log10,rm4_clustered_loess.log10)
rm_new.log10.mat = collateSamples(samples.loess)

plotDendrograms(rm_new.log10.mat,"rm_log10")

