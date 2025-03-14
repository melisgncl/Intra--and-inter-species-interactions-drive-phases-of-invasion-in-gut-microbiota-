## R code
## Louis Gauthier & Melis Gencel, August 2020

#######################    READ ME!   #####################################

###	Here we look at the similarity of clusters across different mice from the same cohorts.
###	This is one way to assess reproducibility of outcomes and dynamics within a cohort.
### It's also the cornerstone for the inference of any genetic effects on the dynamics.

#######################    ^^^^^^^^   #####################################

source("/Users/melis/Dropbox/mouse_gut/github/src//visualization/0_config/0_config.R")

library(readr)
library(hash)
library(phylogram)
library(gplots)

z_palette <- colorRampPalette(c("red", "white", "blue"))(n = 100)
overlap_palette <- colorRampPalette(c("white","blue"))(n = 100)

##########################
#### distance metrics ####
##########################

overlapCoef <- function(s1,s2) {
  i = length(intersect(s1,s2))
  u = min(length(s1),length(s2))
  return(i/u)
}


#######################################
#### estimated Z-score calculation ####
#######################################

# computes  overlap index of simulated samples from sets X and Y
sampleSimul <- function(s1,s2,X,Y){
  xs = sample(x=X,size=s1)
  ys = sample(x=Y,size=s2)
  return(length(intersect(xs,ys))/min(length(xs),length(ys)))
}

# computes zscore for value using distribution
zScore <- function(simul.dist,x.value){
  mu = mean(simul.dist)
  sd = sd(simul.dist)
  return((x.value - mu)/sd)
}

### args
### s1: list of barcodes from cluster 1
### s2: list of barcodes from cluster 2
### X: library of barcodes matching cluster 1
### Y: library of barcodes matching cluster 2 
empiricalZ <- function(s1, s2, X, Y, n_trials){
  # calculate size of intersection
  x.value = length(intersect(s1,s2))/min(length(s1),length(s2))
  simul.dist = replicate(n_trials,sampleSimul(length(s1),length(s2),X,Y))
  z = zScore(simul.dist,x.value)
  return(z)
}

# calculates the z-score matrix for a list of libraries libraries whole pooled library for the mouse whole hash table shows the which barcode in which cluster
computeZmatrix <- function(libraries,hash.table){
  nams <- names(libraries)
  
  # dimensions of matrix
  nrow = sum(sapply(hash.table,function(x) length(x)))
  ncol = nrow
  mat <- matrix(, nrow = nrow, ncol = ncol)
  
  row_ref=1
  col_ref=1
  for(i in seq_along(libraries)){
    for(j in seq_along(libraries)){
      col_shift = 0
      for(x in hash.table[[nams[i]]]){
        row_shift = 0
        for(y in hash.table[[nams[j]]]){
          mat[row_ref+col_shift,col_ref+row_shift] <- empiricalZ(x,y,libraries[[i]],libraries[[j]],1000)
          row_shift = row_shift + 1
        }
        col_shift = col_shift + 1
      }
      col_ref = col_ref + length(hash.table[[nams[j]]])
    }
    row_ref = row_ref + length(hash.table[[nams[i]]])
    col_ref = 1
  }
  return(mat)
}

# calculates the overlap coefficient for a list of libraries
computeOverlap <- function(libraries,hash.table){
  nams <- names(libraries)
  
  # dimensions of matrix
  nrow = sum(sapply(hash.table,function(x) length(x)))
  ncol = nrow
  mat <- matrix(, nrow = nrow, ncol = ncol)
  
  row_ref=1
  col_ref=1
  for(i in seq_along(libraries)){
    for(j in seq_along(libraries)){
      col_shift = 0
      for(x in hash.table[[nams[i]]]){
        row_shift = 0
        for(y in hash.table[[nams[j]]]){
          mat[row_ref+col_shift,col_ref+row_shift] <-  overlapCoef(x,y)
          row_shift = row_shift + 1
        }
        col_shift = col_shift + 1
      }
      col_ref = col_ref + length(hash.table[[nams[j]]])
    }
    row_ref = row_ref + length(hash.table[[nams[i]]])
    col_ref = 1
  }
  return(mat)
}

## split clusters
splitClusters <- function(clusters){
  composition = split(clusters,f = clusters$cluster)
  composition = lapply(composition, "[", , "Center")
  
  # use 'names' to iterate through list of lists
  lapply(names(composition), function(x) length(composition[[x]]))
  return(composition)
}

#### IMPORT INDEX
##all barcode pool from cohorts with also gavaged library###
Gavage.index <- read_csv(paste0("data/",data_folder,"Gavage_index_overlap.csv"))
mergeFormat <- function(series,sample.name){
	center = unique(series[,c("ID","cluster")])
	center = merge(center,Gavage.index,by.x = "ID",by.y = paste(sample.name,"ID",sep="."))[,2:3]
	return(center)
}

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


########################
#### Germ-free mice ####
########################









#### IMPORT DATA
gf1_clustered_series <- read_csv(paste0(outdir, "/clustering/gf1/gf1_clustered_series_log10.csv"))
gf2_clustered_series <- read_csv(paste0(outdir, "/clustering/gf2/gf2_clustered_series_log10.csv"))
gf3_clustered_series <- read_csv(paste0(outdir, "/clustering/gf3/gf3_clustered_series_log10.csv"))
gf4_clustered_series <- read_csv(paste0(outdir, "/clustering/gf4/gf4_clustered_series_log10.csv"))

# get barcode Center (universal)
gf1.center = mergeFormat(gf1_clustered_series,"gf1")
gf2.center = mergeFormat(gf2_clustered_series,"gf2")
gf3.center = mergeFormat(gf3_clustered_series,"gf3")
gf4.center = mergeFormat(gf4_clustered_series,"gf4")

# split to create lists of barcode per cluster
gf1.composition = splitClusters(gf1.center)
gf2.composition = splitClusters(gf2.center)
gf3.composition = splitClusters(gf3.center)
gf4.composition = splitClusters(gf4.center)

# create hash table of clusters
gf.hash <- hash()
gf.hash[["gf1"]] = gf1.composition
gf.hash[["gf2"]] = gf2.composition
gf.hash[["gf3"]] = gf3.composition
gf.hash[["gf4"]] = gf4.composition

# get "library" of barcodes per sample
gf1.ids = as.data.frame(unique(gf1_clustered_series$ID))
gf2.ids = as.data.frame(unique(gf2_clustered_series$ID))
gf3.ids = as.data.frame(unique(gf3_clustered_series$ID))
gf4.ids = as.data.frame(unique(gf4_clustered_series$ID))

colnames(gf1.ids)="id"
colnames(gf2.ids)="id"
colnames(gf3.ids)="id"
colnames(gf4.ids)="id"

gf1.library = merge(gf1.ids,Gavage.index,by.x = "id",by.y = "gf1.ID")[,2]
gf2.library = merge(gf2.ids,Gavage.index,by.x = "id",by.y = "gf2.ID")[,2]
gf3.library = merge(gf3.ids,Gavage.index,by.x = "id",by.y = "gf3.ID")[,2]
gf4.library = merge(gf4.ids,Gavage.index,by.x = "id",by.y = "gf4.ID")[,2]

# list of libraries to loop through
gf.libraries = list(gf1.library,gf2.library,gf3.library,gf4.library)
names(gf.libraries) = c("gf1","gf2","gf3","gf4")





rm1_clustered_series <- read_csv(paste0(outdir, "/clustering/rm1/rm1_clustered_series_log10.csv"))
rm2_clustered_series <- read_csv(paste0(outdir, "/clustering/rm2/rm2_clustered_series_log10.csv"))
rm3_clustered_series <- read_csv(paste0(outdir, "/clustering/rm3/rm3_clustered_series_log10.csv"))
rm4_clustered_series <- read_csv(paste0(outdir, "/clustering/rm4/rm4_clustered_series_log10.csv"))

# get barcode Center (universal)
rm1.center = mergeFormat(rm1_clustered_series,"rm1")
rm2.center = mergeFormat(rm2_clustered_series,"rm2")
rm3.center = mergeFormat(rm3_clustered_series,"rm3")
rm4.center = mergeFormat(rm4_clustered_series,"rm4")

# split to create lists of barcode per cluster
rm1.composition = splitClusters(rm1.center)
rm2.composition = splitClusters(rm2.center)
rm3.composition = splitClusters(rm3.center)
rm4.composition = splitClusters(rm4.center)

# create hash table of clusters
rm.hash <- hash()
rm.hash[["rm1"]] = rm1.composition
rm.hash[["rm2"]] = rm2.composition
rm.hash[["rm3"]] = rm3.composition
rm.hash[["rm4"]] = rm4.composition

# get "library" of barcodes per sample
rm1.ids = as.data.frame(unique(rm1_clustered_series$ID))
rm2.ids = as.data.frame(unique(rm2_clustered_series$ID))
rm3.ids = as.data.frame(unique(rm3_clustered_series$ID))
rm4.ids = as.data.frame(unique(rm4_clustered_series$ID))

colnames(rm1.ids)="id"
colnames(rm2.ids)="id"
colnames(rm3.ids)="id"
colnames(rm4.ids)="id"

rm1.library = merge(rm1.ids,Gavage.index,by.x = "id",by.y = "rm1.ID")[,2]
rm2.library = merge(rm2.ids,Gavage.index,by.x = "id",by.y = "rm2.ID")[,2]
rm3.library = merge(rm3.ids,Gavage.index,by.x = "id",by.y = "rm3.ID")[,2]
rm4.library = merge(rm4.ids,Gavage.index,by.x = "id",by.y = "rm4.ID")[,2]

# list of libraries to loop through
rm.libraries = list(rm1.library,rm2.library,rm3.library,rm4.library)
names(rm.libraries) = c("rm1","rm2","rm3","rm4")


get_overlap_graphs=function(libraries,hash,sample,cor_name){
  library(tidyverse)
  z.matrix = computeZmatrix(libraries,hash)
  ###double tailed pvalues###
  pvalue2sided=-log10(2*pnorm(-abs(z.matrix)))  %>% replace(., is.infinite(.), 0)
  #pvalue2sided_2=-log10(2*pnorm(z.matrix, lower.tail=FALSE))
  
  overlap.matrix = computeOverlap(libraries,hash)
  mat = read_csv(paste0(outdir, "/clusters_of_clusters/",cor_name,"_log10.matrix"))
  rownames(mat) = colnames(mat)
  clust <- hclust(as.dist(mat),method="average")
  as.dendrogram(clust) -> dend
  labels_colors(dend) <- cluster.colors[colnames(mat)][order.dendrogram(dend)]
  col_labels<-cluster.colors[colnames(mat)]
  colnames(z.matrix) = colnames(mat)
  rownames(z.matrix) = colnames(mat)
  colnames(overlap.matrix) = colnames(mat)
  rownames(overlap.matrix) = colnames(mat)
  colnames(pvalue2sided) = colnames(mat)
  rownames(pvalue2sided) = colnames(mat)
  #write_csv(as.data.frame(z.matrix),file = paste0(outdir, "/overlap_of_clusters/",sample,"_z.matrix"),col_names = TRUE)
  write_csv(as.data.frame(overlap.matrix ),file = paste0(outdir, "/overlap_of_clusters/",sample,"_overlap_.matrix"),col_names = TRUE)
  #write_csv(as.data.frame(pvalue2sided),file = paste0(outdir, "/overlap_of_clusters/",sample,"_p_.matrix"),col_names = TRUE)
  
  ####heatmap_plots##
  cairo_ps(paste0(outdir, "/overlap_of_clusters/",sample,"_clusters_matrix_jaccard.eps"),width = 5.5,height = 5)
  par(mar = c(2,2,2,2))
  heatmap.2(overlap.matrix,Rowv = dend,Colv = dend,col=overlap_palette,density.info = "none",trace = "none",
            key.xlab="Jaccard index",cexRow = 0.5,cexCol = 0.5,colRow = col_labels, colCol = col_labels)
  dev.off()
  ##### overlap versus correlation and color them according to p values from z matrix
  df_z=data.frame(Var1=t(combn(colnames(z.matrix),2)),Var2=z.matrix[lower.tri(z.matrix)])
  df_z$id <- paste(df_z$Var1.1,df_z$Var1.2 ,sep="_")
  names(df_z)[3]="zmatrix"
  
  df_p=data.frame(Var1=t(combn(colnames(pvalue2sided),2)),Var2=pvalue2sided[lower.tri(pvalue2sided)])
  df_p$id <- paste(df_p$Var1.1,df_p$Var1.2 ,sep="_")
  names(df_p)[3]="pmatrix"
  
  df_o=data.frame(Var1=t(combn(colnames(overlap.matrix),2)),Var2=overlap.matrix[lower.tri(overlap.matrix)])
  df_o$id <- paste(df_o$Var1.1,df_o$Var1.2 ,sep="_")
  names(df_o)[3]="omatrix"
  
  mat=(mat-1)*-1
  df_c=data.frame(Var1=t(combn(colnames(mat),2)),Var2=mat[lower.tri(mat)])
  df_c$id <- paste(df_c$Var1.1,df_c$Var1.2 ,sep="_")
  names(df_c)[3]="cmatrix"
  
  overlap_statics=list(df_z,df_p,df_o,df_c) %>% reduce(full_join, by = "id")
  write_csv(overlap_statics,paste0(outdir, "/overlap_of_clusters/",sample,"_statistic_data.csv"))
  
  grouped_df_p = overlap_statics%>% group_by(id) %>% filter(pmatrix >= -log10(0.05)) %>% ungroup()
  grouped_df_p_not = overlap_statics%>% group_by(id) %>% filter(pmatrix < -log10(0.05)) %>% ungroup()
  
  
  
  p=ggplot()+ geom_point(aes(omatrix,cmatrix),colour="#6a6a6a",data=grouped_df_p_not,size=4,shape=1,stroke = 1.2) + 
    geom_point(aes(omatrix,cmatrix,size=pmatrix),color="blue", data=grouped_df_p,shape=1,stroke = 1.4)+ theme_Publication() +
    xlab("Overlap coefficent") + ylab("Cluster Correlation") + theme( aspect.ratio=1)+xlim(0,1) +
    scale_size_continuous(range = c(4, 10)) +theme(legend.position = "none") 
  
  
  
  ggsave(p,filename=paste0(outdir, "/overlap_of_clusters/",sample,"_overlap-score_Correlation_value.eps"),width = 8.25,height = 8.25 )
}



###########################
#### Compute distances ####
###########################

get_overlap_graphs(rm.libraries,rm.hash,"rm1-4","rm")


###gf###
get_overlap_graphs(gf.libraries,gf.hash,"gf1-4","gf")




