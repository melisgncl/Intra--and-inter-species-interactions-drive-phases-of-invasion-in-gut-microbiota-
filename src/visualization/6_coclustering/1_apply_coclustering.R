## R code
##Melis Gencel

#######################    READ ME!   #####################################

###	Here we co-cluster the community dynamics(16S) with clone dynamics(E.coli lineages)
##to understand the  heteregenous interaction between intra and interspecies.

### we used z-normalized  series remove amplitude of the dynamics 
#then perform hierhical  clustering  based on shape based distance 
#as a method to efficiently group together subsequences characterized by their shape similarity 
#despite differences in amplitude, duration and alignment(similar to k-shaped clustering)
#######################    ^^^^^^^^   #####################################



source("/Users/melis/Dropbox/mouse_gut/github/src//visualization/0_config/0_config.R")

library(svglite)
library(gridExtra)
quartzFonts(arial = c("Arial", "Arial Black", "Arial Italic", "Arial Black Italic"))
par(family = "Arial")

#cluster_cols = c("#DD0E1C","#189529","#E0640B","#1539BC","#7F13A0","#43756F")

countZeroes = function(col) {
  sum(col==0)<7
}

##melt and take the upper triangle of the distance matrix
melt_dist <- function(dist, order = NULL, dist_name = 'dist') { ###melting distance matrix)
  if(!is.null(order)){
    dist <- dist[order, order]
  } else {
    order <- row.names(dist)
  }
  diag(dist) <- NA
  dist[upper.tri(dist)] <- NA
  dist_df <- as.data.frame(dist)
  dist_df$iso1 <- row.names(dist)
  dist_df <- dist_df %>%
    tidyr::gather_(key = "iso2", value = lazyeval::interp("dist_name", dist_name = as.name(dist_name)), order, na.rm = T)
  return(dist_df)
}


transformTaxa = function(taxa.data,cluster.data){
  # exclude taxa with too many zero-value entries
  taxa.data = taxa.data[,apply(taxa.data,2,countZeroes)]
  
  ##### log transform data
  # deal with zero values (0 = resolution limit)
  # log(x+0.000001)
  taxa.series = tslist(t(log10(taxa.data+0.000001)))
  
  # interpolate taxa series to length of barcode cluster series
  taxa.series <- reinterpolate(taxa.series, new.length = nrow(cluster.data))
  
  return(taxa.series)
}

###### ALL-AGAINST-ALL COMPARISON AND CLUSTERING
### shape-based distance analysis
runcoclustering = function(series, sample_name, nclust,cluster.series,taxa.series){
   ### depending of the silluotte analysis choose the best cluster number if necesaary
  
  cfg <- compare_clusterings_configs(
    "h",
    k = 2:(length(names(series))-1),
    preprocs = pdc_configs(
      type = "preproc",
      zscore = list()),
    controls = list(
      hierarchical = hierarchical_control(method = "average", symmetric = FALSE)
    ),
    distances = pdc_configs("d", sbd = list())
  )
  
  evaluator <- cvi_evaluators("Sil")
  comparison <- compare_clusterings(series, "h", cfg, 
                                    score.clus = evaluator$score,
                                    pick.clus = evaluator$pick,
                                    seed=456L)

  
  hc_sbd <- repeat_clustering(series, comparison, comparison$pick$config_id)
  png(paste(outdir,"/coclustering/",sample_name,"_clusters.png",sep=""),res = 300,width = 6,height = 4, units="in")
  plot(hc_sbd,type="sc")
  dev.off()

  ##save distance matrix cophenetic(hc_sbd) to calculate the mixing index later
  hc_dist = hc_sbd %>% cophenetic() %>% as.matrix() %>% 
    melt_dist() %>% mutate(id=paste(iso1,iso2,sep="_"))
  write_tsv(hc_dist, paste(outdir,"/coclustering/",sample_name,"_distance.tsv",sep=""),col_names=TRUE)
  
  sbd.clust <- data.frame(dtwclust = hc_sbd@cluster, series = hc_sbd$labels)
  sbd.clust = sbd.clust[order(sbd.clust$dtwclust),]
  write_tsv(sbd.clust, paste(outdir,"/coclustering/",sample_name,"_clusters.tsv",sep=""),col_names=TRUE)
  as.dendrogram(hc_sbd) -> dend
  
  ##cluster tree with cluster colors
  as.dendrogram(hc_sbd) -> dend
   n_sbd=names(hc_sbd@cluster)
  n_sbd=n_sbd[!(n_sbd %in% c("Enterobacteriaceae","Moraxellaceae","Desulfovibrionaceae","Lachnospiraceae","Ruminococcaceae",
                             "Peptostreptococcaceae","Peptococcaceae","Clostridiales_vadinBB60_group",
                             "Clostridiaceae_1","Paenibacillaceae","Lactobacillaceae","Muribaculaceae",
                             "Bacteroidaceae","Rikenellaceae","Prevotellaceae","Marinifilaceae",
                             "Tannerellaceae","Anaeroplasmataceae", "Deferribacteraceae", "Saccharimonadaceae",
                             "Xiphinematobacteraceae","Methanobacteriaceae","Erysipelotrichaceae","other","Acholeplasmataceae","Oscillospiraceae", "Akkermansiaceae","[Eubacterium] coprostanoligenes group",
                             "Enterococcaceae","Bacillaceae","Corynebacteriaceae","Microbacteriaceae","Erwiniaceae","Veillonellaceae","Sutterellaceae"))]
  names(cluster.colors)=n_sbd
  cluster.colors=cluster.colors[!is.na(names(cluster.colors))]
  labels_colors(dend)=cluster.colors[names(cluster.colors)][order.dendrogram(dend)]
  cairo_ps(paste0(outdir,"/coclustering/",sample_name,"_sbd_cluster_color_dendrograms.eps"),width = 8,height = 8)
  par(mar = c(4,2,2,8))
  plot(dend, horiz = T,main="Distance: SBD",xlab="Height") 
  dev.off()
}


####rm##
rm1_16S_taxa <- read_delim("data/16S/taxa_long_format/rm1_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm2_16S_taxa <- read_delim("data/16S/taxa_long_format/rm2_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm3_16S_taxa <- read_delim("data/16S/taxa_long_format//rm3_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm4_16S_taxa <- read_delim("data/16S/taxa_long_format//rm4_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)


outdir="reports/figures"



rm1_clustered_loess_log10 <- read_csv(paste0(outdir,"/clustering/rm1/rm1_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))
rm2_clustered_loess_log10 <- read_csv(paste0(outdir,"/clustering/rm2/rm2_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))
rm3_clustered_loess_log10 <- read_csv(paste0(outdir,"/clustering/rm3/rm3_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))
rm4_clustered_loess_log10 <- read_csv(paste0(outdir,"/clustering/rm4/rm4_clustered_loess_log10.csv"), 
                                      col_types = cols(time = col_skip()))


taxa.series = transformTaxa(rm1_16S_taxa,rm1_clustered_loess_log10)
cluster.series = tslist(t(rm1_clustered_loess_log10))
series = append(cluster.series,taxa.series)
runcoclustering(series,"rm1",4,cluster.series,taxa.series)

taxa.series = transformTaxa(rm2_16S_taxa,rm2_clustered_loess_log10)
cluster.series = tslist(t(rm2_clustered_loess_log10))
series = append(cluster.series,taxa.series)
runcoclustering(series,"rm2",4,cluster.series,taxa.series)

taxa.series = transformTaxa(rm3_16S_taxa,rm3_clustered_loess_log10)
cluster.series = tslist(t(rm3_clustered_loess_log10))
series = append(cluster.series,taxa.series)
runcoclustering(series,"rm3",4,cluster.series,taxa.series)

taxa.series = transformTaxa(rm4_16S_taxa,rm4_clustered_loess_log10)
cluster.series = tslist(t(rm4_clustered_loess_log10))
series = append(cluster.series,taxa.series)
runcoclustering (series,"rm4",4,cluster.series,taxa.series)

