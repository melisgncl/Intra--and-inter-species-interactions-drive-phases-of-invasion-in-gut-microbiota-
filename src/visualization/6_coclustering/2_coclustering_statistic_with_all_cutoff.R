


####first name of files 16s the second part of the files barcode dynamics
args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

source("/Users/melis/Dropbox/mouse_gut/github/src/visualization/0_config/0_config.R")



library(svglite)
library(gridExtra)
library(grid)
library("sm")
library(phylogram)
quartzFonts(arial = c("Arial", "Arial Black", "Arial Italic", "Arial Black Italic"))
par(family = "Arial")


taxa.file = args[1]
series.file = args[2]
sample = args[3]
condition = args[4]
cutoff=args[5]


`%notin%` <- Negate(`%in%`)



s_list=identity_list=c("Enterobacteriaceae",
                       "Moraxellaceae",
                       "Desulfovibrionaceae",
                       "Lachnospiraceae",
                       "Ruminococcaceae",
                       "Peptostreptococcaceae",
                       "Clostridiaceae_1",
                       "Clostridiales_vadinBB60_group",
                       "Peptococcaceae",
                       "Paenibacillaceae",
                       "Lactobacillaceae",
                       "Erysipelotrichaceae",
                       "Muribaculaceae",
                       "Bacteroidaceae",
                       "Rikenellaceae",
                       "Prevotellaceae",
                       "Marinifilaceae",
                       "Tannerellaceae",
                       "Anaeroplasmataceae",
                       "Deferribacteraceae",
                       "Saccharimonadaceae",
                       "Xiphinematobacteraceae",
                       "Methanobacteriaceae",
                       "other","Sutterellaceae","Acholeplasmataceae","Oscillospiraceae",
                       "Akkermansiaceae","[Eubacterium] coprostanoligenes group",
                       "Enterococcaceae" ,"Bacillaceae"  ,"Corynebacteriaceae",
                       "Microbacteriaceae" ,"Erwiniaceae","Veillonellaceae","Comamonadaceae","Xanthomonadaceae",
                       "Leptotrichiaceae","Pseudomonadaceae","Sphingomonadaceae","Weeksellaceae")

countZeroes = function(col) {
  sum(col==0)<7
}

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
runcoclustering = function(series, sample_name, nclust,cluster.series,taxa.series,condition,cutoff){
  
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
  
  ##save distance matrix cophenetic and name the pairs as specie-specie(s) , clone-specie(m) or clone-clone(c)
  hc_dist = hc_sbd %>% cophenetic() %>% as.matrix() %>% 
    melt_dist() %>% mutate(id=paste(iso1,iso2,sep="_")) %>% mutate(pair=case_when(iso1 %in% s_list & iso2 %in% s_list ~ "S",
                                                                                  iso1 %in% c_list & iso2 %in% c_list ~ "C",
                                                                                  iso1 %in% s_list & iso2 %in% c_list ~ "M"))
  hc_dist$Condition=condition
  clone=subset(hc_dist,hc_dist$pair=="C")
  species=subset(hc_dist,hc_dist$pair=="S")
  pairwise=subset(hc_dist,hc_dist$pair =="M")
  ###calculate the mxing index from their cumulative disturbution##
  ks_stat=ks.test(clone$dist,pairwise$dist)
  stat_df=c(ks_stat$statistic,ks_stat$p.value)
  names(stat_df)[1]="ks_stat"
  names(stat_df)[2]="p.value"
  stat_df$Condition=condition
  stat_df$sample_name=sample_name
  stat_df$cutoff=cutoff
  
  
  ##ks distance and distance matrix that comes from dendogram
  write_tsv(as.data.frame(stat_df),paste("reports/figures//mixing_index/stat_csv/",sample_name,"_",cutoff,"_stat.tsv",sep=""),col_names=TRUE)
  #write_tsv(hc_dist, paste("reports/figures/mixing_index/stat_csv/",sample_name,"_",cutoff ,"_distance.tsv",sep=""),col_names=TRUE)
  
  as.dendrogram(hc_sbd) -> dend
  
  n_sbd=names(hc_sbd@cluster)
  n_sbd=n_sbd[!(n_sbd %in% c("Enterobacteriaceae","Moraxellaceae","Desulfovibrionaceae","Lachnospiraceae","Ruminococcaceae",
                             "Peptostreptococcaceae","Peptococcaceae","Clostridiales_vadinBB60_group",
                             "Clostridiaceae_1","Paenibacillaceae","Lactobacillaceae","Muribaculaceae",
                             "Bacteroidaceae","Rikenellaceae","Prevotellaceae","Marinifilaceae",
                             "Tannerellaceae","Anaeroplasmataceae", "Deferribacteraceae", "Saccharimonadaceae",
                             "Xiphinematobacteraceae","Methanobacteriaceae","Erysipelotrichaceae","other","Acholeplasmataceae","Oscillospiraceae", "Akkermansiaceae","[Eubacterium] coprostanoligenes group",
                             "Enterococcaceae","Bacillaceae","Corynebacteriaceae","Microbacteriaceae","Erwiniaceae","Veillonellaceae","Sutterellaceae",
                             "Comamonadaceae","Xanthomonadaceae","Leptotrichiaceae","Pseudomonadaceae","Sphingomonadaceae","Weeksellaceae"))]
  
  
  
  cluster.colors=c("#3cb44b","#4363d8","#e6194B","#e8ca00","#911eb4","#f58231","#22766dff","#42d4f4","#f032e6","#9A6324",
                   "#2F4C39", "#1D3F6E","#94170f","#665679","#F17829","#97A69C","#606EA9","#A9606E","#A99060","#F8F1AE",
                   "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8","#003D18","#003D18","#82a479","#c74b0e",
                   "#77b5fe","#ccff00")
  names(cluster.colors)=n_sbd
  cluster.colors=cluster.colors[!is.na(names(cluster.colors))]
  labels_colors(dend)=cluster.colors[names(cluster.colors)][order.dendrogram(dend)]
  
  ###plot the dendograms##
  #cairo_ps(paste0("reports/figures/mixing_index/",condition,"/",sample_name,"/",cutoff, "_sbd_cluster_color_dendrograms.eps"),width = 6,height = 4)
  #par(mar = c(4,2,2,8))
  #plot(dend, horiz = T,main="Distance: SBD",xlab="Height") 
  #dev.off()
}


##taxa_data
taxa.data=read_delim(taxa.file, 
                     "\t", escape_double = FALSE)

taxa.data$Time=NULL

nrow(taxa.data)

##clone_data
cluster.data=read_csv(series.file)

cluster.data=cluster.data[cluster.data$time %in% (seq(min(cluster.data$time),nrow(taxa.data)*10)) ,]

cluster.data$time=NULL

c_list=colnames(cluster.data)

taxa.series = transformTaxa(taxa.data,cluster.data)
cluster.series = tslist(t(cluster.data))
series = append(cluster.series,taxa.series)
runcoclustering(series,sample,4,cluster.series,taxa.series,condition,cutoff)




