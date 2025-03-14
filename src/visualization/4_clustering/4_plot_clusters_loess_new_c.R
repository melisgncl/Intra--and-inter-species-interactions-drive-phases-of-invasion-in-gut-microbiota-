#!/usr/bin/env Rscript

## R script
## Melis Gencel & Louis Gauthier

#######################    READ ME!   #####################################

### Here we plot the clusters ordered by  their last day average frequecny
###	We also compute and extract the loess fit for each cluster.
### We with log-transformed frequencies.

#######################    ^^^^^^^^   #####################################

# ex usage: Rscript src/visualization/4_clustering/4_plot_clusters.R data/clustering/GFM1-4/clusters_GFM2_average data/clustering/GFM1-4/GFM2_filtered.csv reports/figures/clustering/GFM2
source("/Users/melis/Dropbox/mouse_gut/github/src/visualization/0_config/0_config.R")
args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")


library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(hash)
library(broom)
library(dplyr)

#library(extrafont)
#extrafont::font_import(pattern="Arial",prompt = FALSE)
#extrafont::loadfonts()

clusters.file = args[1]
series.file = args[2]
sample = args[3]
type = args[4]






##############

plotClusterLog10 <- function(df, cluster, sample, color) {
  span <- 0.2  # Initial value of span
  max_span <- 0.8  # Maximum value of span
  success <- FALSE  # Flag to track successful plot generation
  
  while (span <= max_span && !success) {
    p <- ggplot(df, aes(x = variable, y = log10(value + 0.0000001))) +
      geom_line(aes(group = ID), color = color, alpha = 0.4) +
      theme_Publication(base_size = 18) +
      ylim(-7, 0) +
      labs(x = "Time (post-gavage)", y = "log10(Barcode frequency)") +
      guides(color = FALSE) +
      scale_x_continuous(breaks = effective.breaks, labels = effective.labels, limits = effective.limits) +
      coord_cartesian(expand = FALSE) +
      stat_smooth(aes(outfit = fit1 <<- ..y..), n = loess.range, se = FALSE, size = 1, color = "black", span = span, method = "loess") +
      annotate("text", y = -0.75, x = 3, label = paste("n", length(unique(df$ID)), sep = " = "), hjust = 0, size = 5) +
      ggtitle(paste("Cluster", cluster, sep = " "))
    
    # Try saving the plot, if successful, set success flag to True
    tryCatch(
      {
        ggsave(filename = paste(sample, "_cluster", cluster, "_log10.eps", sep = ""), width = 5.5, height = 4, device = cairo_ps)
        success <- TRUE
      },
      error = function(e) {
        span <- span + 0.01  # Increase the value of span
      }
    )
  }
  
  if (success) {
    return(list(p, fit1))
  } else {
    stop("No successful plot was generated.")
  }
}


####hash.tables###
breaks.hash <- hash()
breaks.hash[["gf"]] = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
breaks.hash[["rm"]] =  c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
breaks.hash[["im"]] =  c(1,2,3,4,5,6,7,8,9)

# create hash table of breaks
labels.hash <- hash()
labels.hash[["gf"]] = c("3h","","12h","","2d","","","","6d","","","","10d","","","","14d")
labels.hash[["rm"]] = c("3h","","12h","","2d","","","","6d","","","","10d","","","","14d","")
labels.hash[["im"]] = c("3h","","12h","","2d","","","","6d")

# create hash table of breaks
limits.hash <- hash()
limits.hash[["gf"]] = c(1,17)
limits.hash[["rm"]] = c(1,18)
limits.hash[["im"]] = c(1,9)

##threshold.hash
cut.hash <- hash()
cut.hash[["gf"]] = 8
cut.hash[["rm"]] = 8
cut.hash[["im"]] = 5

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]
cut=cut.hash[[type]]






clusters <- read_csv(file = clusters.file, col_names = FALSE)[1]
colnames(clusters)=c("cluster")
nRank = nrow(clusters)
clusters$rank = seq(1:nRank)


series.filtered <- read_csv(file = series.file, col_names = TRUE)
series.filtered$points = NULL
series.filtered$rank=seq(1:nRank)

grouped.clusters = merge(series.filtered,clusters,by.x = "rank",by.y = "rank")
series.reshaped = reshape2::melt(grouped.clusters,id.vars = c('ID','cluster','rank','mean'))
series.reshaped$variable = as.numeric(as.character(series.reshaped$variable))
print(colnames(series.reshaped))
####filter clusters that has less than 8 members###
series.reshaped.1=series.reshaped %>%  group_by(cluster) %>% dplyr::filter(length(unique(ID)) >= cut)


if(nrow(series.reshaped.1) != nrow(series.reshaped)){
  
  #cat("You might be ignoring dominant barcodes, please enter a minimum mean frequency that must be reached by at least one of the members for the cluster to be considered: ")
  minimum_freq <-1e-03
  
  series.reshaped.2 = series.reshaped %>%  dplyr::group_by(cluster) %>% dplyr::filter(length(unique(ID)) < cut) %>% dplyr::mutate(mean_freq = mean(value)) %>% dplyr::filter(mean_freq >= minimum_freq)
  series.reshaped.2$mean_freq = NULL
  
  series.reshaped = rbind(series.reshaped.1, series.reshaped.2)
  
}





series_order=subset(series.reshaped,series.reshaped$variable==max(unique(series.reshaped$variable)))
####order clusters accoridng to their average  ending frequency###
series_order=series_order %>% 
  group_by(cluster) %>% 
  summarise(average = mean(value)) 

series_order=series_order[order(series_order$average,decreasing = TRUE), ] 
series_order$cluster_ranked=rownames(series_order)
series.reshaped=merge(series.reshaped,series_order,by="cluster")
series.reshaped$cluster=NULL
names(series.reshaped)[7]="cluster"

clusters.ranked=series_order$cluster_ranked

cluster.colors=c("#3cb44b","#4363d8","#e6194B","#e8ca00","#911eb4","#f58231","#22766dff","#42d4f4","#f032e6","#9A6324",
                 "#2F4C39", "#1D3F6E","#94170f","#665679","#F17829","#97A69C","#606EA9","#A9606E","#A99060","#F8F1AE",
                 "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8","#003D18","#003D18","#82a479","#c74b0e",
                 "#77b5fe","#ccff00")

max.range = max(series.reshaped$variable)-min(series.reshaped$variable)
loess.range = (max.range*10)+1


###########################
### log10-transformed loess
###

plotList = list()

fit.list <- list()

for(i in seq_along(unique(series.reshaped$cluster))){
  l = plotClusterLog10(series.reshaped[series.reshaped$cluster==clusters.ranked[i],], clusters.ranked[i], sample, cluster.colors[i])
  
  p = l[1]
  fit = l[2]
  plotList = c(plotList,p)
  fit.list = c(fit.list,fit)
}

top10.log =cowplot::plot_grid(plotlist=plotList, align = "hv",axis=c("tblr"))
ggsave(top10.log,filename = paste(sample, "_clusters_log10.eps", sep=""),width = 18,height=18,device = cairo_ps)
write_csv(series.reshaped,file = paste(sample, "_clustered_series_log10.csv",sep=""),col_names = TRUE)

clusters.loess = as.data.frame(fit.list[1])
colnames(clusters.loess)=paste("C",clusters.ranked[1],sep="")

for(fit in seq_along(fit.list[c(-1)])){
  
  clusters.loess = cbind(clusters.loess,as.data.frame(fit.list[fit+1]))
  colnames(clusters.loess)[ncol(clusters.loess)] = paste("C",clusters.ranked[fit+1],sep="")
}

clusters.loess$time= seq(from=min(series.reshaped$variable)*10, to=(min(series.reshaped$variable)*10 +nrow(clusters.loess) -1))

cluster.df = reshape2::melt(clusters.loess,id="time") ##arranged time according to experiment


loess.plot = ggplot(cluster.df) + geom_line(aes(time/10,value,group=variable,color=variable),size=1) + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  theme_Publication() + scale_color_manual(values = cluster.colors,name="cluster") + ylab("fit") + xlab("Time (post-gavage)") + ylim(-7,-1)
ggsave(loess.plot,filename = paste(sample, "_loess_clusters_log10.eps", sep=""),width = 8.25,height = 6,device = cairo_ps)

loess.plot2= loess.plot+ theme(legend.position = "none") + coord_cartesian(expand = FALSE)
ggsave(loess.plot2,filename = paste(sample, "_loess_clusters_log10_without_legend.eps", sep=""),width = 8.25,height = 6,device = cairo_ps)


write_csv(clusters.loess,file = paste(sample, "_clustered_loess_log10.csv",sep=""),col_names = TRUE)
