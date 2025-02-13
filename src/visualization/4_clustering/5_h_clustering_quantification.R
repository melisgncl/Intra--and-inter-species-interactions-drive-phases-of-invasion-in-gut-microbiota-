# Melis Gencel
#######hierarchical clustering quantification####
## Calculate Euclidean distances between clusters depending on cutoff ##
source("/Users/melis/Dropbox/mouse_gut/github/src/visualization/0_config/0_config.R")
library(TSdist)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Function to melt distance matrix
melt_dist <- function(dist, order = NULL, dist_name = 'dist') {
  if (!is.null(order)) {
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

# Function for hierarchical clustering statistics
hc_statistic <- function(sample, treshold, outdir) {
  mydir <- paste0(outdir, "/clustering_control/", sample)
  myfiles <- list.files(path = mydir, pattern = "*clustered_loess_log10.csv", full.names = TRUE, recursive = TRUE)
  df <- lapply(myfiles, function(x) {
    DF <- read_csv(x)
    if (ncol(DF) <= 2) {##remove ones with one cluster
      return(NULL)
    }
    DF$time <- NULL
    subject <- sapply(strsplit(x, paste0("/", sample)), `[`, 2)
    subject <- sapply(strsplit(subject, "/"), `[`, 2)
    cutoff_value <- as.numeric(subject)
    print(paste("Parsed cutoff value from file", x, ":", cutoff_value))  # Diagnostic print statement
    DF$cutoff <- cutoff_value
    return(DF)
  })
  
  # Remove NULL elements (data frames with only two columns)
  df <- df[!sapply(df, is.null)]
  
  if (length(df) == 0) {
    cat("No valid data frames found in directory:", mydir, "\n")
    return(NULL)
  }

  cutoff <- unique(as.numeric(unlist(do.call(rbind, lapply(df, "[", , "cutoff")))))
  
  distance_pairwise <- sapply(df, function(x) dist(t(x[1:(length(x) - 1)]), method = "TSDistances", distance = "euclidean"))
  
  distance_pairwise <- lapply(distance_pairwise, function(x) melt_dist(as.matrix(x)))
  
 
  
  tf <- do.call(rbind, (map2(distance_pairwise, cutoff, ~cbind(.x, cutoff = .y))))
  
  tf <- tf %>% mutate(cluster = as.numeric(gsub("C", "", tf$iso1))) %>% group_by(cutoff) %>% mutate(cluster = max(cluster), id = paste(iso1, iso2, sep = "_"))
  
  tf2 <- tf %>% group_by(cutoff) %>% summarise(dist_small = min(dist), cluster = max(cluster)) %>% filter(dist_small <= 100)
  
  #p1 <- ggplot(tf2, aes(dist_small, cluster, color = cutoff)) + geom_violin(size = 2.5) +
    #theme_Publication()
  
  scale <- max(tf2$dist_small) / max(tf2$cluster)
  
  p2 <- ggplot(tf2, aes(cutoff, dist_small)) + geom_line(color = "black", size = 2.5) +
    theme_Publication() +
    geom_line(aes(cutoff, as.integer(cluster) * scale), size = 2, color = "darkblue") + 
    scale_y_continuous(sec.axis = sec_axis(~./scale, name = "Cluster number")) +
    scale_x_reverse() + 
    xlab("Threshold") + 
    geom_vline(xintercept = treshold, linetype = "solid", color = "darkred")
  
  ggsave(p2, filename = paste0(outdir, "/clustering_control/", sample, "_small_dist.eps"), width = 8.25, height = 6, limitsize = FALSE)
  
  return(p2)
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the output directory as a command-line argument.")
}




#outdir <- args[1]
outdir="reports/figures/"
#outdir="reports/regex_Delcorrect//"



rm1 <- hc_statistic("rm1", 0.85, outdir)
rm2 <- hc_statistic("rm2", 0.54, outdir)
rm3 <- hc_statistic("rm3", 0.63, outdir)
rm4 <- hc_statistic("rm4", 0.43, outdir)

gf1 <- hc_statistic("gf1", 0.22, outdir)
gf2 <- hc_statistic("gf2", 0.45, outdir)
gf3 <- hc_statistic("gf3", 0.38, outdir)
gf4 <- hc_statistic("gf4", 0.34, outdir)

cowplot::plot_grid(rm1, rm2, rm3, rm4, gf1, gf2, gf3, gf4, align = "hv", rel_widths = c(1, 1, 1, 1, 1, 1, 1), rel_heights = c(1, 1, 1, 1, 1, 1, 1), axis = c("tblr"), nrow = 2) %>%
  ggsave(filename = paste0(outdir, "/clustering_control/all_stat_graph_smallest_distance.eps"), width = 30, height = 12)


