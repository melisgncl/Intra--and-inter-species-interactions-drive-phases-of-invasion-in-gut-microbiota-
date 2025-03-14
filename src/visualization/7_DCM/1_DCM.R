###  Melis Gencel
### DCM analysis
#######################    READ ME!   #####################################

##### all the function that is necessary for DCM ###

#######################    ^^^^^^^^   #####################################

##source main function from 0_config

source("/Users/melis/Dropbox/mouse_gut/github/src//visualization/0_config/0_config.R")
library(tidyr)
library(ggforce)
library(kernlab)
library(zoo)
library(changepoint)
library(ggplot2)
library(stats)
library(changepoint.geo)




set.seed(170)

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

countZeroes = function(col,n) {
  ##7 for rm and im while 4 for nc
  sum(col==0)<n
} 

##use if you have time points are missed in between
interpolate_nas <- function(df) {
  # Replace all zeros with NA
  df[df == 0] <- NA
  
  # Apply linear interpolation to fill NA values
  df <- apply(df, 2, function(x) {
    # Perform interpolation
    interpolated_x = na.approx(x, na.rm = FALSE)
    
    # Check if the last value is NA and replace it if necessary
    if (is.na(interpolated_x[length(interpolated_x)])) {
      interpolated_x[length(interpolated_x)] <- interpolated_x[length(interpolated_x) - 1]
    }
    
    return(interpolated_x)
  })
  
  # Convert any remaining NAs to zeros
  df[is.na(df)] <- 0
  
  return(df)
}

##filter the values if the freq i lower than 10^-4
filter_columns_by_mean <- function(data) {
  # Calculate the mean of each column
  means <- colMeans(data)
  
  # Filter columns where the mean is less than 10^-4
  filtered_data <- data[, means >= 5e-3]
  
  # Return the filtered data frame
  return(filtered_data)
}



####linear interpolation family time series like their clonal clusters
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

### calculate the  ż for invaded cohorts  with clonal clusters and family compostions
get_time_derivative_invaded = function(taxa.data,sample,span,n){
  ##cluster data load depending on chotost 
  clusters.loess=read_csv(file= paste(outdir,"/clustering/",sample,"/",sample,"_clustered_loess_log10.csv",sep=""))
  clusters=clusters.loess
  clusters$time=NULL
  # Remove zero counts from taxa data and scale by log10
  ##'n' specifies the maximum allowable number of zeros in a column for it to be included in the analysis
  taxa.data<- taxa.data[, apply(taxa.data, 2, countZeroes,n)]
  ###remove the data points from taxa data if they are not in cluster.loeesd data
  taxa.data=taxa.data[c((min(clusters.loess$time)/10):(max(clusters.loess$time)/10)), ]
  ##### log transform data
  # deal with zero values (0 = resolution limit)
  # log(x+0.000001)
  taxa.data = tslist(t(log10(taxa.data+0.000001)))
  ####Interpolate points depending on real time data points
  Time=seq(from=min(clusters.loess$time)/10, to=max(clusters.loess$time)/10)
  xx <- seq(from=min(clusters.loess$time)/10, to=max(clusters.loess$time)/10, length.out=(length(clusters.loess$time)))
  # Apply loess smoothing
  loess_models_taxa=sapply(taxa.data, function(x) predict(loess(x~Time,span=span),xx,se = FALSE))
  complete_loess_series <- append(clusters,tslist(t(loess_models_taxa)))
  # Calculate time derivative using central differences
  derivative_series=sapply(complete_loess_series, function(y)  tslist(t((diff(y,lag=2)/diff(xx,lag=2)))))
  return(list(complete_loess_series,derivative_series,clusters.loess$time))
}


### calculate the  ż for lineages
get_time_derivative_lineage=function(sample){
  complete_loess_series=read_csv(file= paste(outdir,"/clustering/",sample,"/",sample,"_clustered_loess_log10.csv",sep=""))
  xx <- seq(from=min(complete_loess_series$time)/10, to=max(complete_loess_series$time)/10, length.out=(length(complete_loess_series$time)))
  # Calculate time derivative using central differences
  complete_loess_series$time=NULL
  derivative_series=sapply(complete_loess_series, function(y)  tslist(t((diff(y,lag=2)/diff(xx,lag=2)))))
  options(scipen=0)
  return(list(tslist(t(complete_loess_series)),derivative_series))
}

### calculate the  ż for taxa only
get_time_derivative_taxa = function(taxa.data,span,n){
  ###modeled taxa##polynomial interpolation
  taxa.data = taxa.data[,apply(taxa.data,2,countZeroes,n)]
  #Use this if there are missing time points in your data sequence.
  taxa.data=interpolate_nas(taxa.data)
  
  ##we know all 16s has their first points
  Time=seq(from=1, to=nrow(taxa.data))
  ##just to match with cluster.loees time otherwise you can increase the length as much you can 
  xx <- seq(from=1, to=nrow(taxa.data), length.out=nrow(taxa.data)*10-9)
  ###we use all the time points from 16S
  taxa.data = tslist(t(log10(taxa.data+0.000001)))
  
  complete_loess_series  = tslist(t(sapply(taxa.data, function(x) predict(loess(x~Time,span=span),xx,se = FALSE))))
  derivative_series=sapply(complete_loess_series  , function(y)  tslist(t((diff(y,lag=2)/diff(xx,lag=2)))))
  return(list(complete_loess_series,derivative_series))

}

###calculate time dependent jacobian

####to get the time-progressive jacobians and its eigen_values and vectors##
calculate_time_dependent_jacobian=function(complete_loess_series,derivative_series,window,sample){
  options(scipen=0)
  total <- nrow(as.data.table(complete_loess_series))
  # Generate sequence points from the first time point to the last, adjusting for the window size
  # Windows are defined to span from one day to two days so we hardcode for experimental part
  spots <- seq(from=10, to=(total), by=window)
  spots[length(spots)]=nrow(as.data.table(derivative_series))
  print(spots)
  # Adjust the window by reducing 2 steps to center it appropriately
  # Calculate the community interaction matrix for each window
  
  jacobian=list()  
  for(k in 1:length(spots)){
    jacobian[[k]] <- sapply(1:ncol(as.data.table(derivative_series)), function(i) {
      sapply(1:ncol(as.data.table(complete_loess_series)), function(j){
        resTmp <- cov(derivative_series[[i]][(1:spots[k])],complete_loess_series[[j]][(1:spots[k])]) 
        
      })
    })
  }
 
   # Assign row and column names from derivative_series to each matrix in results
  for(i in seq_along(jacobian)) {
    rownames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
    colnames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
  }
  
  
  return(jacobian)
  
  }
  

get_eigen_values=function(jacobian,sample,time_increase=FALSE){
  
  # Calculate eigenvalues and eigenvectors for each jacobian matrix
  eigen_values=sapply(jacobian, function(y) eigen(y, symmetric= FALSE)$values)
  # Extract real and imaginary parts of eigenvalues
  real_eigen_value=sapply(as.data.table(eigen_values), function(y) Re(y))
  im_eigen_value=sapply(as.data.table(eigen_values), function(y) Im(y))
 
  eigen_vector=lapply(jacobian, function(y) eigen(y, symmetric= FALSE)$vector)
  
  #Calculate the magnitude of eigenvalues for sorting
  conjugate= sqrt(real_eigen_value^2+im_eigen_value^2)
  
  # Organize data into data frames for further manipulation
  real_eigen_value=melt(real_eigen_value)  %>% mutate(eigen="R")
  im_eigen_value=melt(im_eigen_value)  %>% mutate(eigen="I")
  conjugate_value=melt(conjugate)  %>% mutate(eigen="C")
  
  ordered_eigen_values =rbind(real_eigen_value,im_eigen_value,conjugate_value) %>% spread(eigen,value)  %>% 
    mutate(Var2=sub("V", "", Var2) ,Var1=as.factor(Var1),Var2=(as.numeric(Var2)))  %>% group_by(Var2) %>%
    arrange(desc(C),.by_group = TRUE) %>%
    mutate(Var3=as.factor(row_number()))  
  
  ###prepare the data for dimension reduction by using time dependent eigenvalue 
  time_dependent_eigenspace= ordered_eigen_values %>% mutate(sample=sample) %>% 
                            select(-C,-Var1)  %>% pivot_longer(cols = -c(Var2,Var3,sample), names_to = "variable", values_to = "value") %>% 
                            group_by(Var2) %>% mutate(Var4=paste0(variable,"_",Var3)) %>%  select(-Var3,-variable) %>% spread(Var4,value) 
  ##if the first time not exist
  if (time_increase) {
    time_dependent_eigenspace$Var2 <- time_dependent_eigenspace$Var2 +1
  }
  
  
  return(list(ordered_eigen_values,time_dependent_eigenspace))
  
  }


###plot the eigenspace evolution###
plot_eigen_value_evolution=function(all_eigen_value,sample){
  eigen_space= sort(as.numeric(unique(all_eigen_value$Var2)))
  plots_eigen_value=list()
  for(i in seq_along(eigen_space)){
    plots_eigen_value[[i]] = ggplot(all_eigen_value[all_eigen_value$Var2 %in% c(1:eigen_space[i]) ,],aes(R,I,color=Var3))+
      geom_point(size=5,alpha=0.3) + 
      geom_path(size=1.2, aes(color=Var3), arrow = arrow(angle = 15, type = "closed"), alpha=0.5) +
      geom_point(aes(R,I,color=Var3),size=4,data=all_eigen_value[all_eigen_value$Var2==eigen_space[i],]) +
      theme_Publication() +xlab(~ paste("Re", "(",lambda [i],")"))+ ylab(~ paste("lm", "(",lambda [i],")")) +
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0) + xlim(min(all_eigen_value$R),max(all_eigen_value$R)) +
      ylim(min(all_eigen_value$I),max(all_eigen_value$I)) + guides(color=FALSE) +scale_color_manual(values = eigen.colors)
    ##save the all plots
    #ggsave(plots_eigen_value[[i]],filename = paste(outdir,"/DCM/eigen_values/",sample, "_eigen_value_",eigen_space[i],"_time_interval.eps", sep=""),width =10.5,height =7.5, limitsize = FALSE,device = cairo_ps )
  }
  ##save them as a panel
  cowplot::plot_grid(plotlist = plots_eigen_value,align = "hv",ncol = 10) %>%  
    ggsave(filename = paste(outdir,"/DCM2/eigen_values_panel/",sample, "_all_eigen_value_time_interval.eps", sep=""),width =50,height =20, limitsize = FALSE,device = cairo_ps)
}

plot_kpca_for_sample <- function(df, sigma, sample) {
  set.seed(172)
  plotlist <- list()
  kpca_all <- list()
  
  # Check if df is a data frame (including grouped_df) and convert to list if necessary
  if (is.data.frame(df)) {
    df <- list(df)
  } else if (!is.list(df)) {
    stop("'df' must be either a data frame or a list of data frames")
  }
  
  # Ensure sigma is a vector; if not, convert it into a vector
  if (!is.vector(sigma)) {
    sigma <- as.vector(sigma)
  }
  
  # Validate that the length of sigma matches the number of data frames
  if (length(sigma) != length(df)) {
    stop("The length of 'sigma' must match the number of data frames in 'df'")
  }
  
  for (i in seq_along(df)) {
    df_m <- df[[i]]
    sigma_c <- sigma[i]
    print(sigma_c)
    
    # Exclude columns "Var2" and "sample" to retain only numeric columns
    df_numeric <- df_m[, !names(df_m) %in% c("Var2", "sample")]
    
    # Center the dataset for PCA without scaling
    df_numeric <- as.data.frame(scale(df_numeric, center = TRUE, scale = FALSE))
    
    # Perform Kernel PCA
    kpca_result <- kpca(~., data = df_numeric, kernel = "rbfdot", kpar = list(sigma = sigma_c), features = 2)
    
    # Combine the KPCA results with the original "Var2" and "sample" columns
    final <- cbind(kpca_result@rotated, df_m[, c("Var2", "sample")])
    colnames(final) <- c('kpca1', 'kpca2', 'phase', 'mouse')
    
    # Create the plot
    plotlist[[i]] <- ggplot(final) +
      geom_point(size = 4, aes(color = phase, x = phase, y = kpca1)) +
      geom_line(size = 1.2, aes(x = phase, y = kpca1), color = "red") +
      geom_point(size = 4, aes(color = phase, x = phase, y = kpca2)) +
      geom_line(size = 1.2, aes(x = phase, y = kpca2), color = "blue") +
      theme_Publication()
    
    kpca_all[[i]] <- final
  }
  
  # Combine all plots into a grid and save the output
  cowplot::plot_grid(plotlist = plotlist, align = "hv", ncol = 4)
  ggsave(filename = paste0(outdir, "/DCM2/pca_", sample, ".eps"), width = 40, height = 8, limitsize = FALSE)
  
  return(list(plotlist, kpca_all))
}

find_significant_changepoints_penalties <- function(data, penalties = c("MBIC", "SIC", "BIC", "Hannan-Quinn"),
                                                    pen_value = 10, 
                                                    test_stat = "Empirical", 
                                                    nquantiles_range = 2:ceiling(nrow(data)), 
                                                    freq_threshold = 0,       flip_kpca1 = FALSE, 
                                                    flip_kpca2 = FALSE) {
  
  library(fitdistrplus)
  library(changepoint.geo)
  
  
  # Optionally flip the PCA components
  if (flip_kpca1) {
    data$kpca1 <- -data$kpca1
  }
  if (flip_kpca2) {
    data$kpca2 <- -data$kpca2
  }
  
  # Initialize a list to store results for each penalty
  penalty_results <- list()
  
  # Loop through each penalty type
  for (penalty in penalties) {
    changepoints_list <- lapply(nquantiles_range, function(nq) {
      geo_result <- geomcp(
        data[, c("kpca1", "kpca2")],
        penalty = penalty,
        pen.value = pen_value,
        test.stat = test_stat,
        nquantiles = nq
      )
      unique(sort(c(geo_result@ang.cpts, geo_result@dist.cpts)))
    })
    
    all_changepoints <- unlist(changepoints_list)
    changepoint_freq <- table(all_changepoints)
    freq_df <- data.frame(
      Changepoint = as.numeric(names(changepoint_freq)),
      Frequency = as.numeric(changepoint_freq)
    )
    freq_df$Score <- freq_df$Frequency / max(freq_df$Frequency)
    significant_points <- freq_df[freq_df$Frequency >= freq_threshold, ]
    ranked_changepoints <- significant_points[order(-significant_points$Frequency), ]
    penalty_results[[penalty]] <- ranked_changepoints
  }
  
  pooled_results <- do.call(rbind, lapply(names(penalty_results), function(penalty) {
    cbind(penalty_results[[penalty]], Penalty = penalty)
  }))
  overall_results <- aggregate(
    cbind(Frequency, Score) ~ Changepoint, 
    data = pooled_results, 
    FUN = mean
  )
  overall_results <- overall_results[order(-overall_results$Frequency), ]
  
  # Percentile-Based Cutoff
  percentile_cutoff <- quantile(overall_results$Frequency, 0.8)
  percentile_significant <- subset(overall_results, Frequency >= percentile_cutoff)
  
  # Z-Score Threshold
  overall_results$z_score <- (overall_results$Frequency - mean(overall_results$Frequency)) / sd(overall_results$Frequency)
  z_score_cutoff <- 1
  z_score_significant <- subset(overall_results, z_score > z_score_cutoff)
  
  # Distribution Fitting
  fit <- fitdist(overall_results$Frequency, "norm")
  dist_cutoff <- qnorm(0.8, mean = fit$estimate["mean"], sd = fit$estimate["sd"])
  dist_significant <- subset(overall_results, Frequency >= dist_cutoff)
  
  # Manual Elbow Method
  find_elbow <- function(scores) {
    x <- 1:length(scores)
    y <- sort(scores, decreasing = TRUE)
    dists <- sapply(1:length(x), function(i) {
      a <- c(x[1], y[1])
      b <- c(x[length(x)], y[length(y)])
      c <- c(x[i], y[i])
      dist_ab <- sqrt(sum((a - b)^2))
      dist_ac <- sqrt(sum((a - c)^2))
      dist_cb <- sqrt(sum((c - b)^2))
      area <- abs((a[1] * (b[2] - c[2]) + b[1] * (c[2] - a[2]) + c[1] * (a[2] - b[2])) / 2)
      height <- (2 * area) / dist_ab
      return(height)
    })
    return(which.max(dists))
  }
  sorted_frequency <- sort(overall_results$Score, decreasing = TRUE)
  elbow_point <- find_elbow(sorted_frequency)
  elbow_cutoff <- sorted_frequency[elbow_point]
  elbow_significant <- subset(overall_results, Frequency >= elbow_cutoff)
  
  plot(overall_results$Score) 
  abline(v = elbow_point, col = "red", lty = 2, lwd = 2)
  points(elbow_point, elbow_cutoff, col = "red", pch = 19, cex = 1.5)
  
  
  return(list(
    PenaltyResults = penalty_results,
    OverallResults = overall_results,
    PercentileSignificant = percentile_significant,
    ZScoreSignificant = z_score_significant,
    DistributionSignificant = dist_significant,
    ElbowSignificant = elbow_significant
  ))
  
  
}








