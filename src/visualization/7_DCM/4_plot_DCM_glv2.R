###  Melis Gencel
### DCM analysis
#######################    READ ME!   #####################################

##### In here we analyzesimulated glv data

#######################    ^^^^^^^^   #####################################

source("/Users/melis/Dropbox/mouse_gut/github/src/visualization/0_config/0_config.R")
#################################################im########################################################
#################################################im########################################################
###functions for GLV simulation
library(deSolve) # integrate ODEs
library(tidyverse) # plotting and wrangling
GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-8] <- 0 # prevent numerical problems
    dxdt <- x * (r + A %*% x)
    return(list(dxdt))
  })
}
# function to plot output
plot_ODE_output <- function(out,n_species){
  out <- as.data.frame(out)
  #out <- out/ rowSums(out)
  colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
  out <- as_tibble(out) %>% gather(species, density, -time)  %>% group_by(time) #%>%
  #mutate(normalized_density = density / sum(density))
  pl <- ggplot(data = out) + 
    aes(x = time, y = density, colour = species) + 
    geom_line(size=1.2) +theme_Publication() 
  show(pl)
  #ggsave(pl,filename=paste0("reports/figures/GLV/",n_species,"_species_dynamics.eps"))
  return(list(out,pl))
}
integrate_GLV <- function(r, A, x0, maxtime = 100, steptime = 0.01,n_species){
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- as.data.frame(ode(y = x0, times = times, 
                           func = GLV, parms = parameters, maxsteps = 10^9,
                           method = "ode45"))
  print(out)
  out <- plot_ODE_output(out,n_species)
  out[[1]]=spread(out[[1]], species, density)
  return(list(out[[1]],out[[2]]))
}
get_time_derivative_glv=function(clusters.loess){
  time=clusters.loess$time
  clusters.loess$time=NULL
  Time=seq(1:nrow(clusters.loess))
  series_whole_loess=clusters.loess
  deriv_series_loess=sapply(series_whole_loess, function(y)  tslist(t((diff(y,lag=2)/diff(time,lag=2)))))
  options(scipen=0)
  return(list(tslist(t(series_whole_loess)),deriv_series_loess,Time))
  
}

get_eigen_values=function(jacobian,sample,time_increase=FALSE){
  
  
  eigen_values=sapply(jacobian, function(y) eigen(y, symmetric= FALSE)$values)
  # Extract real and imaginary parts of eigenvalues
  real_eigen_value=sapply(as.data.table(eigen_values), function(y) Re(y))
  im_eigen_value=sapply(as.data.table(eigen_values), function(y) Im(y))
  
  eigen_vector=lapply(jacobian, function(y) eigen(y, symmetric= FALSE)$vector)
  
  #Calculate the magnitude of eigenvalues for sorting
  conjugate= sqrt(real_eigen_value^2+im_eigen_value^2)
  
  # Organize data into data frames for further manipulation
  real_eigen_value=reshape2::melt(real_eigen_value)  %>% mutate(eigen="R")
  im_eigen_value=reshape2::melt(im_eigen_value)  %>% mutate(eigen="I")
  conjugate_value=reshape2::melt(conjugate)  %>% mutate(eigen="C")
  
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
calculate_time_dependent_jacobian=function(complete_loess_series,derivative_series,window,start=1){
  options(scipen=0)
  total <- nrow(as.data.table(derivative_series))
  # Generate sequence points from the first time point to the last, adjusting for the window size
  # Windows are defined to span from one day to two days
  spots <- seq(from=start, to=(total), by=window)
  if (spots[length(spots)] != total) {
    spots <- c(spots, total)
  }
  print(spots)
  # Adjust the window by reducing 2 steps to center it appropriately
  # Calculate the community interaction matrix for each window
  
  jacobian=list()  
  for(k in 2:length(spots)){
    jacobian[[k]] <- sapply(1:ncol(as.data.table(derivative_series)), function(i) {
      sapply(1:ncol(as.data.table(complete_loess_series)), function(j){
        resTmp <- cov(derivative_series[[i]][(spots[1]:spots[k])],complete_loess_series[[j]][(spots[1]:spots[k])]) 
        
      })
    })
  }
  jacobian <- jacobian[-1]
  # Assign row and column names from derivative_series to each matrix in results
  for(i in seq_along(jacobian)) {
    rownames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
    colnames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
  }
  
  
  return(jacobian)
  
}
plot_kpca_for_sample <- function(df, sigma, sample) {
  set.seed(172)
  plotlist <- list()
  kpca_all <- list()
  library(kernlab)
  
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


plot_correlation_and_pvalues=function(df,A,steptime){
  library(reshape2)
  melted_list <- lapply(df, melt)
  # Add a column for list number
  melted_list <- lapply(seq_along(melted_list), function(i) {
    melted_data <- melted_list[[i]]
    melted_data$list_number <- i
    melted_data
  })
  # Combine all melted dataframes into one
  combined_df <- do.call(rbind, melted_list)
  # Rename the columns for clarity
  colnames(combined_df) <- c("V1", "V2", "Value","List_Number")
  combined_df$sp=paste0(combined_df$V2,"_",combined_df$V1)
  #####M calculated by (Dx)A
  #M=melt(t(diag(res_2_2$average) %*%  A_2))
  M=reshape2::melt(t(A))
  M$sp=paste0(M$Var2,"_",M$Var1)
  M=M[c(3,4)]
  ####corellation
  corr <- combined_df %>%
    group_by(List_Number) %>%
    summarize(correlation = cor.test(M$value, Value)$estimate)
  corr$steptime=steptime
  pvalue= combined_df %>%
    group_by(List_Number) %>%
    summarize(pvalue = cor.test(M$value, Value)$p.value)
  pvalue$steptime=steptime
  ggplot()+geom_line(data=corr,aes(List_Number,correlation)) +theme_Publication() +xlab("Progressive Time interval")
  ggplot()+geom_line(data=pvalue,aes(List_Number,-log10(pvalue))) +theme_Publication() +xlab("Progressive Time interval")
  
  return(list(corr,pvalue,M,combined_df))
  
}
matrix_correlation <- function(A, B) {
  A_flat <- as.vector(A)
  B_flat <- as.vector(B)
  return(cor(A_flat, B_flat))
}

perturb_and_solve2 <- function(A_3, r_3, seed, max_attempts = 1000, tolerance = 1e-6, target_range = c(4, 25),
                               correlation_range = c(0, 0.6)) {
  set.seed(seed)  # To make the perturbation reproducible
  attempts <- 0
  feasible_solution_found <- FALSE
  
  while (attempts < max_attempts && !feasible_solution_found) {
    attempts <- attempts + 1
    
    # Perturb A_3 slightly
    A_3_perturbed <- A_3 + matrix(runif(nrow(A_3) * ncol(A_3), -0.1, 0.1), nrow = nrow(A_3), ncol = ncol(A_3))
    
    # Perturb r_3 slightly but ensure non-negative values
    r_3_perturbed <- pmax(r_3 + runif(length(r_3), -0.5, 0.5), 0)  # pmax ensures r_3_perturbed >= 0
    
    # Calculate the correlation between the original and perturbed matrix
    corr_value <- matrix_correlation(A_3, A_3_perturbed)
    
    # Ensure correlation is within the specified range
    if (corr_value >= correlation_range[1] && corr_value <= correlation_range[2]) {
      # Check if A_3_perturbed is invertible
      if (abs(det(A_3_perturbed)) > tolerance) {
        # Solve the system
        solution <- solve(A_3_perturbed, -r_3_perturbed)
        
        # Check if the solution is feasible (all values within the target range)
        if (all(solution > target_range[1] & solution < target_range[2])) {
          feasible_solution_found <- TRUE
          cat("Feasible solution found after", attempts, "attempts.\n")
          cat("Correlation between original and perturbed matrix:", round(corr_value, 3), "\n")
          print("Perturbed A_3:")
          print(round(A_3_perturbed, 3))
          print("Perturbed r_3 (non-negative):")
          print(round(r_3_perturbed, 3))
          print("Solution:")
          print(solution)
          return(list(A_3 = A_3_perturbed, r_3 = r_3_perturbed, solution = solution, correlation = corr_value))
        }
      }
    }
  }
  
  if (!feasible_solution_found) {
    cat("No feasible solution found after", max_attempts, "attempts.\n")
  }
  
  return(NULL)
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

############### 5 species  perturbed sytems ############################
##################DATA#################
set.seed(3) # for reproducibility
r_3 =  round(c(0.01, 0.001, 0.944589198242609, 0.00099999999999989, 1.17206428893925),3)
A_3 =matrix(c(
  0, 0.014590970941945, -0.00904270557385977, -0.00539161707268522, -0.00234002641389577,
  -0.0113493483829496, -0.0335522256500586, 0.0764793595400384, -0.00843814532197062, -0.0218364969836796,
  -0.00374025124828574, -0.0654365431717272, 0, -0.000927135342128989, -0.0132018969945225,
  -0.0157245111849741, -0.0159667586886125, -0.0188617119402023, -0.0386674697325392, 0.102726096554821,
  -0.000780179936723004, -0.00226628907390601, 0.00138916030935066, -0.105880681021443, 0
), nrow=5, byrow=TRUE)

A_3 <- round(A_3, 3)
x0_3=rep(15,5)

##solve the sytem
res_3 <- integrate_GLV(r_3, A_3, x0_3,20,0.1)

###now perturb the Aij 
result <- perturb_and_solve2(A_3, r_3,399,correlation_range = c(0, 0.6)) 

A_3_new=result[[1]]
r_3_new=result[[2]]

final_abundances <- as.numeric(tail(res_3[[1]], 1)[-1] [1, ])  # Extract final abundances, excluding the time column
res_3_new <- integrate_GLV(r_3_new, A_3_new, final_abundances, maxtime = 30, steptime = 0.1, n_species = 5)
#res_3=res_3_new

out_1 <- res_3[[1]]  # Initial simulation results
out_2 <- res_3_new[[1]][-1,]  # New simulation results

# Adjust the time in the new results to continue from where the first simulation ended
out_2$time <- out_2$time + max(out_1$time)

# Combine the two time series
combined_out <- rbind(out_1, out_2)

# Print combined output
print(combined_out)

###start the time from 1
combined_out$time=combined_out$time+1
p1=plot_ODE_output(combined_out,3)[[2]]

d_3=get_time_derivative_glv(combined_out)
j3=calculate_time_dependent_jacobian(d_3[[1]],d_3[[2]],10)
e3=get_eigen_values(j3,"three_sp")


####run pca on eigen and find phases
r_kpca=plot_kpca_for_sample(e3[[2]],2,"5_species")[[2]]
changepoints=find_significant_changepoints_penalties(r_kpca[[1]][-4],pen_value = 30)$OverallResults
changepoints=changepoints[changepoints$Score>=0.8,]$Changepoint



