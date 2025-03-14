# Load necessary package for plotting
library(gplots) # for heatmap.2
set.seed(456L)
overlap_palette <- colorRampPalette(c("white","blue"))(n = 100)


rename_columns <- function(series, prefix) {
  names(series) <- paste(prefix, names(series), sep = ".")
  return(series)
}
taxa.series = transformTaxa(rm1_16S_taxa,rm1_clustered_loess_log10)
cluster.series = tslist(t(rm1_clustered_loess_log10))
series = append(cluster.series,taxa.series)

rm1_series=rename_columns(series,"rm1")

taxa.series = transformTaxa(rm2_16S_taxa,rm2_clustered_loess_log10)
cluster.series = tslist(t(rm2_clustered_loess_log10))
series = append(cluster.series,taxa.series)

rm2_series=rename_columns(series,"rm2")

taxa.series = transformTaxa(rm3_16S_taxa,rm3_clustered_loess_log10)
cluster.series = tslist(t(rm3_clustered_loess_log10))
series = append(cluster.series,taxa.series)

rm3_series=rename_columns(series,"rm3")

taxa.series = transformTaxa(rm4_16S_taxa,rm4_clustered_loess_log10)
cluster.series = tslist(t(rm4_clustered_loess_log10))
series = append(cluster.series,taxa.series)

rm4_series=rename_columns(series,"rm4")

combined_series <- append(rm1_series, rm2_series)
combined_series <- append(combined_series, rm3_series)
combined_series <- append(combined_series,rm4_series)


sbD <- proxy::dist(combined_series, method = "SBD", znorm = TRUE)

# Output the resulting SBD distance matrix
print(sbD)
hc <- hclust(sbD, method = "average")
dendrogram <- as.dendrogram(hc)

par(mar = c(2,2,2,2))
setEPS()
postscript("reports/barcodeCounter_Delcorrect/nature_com/all_sbd.eps",width = 15,height = 10)
heatmap.2(as.matrix(sbD),Rowv = dendrogram,Colv = dendrogram,
          col = colorRampPalette(c("navy", "white", "firebrick"))(50),
          density.info = "none",trace = "none",
          key.xlab="SBD Value",
          cexRow = 0.2,cexCol = 0.2)
dev.off()

overlap.matrix=read_csv("reports/barcodeCounter_Delcorrect/overlap_of_clusters/rm1-4_overlap_.matrix",col_names = TRUE)
rownames(overlap.matrix) = colnames(overlap.matrix)
overlap_columns <- colnames(overlap.matrix)
sbd_columns <- colnames(as.matrix(sbD))
# Find the columns that are in the SBD matrix but not in the overlap matrix
missing_columns <- setdiff(sbd_columns, overlap_columns)
new_columns <- matrix(0, nrow = nrow(overlap.matrix), ncol = length(missing_columns))
colnames(new_columns) <- missing_columns
new_columns_df <- as_tibble(new_columns)
# Bind the new columns with zeros to the overlap.matrix
updated_overlap_matrix <- bind_cols(overlap.matrix, new_columns_df)
for (missing_var in missing_columns) {
  zero_row <- setNames(rep(0, ncol(updated_overlap_matrix)), colnames(updated_overlap_matrix))
  updated_overlap_matrix <- add_row(updated_overlap_matrix, !!!zero_row)
  rownames(updated_overlap_matrix)[nrow(updated_overlap_matrix)] <- missing_var
}

rownames(updated_overlap_matrix) = colnames(updated_overlap_matrix)
updated_overlap_matrix <- updated_overlap_matrix[rownames(as.matrix(sbD)), colnames(as.matrix(sbD))]
rownames(updated_overlap_matrix) = colnames(updated_overlap_matrix)


par(mar = c(2,2,2,2))
setEPS()
postscript("reports/barcodeCounter_Delcorrect/nature_com/all_sbd_overlap.eps",width = 15,height = 10)
heatmap.2(as.matrix(updated_overlap_matrix),
          Rowv = dendrogram,Colv = dendrogram,
          col=overlap_palette,
          density.info = "none",trace = "none",
          key.xlab="Jaccard index",cexRow = 0.2,cexCol = 0.2)

dev.off()


cluster_assignments <- cutree(hc, k=4)

# Assign cluster colors for visualization purposes
dendrogram <- color_branches(dendrogram,  k=4)
plot(dendrogram)

# Print out the cluster assignments
print(cluster_assignments)

# Create an empty list to store data for each cluster
clustered_matrices <- list()

# Convert the sbD distance object into a regular matrix

# Convert sbD distance object into a regular matrix
sbD_matrix <- as.matrix(sbD)

# Extract row and column order from dendrogram and convert indices to actual names
dend_order <- order.dendrogram(dendrogram)
ordered_names <- rownames(sbD_matrix)[dend_order]

# Debugging: Print ordered names
cat("Ordered Names from Dendrogram:\n", ordered_names, "\n")

# Create an empty list to store sub-matrices for each combination of row and column clusters
cluster_combination_matrices <- list()

# Step 1: Extract sub-matrices for each combination of row and column clusters
for (row_cluster in 1:4) {
  for (col_cluster in 1:4) {
    # Get the names of the elements that belong to row_cluster and col_cluster
    row_elements <- names(cluster_assignments[cluster_assignments == row_cluster])
    col_elements <- names(cluster_assignments[cluster_assignments == col_cluster])
    
    # Only proceed if both row_elements and col_elements are not empty
    if (length(row_elements) > 0 && length(col_elements) > 0) {
      # Extract the sub-matrix for the specific row and column cluster combination
      sub_matrix <- sbD_matrix[row_elements, col_elements, drop = FALSE]
      # Save the extracted sub-matrix
      cluster_combination_matrices[[paste("row", row_cluster, "col", col_cluster, sep = "_")]] <- sub_matrix
    } else {
      # Save an empty sub-matrix if the cluster is empty
      cluster_combination_matrices[[paste("row", row_cluster, "col", col_cluster, sep = "_")]] <- matrix(nrow = 0, ncol = 0)
    }
  }
}



# Step 2: Reorder each sub-matrix based on the original dendrogram order
# Reorder each of the sub-matrices based on the ordered names
for (name in names(cluster_combination_matrices)) {
  cluster_matrix <- cluster_combination_matrices[[name]]
  
  # Only reorder if the matrix is non-empty and has more than one row/column
  if (!is.null(cluster_matrix) && nrow(cluster_matrix) > 0 && ncol(cluster_matrix) > 0) {
    # Get the elements that are present in the current sub-matrix and are also in dendrogram order
    ordered_row_elements <- intersect(ordered_names, rownames(cluster_matrix))
    ordered_col_elements <- intersect(ordered_names, colnames(cluster_matrix))
    
    # Reorder the sub-matrix
    cluster_combination_matrices[[name]] <- cluster_matrix[ordered_row_elements, ordered_col_elements, drop = FALSE]
  }
}
library(gplots)

# Step 1: Define a consistent color gradient from blue to white to red
color_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(50)

# Step 2: Define fixed breaks from 0 to 1
# This ensures that the entire range from 0 to 1 is covered consistently
breaks <- seq(0, 1, length.out = 51)  # 50 intervals in the range 0 to 1

# Plotting the heatmaps with a consistent color scale across all matrices
for (name in names(cluster_combination_matrices)) {
  cluster_matrix <- cluster_combination_matrices[[name]]
  
  # Check if the sub-matrix has more than one row and column to be plotted
  if (!is.null(cluster_matrix) && nrow(cluster_matrix) > 1 && ncol(cluster_matrix) > 1) {
    # Set up PNG device for the current cluster combination
    png_filename <- paste(name, "_ordered_heatmap2.png", sep = "")
    png(png_filename, width = 1000, height = 800, units = "px", res = 150)
    
    # Plot the heatmap using heatmap.2() with consistent color scaling
    heatmap.2(
      cluster_matrix,                           # Matrix to plot
      dendrogram = "none",                      # Disable dendrograms for both rows and columns
      Rowv = FALSE,                             # Disable row reordering
      Colv = FALSE,                             # Disable column reordering
      col = color_palette,                      # Consistent color gradient for all heatmaps
      breaks = breaks,                          # Fixed color breaks from 0 to 1
      trace = "none",                           # Disable trace lines inside the heatmap
      main = paste("Heatmap of Cluster Combination -", name), # Main title for the heatmap
      margins = c(10, 10),                      # Adjust margins for better visualization
      xlab = "Cluster Elements (Columns)",      # Label for x-axis
      ylab = "Cluster Elements (Rows)",         # Label for y-axis
      key = TRUE,                               # Add a color key for reference
      density.info = "none"                     # Disable density plot for the color key
    )
    
    # Close the PNG device
    dev.off()
  } else {
    # If the sub-matrix has less than two elements, skip the heatmap
    cat(name, "has less than two elements. Skipping heatmap.\n")
  }
}


library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Step 1: Melt all matrices in the list and retain cluster information
combined_data <- data.frame()

# List of clusters to include and their new names
# selected_clusters <- c("row_1_col_1", "row_1_col_3", "row_1_col_2", "row_3_col_3", "row_3_col_2", "row_2_col_2")
# new_names <- c("1", "2", "3", "4", "5", "6")

selected_clusters <- c("row_1_col_1", "row_1_col_3", "row_3_col_3")
new_names <- c("1", "3", "2")

library(ggbeeswarm)

library(dplyr)
library(reshape2)

combined_data <- data.frame()  # Initialize combined data frame

# Loop through the cluster combination matrices
for (name in names(cluster_combination_matrices)) {
  if (name %in% selected_clusters) {
    cluster_matrix <- cluster_combination_matrices[[name]]
    
    # Proceed only if the matrix has more than one row and column
    if (!is.null(cluster_matrix) && nrow(cluster_matrix) > 1 && ncol(cluster_matrix) > 1) {
      
      # Melt the matrix into long format
      melted_matrix <- melt(cluster_matrix)
      
      # Rename columns for convenience
      colnames(melted_matrix) <- c("Row", "Column", "Value")
      
      # Step 2: Remove duplicate points within the cluster due to symmetry
      # Create a consistent identifier by always putting the smaller value first
      melted_matrix <- melted_matrix %>%
        mutate(UniquePair = paste0(pmin(as.character(Row), as.character(Column)), "_", 
                                   pmax(as.character(Row), as.character(Column)))) %>%
        distinct(UniquePair, .keep_all = TRUE)
      
      # Remove rows where UniquePair has the form "X_X" (e.g., "rm2.Acholeplasmataceae_rm2.Acholeplasmataceae")
      melted_matrix <- melted_matrix %>%
        filter(!grepl("^(.*)_\\1$", UniquePair))
      
      # Add a column for description
      melted_matrix <- melted_matrix %>%
        mutate(
          Description = case_when(
            grepl("^rm[0-9]+\\.C[0-9]+_rm[0-9]+\\.C[0-9]+$", UniquePair) ~ "clone-clone",
            grepl("^rm[0-9]+\\.C[0-9]+_rm[0-9]+\\.[A-Za-z]+$", UniquePair) |
              grepl("^rm[0-9]+\\.[A-Za-z]+_rm[0-9]+\\.C[0-9]+$", UniquePair) ~ "clone-species",
            grepl("^rm[0-9]+\\.[A-Za-z]+_rm[0-9]+\\.[A-Za-z]+$", UniquePair) ~ "species-species",
            TRUE ~ "unknown"
          )
        )
      
      # Rename the cluster with the new given sequence
      new_name <- new_names[which(selected_clusters == name)]
      melted_matrix$Cluster <- new_name
      
      # Combine with the rest of the data
      combined_data <- rbind(combined_data, melted_matrix)
    }
  }
}



# Display combined data for debugging purposes
print(combined_data)

# Define comparisons
comparisons <- list(
  c("2", "3"),
  c("1", "3"))
  

# Step 3: Create a side-by-side violin plot for selected clusters
violin_plot=ggplot(combined_data, aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 0.2) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 

# Step 5: Save the plot to a file
ggsave("selected_violin_plot_with_pvalues.png", plot = violin_plot, width = 14, height = 8, units = "in")


####separete clone-clone and clone-species

violin_plot=ggplot(combined_data[combined_data$Description=='clone-clone',], aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 0.2) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 

ggsave("selected_violin_plot_with_pvalues_clone.png", plot = violin_plot, width = 14, height = 8, units = "in")


violin_plot=ggplot(combined_data[combined_data$Description=='clone-species',], aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 0.2) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 

ggsave("selected_violin_plot_with_pvalues_clone_species.png", plot = violin_plot, width = 14, height = 8, units = "in")





filtered_combined_data <- combined_data %>%
  filter(Description == 'clone-clone') %>%
  separate(UniquePair, into = c("First", "Second"), sep = "_") %>%
  filter(gsub("\\..*", "", First) != gsub("\\..*", "", Second)) %>%
  unite("UniquePair", First, Second, sep = "_")


violin_plot=ggplot(filtered_combined_data, aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 0.2) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 

ggsave("selected_violin_plot_with_pvalues_clone_clone_filetred.png", plot = violin_plot, width = 14, height = 8, units = "in")










updated_overlap_matrix=as.matrix(updated_overlap_matrix)
o_cluster_combination_matrices <- list()

# Step 1: Extract sub-matrices for each combination of row and column clusters
for (row_cluster in 1:4) {
  for (col_cluster in 1:4) {
    # Get the names of the elements that belong to row_cluster and col_cluster
    row_elements <- names(cluster_assignments[cluster_assignments == row_cluster])
    col_elements <- names(cluster_assignments[cluster_assignments == col_cluster])
    
    # Only proceed if both row_elements and col_elements are not empty
    if (length(row_elements) > 0 && length(col_elements) > 0) {
      # Extract the sub-matrix for the specific row and column cluster combination
      sub_matrix <- updated_overlap_matrix[row_elements, col_elements, drop = FALSE]
      # Save the extracted sub-matrix
      o_cluster_combination_matrices[[paste("row", row_cluster, "col", col_cluster, sep = "_")]] <- sub_matrix
    } else {
      # Save an empty sub-matrix if the cluster is empty
      o_cluster_combination_matrices[[paste("row", row_cluster, "col", col_cluster, sep = "_")]] <- matrix(nrow = 0, ncol = 0)
    }
  }
}



# Step 2: Reorder each sub-matrix based on the original dendrogram order
# Reorder each of the sub-matrices based on the ordered names
for (name in names(o_cluster_combination_matrices)) {
  cluster_matrix <- o_cluster_combination_matrices[[name]]
  
  # Only reorder if the matrix is non-empty and has more than one row/column
  if (!is.null(cluster_matrix) && nrow(cluster_matrix) > 0 && ncol(cluster_matrix) > 0) {
    # Get the elements that are present in the current sub-matrix and are also in dendrogram order
    ordered_row_elements <- intersect(ordered_names, rownames(cluster_matrix))
    ordered_col_elements <- intersect(ordered_names, colnames(cluster_matrix))
    
    # Reorder the sub-matrix
    o_cluster_combination_matrices[[name]] <- cluster_matrix[ordered_row_elements, ordered_col_elements, drop = FALSE]
  }

}



# Step 1: Melt all matrices in the list and retain cluster information
combined_data <- data.frame()

# List of clusters to include and their new names
#selected_clusters <- c("row_1_col_1", "row_1_col_3", "row_1_col_2", "row_3_col_3", "row_3_col_2", "row_2_col_2")
#new_names <- c("1", "2", "3", "4", "5", "6")

selected_clusters <- c("row_1_col_1", "row_1_col_3", "row_3_col_3")
new_names <- c("1", "3", "2")


# Loop through the cluster combination matrices
for (name in names(o_cluster_combination_matrices)) {
  if (name %in% selected_clusters) {
    cluster_matrix <-o_cluster_combination_matrices[[name]]
    
    # Proceed only if the matrix has more than one row and column
    if (!is.null(cluster_matrix) && nrow(cluster_matrix) > 1 && ncol(cluster_matrix) > 1) {
      
      # Melt the matrix into long format
      melted_matrix <- melt(cluster_matrix)
      
      # Rename columns for convenience
      colnames(melted_matrix) <- c("Row", "Column", "Value")
      
      # Step 2: Remove duplicate points within the cluster due to symmetry
      # Create a consistent identifier by always putting the smaller value first
      melted_matrix <- melted_matrix %>%
        mutate(UniquePair = paste0(pmin(as.character(Row), as.character(Column)), "_", 
                                   pmax(as.character(Row), as.character(Column)))) %>%
        distinct(UniquePair, .keep_all = TRUE)
      
      # Remove rows where UniquePair has the form "X_X" (e.g., "rm2.Acholeplasmataceae_rm2.Acholeplasmataceae")
      melted_matrix <- melted_matrix %>%
        filter(!grepl("^(.*)_\\1$", UniquePair))
      
      # Add a column for description
      melted_matrix <- melted_matrix %>%
        mutate(
          Description = case_when(
            grepl("^rm[0-9]+\\.C[0-9]+_rm[0-9]+\\.C[0-9]+$", UniquePair) ~ "clone-clone",
            grepl("^rm[0-9]+\\.C[0-9]+_rm[0-9]+\\.[A-Za-z]+$", UniquePair) |
              grepl("^rm[0-9]+\\.[A-Za-z]+_rm[0-9]+\\.C[0-9]+$", UniquePair) ~ "clone-species",
            grepl("^rm[0-9]+\\.[A-Za-z]+_rm[0-9]+\\.[A-Za-z]+$", UniquePair) ~ "species-species",
            TRUE ~ "unknown"
          )
        )
      
      # Rename the cluster with the new given sequence
      new_name <- new_names[which(selected_clusters == name)]
      melted_matrix$Cluster <- new_name
      
      # Combine with the rest of the data
      combined_data <- rbind(combined_data, melted_matrix)
    }
  }
}

library(dplyr)

# Step 1: Define a vector of allowed patterns
allowed_patterns <- paste0("C", 1:20) # This generates "C1", "C2", "C3", ..., "C20". Adjust the range as needed.







# Define comparisons
comparisons <- list(
  c("2", "3"),
  c("1", "3"))

violin_plot=ggplot(combined_data[combined_data$Description=='clone-clone',], aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 0.2) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 

# Step 5: Save the plot to a file
ggsave("selected_violin_plot_with_ovelap.png", plot = violin_plot, width = 14, height = 8, units = "in")

###remove same mouse clones pairs from here
filtered_combined_data <- combined_data %>%
  filter(Description == 'clone-clone') %>%
  separate(UniquePair, into = c("First", "Second"), sep = "_") %>%
  filter(gsub("\\..*", "", First) != gsub("\\..*", "", Second)) %>%
  unite("UniquePair", First, Second, sep = "_")


violin_plot=ggplot(filtered_combined_data, aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 0.2) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 

# Step 5: Save the plot to a file
ggsave("selected_violin_plot_with_ovelap_filtered.png", plot = violin_plot, width = 14, height = 8, units = "in")



####check their statistic according to zscore
df=read_csv("reports/barcodeCounter_Delcorrect/overlap_of_clusters/rm1-4_statistic_data.csv")
subset_df=df[df$id %in% filtered_combined_data$UniquePair, ]
# Assume subset_df contains a column 'p_value' with the pre-calculated p-values
# Add a color column based on the p-values
subset_df$color <- ifelse(subset_df$pmatrix >= -log10(0.05), "blue", "grey")
# Select only the necessary columns from subset_df
subset_df_selected <- subset_df %>%
  select(id, pmatrix, color)

# Merge the filtered_combined_data with the selected columns
merged_df <- filtered_combined_data %>%
  left_join(subset_df_selected, by = c("UniquePair" = "id"))



# Create the violin plot with color-coded points based on p-value
violin_plot <- ggplot(merged_df, aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE) +
  geom_quasirandom(aes(color = color), size = 1) + # Color points based on p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") +
  scale_color_manual(values = c("grey" = "grey", "blue" = "blue")) +
  labs(x = "Cluster", y = "Value") +
  theme(legend.position = "none")

# Display the plot
print(violin_plot)

ggsave("selected_violin_plot_with_ovelap_filtered_pvalues.png", plot = violin_plot, width = 14, height = 8, units = "in")


library(dplyr)

# Save data for each cluster as a separate CSV file
clusters <- c("1", "4")  # Specify the clusters you want to save

for (cluster in clusters) {
  # Filter the data for the specific cluster
  cluster_data <- merged_df %>% filter(Cluster == cluster)
  
  # Save the filtered data as a CSV file
  write.csv(cluster_data, file = paste0("cluster_", cluster, "_data.csv"), row.names = FALSE)
}





######alll of them sbd #######

violin_plot1=ggplot(combined_data, aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 0.5) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 

violin_plot2=ggplot(filtered_combined_data, aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_quasirandom(size = 1) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") 



violin_plot3 <- ggplot(merged_df, aes(x = Cluster, y = Value)) + 
  theme_Publication() + 
  geom_violin(trim = FALSE) +
  geom_quasirandom(aes(color = color), size = 1) + # Color points based on p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(comparisons = comparisons, 
              test = "wilcox.test",
              map_signif_level = TRUE, 
              step_increase = 0.1, 
              textsize = 8, 
              tip_length = 0, 
              vjust = 0.5,
              color = "black") +
  scale_color_manual(values = c("grey" = "grey", "blue" = "blue")) +
  labs(x = "Cluster", y = "Value") +
  theme(legend.position = "none")

# Display the plot
print(violin_plot)


# Load necessary libraries
library(ggplot2)
library(ggsignif)
library(ggbeeswarm)
library(cowplot)

# Assuming violin_plot1, violin_plot2, and violin_plot3 are already defined
# Combine the plots side by side using cowplot
combined_plot <- plot_grid(
  violin_plot1, violin_plot2, violin_plot3,
  labels = c("A", "B", "C"),  # Optional: Adding labels to each plot
  ncol = 3                     # Arrange plots in a single row with 3 columns
)

# Save the combined plot as an EPS file
ggsave("combined_violin_plots.eps", combined_plot, width = 18, height = 6)

