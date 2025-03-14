## CODE BLOCK
## Louis Gauthier, August 2020

#######################    READ ME!   #####################################

##### Merge barcode list with clusters to fetch the top N barcodes per sample

#######################    ^^^^^^^^   #####################################

n_intersect=1000

###############
# read pooled cluster tables

gf1_cluster <- read_csv(paste0("data/processed_samples/GFM1/GFM1_cluster.csv"))
gf2_cluster <- read_csv(paste0("data/processed_samples/GFM2/GFM2_cluster.csv"))
gf3_cluster <- read_csv(paste0("data/processed_samples/GFM3/GFM3_cluster.csv"))
gf4_cluster <- read_csv(paste0("data/processed_samples/GFM4/GFM4_cluster.csv"))

rm1_cluster <- read_csv(paste0("data/processed_samples/M1/M1_cluster.csv"))
rm2_cluster <- read_csv(paste0("data/processed_samples/M2/M2_cluster.csv"))
rm3_cluster <- read_csv(paste0("data/processed_samples/M3/M3_cluster.csv"))
rm4_cluster <- read_csv(paste0("data/processed_samples/M4/M4_cluster.csv"))

im1_cluster <- read_csv(paste0("data/processed_samples/IM1/IM1_cluster.csv"))
im2_cluster <- read_csv(paste0("data/processed_samples/IM2/IM2_cluster.csv"))
im3_cluster <- read_csv(paste0("data/processed_samples/IM3/IM3_cluster.csv"))
im4_cluster <- read_csv(paste0("data/processed_samples/IM4/IM4_cluster.csv"))



merged_cluster <- gf1_cluster %>%
  full_join(gf2_cluster, by = "Center", suffix = c(".gf1", ".gf2")) %>%
  full_join(gf3_cluster, by = "Center", suffix = c("", ".gf3")) %>%
  full_join(gf4_cluster, by = "Center", suffix = c("", ".gf4")) %>%
  full_join(rm1_cluster, by = "Center", suffix = c("", ".rm1")) %>%
  full_join(rm2_cluster, by = "Center", suffix = c("", ".rm2")) %>%
  full_join(rm3_cluster, by = "Center", suffix = c("", ".rm3")) %>%
  full_join(rm4_cluster, by = "Center", suffix = c("", ".rm4")) %>%
  full_join(im1_cluster, by = "Center", suffix = c("", ".im1")) %>%
  full_join(im2_cluster, by = "Center", suffix = c("", ".im2")) %>%
  full_join(im3_cluster, by = "Center", suffix = c("", ".im3")) %>%
  full_join(im4_cluster, by = "Center", suffix = c("", ".im4"))

final_cluster <- merged_cluster %>%
  select(Center,
         Cluster.ID.gf1, Cluster.ID.gf2, Cluster.ID, Cluster.ID.gf4,
         Cluster.ID.rm1, Cluster.ID.rm2, Cluster.ID.rm3, Cluster.ID.rm4,
         Cluster.ID.im1, Cluster.ID.im2, Cluster.ID.im3, Cluster.ID.im4) %>%
  rename(gf1.ID = Cluster.ID.gf1, gf2.ID = Cluster.ID.gf2, gf3.ID = Cluster.ID,
         gf4.ID = Cluster.ID.gf4, rm1.ID = Cluster.ID.rm1, rm2.ID = Cluster.ID.rm2,
         rm3.ID = Cluster.ID.rm3, rm4.ID = Cluster.ID.rm4, im1.ID = Cluster.ID.im1,
         im2.ID = Cluster.ID.im2, im3.ID = Cluster.ID.im3, im4.ID = Cluster.ID.im4)

# View the final merged data frame
print(final_cluster)

write_csv(final_cluster,paste0("data/processed_samples/Gavage_index_overlap.csv"))






fetchTop <- function(reshaped_df, sample_cluster) {
	df.top = unique(reshaped_df[,1:4])
	df.top.final = df.top[order(-df.top$final),]
	df.top.max = df.top[order(-df.top$max),]

	df.top.final = df.top.final[1:n_intersect,]
	df.top.max = df.top.max[1:n_intersect,]

	# match back to the barcode cluster sequence
	df.top.final = merge(df.top.final,sample_cluster,by.x = "ID",by.y = "Cluster.ID")
	df.top.final$Cluster.Score=NULL
	df.top.final$time_point_1=NULL

	df.top.final = df.top.final[order(-df.top.final$final),]

	df.top.max = merge(df.top.max,sample_cluster,by.x = "ID",by.y = "Cluster.ID")
	df.top.max$Cluster.Score=NULL
	df.top.max$time_point_1=NULL

	df.top.max = df.top.max[order(-df.top.max$max),]

	return(list(df.top.final,df.top.max))
}

###############
# get top N IDs based on max frequency


###gf###
x = fetchTop(gf1.df,gf1_cluster)
gf1.top.final = as.data.frame(x[1])
gf1.top.max = as.data.frame(x[2])

x = fetchTop(gf2.df,gf2_cluster)
gf2.top.final = as.data.frame(x[1])
gf2.top.max = as.data.frame(x[2])

x = fetchTop(gf3.df,gf3_cluster)
gf3.top.final = as.data.frame(x[1])
gf3.top.max = as.data.frame(x[2])

x = fetchTop(gf4.df,gf4_cluster)
gf4.top.final = as.data.frame(x[1])
gf4.top.max = as.data.frame(x[2])

###rm###
x= fetchTop(rm1.df,rm1_cluster)
rm1.top.final = as.data.frame(x[1])
rm1.top.max = as.data.frame(x[2])

x = fetchTop(rm2.df,rm2_cluster)
rm2.top.final = as.data.frame(x[1])
rm2.top.max = as.data.frame(x[2])

x = fetchTop(rm3.df,rm3_cluster)
rm3.top.final = as.data.frame(x[1])
rm3.top.max = as.data.frame(x[2])

x = fetchTop(rm4.df,rm4_cluster)
rm4.top.final = as.data.frame(x[1])
rm4.top.max = as.data.frame(x[2])

###im###
x= fetchTop(im1.df,im1_cluster)
im1.top.final = as.data.frame(x[1])
im1.top.max = as.data.frame(x[2])

x = fetchTop(im2.df,im2_cluster)
im2.top.final = as.data.frame(x[1])
im2.top.max = as.data.frame(x[2])

x = fetchTop(im3.df,im3_cluster)
im3.top.final = as.data.frame(x[1])
im3.top.max = as.data.frame(x[2])

x = fetchTop(im4.df,im4_cluster)
im4.top.final = as.data.frame(x[1])
im4.top.max = as.data.frame(x[2])

