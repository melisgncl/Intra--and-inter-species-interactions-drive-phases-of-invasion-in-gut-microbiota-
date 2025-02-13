#!/usr/bin/env Rscript

#######################    READ ME!   #####################################

#####	This Rscript filters a dataframe of barcode frequencies based on some minimum number of non-zero time points


#######################    ^^^^^^^^   #####################################

## R script 
## Louis Gauthier & Melis Gencel



# args[1] input file
# args[2] output file
# args[3] flag to remove one additional column (in case of low CFU)
# args[4] time point threshold (keep barcodes with at least x non-zero time points)

# ex usage: Rscript src/visualization/5_clustering/1_filter_data.R data/processed_samples/GFM1-4/Sample_GFM4_clustering.txt data/clustering/GFM1-4/GFM4_filtered.csv 0 14

args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

library(readr)
library(reshape2)

filename = args[1]
outfile = args[2]
flag = as.integer(args[3])
time.cut = as.integer(args[4])
mean.cut=as.numeric(args[5])

df <- read_tsv(file = filename, col_names = TRUE)

sample = reshape2::dcast(df, ID ~ Time, value.var = 'Reads')
sample$`0`=NULL
if(flag){
  sample$`1`=NULL
}
m = as.matrix(sample)
mat = as.data.frame(sweep(m,2,colSums(m,na.rm = TRUE),`/`))
mat$ID = sample$ID
mat[,"mean"] = apply(mat[,-1],1, mean,na.rm=TRUE)

z = is.na.data.frame(mat)
mat[z]=0
mat$points = apply(mat[,-c(1,ncol(mat))], 1, function(c)sum(c!=0))

sample.clustering = mat[mat$points>time.cut & mat$mean> mean.cut,]

# indicate how many trajectories pass filter
nrow(sample.clustering)

write_csv(sample.clustering,file = outfile,col_names = TRUE)
