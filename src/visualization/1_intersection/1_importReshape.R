## CODE BLOCK
## Louis Gauthier & Melis Gencel

#######################    READ ME!   #####################################

# Import each sample's clustering data and reshape in convenient format

#######################    ^^^^^^^^   #####################################

### 1. import samples

###innate (im)


Sample_im1 <- read_delim(paste0("data/processed_samples/IM1/Sample_IM1_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_im2 <- read_delim(paste0("data/processed_samples/IM2/Sample_IM2_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_im3 <- read_delim(paste0("data/processed_samples/IM3/Sample_IM3_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_im4 <- read_delim(paste0("data/processed_samples/IM4/Sample_IM4_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)


####reduced microbiota (rm)

Sample_rm1 <- read_delim(paste0("data/processed_samples/M1/Sample_M1_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_rm2 <- read_delim(paste0("data/processed_samples/M2/Sample_M2_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_rm3 <- read_delim(paste0("data/processed_samples/M3/Sample_M3_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_rm4 <- read_delim(paste0("data/processed_samples/M4/Sample_M4_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)



######gf###

Sample_gf1 <- read_delim(paste0("data/processed_samples/GFM1/Sample_GFM1_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_gf2 <- read_delim(paste0("data/processed_samples/GFM2/Sample_GFM2_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_gf3 <- read_delim(paste0("data/processed_samples/GFM3/Sample_GFM3_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Sample_gf4 <- read_delim(paste0("data/processed_samples/GFM4/Sample_GFM4_clustering.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE)





### 2. define function to format data for plotting

reshapeDF <- function(raw_df) {
  library(reshape2)
  library(dplyr)
  
  sample = raw_df
  
  test = reshape2::dcast(sample, ID ~ Time, value.var = 'Reads')
  m = as.matrix(test[,-1])
  mat = as.data.frame(sweep(m,2,colSums(m,na.rm = TRUE),`/`))
  mat$ID = test$ID
  mat[,"mean"] = apply(mat[,-ncol(mat)],1, mean,na.rm=TRUE)
  mat[,"max"] = apply(mat[,-c(ncol(mat),ncol(mat)-1)],1, max,na.rm=TRUE)
  mat$start = mat[,1]
  mat$final = mat[,ncol(mat)-4]
  
  
  df = reshape2::melt(mat,id.vars = c('ID','max','start','final','mean'))
  df$variable = as.numeric(levels(df$variable))[df$variable]
  
  df$ID = as.factor(df$ID)
  
  if (any(is.na(df$value))) {
    df[is.na(df$value), ]$value <- 0
  }
  
  tf = df[order(df$max),]
  
  return(tf)
}

### 3. reshape data


###gf###
gf1.df = reshapeDF(Sample_gf1)
gf2.df = reshapeDF(Sample_gf2)
gf3.df = reshapeDF(Sample_gf3)
gf4.df = reshapeDF(Sample_gf4)

###rm###
rm1.df = reshapeDF(Sample_rm1)
rm2.df = reshapeDF(Sample_rm2)
rm3.df = reshapeDF(Sample_rm3)
rm4.df = reshapeDF(Sample_rm4)

###im###
im1.df = reshapeDF(Sample_im1)
im2.df = reshapeDF(Sample_im2)
im3.df = reshapeDF(Sample_im3)
im4.df = reshapeDF(Sample_im4)

