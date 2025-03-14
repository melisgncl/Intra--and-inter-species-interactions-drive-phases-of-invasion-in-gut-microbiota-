###  Melis Gencel
### DCM analysis
#######################    READ ME!   #####################################

##### In here we analyze  nc and im which are not colonized

#######################    ^^^^^^^^   #####################################

source("/Users/melis/Dropbox/mouse_gut/github/src/visualization/7_DCM/1_DCM.R")
#################################################im########################################################
#################################################im########################################################

##### im with E.coli
plot_eigen_and_pca=function(taxa,span,n,sample,exclude_Enterobacteriaceae=FALSE){
    ###if you want to remove e.coli
    if (exclude_Enterobacteriaceae) {
      taxa <- select(as.data.frame(taxa), -Enterobacteriaceae)
    }     
    d=get_time_derivative_taxa(taxa,span,n)
    j=calculate_time_dependent_jacobian(d[[1]],d[[2]],10,sample)
    e=get_eigen_values(j,sample)
    ##plot_eigen
    plot_eigen_value_evolution(e[[1]],sample)
    return(e)       
      
}

####im with Enterobacteriaceae

im1_16S_taxa <- read_delim("data/16S/taxa_long_format/im1_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
im2_16S_taxa <- read_delim("data/16S/taxa_long_format//im2_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
im3_16S_taxa <- read_delim("data/16S/taxa_long_format//im3_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
im4_16S_taxa <- read_delim("data/16S/taxa_long_format//im4_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)

e_im1=plot_eigen_and_pca(im1_16S_taxa,0.25,7,"im1")
e_im2=plot_eigen_and_pca(im2_16S_taxa,0.25,7,"im2")
e_im3=plot_eigen_and_pca(im3_16S_taxa,0.25,7,"im3")
e_im4=plot_eigen_and_pca(im4_16S_taxa,0.25,7,"im4")

im=list(e_im1[[2]],e_im2[[2]],e_im3[[2]],e_im4[[2]])

im_eigen_kpca_with_Enterobacteriaceae=plot_kpca_for_sample(im,c(2,2,2,2),"im")[[2]]


##without Enterobacteriaceae

e_im1_w=plot_eigen_and_pca(im1_16S_taxa,0.25,7,"im1",exclude_Enterobacteriaceae=TRUE)
e_im2_w=plot_eigen_and_pca(im2_16S_taxa,0.25,7,"im2",exclude_Enterobacteriaceae=TRUE)
e_im3_w=plot_eigen_and_pca(im3_16S_taxa,0.25,7,"im3",exclude_Enterobacteriaceae=TRUE)
e_im4_w=plot_eigen_and_pca(im4_16S_taxa,0.25,7,"im4",exclude_Enterobacteriaceae=TRUE)

im=list(e_im1_w[[2]],e_im2_w[[2]],e_im3_w[[2]],e_im4_w[[2]])

im_eigen_kpca_without_Enterobacteriaceae=plot_kpca_for_sample(im,c(2,2,2,2),"im_without_E")[[2]]


###############################################NC###################

nc1_16S_taxa <- read_delim("data/16S/taxa_long_format/nc1_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)

nc2_16S_taxa <- read_delim("data/16S/taxa_long_format/nc2_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)

nc3_16S_taxa <- read_delim("data/16S/taxa_long_format/nc3_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)

nc4_16S_taxa <- read_delim("data/16S/taxa_long_format/nc4_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)

###interpolate nc3 nc 7 and nc8 coloumns
# Function to interpolate between two rows and return the new row
interpolate_rows <- function(row1, row2) {
  return((row1 + row2) / 2)
}

# Interpolate between rows 3 and 4
new_row_7 <- interpolate_rows(nc3_16S_taxa[6, ], nc3_16S_taxa [7, ])
# Interpolate between (original) rows 4 and 5
new_row_8 <- interpolate_rows(new_row_7, nc3_16S_taxa [7, ])
nc3_16S_taxa <- rbind(nc3_16S_taxa[1:6, ], new_row_7, new_row_8, nc3_16S_taxa[7:8,])


nc1_16S_taxa=filter_columns_by_mean(nc1_16S_taxa)
nc2_16S_taxa=filter_columns_by_mean(nc2_16S_taxa)
nc3_16S_taxa=filter_columns_by_mean(nc3_16S_taxa)
nc4_16S_taxa=filter_columns_by_mean(nc4_16S_taxa)

e_nc1=plot_eigen_and_pca(nc1_16S_taxa,0.2,4,"nc1")
e_nc2=plot_eigen_and_pca(nc2_16S_taxa,0.2,4,"nc2")
e_nc3=plot_eigen_and_pca(nc3_16S_taxa,0.2,4,"nc3")
e_nc4=plot_eigen_and_pca(nc4_16S_taxa,0.2,4,"nc4")

nc=list(e_nc1[[2]],e_nc2[[2]],e_nc3[[2]],e_nc4[[2]])

nc_eigen_kpca=plot_kpca_for_sample(nc,c(2,2,2,2),"nc")[[2]]

