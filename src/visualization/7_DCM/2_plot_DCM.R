###  Melis Gencel
### DCM analysis
#######################    READ ME!   #####################################

##### In here we analyze gavaged cohorts rm and gf later we anlayze the nc and im which are not colonized

#######################    ^^^^^^^^   #####################################

source("/Users/melis/Dropbox/mouse_gut/github/src//visualization/7_DCM/1_DCM.R")
#################################################im########################################################
#################################################im########################################################

###rm cohort###
outdir="reports/figures/"


rm1_16S_taxa <- read_delim("data/16S/taxa_long_format/rm1_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm2_16S_taxa <- read_delim("data/16S/taxa_long_format//rm2_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm3_16S_taxa <- read_delim("data/16S/taxa_long_format//rm3_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm4_16S_taxa <- read_delim("data/16S/taxa_long_format//rm4_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)



##calculate z dot##
d_rm1=get_time_derivative_invaded(rm1_16S_taxa,"rm1",0.2,7)
d_rm2=get_time_derivative_invaded(rm2_16S_taxa,"rm2",0.2,7)
d_rm3=get_time_derivative_invaded(rm3_16S_taxa,"rm3",0.2,7)
d_rm4=get_time_derivative_invaded(rm4_16S_taxa,"rm4",0.2,7)


###calculate time dependent jacobian
j_rm1=calculate_time_dependent_jacobian(d_rm1[[1]],d_rm1[[2]],10,"rm1")
j_rm2=calculate_time_dependent_jacobian(d_rm2[[1]],d_rm2[[2]],10,"rm2")
j_rm3=calculate_time_dependent_jacobian(d_rm3[[1]],d_rm3[[2]],10,"rm3")
j_rm4=calculate_time_dependent_jacobian(d_rm4[[1]],d_rm4[[2]],10,"rm4")

##get_eigen_values

e_rm1=get_eigen_values(j_rm1,"rm1")
e_rm2=get_eigen_values(j_rm2,"rm2",time_increase = TRUE)
e_rm3=get_eigen_values(j_rm3,"rm3",time_increase = TRUE)
e_rm4=get_eigen_values(j_rm4,"rm4")

##plot eigenvalues##
plot_eigen_value_evolution(e_rm1[[1]],"rm1")
plot_eigen_value_evolution(e_rm2[[1]],"rm2")
plot_eigen_value_evolution(e_rm3[[1]],"rm3")
plot_eigen_value_evolution(e_rm4[[1]],"rm4")

####determine alternative phases using PCA
rm=list(e_rm1[[2]],e_rm2[[2]],e_rm3[[2]],e_rm4[[2]])
rm_eigen_kpca=plot_kpca_for_sample(rm,c(0.5,0.5,0.5,1),"rm")[[2]]

find_significant_changepoints_penalties(rm_eigen_kpca[[1]],pen_value = 30,
                                        flip_kpca2 = TRUE,flip_kpca1 = FALSE)

find_significant_changepoints_penalties(rm_eigen_kpca[[2]],pen_value = 30,
                                        flip_kpca2 = TRUE,flip_kpca1 = TRUE)$OverallResults

find_significant_changepoints_penalties(rm_eigen_kpca[[3]],pen_value = 30,
                                        flip_kpca2 = FALSE,flip_kpca1 = TRUE)$OverallResults

find_significant_changepoints_penalties(rm_eigen_kpca[[4]],pen_value = 30,
                                        flip_kpca2 = TRUE,flip_kpca1 = TRUE)$OverallResults

###########################GF####################################################


d_gf1=get_time_derivative_lineage("gf1")
d_gf2=get_time_derivative_lineage("gf2")
d_gf3=get_time_derivative_lineage("gf3")
d_gf4=get_time_derivative_lineage("gf4")

###calculate time dependent jacobian
j_gf1=calculate_time_dependent_jacobian(d_gf1[[1]],d_gf1[[2]],10,"gf1")
j_gf2=calculate_time_dependent_jacobian(d_gf2[[1]],d_gf2[[2]],10,"gf2")
j_gf3=calculate_time_dependent_jacobian(d_gf3[[1]],d_gf3[[2]],10,"gf3")
j_gf4=calculate_time_dependent_jacobian(d_gf4[[1]],d_gf4[[2]],10,"gf4")

##get_eigen_values

e_gf1=get_eigen_values(j_gf1,"gf1")
e_gf2=get_eigen_values(j_gf2,"gf2")
e_gf3=get_eigen_values(j_gf3,"gf3")
e_gf4=get_eigen_values(j_gf4,"gf4")

##plot eigenvalues##
plot_eigen_value_evolution(e_gf1[[1]],"gf1")
plot_eigen_value_evolution(e_gf2[[1]],"gf2")
plot_eigen_value_evolution(e_gf3[[1]],"gf3")
plot_eigen_value_evolution(e_gf4[[1]],"gf4")

####detegfine alternative phases using PCA
gf=list(e_gf1[[2]],e_gf2[[2]],e_gf3[[2]],e_gf4[[2]])
gf_eigen_kpca=plot_kpca_for_sample(gf,c(1,1,1,1),"gf")[[2]]

find_significant_changepoints_penalties(gf_eigen_kpca[[1]],pen_value = 30,
                                        flip_kpca2 = FALSE,flip_kpca1 = TRUE)$OverallResults
find_significant_changepoints_penalties(gf_eigen_kpca[[2]],pen_value = 30,
                                        flip_kpca2 = TRUE,flip_kpca1 = TRUE)$OverallResults
find_significant_changepoints_penalties(gf_eigen_kpca[[3]],pen_value = 30,
                                        flip_kpca2 = FALSE,flip_kpca1 = FALSE)$OverallResults
find_significant_changepoints_penalties(gf_eigen_kpca[[4]],pen_value = 30,
                                        flip_kpca2 = TRUE,flip_kpca1 = FALSE)$OverallResults
