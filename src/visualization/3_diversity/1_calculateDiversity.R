## CODE BLOCK
## Louis Gauthier & Melis Gencel

#######################    READ ME!   #####################################

##### Compute barcode diversity for all samples
##### Create dataframes with diversities per cohort

#######################    ^^^^^^^^   #####################################

# CALCULATE DIVERSITY (q = 0): number of lineages with nonzero frequency (species richness)
calculate_q_0 <- function(mat) {
  matbool = mat
  matbool[] = TRUE
  matbool[mat==0] <- FALSE
  q_0 = as.data.frame(colSums(matbool))
  colnames(q_0)="q_0"
  return(q_0)
}

# CALCULATE DIVERSITY (q = 1): Shannon diversity
calculate_q_1 <- function(mat) {
  matbool = mat
  matbool[] = TRUE
  matbool[mat==0] <- FALSE
  q_1 = as.data.frame(exp(sapply(mat, function(x) entropy.empirical(x,unit = "log"))))
  colnames(q_1)="q_1"
  return(q_1)
}

# CALCULATE DIVERSITY (q = infinity): reciprocal of the maximum lineage frequency
calculate_q_inf <- function(mat) {
  colMax <- function(data) sapply(data, max, na.rm = TRUE)
  q_inf = as.data.frame(1/colMax(mat))
  colnames(q_inf)="q_inf"
  return(q_inf)
}

# GENERATE DIVERSITY DATAFRAME FOR ONE SIMULATION
calculate_diversity <- function(generations,maxgen=12,step=1,sample=FALSE){
  
  m = as.matrix(generations[,-1])
  
  if(sample){
    m =floor(m*0.8)
  }
  
  mat = as.data.frame(sweep(m,2,colSums(m,na.rm = TRUE),`/`))
  mat$ID = generations$X1
  
  
  q0 = calculate_q_0(mat)
  q1 = calculate_q_1(mat)
  qinf = calculate_q_inf(mat)
  qall = cbind(q0,q1,qinf)
  qall$Generations = seq(from=0,to=maxgen,by=step)
  return(qall)
}
#################

format_sample <- function(sample){  
  casted = reshape2::dcast(sample, ID ~ Time, value.var = 'Reads')

  casted[is.na(casted)] <- 0
  return(casted)
}








