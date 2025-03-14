## CODE BLOCK
## Louis Gauthier &Melis Gencel
###this config for article figures and main functions

#######################    READ ME!   #####################################

# Here we define global variables and functions that are used by one or multiple code blocks 

########################## 0. SET WORKDIR  ################################


setwd("/Users/melis///Dropbox/mouse_gut/github/")
########################## 1. PACKAGES  ###################################

data_folder="figures"
outdir=paste0("reports/",data_folder)

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}


library(readr)
library(ggplot2)
library(scales)
library(ggpubr)
library(entropy)
library(reshape2)
library(data.table)
library(dplyr)
library(ggnewscale)
library(UpSetR)
library(grid)
library(cowplot)
library(dtwclust)
library(dendextend)
library(ggdendro)
library(gplots)
library(viridis)
library(hash)
library(gridGraphics)

########################## 2. PLOT THEMES #################################

theme_Publication <- function(base_size=32, base_family="Helvetica",aspect.ratio = 0.75) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill="#FCFCFC"),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size=1.5),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            aspect.ratio=aspect.ratio
    ))
}

theme_Publication_noYaxis <- function(base_size=32, base_family="Helvetica",aspect.ratio = 0.75) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill="#FCFCFC"),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size=1.5),
            axis.title = element_text(size = rel(1)),
            axis.title.x = element_text(vjust = -0.2),
            axis.title.y = element_blank(),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            aspect.ratio=aspect.ratio
    ))
}

theme_Publication_bottomLegend <- function(base_size=32, base_family="Helvetica",aspect.ratio = 0.75) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill="#FCFCFC"),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size=1.5),
            axis.title = element_text(size = rel(1)),
            axis.title.x = element_text(vjust = -0.2),
            axis.title.y = element_blank(),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            legend.position=c(0.35,0.15),
            legend.direction="horizontal",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            aspect.ratio=aspect.ratio
    ))
}


theme_Matrix <- function(base_size=32, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill="#FCFCFC"),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size=1.5),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle = 45,hjust=1),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            legend.title = element_text(angle = -90),
            legend.title.align = 0.5,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm")
    ))
}


########################## 3.1 PLOT BREAKS & LABELS ########################
breaks.hash <- hash()
breaks.hash[["gf"]] = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
breaks.hash[["rm"]] =  c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
breaks.hash[["im"]] =  c(1,2,3,4,5,6,7,8,9)
breaks.hash[["nc"]] =  c(1,2,3,4,5,6,7,8,9,10)

# create hash table of breaks
labels.hash <- hash()
labels.hash[["gf"]] = c("3h","","12h","","2d","","","","6d","","","","10d","","","","14d")
labels.hash[["rm"]] = c("3h","","12h","","2d","","","","6d","","","","10d","","","","14d","")
labels.hash[["im"]] = c("3h","","12h","","2d","","","","6d")
labels.hash[["nc"]] = c("1d","","3d","","5d","","7d","","9d","")
# create hash table of breaks
limits.hash <- hash()
limits.hash[["gf"]] = c(1,17)
limits.hash[["rm"]] = c(1,18)
limits.hash[["im"]] = c(1,9)
limits.hash[["nc"]] = c(1,10)


########################## 3.2 COLOR SCHEMES  ############################
pal1=c("#045a8d","#2b8cbe","#74a9cf","#bdc9e1") ##blue ##innate##
pal2 = c("#E8A20C","#FFAF19","#E2522F","#CB6015")##orange hue## ## rm##
pal3=c("#c2e699","#78c679","#31a354","#19804F") #	germ-free			green hue 
pal4 = c("#7a0177","#c51b8a","#9138A7","#552586")##purple ##nc##

# color scheme for cohorts (group color)
c1 = "#bdc9e1"  		
c2="#CB6015"       
c3="#19804F"


eigen.colors=c('#9b2335','#0F4C81','#e15d44', '#6667AB','#aebf9a','#c5a2e8','#e1bc44',
               '#884c5e', '#218EAE', '#009B77', '#fabebe', 
               '#008080', '#e6beff', '#9a6324', '#C96115', '#800000', '#808000', '#ffd8b1', 
               '#00756d',"#B0BA97", "#E940E2", "#69A1AD", "#E447AD", "#84F046" ,"#ADEAA5" ,
               "#E2ECBC", "#BC31EB", "#E29B3F" ,"#904D7E", "#C756D7" ,"#B5DFE9", "#D8537A" ,"#DE898C", "#60ABE4",
               "#9CE2CA", "#8A764C" ,"#E582CF", "#E7BBAC", "#ABE36B")


cluster.colors=c("#3cb44b","#4363d8","#e6194B","#e8ca00","#911eb4","#f58231","#22766dff","#42d4f4","#f032e6","#9A6324",
                 "#2F4C39", "#1D3F6E","#94170f","#665679","#F17829","#97A69C","#606EA9","#A9606E","#A99060","#F8F1AE",
                 "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8","#003D18","#003D18","#82a479","#c74b0e",
                 "#77b5fe","#ccff00")

family.colors = c("Enterobacteriaceae" = "#003D18",
                  "Moraxellaceae" = "#006D2C",
                  "Desulfovibrionaceae" = "#00A341",
                  "Lachnospiraceae" = "#063D74",
                  "Ruminococcaceae" = "#08519C",
                  "Peptostreptococcaceae" = "#0A66C2",
                  "Clostridiaceae_1" = "#77B8F8",
                  "Clostridiales_vadinBB60_group" = "#51A3F6",
                  "Peptococcaceae" = "#1685F3",
                  "Paenibacillaceae" = "#561AA3",
                  "Lactobacillaceae" = "#7E1CFF",
                  "Erysipelotrichaceae" = "#a6007a",
                  "Muribaculaceae" = "#9D0211",
                  "Bacteroidaceae" = "#BD0026",
                  "Rikenellaceae" = "#F03B20",
                  "Prevotellaceae" = "#FD8D3C",
                  "Marinifilaceae" = "#FEB24C",
                  "Tannerellaceae" = "#FED976",
                  "Anaeroplasmataceae" = "#633228",
                  "Deferribacteraceae" = "#ff028d",
                  "Saccharimonadaceae" = "#A18E66",
                  "Xiphinematobacteraceae" = "#CEA2FD",
                  "Methanobacteriaceae" = "#a6a600",
                  "other" = "black","Sutterellaceae"= "#5F021F","Acholeplasmataceae"="#b21561","Oscillospiraceae"="#e17d68",
                  "Akkermansiaceae"="#30bc9e","[Eubacterium] coprostanoligenes group"="#aa862f",
                  "Enterococcaceae" ="#a9cea8","Bacillaceae"  ="#8f7e83","Corynebacteriaceae"="#1f194d",
                  "Microbacteriaceae" ="#8bbada","Erwiniaceae"="#d4e125","Veillonellaceae"="#2a01ae","Comamonadaceae"="#b087b7","Xanthomonadaceae"="#640b66",
                  "Leptotrichiaceae"="#be6a5b","Pseudomonadaceae"="#29e2b2","Sphingomonadaceae"="#fffc78","Weeksellaceae"="#527d54")





s_list=identity_list=c("Enterobacteriaceae",
                       "Moraxellaceae",
                       "Desulfovibrionaceae",
                       "Lachnospiraceae",
                       "Ruminococcaceae",
                       "Peptostreptococcaceae",
                       "Clostridiaceae_1",
                       "Clostridiales_vadinBB60_group",
                       "Peptococcaceae",
                       "Paenibacillaceae",
                       "Lactobacillaceae",
                       "Erysipelotrichaceae",
                       "Muribaculaceae",
                       "Bacteroidaceae",
                       "Rikenellaceae",
                       "Prevotellaceae",
                       "Marinifilaceae",
                       "Tannerellaceae",
                       "Anaeroplasmataceae",
                       "Deferribacteraceae",
                       "Saccharimonadaceae",
                       "Xiphinematobacteraceae",
                       "Methanobacteriaceae",
                       "other","Sutterellaceae","Acholeplasmataceae","Oscillospiraceae",
                       "Akkermansiaceae","[Eubacterium] coprostanoligenes group",
                       "Enterococcaceae" ,"Bacillaceae"  ,"Corynebacteriaceae",
                       "Microbacteriaceae" ,"Erwiniaceae","Veillonellaceae","Comamonadaceae","Xanthomonadaceae",
                       "Leptotrichiaceae","Pseudomonadaceae","Sphingomonadaceae","Weeksellaceae")

c_list=c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")

bottom_genera=read_csv("src/visualization/0_config/bottom_genera_color.csv")
family_bottom=as.character(bottom_genera$hex)
names(family_bottom)=as.character(bottom_genera$all.main.families)
colors_16=c(family.colors,family_bottom)


########################## 4.1 GLOBAL METHODS  ############################

reformatTime <- function(df){
  sample = df
  tryCatch({
    try(sample[sample$Time==0,]$Time <- 0)
    try(sample[sample$Time==3,]$Time <- 1)
    try(sample[sample$Time==6,]$Time <- 2)
    try(sample[sample$Time==12,]$Time <- 3)
    try(sample[sample$Time==24,]$Time <- 4)
    try(sample[sample$Time==384,]$Time <- 19)
    try(sample[sample$Time==360,]$Time <- 18)
    try(sample[sample$Time==336,]$Time <- 17)
    try(sample[sample$Time==312,]$Time <- 16)
    try(sample[sample$Time==288,]$Time <- 15)
    try(sample[sample$Time==264,]$Time <- 14)
    try(sample[sample$Time==240,]$Time <- 13)
    try(sample[sample$Time==216,]$Time <- 12)
    try(sample[sample$Time==192,]$Time <- 11)
    try(sample[sample$Time==168,]$Time <- 10)
    try(sample[sample$Time==144,]$Time <- 9)
    try(sample[sample$Time==120,]$Time <- 8)
    try(sample[sample$Time==96,]$Time <- 7)
    try(sample[sample$Time==72,]$Time <- 6)
    try(sample[sample$Time==48,]$Time <- 5)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(sample)
}

#########define function to format data for plotting######

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
  df[is.na(df$value),]$value = 0
  
  tf = df[order(df$max),]
  
  return(tf)
}




#####to grab the heatmaps####
grab_grob <- function(){
  grid.echo()
  grid.grab()
}


##barcode coloring
all.top.max=read_csv("src/visualization/0_config/all_top_max2.csv")

###topcolors#
top_colors2=read_csv("src/visualization/0_config/top_colors3.csv")
## we create a very long list of colors for low-frequency barcodes. 
## repetition of colors isn't an issue, we just want each barcodes to be colored instead of assigning gray values
long.color.list = rep(top_colors2$hex,50)
long.color.list.random = sample(long.color.list)
