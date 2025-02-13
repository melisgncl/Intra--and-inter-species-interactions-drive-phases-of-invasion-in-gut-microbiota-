### Louis Gauthier & Melis Gencel

#######################    READ ME!   #####################################

### read 16S data and create color profiles for families

#######################    ^^^^^^^^   #####################################

### 1. import   16S data

source("~/Desktop/mouse_barcoding_last/src/visualization/0_config/0_config.R")
library(tidyr)

im1_16S_abundances <- read_csv("data/16S/formatted/im1_16S_abundances.csv")
im2_16S_abundances <- read_csv("data/16S/formatted/im2_16S_abundances.csv")
im3_16S_abundances <- read_csv("data/16S/formatted/im3_16S_abundances.csv")
im4_16S_abundances <- read_csv("data/16S/formatted/im4_16S_abundances.csv")

rm1_16S_abundances <- read_csv("data/16S/formatted/rm1_16S_abundances.csv")
rm2_16S_abundances <- read_csv("data/16S/formatted/rm2_16S_abundances.csv")
rm3_16S_abundances <- read_csv("data/16S/formatted/rm3_16S_abundances.csv")
rm4_16S_abundances <- read_csv("data/16S/formatted/rm4_16S_abundances.csv")

###CONTROL_SAMPLE####
###control group
nc1_16S_abundances <- read_csv("data/16S/formatted//nc1_16S_abundances.csv")
nc2_16S_abundances <- read_csv("data/16S/formatted//nc2_16S_abundances.csv")
nc3_16S_abundances <- read_csv("data/16S/formatted//nc3_16S_abundances.csv")
nc4_16S_abundances <- read_csv("data/16S/formatted//nc4_16S_abundances.csv")



##to assign colors and recode factor levels so that all taxa share the same ordering & color across samples for the future
## these are the HARDCODED taxa we want to order and color
####color codes##
genus.colors = c("Escherichia/Shigella" = "#003D18",
                 "Desulfovibrio" = "#00A341",
                 "Paenibacillus" = "#561AA3",
                 "Lachnospiraceae_NK4A136_group" = "#063D74",
                 "Shuttleworthia" = "#060874",
                 "Acetatifactor" = "#0A66C2",
                 "Lachnospiraceae_UCG-001" = "#1685F3",
                 "Clostridium_sensu_stricto_1" = "#51A3F6",
                 "Ruminiclostridium_9" = "#08519C",
                 "Romboutsia" = "#77B8F8",
                 "Bacteroides" = "#BD0026",
                 "Prevotellaceae_UCG-001" = "#FD8D3C",
                 "Alloprevotella" = "#FEB24C",
                 "Alistipes" = "#FED976",
                 "Methanobacterium" = "#a6a600",
                 "other" = "black","Parasutterella"="#571150","Ruminococcaceae_UCG-005"="#ab1e96",
                 "Lachnoclostridium"="#8bbada", "Ruminiclostridium_6"="#48cb49","Lactobacillus"="#2ec7f4",
                 "Akkermansia" ="#a37110","Anaeroplasma"="#fee24b","Veillonella"="#2a01ae","Oscillibacter"="#e17d68","")


major.genus = c("Escherichia/Shigella","Desulfovibrio","Paenibacillus","Lachnospiraceae_NK4A136_group","Acetatifactor","Lachnospiraceae_UCG-001",
                "Clostridium_sensu_stricto_1","Bacteroides","Ruminiclostridium_9","Shuttleworthia","Prevotellaceae_UCG-001",
                "Alloprevotella","Alistipes","Romboutsia","Methanobacterium","other","Parasutterella","Ruminiclostridium_6","Lachnoclostridium","Ruminococcaceae_UCG-005",
                "Lactobacillus","Akkermansia","Anaeroplasma" ,"Veillonella","Oscillibacter")





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

####these are the families on low frequency we assign them random colors
bottom_genera=read_csv("src/visualization/0_config/bottom_genera_color.csv")
family_bottom=as.character(bottom_genera$hex)
names(family_bottom)=as.character(bottom_genera$all.main.families)
family.colors=c(family.colors,family_bottom)

major.family=c("Enterobacteriaceae","Moraxellaceae","Desulfovibrionaceae","Lachnospiraceae","Ruminococcaceae",
               "Peptostreptococcaceae","Peptococcaceae","Clostridiales_vadinBB60_group",
               "Clostridiaceae_1","Paenibacillaceae","Lactobacillaceae","Muribaculaceae",
               "Bacteroidaceae","Rikenellaceae","Prevotellaceae","Marinifilaceae",
               "Tannerellaceae","Anaeroplasmataceae", "Deferribacteraceae", "Saccharimonadaceae",
               "Xiphinematobacteraceae","Methanobacteriaceae","Erysipelotrichaceae","other","Acholeplasmataceae","Oscillospiraceae", "Akkermansiaceae","[Eubacterium] coprostanoligenes group",
               "Enterococcaceae","Bacillaceae","Corynebacteriaceae","Microbacteriaceae","Erwiniaceae","Veillonellaceae","Sutterellaceae","Comamonadaceae","Xanthomonadaceae",
               "Leptotrichiaceae","Pseudomonadaceae","Sphingomonadaceae","Weeksellaceae")

phylum.colors = c("Bacteroidetes" = "#BD0026",
                  "Firmicutes" = "#0A66C2",
                  "Proteobacteria" = "#00A341",
                  "Deferribacteres" = "#ff028d",
                  "Patescibacteria" = "#A18E66",
                  "Tenericutes" = "#633228")

major.phylum = c("Proteobacteria","Firmicutes","Bacteroidetes","Patescibacteria","Tenericutes","Deferribacteres")






