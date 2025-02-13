### Louis Gauthier & Melis Gencel

#######################    READ ME!   #####################################

##### Calculate relative abundances at the level of family#######
####Plot diveristy and freqeuncies before filtering them for the further analysis##

#######################    ^^^^^^^^   #####################################

### 1. calculate relative abundance per time point

im1_16S_abundances = data.table(im1_16S_abundances)
im2_16S_abundances = data.table(im2_16S_abundances)
im3_16S_abundances = data.table(im3_16S_abundances)
im4_16S_abundances = data.table(im4_16S_abundances)

im1_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
im2_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
im3_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
im4_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]

rm1_16S_abundances = data.table(rm1_16S_abundances)
rm2_16S_abundances = data.table(rm2_16S_abundances)
rm3_16S_abundances = data.table(rm3_16S_abundances)
rm4_16S_abundances = data.table(rm4_16S_abundances)

rm1_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
rm2_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
rm3_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
rm4_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]



##control

nc1_16S_abundances = data.table(nc1_16S_abundances)
nc2_16S_abundances = data.table(nc2_16S_abundances)
nc3_16S_abundances = data.table(nc3_16S_abundances)
nc4_16S_abundances = data.table(nc4_16S_abundances)


nc1_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
nc2_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
nc3_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]
nc4_16S_abundances[, Abundance.relative := Abundance/sum(Abundance), by="Time"]




### 3. group by time+family and sum abundances

im1_16S.family = im1_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
im2_16S.family = im2_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
im3_16S.family = im3_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
im4_16S.family = im4_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))


rm1_16S.family = rm1_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
rm2_16S.family = rm2_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
rm3_16S.family = rm3_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
rm4_16S.family = rm4_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))



nc1_16S.family = nc1_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
nc2_16S.family = nc2_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
nc3_16S.family = nc3_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))
nc4_16S.family = nc4_16S_abundances %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.relative),Abundance=sum(Abundance))







####PLOT the diveristy and frequencies at family level without filtering##
source("~/Desktop/mouse_barcoding_last/src/visualization/5_16S/2_family and_species_diveristy.R")

source("~/Desktop/mouse_barcoding_last/src/visualization/5_16S/4_family_plot.R")






### 4. organize by group, filtering by minimum mean frequency for the rest of the analysis

combined.im  = rbind(im1_16S.family,im2_16S.family,im3_16S.family,im4_16S.family)

im.mean.abundances = aggregate(. ~ Family, combined.im[,-1], mean)
im.mean.abundances = im.mean.abundances[order(im.mean.abundances$Abundance.family,decreasing = TRUE),]
im.main.families = im.mean.abundances[im.mean.abundances$Abundance.family>0.001,]$Family


combined.rm = rbind(rm1_16S.family,rm2_16S.family,rm3_16S.family,rm4_16S.family)
rm.mean.abundances = aggregate(. ~ Family, combined.rm[,-1], mean)
rm.mean.abundances = rm.mean.abundances[order(rm.mean.abundances$Abundance.family,decreasing = TRUE),]
rm.main.families = rm.mean.abundances[rm.mean.abundances$Abundance.family>0.001,]$Family


combined.nc = rbind(nc1_16S.family,nc2_16S.family,nc3_16S.family,nc4_16S.family)
nc.mean.abundances = aggregate(. ~ Family, combined.nc[,-1], mean)
nc.mean.abundances = nc.mean.abundances[order(nc.mean.abundances$Abundance.family,decreasing = TRUE),]
nc.main.families = nc.mean.abundances[nc.mean.abundances$Abundance.family>0.001,]$Family


### these are the families that should get named colors
### the remaining families will be set to "other"
## finally we regroup the abundances of "other" as a single value

all.families = union(im.mean.abundances$Family,c(rm.mean.abundances$Family,nc.mean.abundances$Family))



im1_16S.family[!im1_16S.family$Family %in% im.main.families,]$Family="other"
im2_16S.family[!im2_16S.family$Family %in% im.main.families,]$Family="other"
im3_16S.family[!im3_16S.family$Family %in% im.main.families,]$Family="other"
im4_16S.family[!im4_16S.family$Family %in% im.main.families,]$Family="other"

im1_16S.family = im1_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
im2_16S.family = im2_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
im3_16S.family = im3_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
im4_16S.family = im4_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))

rm1_16S.family[!rm1_16S.family$Family %in% rm.main.families,]$Family="other"
rm2_16S.family[!rm2_16S.family$Family %in% rm.main.families,]$Family="other"
rm3_16S.family[!rm3_16S.family$Family %in% rm.main.families,]$Family="other"
rm4_16S.family[!rm4_16S.family$Family %in% rm.main.families,]$Family="other"

rm1_16S.family = rm1_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
rm2_16S.family = rm2_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
rm3_16S.family = rm3_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
rm4_16S.family = rm4_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))



nc1_16S.family[!nc1_16S.family$Family %in% nc.main.families,]$Family="other"
nc2_16S.family[!nc2_16S.family$Family %in% nc.main.families,]$Family="other"
nc3_16S.family[!nc3_16S.family$Family %in% nc.main.families,]$Family="other"
nc4_16S.family[!nc4_16S.family$Family %in% nc.main.families,]$Family="other"

nc1_16S.family = nc1_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
nc2_16S.family = nc2_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
nc3_16S.family = nc3_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))
nc4_16S.family = nc4_16S.family %>% group_by(Time,Family) %>% summarise(Abundance.family = sum(Abundance.family),Abundance=sum(Abundance))


####### FOR Co-clustering ANALYSIS ########

im1_16S.df = im1_16S.family[im1_16S.family$Family %in% im.main.families,]
im2_16S.df = im2_16S.family[im2_16S.family$Family %in% im.main.families,]
im3_16S.df = im3_16S.family[im3_16S.family$Family %in% im.main.families,]
im4_16S.df = im4_16S.family[im4_16S.family$Family %in% im.main.families,]


im1_16S.df$Time=as.integer(im1_16S.df$Time)
im2_16S.df$Time=as.integer(im2_16S.df$Time)
im3_16S.df$Time=as.integer(im3_16S.df$Time)
im4_16S.df$Time=as.integer(im4_16S.df$Time)

# cast in foimat suitable for co-clustering
im1_16S.long = reshape2::dcast(im1_16S.df, Time ~ Family, value.var = 'Abundance.family')
im2_16S.long = reshape2::dcast(im2_16S.df, Time ~ Family, value.var = 'Abundance.family')
im3_16S.long = reshape2::dcast(im3_16S.df, Time ~ Family, value.var = 'Abundance.family')
im4_16S.long = reshape2::dcast(im4_16S.df, Time ~ Family, value.var = 'Abundance.family')



write_tsv(im1_16S.long,file = "data/16S/taxa_long_format/im1_16S_taxa.tsv",col_names = TRUE)
write_tsv(im2_16S.long,file = "data/16S/taxa_long_format/im2_16S_taxa.tsv",col_names = TRUE)
write_tsv(im3_16S.long,file = "data/16S/taxa_long_format/im3_16S_taxa.tsv",col_names = TRUE)
write_tsv(im4_16S.long,file = "data/16S/taxa_long_format/im4_16S_taxa.tsv",col_names = TRUE)




rm1_16S.df = rm1_16S.family[rm1_16S.family$Family %in% rm.main.families,]
rm2_16S.df = rm2_16S.family[rm2_16S.family$Family %in% rm.main.families,]
rm3_16S.df = rm3_16S.family[rm3_16S.family$Family %in% rm.main.families,]
rm4_16S.df = rm4_16S.family[rm4_16S.family$Family %in% rm.main.families,]


rm1_16S.df$Time=as.integer(rm1_16S.df$Time)
rm2_16S.df$Time=as.integer(rm2_16S.df$Time)
rm3_16S.df$Time=as.integer(rm3_16S.df$Time)
rm4_16S.df$Time=as.integer(rm4_16S.df$Time)




# cast in format suitable for co-clustering
rm1_16S.long = reshape2::dcast(rm1_16S.df, Time ~ Family, value.var = 'Abundance.family')
rm2_16S.long = reshape2::dcast(rm2_16S.df, Time ~ Family, value.var = 'Abundance.family')
rm3_16S.long = reshape2::dcast(rm3_16S.df, Time ~ Family, value.var = 'Abundance.family')
rm4_16S.long = reshape2::dcast(rm4_16S.df, Time ~ Family, value.var = 'Abundance.family')
##The only 13th time point is interpolated for the sake of analysis 
rm1_13=names(rm1_16S.long)
rm1_13=(rm1_16S.long [12,]+rm1_16S.long[13,])/2
rm1_16S.long=rbind(rm1_16S.long,rm1_13)
rm1_16S.long<- rm1_16S.long[order(rm1_16S.long$Time),] 


write_tsv(rm1_16S.long,file = "data/16S/taxa_long_format/rm1_16S_taxa.tsv",col_names = TRUE)
write_tsv(rm2_16S.long,file = "data/16S/taxa_long_format/rm2_16S_taxa.tsv",col_names = TRUE)
write_tsv(rm3_16S.long,file = "data/16S/taxa_long_format/rm3_16S_taxa.tsv",col_names = TRUE)
write_tsv(rm4_16S.long,file = "data/16S/taxa_long_format/rm4_16S_taxa.tsv",col_names = TRUE)



#####control###
nc1_16S.df = nc1_16S.family[nc1_16S.family$Family %in% nc.main.families,]
nc2_16S.df = nc2_16S.family[nc2_16S.family$Family %in% nc.main.families,]
nc3_16S.df = nc3_16S.family[nc3_16S.family$Family %in% nc.main.families,]
nc4_16S.df = nc4_16S.family[nc4_16S.family$Family %in% nc.main.families,]

nc1_16S.long = reshape2::dcast(nc1_16S.df, Time ~ Family, value.var = 'Abundance.family')
nc2_16S.long = reshape2::dcast(nc2_16S.df, Time ~ Family, value.var = 'Abundance.family')
nc3_16S.long = reshape2::dcast(nc3_16S.df, Time ~ Family, value.var = 'Abundance.family')
nc4_16S.long = reshape2::dcast(nc4_16S.df, Time ~ Family, value.var = 'Abundance.family')

write_tsv(nc1_16S.long,file = "data/16S/taxa_long_format/nc1_16S_taxa.tsv",col_names = TRUE)
write_tsv(nc2_16S.long,file = "data/16S/taxa_long_format/nc2_16S_taxa.tsv",col_names = TRUE)
write_tsv(nc3_16S.long,file = "data/16S/taxa_long_format/nc3_16S_taxa.tsv",col_names = TRUE)
write_tsv(nc4_16S.long,file = "data/16S/taxa_long_format/nc4_16S_taxa.tsv",col_names = TRUE)







