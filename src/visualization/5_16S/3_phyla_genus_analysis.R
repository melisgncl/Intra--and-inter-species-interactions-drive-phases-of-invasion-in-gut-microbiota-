### Louis Gauthier & Melis Gencel
### read 16S data and plot abundances of phyla and genus levels 

#######################    READ ME!   #####################################

##### Calculate relative abundances at the level of phyla and genus

#######################    ^^^^^^^^   #####################################

### 1. group by time+phyla and sum abundances


im1_16S.phyla = im1_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
im2_16S.phyla = im2_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
im3_16S.phyla = im3_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
im4_16S.phyla = im4_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))


rm1_16S.phyla = rm1_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
rm2_16S.phyla = rm2_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
rm3_16S.phyla = rm3_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
rm4_16S.phyla = rm4_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))


nc1_16S.phyla = nc1_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
nc2_16S.phyla = nc2_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
nc3_16S.phyla = nc3_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))
nc4_16S.phyla = nc4_16S_abundances %>% group_by(Time,Phylum) %>% summarise(Abundance.phyla = sum(Abundance.relative))


#### 2. group by time+genus and sum abundances

im1_16S.genus = im1_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
im2_16S.genus = im2_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
im3_16S.genus = im3_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
im4_16S.genus = im4_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))


rm1_16S.genus = rm1_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
rm2_16S.genus = rm2_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
rm3_16S.genus = rm3_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
rm4_16S.genus = rm4_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))

nc1_16S.genus = nc1_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
nc2_16S.genus = nc2_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
nc3_16S.genus = nc3_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))
nc4_16S.genus = nc4_16S_abundances %>% group_by(Time,Genus) %>% summarise(Abundance.genus = sum(Abundance.relative))



###3 order them according to major groups### abundance higher than 0.05 are top.genera

#########///FAMILY///#############
im1_16S.family$Family = ordered(im1_16S.family$Family,levels=c(major.family))
im2_16S.family$Family = ordered(im2_16S.family$Family,levels=c(major.family))
im3_16S.family$Family = ordered(im3_16S.family$Family,levels=c(major.family))
im4_16S.family$Family = ordered(im4_16S.family$Family,levels=c(major.family))

rm1_16S.family$Family = ordered(rm1_16S.family$Family,levels=c(major.family))
rm2_16S.family$Family = ordered(rm2_16S.family$Family,levels=c(major.family))
rm3_16S.family$Family = ordered(rm3_16S.family$Family,levels=c(major.family))
rm4_16S.family$Family = ordered(rm4_16S.family$Family,levels=c(major.family))

nc1_16S.family$Family = ordered(nc1_16S.family$Family,levels=c(major.family))
nc2_16S.family$Family = ordered(nc2_16S.family$Family,levels=c(major.family))
nc3_16S.family$Family = ordered(nc3_16S.family$Family,levels=c(major.family))
nc4_16S.family$Family = ordered(nc4_16S.family$Family,levels=c(major.family))



#########///PHYLA///#############
im1_16S.phyla$Phylum = ordered(im1_16S.phyla$Phylum,levels=c(major.phylum))
im2_16S.phyla$Phylum = ordered(im2_16S.phyla$Phylum,levels=c(major.phylum))
im3_16S.phyla$Phylum = ordered(im3_16S.phyla$Phylum,levels=c(major.phylum))
im4_16S.phyla$Phylum = ordered(im4_16S.phyla$Phylum,levels=c(major.phylum))

rm1_16S.phyla$Phylum = ordered(rm1_16S.phyla$Phylum,levels=c(major.phylum))
rm2_16S.phyla$Phylum = ordered(rm2_16S.phyla$Phylum,levels=c(major.phylum))
rm3_16S.phyla$Phylum = ordered(rm3_16S.phyla$Phylum,levels=c(major.phylum))
rm4_16S.phyla$Phylum = ordered(rm4_16S.phyla$Phylum,levels=c(major.phylum))

nc1_16S.phyla$Phylum = ordered(nc1_16S.phyla$Phylum,levels=c(major.phylum))
nc2_16S.phyla$Phylum = ordered(nc2_16S.phyla$Phylum,levels=c(major.phylum))
nc3_16S.phyla$Phylum = ordered(nc3_16S.phyla$Phylum,levels=c(major.phylum))
nc4_16S.phyla$Phylum = ordered(nc4_16S.phyla$Phylum,levels=c(major.phylum))



#########///GENUS///#############
im1_16S.genus$Genus = ordered(im1_16S.genus$Genus,levels=c(major.genus))
im2_16S.genus$Genus = ordered(im2_16S.genus$Genus,levels=c(major.genus))
im3_16S.genus$Genus = ordered(im3_16S.genus$Genus,levels=c(major.genus))
im4_16S.genus$Genus = ordered(im4_16S.genus$Genus,levels=c(major.genus))

rm1_16S.genus$Genus = ordered(rm1_16S.genus$Genus,levels=c(major.genus))
rm2_16S.genus$Genus = ordered(rm2_16S.genus$Genus,levels=c(major.genus))
rm3_16S.genus$Genus = ordered(rm3_16S.genus$Genus,levels=c(major.genus))
rm4_16S.genus$Genus = ordered(rm4_16S.genus$Genus,levels=c(major.genus))


nc1_16S.genus$Genus = ordered(nc1_16S.genus$Genus,levels=c(major.genus))
nc2_16S.genus$Genus = ordered(nc2_16S.genus$Genus,levels=c(major.genus))
nc3_16S.genus$Genus = ordered(nc3_16S.genus$Genus,levels=c(major.genus))
nc4_16S.genus$Genus = ordered(nc4_16S.genus$Genus,levels=c(major.genus))


