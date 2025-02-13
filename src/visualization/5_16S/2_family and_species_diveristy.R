### Louis Gauthier & Melis Gencel
#######################    READ ME!   #####################################

#Compute and plot 16S diversities per sample and cohort
#Diversity function comes from 3_diversity###
#calulate diveristy at the level of family and ASV
#######################    ^^^^^^^^   #####################################

##### Compute and plot 16S diversities per sample and cohort

#######################    ^^^^^^^^   #####################################
 ##source the diveristy function from 3_diveristy###
source("~/Desktop/mouse_barcoding_last/src/visualization/3_diversity/1_calculateDiversity.R")



format_OTU <- function(sample){  
  casted = reshape2::dcast(sample, Family ~ Time, value.var = 'Abundance')
  casted[is.na(casted)] <- 0
  return(casted)
}

#####################################################

div.im1 = calculate_diversity(format_OTU(im1_16S.family),maxgen = 8,step = 1)
div.im1$Generations = as.double(row.names(div.im1))
div.im2 = calculate_diversity(format_OTU(im2_16S.family),maxgen = 8,step = 1)
div.im2$Generations = as.double(row.names(div.im2))
div.im3 = calculate_diversity(format_OTU(im3_16S.family),maxgen = 8,step = 1)
div.im3$Generations = as.double(row.names(div.im3))
div.im4 = calculate_diversity(format_OTU(im4_16S.family),maxgen = 8,step = 1)
div.im4$Generations = as.double(row.names(div.im4))

div.im1$Sample = "im1"
div.im2$Sample = "im2"
div.im3$Sample = "im3"
div.im4$Sample = "im4"
div.im1to4 = rbind(div.im1,div.im2,div.im3,div.im4)

type="im"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]


im1to4.q0 = ggplot(div.im1to4) + geom_line(aes(Generations,log10(q_0),color=Sample),size=1.5) + 
  theme_Publication() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) + 
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal1) + guides(color = guide_legend(override.aes = list(size=8,shape=15))) +
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)}) +ylab("")
p = im1to4.q0 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits)+
coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/im_q0_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

im1to4.q1 = ggplot(div.im1to4) + geom_line(aes(Generations,log10(q_1),color=Sample),size=1.5) + 
  theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal1) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))+
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)})
p = im1to4.q1 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/im_q1_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

im1to4.qinf = ggplot(div.im1to4) + geom_line(aes(Generations,log10(q_inf),color=Sample),size=1.5)  + 
  theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal1) + guides(color = guide_legend(override.aes = list(size=8,shape=15))) +
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)})
p = im1to4.qinf + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/im_qinf_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

im1to4.plot = ggarrange(im1to4.q0,im1to4.q1,im1to4.qinf, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
ggsave(im1to4.plot,filename = "reports/figures/16S/im_family_diversity.eps",fonts=c("serif"),width = 15,height = 6,)
#####################################################
##The only 13th time point is removed for the sake of analysis 
rm1_16S.family_wo13=rm1_16S.family[!(rm1_16S.family$Time==13),]
  


div.rm1 = calculate_diversity(format_OTU(rm1_16S.family_wo13),maxgen = 16,step = 1)
div.rm1$Generations = as.double(row.names(div.rm1))
div.rm2 = calculate_diversity(format_OTU(rm2_16S.family),maxgen = 17,step = 1)
div.rm2$Generations = as.double(row.names(div.rm2))
div.rm3 = calculate_diversity(format_OTU(rm3_16S.family),maxgen = 17,step = 1)
div.rm3$Generations = as.double(row.names(div.rm3))
div.rm4 = calculate_diversity(format_OTU(rm4_16S.family),maxgen = 17,step = 1)
div.rm4$Generations = as.double(row.names(div.rm4))

div.rm1$Sample = "rm1"
div.rm2$Sample = "rm2"
div.rm3$Sample = "rm3"
div.rm4$Sample = "rm4"
div.rm1to4 = rbind(div.rm1,div.rm2,div.rm3,div.rm4)

type="rm"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]


rm1to4.q0 = ggplot(div.rm1to4) + geom_line(aes(Generations,log10(q_0),color=Sample),size=1.5) + 
  theme_Publication() + xlab("Trme (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) + 
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal2) + guides(color = guide_legend(override.aes = list(size=8,shape=15))) +
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)}) +ylab("")
p = rm1to4.q0 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits)+
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/rm_q0_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

rm1to4.q1 = ggplot(div.rm1to4) + geom_line(aes(Generations,log10(q_1),color=Sample),size=1.5) + 
  theme_Publication_noYaxis() + xlab("Trme (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal2) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))+
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)})
p = rm1to4.q1 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/rm_q1_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

rm1to4.qinf = ggplot(div.rm1to4) + geom_line(aes(Generations,log10(q_inf),color=Sample),size=1.5)  + 
  theme_Publication_noYaxis() + xlab("Trme (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal2) + guides(color = guide_legend(override.aes = list(size=8,shape=15))) +
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)})
p = rm1to4.qinf + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/rm_qinf_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

rm1to4.plot = ggarrange(rm1to4.q0,rm1to4.q1,rm1to4.qinf, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
ggsave(rm1to4.plot,filename = "reports/figures/16S/rm_family_diversity.eps",fonts=c("serif"),width = 15,height = 6)

#####################################################

div.nc1 = calculate_diversity(format_OTU(nc1_16S.family),maxgen = 9,step = 1)
div.nc1$Generations = as.double(row.names(div.nc1))
div.nc2 = calculate_diversity(format_OTU(nc2_16S.family),maxgen = 9,step = 1)
div.nc2$Generations = as.double(row.names(div.nc2))
div.nc3 = calculate_diversity(format_OTU(nc3_16S.family),maxgen = 7,step = 1)
div.nc3$Generations = as.double(row.names(div.nc3))
div.nc4 = calculate_diversity(format_OTU(nc4_16S.family),maxgen = 9,step = 1)
div.nc4$Generations = as.double(row.names(div.nc4))

div.nc1$Sample = "nc1"
div.nc2$Sample = "nc2"
div.nc3$Sample = "nc3"
div.nc4$Sample = "nc4"
div.nc1to4 = rbind(div.nc1,div.nc2,div.nc3,div.nc4)

type="nc"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]


nc1to4.q0 = ggplot(div.nc1to4) + geom_line(aes(Generations,log10(q_0),color=Sample),size=1.5) + 
  theme_Publication() + xlab("Tnce (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) + 
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal4) + guides(color = guide_legend(override.aes = list(size=8,shape=15))) +
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)}) +ylab("")
p = nc1to4.q0 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits)+
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/nc_q0_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

nc1to4.q1 = ggplot(div.nc1to4) + geom_line(aes(Generations,log10(q_1),color=Sample),size=1.5) + 
  theme_Publication_noYaxis() + xlab("Tnce (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal4) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))+
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)})
p = nc1to4.q1 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/nc_q1_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

nc1to4.qinf = ggplot(div.nc1to4) + geom_line(aes(Generations,log10(q_inf),color=Sample),size=1.5)  + 
  theme_Publication_noYaxis() + xlab("Tnce (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE)  + scale_color_manual(values = pal4) + guides(color = guide_legend(override.aes = list(size=8,shape=15))) +
  scale_y_continuous(breaks=c(0,1),limits=c(0,2),labels=function(n){format(10^n, scientific = TRUE)})
p = nc1to4.qinf + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = "reports/figures/16S/nc_qinf_family_diversity.eps",fonts=c("serif"),width = 8.25,height = 6,device = cairo_ps())

nc1to4.plot = ggarrange(nc1to4.q0,nc1to4.q1,nc1to4.qinf, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
ggsave(nc1to4.plot,filename = "reports/figures/16S/nc_family_diversity.eps",fonts=c("serif"),width = 15,height = 6)




