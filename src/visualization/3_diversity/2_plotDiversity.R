## CODE BLOCK
## Louis Gauthier & Melis Gencel

#######################    READ ME!   #####################################

### Plot barcode diversities for each sample
### Each panel of each cohort is saved individually
### Then we save a multi-panel figure with all diversities and the legend.

#######################    ^^^^^^^^   #####################################

outdir=paste0("reports/figures/")

if (!dir.exists(paste0(outdir, "/diversity"))) {
  dir.create(paste0(outdir, "/diversity"), recursive = TRUE)
}

#####################################################

div.gf1 = calculate_diversity(format_sample(Sample_gf1),maxgen = 14,step = 1)
div.gf1$Generations = as.double(row.names(div.gf1))
div.gf2 = calculate_diversity(format_sample(Sample_gf2),maxgen = 14,step = 1)
div.gf2$Generations = as.double(row.names(div.gf2))
div.gf3 = calculate_diversity(format_sample(Sample_gf3),maxgen = 14,step = 1)
div.gf3$Generations = as.double(row.names(div.gf3))
div.gf4 = calculate_diversity(format_sample(Sample_gf4),maxgen = 14,step = 1)
div.gf4$Generations = as.double(row.names(div.gf4))


div.gf1$Sample = "gf1"
div.gf2$Sample = "gf2"
div.gf3$Sample = "gf3"
div.gf4$Sample = "gf4"
div.gf1to4 = rbind(div.gf1,div.gf2,div.gf3,div.gf4)





#####################################################

div.rm1 = calculate_diversity(format_sample(Sample_rm1),maxgen = 17,step = 1)
div.rm1$Generations = as.double(row.names(div.rm1))
div.rm2 = calculate_diversity(format_sample(Sample_rm2),maxgen = 16,step = 1)
div.rm2$Generations = as.double(row.names(div.rm2))
div.rm3 = calculate_diversity(format_sample(Sample_rm3),maxgen = 16,step = 1)
div.rm3$Generations = as.double(row.names(div.rm3))
div.rm4 = calculate_diversity(format_sample(Sample_rm4),maxgen = 17,step = 1)
div.rm4$Generations = as.double(row.names(div.rm4))


div.rm1$Sample = "rm1"
div.rm2$Sample = "rm2"
div.rm3$Sample = "rm3"
div.rm4$Sample = "rm4"
div.rm1to4 = rbind(div.rm1,div.rm2,div.rm3,div.rm4)



#####################################################



div.im1 = calculate_diversity(format_sample(Sample_im1),maxgen = 5,step = 1)
div.im1$Generations = as.double(row.names(div.im1))
div.im2 = calculate_diversity(format_sample(Sample_im2),maxgen = 5,step = 1)
div.im2$Generations = as.double(row.names(div.im2))
div.im3 = calculate_diversity(format_sample(Sample_im3),maxgen = 7,step = 1)
div.im3$Generations = as.double(row.names(div.im3))
div.im4 = calculate_diversity(format_sample(Sample_im4),maxgen = 7,step = 1)
div.im4$Generations = as.double(row.names(div.im4))


div.im1$Sample = "im1"
div.im2$Sample = "im2"
div.im3$Sample = "im3"
div.im4$Sample = "im4"
div.im1to4 = rbind(div.im1,div.im2,div.im3,div.im4)


#####PLOT####





###im###
type="im"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]


im1to4.q0 = ggplot(div.im1to4) + geom_line(aes(Generations,log10(q_0),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
  theme_Publication() + xlab("Time (post-gavage)") +ylab("diversity")+ scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) + 
  coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal1) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))

p = im1to4.q0 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = paste0(outdir,"/diversity/im_q0_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

im1to4.q1 = ggplot(div.im1to4) + geom_line(aes(Generations,log10(q_1),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
  theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal1) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))
p = im1to4.q1 + guides(color = FALSE) + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = paste0(outdir,"/diversity/im_q1_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

im1to4.qinf = ggplot(div.im1to4) + geom_line(aes(Generations,log10(q_inf),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
  theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal1) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))
p = im1to4.qinf + guides(color = FALSE) + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE) 
ggsave(p,filename = paste0(outdir,"//diversity/im_qinf_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

im1to4.plot = ggarrange(im1to4.q0,im1to4.q1,im1to4.qinf, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
ggsave(im1to4.plot,filename = paste0(outdir,"//diversity/im_diversity_all.eps"),fonts=c("serif"),width = 21,height = 6)


###rm###
type="rm"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]


rm1to4.q0 = ggplot(div.rm1to4) + geom_line(aes(Generations,log10(q_0),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
  theme_Publication() + xlab("Time (post-gavage)") +ylab("diversity")+ scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) + 
  coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal3) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))

p = rm1to4.q0 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = paste0(outdir,"//diversity/rm_q0_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

rm1to4.q1 = ggplot(div.rm1to4) + geom_line(aes(Generations,log10(q_1),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
  theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal3) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))
p = rm1to4.q1 + guides(color = FALSE) + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE)
ggsave(p,filename = paste0(outdir,"//diversity/rm_q1_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

rm1to4.qinf = ggplot(div.rm1to4) + geom_line(aes(Generations,log10(q_inf),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
  theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
  coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal3) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))
p = rm1to4.qinf + guides(color = FALSE) + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
  coord_cartesian(expand = TRUE) 
ggsave(p,filename = paste0(outdir,"//diversity/rm_qinf_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

rm1to4.plot = ggarrange(rm1to4.q0,rm1to4.q1,rm1to4.qinf, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
ggsave(rm1to4.plot,filename = paste0(outdir,"//diversity/rm_diversity_all.eps"),fonts=c("serif"),width = 21,height = 6)








###gf###

type="gf"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]


gf1to4.q0 = ggplot(div.gf1to4) + geom_line(aes(Generations,log10(q_0),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
            theme_Publication() + xlab("Time (post-gavage)") +ylab("diversity")+ scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) + 
            coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal3) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))

p = gf1to4.q0 + guides(color = FALSE) + theme_Publication_noYaxis() + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
            coord_cartesian(expand = TRUE)
ggsave(p,filename = paste0(outdir,"//diversity/gf_q0_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

gf1to4.q1 = ggplot(div.gf1to4) + geom_line(aes(Generations,log10(q_1),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
            theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
            coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal3) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))
p = gf1to4.q1 + guides(color = FALSE) + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
            coord_cartesian(expand = TRUE)
ggsave(p,filename = paste0(outdir,"//diversity/gf_q1_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

gf1to4.qinf = ggplot(div.gf1to4) + geom_line(aes(Generations,log10(q_inf),color=Sample),size=1.5) + geom_hline(yintercept =3,linetype="dashed") + 
              theme_Publication_noYaxis() + xlab("Time (post-gavage)") + scale_x_continuous(breaks = effective.breaks,labels = effective.labels, limits = effective.limits) +
              coord_cartesian(expand = FALSE) + scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6),labels=function(n){format(10^n, scientific = TRUE)}) + scale_color_manual(values = pal3) + guides(color = guide_legend(override.aes = list(size=8,shape=15)))
p = gf1to4.qinf + guides(color = FALSE) + scale_x_continuous(breaks = effective.breaks,labels = effective.labels,limits=effective.limits) + 
            coord_cartesian(expand = TRUE) 
ggsave(p,filename = paste0(outdir,"//diversity/gf_qinf_diversity.eps"),fonts=c("serif"),width = 8.25,height = 6)

gf1to4.plot = ggarrange(gf1to4.q0,gf1to4.q1,gf1to4.qinf, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
ggsave(gf1to4.plot,filename = paste0(outdir,"//diversity/gf_diversity_all.eps"),fonts=c("serif"),width = 21,height = 6)
############################

