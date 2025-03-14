## Louis & Melis Gencel

#######################    READ ME!   #####################################

##### Plot 16S abundance area plots (relative frequency, Muller-like)
##### Plot 16S abundance line plots (linear)
##### Plot 16S abundance line plots (log)

#######################    ^^^^^^^^   #####################################

linearLinePlot <- function(df,mybreaks,mylabels,mylimits,sample.name){
        mean.abundances = aggregate(. ~ Family,df[,c(-1,-4)], mean)
        mean.abundances = mean.abundances[order(mean.abundances$Abundance.family,decreasing = TRUE),]
        families = mean.abundances[mean.abundances$Abundance.family>0,]$Family
        df=df[df$Family %in% families,] 
   p = ggplot(df) + geom_line(aes(x=Time,y=Abundance.family,group=Family,color=Family),size = 1.5) + 
        ylim(0,1) + scale_x_continuous(limits = mylimits, breaks = mybreaks,labels = mylabels) +
        labs(x = "Time (post-gavage)",y="family composition") + theme_Publication_noYaxis() + 
        scale_color_manual(values=family.colors) +
        guides(color = guide_legend(override.aes = list(size=8,shape=15)))
    q = p + theme_Publication() + guides(color=FALSE) + coord_cartesian(expand = FALSE)
    ggsave(q,filename = paste("reports/figures/16S/",sample.name,"_family_line.jpeg",sep=""),fonts=c("serif"),width = 11,height = 8,dpi=300)
    return(p)
}

logLinePlot <- function(df,mybreaks,mylabels,mylimits,sample.name){
        mean.abundances = aggregate(. ~ Family,df[,c(-1,-4)], mean)
        mean.abundances = mean.abundances[order(mean.abundances$Abundance.family,decreasing = TRUE),]
        families = mean.abundances[mean.abundances$Abundance.family>0,]$Family
        df=df[df$Family %in% families,]
  p = ggplot(df) + geom_line(aes(x=Time,y=Abundance.family,group=Family,color=Family),size = 1.5) + 
        scale_x_continuous(limits = mylimits, breaks = mybreaks,labels = mylabels) +
        labs(x = "Time (post-gavage)",y="family composition") + theme_Publication_noYaxis() + 
        scale_color_manual(values=family.colors) + 
        scale_y_log10(limits=c(1e-7,1e0),breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        guides(color = guide_legend(override.aes = list(size=8,shape=15)))
    q = p + theme_Publication() + guides(color=FALSE) + coord_cartesian(expand = FALSE)
    #ggsave(q,filename = paste("reports/figures/16S/",sample.name,"_family_logline.jpeg",sep=""),fonts=c("serif"),width = 11,height = 8,dpi=300)
    return(q)
}

areaPlot <- function(df,mybreaks,mylabels,mylimits,sample.name){
  mean.abundances = aggregate(. ~ Family,df[,c(-1,-4)], mean)
  mean.abundances = mean.abundances[order(mean.abundances$Abundance.family,decreasing = TRUE),]
  families = mean.abundances[mean.abundances$Abundance.family>0,]$Family
  df=df[df$Family %in% families,]  
  p = ggplot(df) + geom_area(aes(x=Time,y=Abundance.family,group=Family,fill=Family)) + 
        ylim(0,1) + scale_x_continuous(limits = mylimits, breaks = mybreaks,labels = mylabels) +
        labs(x = "Time (post-gavage)",y="family composition") + theme_Publication_noYaxis() + 
        scale_fill_manual(values=family.colors) + coord_cartesian(expand = FALSE) +
        guides(color = guide_legend(override.aes = list(size=8,shape=15)))
    q = p + theme_Publication() + guides(fill=FALSE) + coord_cartesian(expand = FALSE)
    ggsave(q,filename = paste("reports/figures/16S/",sample.name,"_family_area.jpeg",sep=""),fonts=c("serif"),width = 11,height = 8,dpi=300)
    return(p)
}

#######################    ^^^^^^^^   #####################################

type="rm"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]




rm1.family.plot = linearLinePlot(rm1_16S.family_wo13,effective.breaks,effective.labels,effective.limits,"rm1")
rm2.family.plot = linearLinePlot(rm2_16S.family,effective.breaks,effective.labels,effective.limits,"rm2")
rm3.family.plot = linearLinePlot(rm3_16S.family,effective.breaks,effective.labels,effective.limits,"rm3")
rm4.family.plot = linearLinePlot(rm4_16S.family,effective.breaks,effective.labels,effective.limits,"rm4")

effectivefamily.lineplot = ggarrange(rm1.family.plot,rm2.family.plot,rm3.family.plot,rm4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.lineplot,filename = "reports/figures/16S/rm_line.png",fonts=c("serif"),width = 20,height = 6,dpi=300)

#######################    ^^^^^^^^   #####################################

rm1.family.plot = logLinePlot(rm1_16S.family,effective.breaks,effective.labels,effective.limits,"rm1")
rm2.family.plot = logLinePlot(rm2_16S.family,effective.breaks,effective.labels,effective.limits,"rm2")
rm3.family.plot = logLinePlot(rm3_16S.family,effective.breaks,effective.labels,effective.limits,"rm3")
rm4.family.plot = logLinePlot(rm4_16S.family,effective.breaks,effective.labels,effective.limits,"rm4")

effectivefamily.logplot = ggarrange(rm1.family.plot,rm2.family.plot,rm3.family.plot,rm4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.logplot,filename = "reports/figures/16S/rm_logline.png",fonts=c("serif"),width = 20,height = 6,dpi=300)

#######################    ^^^^^^^^   #####################################

rm1.family.plot = areaPlot(rm1_16S.family_wo13,effective.breaks,effective.labels,effective.limits,"rm1")
rm2.family.plot = areaPlot(rm2_16S.family,effective.breaks,effective.labels,effective.limits,"rm2")
rm3.family.plot = areaPlot(rm3_16S.family,effective.breaks,effective.labels,effective.limits,"rm3")
rm4.family.plot = areaPlot(rm4_16S.family,effective.breaks,effective.labels,effective.limits,"rm4")

effectivefamily.areaplot = ggarrange(rm1.family.plot,rm2.family.plot,rm3.family.plot,rm4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.areaplot,filename = "reports/figures/16S/rm_area.png",fonts=c("serif"),width = 20,height = 6,dpi=300)
#######################    ^^^^^^^^   #####################################

type="im"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]




im1.family.plot = linearLinePlot(im1_16S.family,effective.breaks,effective.labels,effective.limits,"im1")
im2.family.plot = linearLinePlot(im2_16S.family,effective.breaks,effective.labels,effective.limits,"im2")
im3.family.plot = linearLinePlot(im3_16S.family,effective.breaks,effective.labels,effective.limits,"im3")
im4.family.plot = linearLinePlot(im4_16S.family,effective.breaks,effective.labels,effective.limits,"im4")

effectivefamily.lineplot = ggarrange(im1.family.plot,im2.family.plot,im3.family.plot,im4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.lineplot,filename = "reports/figures/16S/im_line.png",fonts=c("serif"),width = 20,height = 6,dpi=300)

#######################    ^^^^^^^^   #####################################

im1.family.plot = logLinePlot(im1_16S.family,effective.breaks,effective.labels,effective.limits,"im1")
im2.family.plot = logLinePlot(im2_16S.family,effective.breaks,effective.labels,effective.limits,"im2")
im3.family.plot = logLinePlot(im3_16S.family,effective.breaks,effective.labels,effective.limits,"im3")
im4.family.plot = logLinePlot(im4_16S.family,effective.breaks,effective.labels,effective.limits,"im4")

effectivefamily.logplot = ggarrange(im1.family.plot,im2.family.plot,im3.family.plot,im4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.logplot,filename = "reports/figures/16S/im_logline.png",fonts=c("serif"),width = 20,height = 6,dpi=300)

#######################    ^^^^^^^^   #####################################

im1.family.plot = areaPlot(im1_16S.family,effective.breaks,effective.labels,effective.limits,"im1")
im2.family.plot = areaPlot(im2_16S.family,effective.breaks,effective.labels,effective.limits,"im2")
im3.family.plot = areaPlot(im3_16S.family,effective.breaks,effective.labels,effective.limits,"im3")
im4.family.plot = areaPlot(im4_16S.family,effective.breaks,effective.labels,effective.limits,"im4")

effectivefamily.areaplot = ggarrange(im1.family.plot,im2.family.plot,im3.family.plot,im4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.areaplot,filename = "reports/figures/16S/im_area.png",fonts=c("serif"),width = 20,height = 6,dpi=300)


type="nc"

effective.breaks = breaks.hash[[type]]
effective.labels = labels.hash[[type]]
effective.limits = limits.hash[[type]]




nc1.family.plot = linearLinePlot(nc1_16S.family,effective.breaks,effective.labels,effective.limits,"nc1")
nc2.family.plot = linearLinePlot(nc2_16S.family,effective.breaks,effective.labels,effective.limits,"nc2")
nc3.family.plot = linearLinePlot(nc3_16S.family,effective.breaks,effective.labels,effective.limits,"nc3")
nc4.family.plot = linearLinePlot(nc4_16S.family,effective.breaks,effective.labels,effective.limits,"nc4")

effectivefamily.lineplot = ggarrange(nc1.family.plot,nc2.family.plot,nc3.family.plot,nc4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.lineplot,filename = "reports/figures/16S/nc_line.png",fonts=c("serif"),width = 20,height = 6,dpi=300)

#######################    ^^^^^^^^   #####################################

nc1.family.plot = logLinePlot(nc1_16S.family,effective.breaks,effective.labels,effective.limits,"nc1")
nc2.family.plot = logLinePlot(nc2_16S.family,effective.breaks,effective.labels,effective.limits,"nc2")
nc3.family.plot = logLinePlot(nc3_16S.family,effective.breaks,effective.labels,effective.limits,"nc3")
nc4.family.plot = logLinePlot(nc4_16S.family,effective.breaks,effective.labels,effective.limits,"nc4")

effectivefamily.logplot = ggarrange(nc1.family.plot,nc2.family.plot,nc3.family.plot,nc4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.logplot,filename = "reports/figures/16S/nc_logline.png",fonts=c("serif"),width = 20,height = 6,dpi=300)

#######################    ^^^^^^^^   #####################################

nc1.family.plot = areaPlot(nc1_16S.family,effective.breaks,effective.labels,effective.limits,"nc1")
nc2.family.plot = areaPlot(nc2_16S.family,effective.breaks,effective.labels,effective.limits,"nc2")
nc3.family.plot = areaPlot(nc3_16S.family,effective.breaks,effective.labels,effective.limits,"nc3")
nc4.family.plot = areaPlot(nc4_16S.family,effective.breaks,effective.labels,effective.limits,"nc4")

effectivefamily.areaplot = ggarrange(nc1.family.plot,nc2.family.plot,nc3.family.plot,nc4.family.plot, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(effectivefamily.areaplot,filename = "reports/figures/16S/nc_area.png",fonts=c("serif"),width = 20,height = 6,dpi=300)

