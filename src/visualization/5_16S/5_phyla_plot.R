im1.phyla.plot = ggplot(im1_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
im2.phyla.plot = ggplot(im2_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
im3.phyla.plot = ggplot(im3_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
im4.phyla.plot = ggplot(im4_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)

im.phyla.plot = ggarrange(im1.phyla.plot,im2.phyla.plot,im3.phyla.plot,im4.phyla.plot, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
ggsave(im.phyla.plot,filename = "reports/figures/16S/im1_4_phyla.eps",fonts=c("serif"),width = 18,height = 16,device = cairo_ps())


rm1.phyla.plot = ggplot(rm1_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
rm2.phyla.plot = ggplot(rm2_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
rm3.phyla.plot = ggplot(rm3_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
rm4.phyla.plot = ggplot(rm4_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
rm.phyla.plot = ggarrange(rm1.phyla.plot,rm2.phyla.plot,rm3.phyla.plot,rm4.phyla.plot, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
ggsave(rm.phyla.plot,filename = "reports/figures/16S/rm1_4_phyla.eps",fonts=c("serif"),width = 18,height = 16,device = cairo_ps())


nc1.phyla.plot = ggplot(nc1_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
nc2.phyla.plot = ggplot(nc2_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
nc3.phyla.plot = ggplot(nc3_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
nc4.phyla.plot = ggplot(nc4_16S.phyla) + geom_area(aes(x=Time,y=Abundance.phyla,group=Phylum,fill=Phylum)) +
  labs(x = "Time (days)",y="Phyla composition") + coord_cartesian(expand = FALSE) + theme_Publication()  + scale_fill_manual(values=phylum.colors)
nc.phyla.plot = ggarrange(nc1.phyla.plot,nc2.phyla.plot,nc3.phyla.plot,nc4.phyla.plot, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
ggsave(nc.phyla.plot,filename = "reports/figures/16S/nc1_4_phyla.eps",fonts=c("serif"),width = 18,height = 16,device = cairo_ps())