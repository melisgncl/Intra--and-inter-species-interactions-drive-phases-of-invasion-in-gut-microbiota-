## CODE BLOCK
## Louis Gauthier & Melis Gencel

#######################    READ ME!   #####################################

##### Compute and plot within- and between-cohort intersection

#######################    ^^^^^^^^   #####################################

# plot within cohort intersection

outdir=paste0("reports/figures/")

if (!dir.exists(paste0(outdir, "/intersection"))) {
  dir.create(paste0(outdir, "/intersection"), recursive = TRUE)
}







gf1_4.Input = list(gf1 = c(gf1.top.final$Center), gf2 = c(gf2.top.final$Center), gf3 = c(gf3.top.final$Center), gf4 = c(gf4.top.final$Center))
gf1_4.upset = upset(fromList(gf1_4.Input),order.by = "freq", point.size = 3, line.size = 1.1,set_size.show = FALSE,mainbar.y.label = "Barcodes unique to intersection",matrix.dot.alpha = 0.25,empty.intersections = TRUE,text.scale = 1.3,main.bar.color = c3,sets.bar.color = c3,mainbar.y.max = n_intersect)


rm1_4.Input = list(rm1 = c(rm1.top.final$Center), rm2 = c(rm2.top.final$Center), rm3 = c(rm3.top.final$Center), rm4 = c(rm4.top.final$Center))
rm1_4.upset = upset(fromList(rm1_4.Input),order.by = "freq", point.size = 3, line.size = 1.1,set_size.show = FALSE,mainbar.y.label = "Barcodes unique to intersection",matrix.dot.alpha = 0.25,empty.intersections = TRUE,text.scale = 1.3,main.bar.color = c2,sets.bar.color = c2,mainbar.y.max = n_intersect)

im1_4.Input = list(im1 = c(im1.top.final$Center), im2 = c(im2.top.final$Center), im3 = c(im3.top.final$Center), im4 = c(im4.top.final$Center))
im1_4.upset = upset(fromList(im1_4.Input),order.by = "freq", point.size = 3, line.size = 1.1,set_size.show = FALSE,mainbar.y.label = "Barcodes unique to intersection",matrix.dot.alpha = 0.25,empty.intersections = TRUE,text.scale = 1.3,main.bar.color = c1,sets.bar.color = c1,mainbar.y.max = n_intersect)




## plot & save

im1_4.upset
grid.edit('arrange',name='arrange2')
vp1 = grid.grab()
gf1_4.upset
grid.edit('arrange',name='arrange2')
vp3 = grid.grab()
rm1_4.upset
grid.edit('arrange',name='arrange2')
vp2 = grid.grab()
cairo_ps(paste0(outdir,"/intersection/within_cohort_sets_gf_rm_im.eps"),width = 18,height = 10)
plot_grid(vp1,vp2,vp3 ,ncol = 3, labels = c('A', 'B','C'))
dev.off()



###############
# calculate between cohort intersection

gf1_4.top.final = rbind(gf1.top.final,gf2.top.final,gf3.top.final,gf4.top.final)
gf1_4.top.final = unique(gf1_4.top.final$Center)

rm1_4.top.final = rbind(rm1.top.final,rm2.top.final,rm3.top.final,rm4.top.final)
rm1_4.top.final = unique(rm1_4.top.final$Center)

im1_4.top.final = rbind(im1.top.final,im2.top.final,im3.top.final,im4.top.final)
im1_4.top.final = unique(im1_4.top.final$Center)



all_cohorts.Input = list(im1_4 = im1_4.top.final, rm1_4 = rm1_4.top.final,gf1_4 = gf1_4.top.final)
all_cohorts.upset = upset(fromList(all_cohorts.Input),order.by = "freq",nsets = 6 ,point.size = 3, line.size = 1,set_size.show = FALSE,mainbar.y.label = "Barcodes unique to intersection",matrix.dot.alpha = 0.25,empty.intersections = TRUE,text.scale = 1)

length(unique(c(im1_4.top.final,gf1_4.top.final,rm1_4.top.final)))

## plot & save


cairo_ps(paste0(outdir,"/intersection/between_cohorts_sets_gf_wm_im.eps"),width = 18,height = 10)
all_cohorts.upset
dev.off()

