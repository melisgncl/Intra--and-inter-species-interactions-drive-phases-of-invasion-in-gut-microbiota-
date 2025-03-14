## CODE BLOCK
## Louis Gauthier & Melis Gencel

#######################    READ ME!   #####################################

#####   Plot barcode dynamics (Muller-style area plot and log-line plot)

#######################    ^^^^^^^^   #####################################

outdir=paste0("reports/figures/")

if (!dir.exists(paste0(outdir, "/dynamics"))) {
  dir.create(paste0(outdir, "/dynamics"), recursive = TRUE)
}




#here we hardcode the selection such that all barcodes with max frequency > 0.0005 are assigned a hex color
all.top.max=read_csv("src/visualization/0_config/all_top_max2.csv")

###topcolors#
top_colors2=read_csv("src/visualization/0_config/top_colors3.csv")
## we create a very long list of colors for low-frequency barcodes. 
## repetition of colors isn't an issue, we just want each barcodes to be colored instead of assigning gray values
long.color.list = rep(top_colors2$hex,50)
long.color.list.random = sample(long.color.list)


plotDynamics <- function(cluster_df, raw_df, sample,type) {
 
  ####get the breaks and labels from has table##
  effective.breaks = breaks.hash[[type]]
  effective.labels = labels.hash[[type]]
  effective.limits = limits.hash[[type]]
  
  
  cluster.top = merge(cluster_df,all.top.max,by = "Center")

  # Remove unwanted columns if they exist
  columns_to_remove <- c("time_point_1", "max", "Cluster.Score")
  
  for (col in columns_to_remove) {
    if (col %in% names(cluster.top)) {
      cluster.top[[col]] <- NULL
    }
  }
    
    df = raw_df
    df = merge(df, cluster.top, by.x = "ID",by.y = "Cluster.ID",all.x = TRUE)
    
    # for all barcodes without a universal color, assign a gray hex value
    df[is.na(df$hex),]$hex="#cccccc"
    
    # remove unwanted columns
    # Check and remove 'Cluster.Score' column if it exists
    if ("Cluster.Score" %in% names(df)) {
      df$Cluster.Score <- NULL
    }
    
    # Check and remove 'time_point_1' column if it exists
    if ("time_point_1" %in% names(df)) {
      df$time_point_1 <- NULL
    }
    
    # convert ID to factor for grouping
    df$ID = as.factor(df$ID)

    # important: order dataframe by max frequency; this will determine the stacking of areas in the plot 
    tf = df[order(df$max),]

    # subset dataframe with barcodes above a certain max frequency
    grouped_tf = tf %>% group_by(hex) %>% filter(max>0.001) %>% ungroup()

    # set color scale using universal hex colors
    mycolors = grouped_tf$hex
    names(mycolors) = grouped_tf$ID

    # map factor levels to max frequencies
    grouped_tf$ID = factor(grouped_tf$ID, levels = unique(grouped_tf$ID[order(grouped_tf$max)]))
    tf$ID = factor(tf$ID, levels = unique(tf$ID[order(tf$max)]))
    #write_csv(tf,file= paste("reports_article/figures/dynamics/", sample, ".csv", sep=""))
    # plot area
    # here we plot a first geom_area with all the barcodes using the long color list
    # then we add a new fill scale, and plot a second geom area on top of the first using only the subsetted, high-freq barcodes
    g = ggplot(tf) + geom_area(aes(x=variable,y=value,group=ID,fill=ID),data=tf) + scale_fill_manual(values = long.color.list.random) + 
      new_scale_fill() + geom_area(aes(x=variable,y=value,group=ID,fill=ID),data=grouped_tf) +
      scale_x_continuous(breaks = effective.breaks,labels = effective.labels) + theme_Publication() + guides(color = FALSE)
      
    #labs(x = "Time (post-gavage)",y="Barcode density") + coord_cartesian(expand = FALSE)

    ## Save as jpeg due to the high information density of this plot
    ggsave(g,filename = paste(outdir,"/dynamics/", sample, "_area.jpeg", sep=""), type = "cairo",width = 8.25,height = 6)
    

    ##all.line = ggplot() + geom_line(aes(x=variable,y=value,group=ID),data=tf,colour="#CCCCCC",alpha=0.3) +
    ## geom_line(aes(x=variable,y=value,group=ID,colour=ID),data=grouped_tf,size=1) +
    ##theme_Publication() + scale_y_log10(limits=c(1e-7,1e0)) + scale_color_manual(values = mycolors,name="Cluster ID") + labs(x = "Time (post-gavage)",y="Barcode frequency") +
    ##scale_x_continuous(breaks = effective.breaks,labels = effective.labels) +
    ## guides(color = FALSE, shape = guide_legend(order = 1))  + coord_cartesian(expand = FALSE)
    ## Save as jpeg due to the high information density of this plot
    ##ggsave(all.line,filename = paste(outdir,"/dynamics/", sample, "_line.jpeg", sep=""), type = "cairo",width = 8.25,height = 6)
    
}


#######################    ^^^^^^^^   #####################################



plotDynamics(gf1_cluster,gf1.df,"gf1","gf")
plotDynamics(gf2_cluster,gf2.df,"gf2","gf")
plotDynamics(gf3_cluster,gf3.df,"gf3","gf")
plotDynamics(gf4_cluster,gf4.df,"gf4","gf")

#######################    ^^^^^^^^   #####################################

plotDynamics(rm1_cluster,rm1.df,"rm1","rm")
plotDynamics(rm2_cluster,rm2.df,"rm2","rm")
plotDynamics(rm3_cluster,rm3.df,"rm3","rm")
plotDynamics(rm4_cluster,rm4.df,"rm4","rm")

#######################    ^^^^^^^^   #####################################

plotDynamics(im1_cluster,im1.df,"im1","im")
plotDynamics(im2_cluster,im2.df,"im2","im")
plotDynamics(im3_cluster,im3.df,"im3","im")
plotDynamics(im4_cluster,im4.df,"im4","im")

