
create_plots <- function() {
  for(df_name in names(dpmsr_set$data$final)){
    plot_dir <- create_dir(str_c(dpmsr_set$file$output_dir, df_name))
    df <- dpmsr_set$data$final[[df_name]]
    volcano_plot(df, df_name, plot_dir)
    df <- df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(df, df_name, plot_dir)
    box_plot(df, df_name, plot_dir)
    densities_plot(df, df_name, plot_dir)
    MDS_plot(df, df_name, plot_dir)   #str_c(y, "Multidimension Scaling"))
    try(PCA_plot(df, df_name, plot_dir), silent = TRUE)
    try(Cluster_plot(df, df_name, plot_dir), silent = TRUE)
    try(heatmap_plot(df, df_name, plot_dir), silent = TRUE)
    }
}

create_select_plots <- function() {
  for(df_name in names(dpmsr_set$data$final)){
    plot_dir <- str_c(dpmsr_set$file$output_dir, "//", df_name, "//")
    df <- dpmsr_set$data$final[[df_name]]
    protein_qc_plots(df, df_name, plot_dir)
  }
}

#-------------------------------------------------------------------------------
volcano_plot <- function(df, df_name, plot_dir)  #x, comp_name, title, plot_dir)
{
  df[is.na(df)] <- 0
  for(i in 1:as.numeric(dpmsr_set$x$comp_number)) {
      plottitle <- str_c(df_name,"_", dpmsr_set$y$comp_groups$comp_name[i])
      y <- unlist(df[ , dpmsr_set$y$comp_groups$fc2[i]])
      z <- unlist(df[ , dpmsr_set$y$comp_groups$pval[i]])
      x <- data.frame(y)
      x <- cbind(x,z)
      file_name <- str_c(plot_dir, df_name, "_", dpmsr_set$y$comp_groups$comp_name[i], "_volcano.png")
      ggplot(x, aes(log2(y), -log10(z))  ) +
        geom_point(alpha=0.4, size=2, aes(color = z)) +
        xlab("log2 fold change") + ylab("-log10 p-value") +
        scale_colour_gradient(low = "blue", high = "black") +
        ggtitle(plottitle)+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "none")  
      #legend.text=element_text(size=10, vjust=-1), 
      #legend.title = element_text(size=10, vjust = 0.2),
      #legend.title = element_blank(),
      #legend.position=c(.9,0.8))
      ggsave(file_name, width=5, height=4)
  }
}


#Bar plot-------------------------------------------------
bar_plot <- function(df,plot_title,plot_dir) {
  namex <- dpmsr_set$design$Label
  datay <- colSums(df, na.rm = TRUE)
  df2 <- data.frame(namex)
  df2$Total_Intensity <- datay
  colnames(df2) <- c("Sample", "Total_Intensity")
  df2$Sample <- factor(df2$Sample, levels = df2$Sample)
  file_name <- str_c(plot_dir, plot_title, "_barplot.png")
  ymax <- max(datay)
  ggplot(data=df2, aes(x=Sample, y=Total_Intensity)) +
    geom_bar(stat="identity", fill=dpmsr_set$design$colorlist)+ theme_classic() + 
    ggtitle(plot_title) + 
    xlab(NULL)+
    ylab(NULL)+
    #scale_y_discrete(labels = NULL) +
    coord_cartesian(ylim=NULL, expand = TRUE) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(size=5, angle = 90,  color="black"),
          axis.text.y = element_text(size=5,  color="black"),
          ) 
  ggsave(file_name, width=5, height=4)
  return("done")
}


#--plot densities------------------------------------------------------------------
densities_plot <- function(x,y,plot_dir) {
  png(filename=str_c(plot_dir, y, "_density.png"), width = 888, height = 571)  
  plotDensities(log2(x), 
                group = dpmsr_set$y$sample_groups$Group,   
                col = dpmsr_set$y$sample_groups$colorlist, 
                main = y)
  dev.off()
}


#Box plot-------------------------------------------------
box_plot <- function(x,y, plot_dir) {
  png(filename=str_c(plot_dir, y, "_boxplot.png"), width = 888, height = 571)
  data_box <- log2(x)
  data_box[data_box ==-Inf ] <- NA
  boxplot(data_box, 
          col = dpmsr_set$design$colorlist, 
          notch = TRUE, 
          boxwex = 0.8,
          main = c(y),
          axes=TRUE,
          horizontal = TRUE,
          las=1,
          par(mar=c(8,10,4,2)))
  dev.off()
}



#MDS Plot-------------------------------------------------
MDS_plot <- function(x,y,plot_dir) {
  png(filename=str_c(plot_dir, y, "_MDS.png"), width = 888, height = 571)  
  plotMDS(log2(x), 
          col = dpmsr_set$design$colorlist, 
          main = y)
  dev.off()
}


#Heat map-------------------------------------------------
heatmap_plot <- function(y,plottitle,plot_dir)
{
  #y <- y[(info_columns_final+1):ncol(y)] # strip off info columns
  y <- log2(y)
  y <- data.matrix(y)
  ## Row- and column-wise clustering 
  hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
  hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
  ## Tree cutting
  mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
  ## Plot heatmap 
  mycol <- redgreen(75) #colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
  #png(filename=str_c(plot_dir, plottitle, "_heatmap.png"), width = 888, height = 571, units="px")  
  #png(filename=str_c(plot_dir, plottitle, "_heatmap2.png"), width = 8880, height = 5710)  
  jpeg(filename=str_c(plot_dir, plottitle, "_heatmap.jpeg"), units="px", quality=100, width = 1776, height = 1146)  
  #tiff(filename=str_c(plot_dir, plottitle, "_heatmap2.tiff"), units="px", compression = "none", pointsize=20, width = 1776, height = 1146)  
  #tiff(filename=str_c(plot_dir, plottitle, "_heatmap1.tiff"), units="px", compression = "none", width = 888, height = 571)  
  heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, labCol= dpmsr_set$design$Group, 
            scale="row", density.info="none", trace="none", RowSideColors=mycolhc, main = plottitle) 
  dev.off()
}




#-------------------------------------------------------------------------------

#PCA 2D 3D-------------------------------------------------
PCA_plot <- function(x,y, plot_dir) {
  #x <- x[(info_columns_final+1):ncol(x)] # strip off info columns
  require(pca3d)
  require(rgl)
  x_transpose <- t(x)
  x_transpose <-data.frame(x_transpose)
  row.names(x_transpose) <- NULL
  x_transpose <-cbind(dpmsr_set$design$Group, x_transpose)
  x_pca <- prcomp(x_transpose[,-1], scale=TRUE)
  test_this <-x_transpose[,1]
  x_gr <- factor(unlist(test_this))
  summary(x_gr)
  pca3d(x_pca, 
        group=x_gr,
        legend = "right",
        palette = dpmsr_set$y$sample_groups$colorlist, 
        radius = 2,
        title = y)
  snapshotPCA3d(file=str_c(plot_dir, y, "_PCA3d", ".png"))
  
  #2D PCA
  df_out <- as.data.frame(x_pca$x)
  file_name <- str_c(plot_dir, y, "_PCA2D.png")
  ggplot(df_out,aes(x=PC1,y=PC2,color=x_gr )) +
    geom_point(size =3) +
    theme(legend.title=element_blank()) +
    ggtitle(y) +
    scale_color_manual(values = dpmsr_set$y$sample_groups$colorlist) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file_name, width=5, height=4)
  
  #Cluster
  file_name <- str_c(plot_dir, y, "_Cluster.png")
  df<- x_transpose
  df <- na.omit(df)
  #the  distance  measure  to  be  used.   This  must  be  one  of  "euclidean",  "maxi-mum", 
  #"manhattan", "canberra", "binary", "minkowski", "pearson", "spearman"or "kendall".
  distance <- get_dist(df, method="euclidean", stand = TRUE)
  fviz_dist(distance,  show_labels = TRUE, gradient = list(low = "#009933", mid = "white", high = "#FF3366")) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(size=5, angle = 90,  color="black"),
          axis.text.y = element_text(size=5,  color="black"))+
    ggtitle(y) 
  ggsave(file_name, width=5, height=4)
  
  return("done")
}


#Cluster-------------------------------------------------
Cluster_plot <- function(x,y, plot_dir) {
  #x <- x[(info_columns_final+1):ncol(x)] # strip off info columns
  df <- t(x)
  df <-data.frame(df)
  #row.names(df) <- NULL
  #df <- na.omit(df)
  df[] <- lapply(df, as.numeric)
  df <- scale(df)
  #df <-cbind(dpmsr_set$design$Header2, df)
  
  #Cluster
  file_name <- str_c(plot_dir, y, "_Cluster.png")
  #the  distance  measure  to  be  used.   This  must  be  one  of  "euclidean",  "maxi-mum", 
  #"manhattan", "canberra", "binary", "minkowski", "pearson", "spearman"or "kendall".
  distance <- get_dist(df, method="euclidean")
  fviz_dist(distance,  show_labels = TRUE, gradient = list(low = "#009933", mid = "white", high = "#FF3366")) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(size=5, angle = 90,  color="black"),
          axis.text.y = element_text(size=5,  color="black"))+
    ggtitle(y) 
  ggsave(file_name, width=5, height=4)
  
  return("done")
}

#Histogram for total intensity-------------------------------------------------
histogram_plot <- function()
{
  x<-dpmsr_set$data$data_to_norm[(dpmsr_set$y$info_columns+1):ncol(dpmsr_set$data$data_to_norm)]
  title<-as.character(dpmsr_set$x$file_prefix)
  plottitle<-"Total Intensity Histogram"
  testthis1 <- as.matrix(log2(x))
  intensity_cutoff <- dpmsr_set$x$int_cutoff
  
  x_mean <- mean(testthis1, na.rm=TRUE)
  x_stdev <- sd(testthis1, na.rm=TRUE)
  new_cutoff <- x_mean - (0.5*x_stdev)
  
  if (new_cutoff < intensity_cutoff){intensity_cutoff<-new_cutoff}
  if (intensity_cutoff < log2(100000)){intensity_cutoff<-log2(100000)}
  
  file_name <- str_c(dpmsr_set$file$output_dir, title, "_histogram.png")
  testgather <- gather(x)
  testgather <- subset(testgather, testgather$value>0)
  testgather$value <- log2(testgather$value)
  dpmsr_set$x$int_cutoff <<- trunc(2^intensity_cutoff)
  
  ggplot(testgather, aes(value))+
    geom_histogram(bins=100, fill="blue") +
    geom_vline(aes(xintercept = intensity_cutoff), color='red'   ) +
    geom_vline(aes(xintercept = log2(5000000)), color='green') +
    geom_vline(aes(xintercept = x_mean)) +
    geom_vline(aes(xintercept = (x_mean + x_stdev))) +
    geom_vline(aes(xintercept = (x_mean - x_stdev))) +
    ggtitle(plottitle)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") 
  ggsave(str_c(dpmsr_set$file$qc_dir,"Intensity_Histogram.png"), width=5, height=4)
  #ggsave(file_name, width=5, height=4)
}





