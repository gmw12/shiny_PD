

#-----------------------------------------------------------------------------------------
qc_apply <- function(){
  #set up tables for collection of QC Spike Data, and %CV's, and total CV's
  dpmsr_set$data$qc_spike <<- data.frame(dpmsr_set$design$"QC Spike Level")
  colnames(dpmsr_set$data$qc_spike)[1] <<- "SpikeLevel"
  dpmsr_set$data$summary_cv <<- data.frame(dpmsr_set$y$sample_groups$Group)
  colnames(dpmsr_set$data$summary_cv)[1] <<- "Group"
  dpmsr_set$data$total_cv <<- data.frame(dpmsr_set$data$final$impute$Accession )
  colnames(dpmsr_set$data$total_cv)[1] <<- "Accession"

  for(df_name in names(dpmsr_set$data$final)){
    cv_stats(dpmsr_set$data$final[[df_name]], df_name)
    qc_spike(dpmsr_set$data$final[[df_name]], df_name)
  }
  
  colnames(dpmsr_set$data$summary_cv) <<- c("Group",names(dpmsr_set$data$final))
  rownames(dpmsr_set$data$summary_cv) <<- NULL
  dpmsr_set$data$qc_spike_final <<- qc_spike_final(dpmsr_set$data$qc_spike)
  cv_grouped_plot()
  qc_spike_plot()
  if (as.logical(dpmsr_set$x$adh_spike)) {try(adh_spike(), silent=TRUE)}
}


#----------------------------------------------------------------------------------------- 
#collect CV statistics
cv_stats <- function(data_out, title){
  cv_start <- dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number+1
  avg_cv <- colMeans(data_out[cv_start:(cv_start+dpmsr_set$y$group_number-1)], na.rm = TRUE)
  dpmsr_set$data$summary_cv <<- cbind(dpmsr_set$data$summary_cv, avg_cv)
  all_cv <- data_out[cv_start:(cv_start+dpmsr_set$y$group_number-1)]
  column_names <- colnames(all_cv)
  column_names <- str_c(title, " ",column_names)
  colnames(all_cv) <- column_names
  dpmsr_set$data$total_cv <<- cbind(dpmsr_set$data$total_cv, all_cv)
}


#qc spike metrics ---------------------------------
qc_spike <- function(data_in, data_title) {
  qc_spike_id <- unlist(strsplit(as.character(dpmsr_set$x$qc_spike_id), split=","))
  spike_protein <-subset(data_in, Accession %in% qc_spike_id)  
  spike_protein <- spike_protein[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
  spike_total <- colSums(spike_protein)
  
  dpmsr_set$data$qc_spike <<- cbind(dpmsr_set$data$qc_spike, spike_total)
  colnames(dpmsr_set$data$qc_spike)[colnames(dpmsr_set$data$qc_spike) == 'spike_total'] <<- data_title
  dpmsr_set$data$qc_spike <<- round(dpmsr_set$data$qc_spike,0)
 }


#qc spike metrics finalize---------------------------------
qc_spike_final <- function(data_in) {
  data_in <- aggregate(data_in, by=list(Category=data_in$SpikeLevel), FUN=mean)
  data_in$Category <- NULL
  row_count <- nrow(data_in)
  
  for(i in 2:row_count) {
    testme <- round((data_in[i,]/data_in[1,]), 2)
    if (i==2){
      qc_stat <- testme
    }else{
      qc_stat <- rbind(qc_stat, testme)
    }
  }
  data_in <-trunc(round(data_in,0))
  data_in <- data_in %>% mutate_all(as.character)
  qc_stat <- qc_stat %>% mutate_all(as.character)
  data_out <- rbind(data_in, qc_stat)
  
  Simple_Excel(data_out, str_c(dpmsr_set$file$qc_dir, "QC_Spike.xlsx"))
  return(data_out)
}

#-----------------------------------------------------------------------------------------
#save dataframe with ADH peptides for QC meteric
adh_spike <- function(){
  ADH_data <-subset(dpmsr_set$data$impute$impute, Accession %in% dpmsr_set$x$adh_list)
  ADH_data$missings <- rowSums(ADH_data[(dpmsr_set$y$info_columns+1):ncol(ADH_data)] ==0)
  ADH_data <- subset(ADH_data, missings==0)
  ADH_data <- ADH_data[ , -ncol(ADH_data)]
  ADH_annotate <- ADH_data[, c("Accession", "Sequence")]
  ADH_data <- ADH_data[(dpmsr_set$y$info_columns+1):ncol(ADH_data)]
  ADH_data <- cbind(ADH_annotate, ADH_data)
  ADH_data <- ADH_data  %>% group_by(Accession, Sequence) %>% summarise_all(funs(sum))
  ADH_data$cv <- percentCV_gw(ADH_data[(dpmsr_set$y$info_columns+1):ncol(ADH_data)])
  ADH_data$sd <- apply(ADH_data[(dpmsr_set$y$info_columns+1):(ncol(ADH_data)-1)], 1, FUN = function(x) {sd(x)})
  ADH_data$av <- apply(ADH_data[(dpmsr_set$y$info_columns+1):(ncol(ADH_data)-1)], 1, FUN = function(x) {mean(x)})
  ADH_data$id <- seq(1, nrow(ADH_data), by=1)
  ADH_data$color <- distinctColorPalette(nrow(ADH_data))
  ADH_data <- ADH_data[order(-ADH_data$av),]
  ADH_data <- ADH_data[1:15,]
  file_name <- str_c(dpmsr_set$file$qc_dir, "ADH_barplot.png")
  barplot_adh(ADH_data$cv[1:12], ADH_data$Sequence[1:12],"ADH Peptide CV's", file_name)
}

#-----------------------------------------------------------------------------------------
cv_grouped_plot <- function() {
  file_name <- str_c(dpmsr_set$file$qc_dir, "CV_barplot.png")
  
  data_in <- dpmsr_set$data$summary_cv
  data_in <- setNames(data.frame(t(data_in[,-1])), data_in[,1])
  data_in$Norm <- row.names(data_in)
  row.names(data_in) <- NULL
  data_plot <- data_in %>% pivot_longer(-Norm, names_to = "Sample", values_to = "CV")
  # Grouped
  ggplot(data_plot, aes(fill=Sample, y=CV, x=Norm)) + 
    geom_bar(position="dodge", stat="identity")
  ggsave(file_name, width=8, height=4)
}

#-----------------------------------------------------------------------------------------
qc_spike_plot <- function() {
  data_in <- dpmsr_set$data$qc_spike
  data_in$Sample <- dpmsr_set$design$Label
  #data_in <- data_in[order(data_in$SpikeLevel),]
  
  for(df_name in names(dpmsr_set$data$final)){
      file_name <- str_c(dpmsr_set$file$qc_dir,df_name, "_QC_Spike_barplot.png")
      plot_title <- df_name
      name1 <- data_in$Sample
      name2 <- data_in$SpikeLevel
      datay <- data_in[[df_name]]
      df2 <- data.frame(name1)
      df2$SpikeLevel <- name2
      df2$Intensity <- datay
      colnames(df2) <- c("Sample", "Spike_Level", "Intensity")
      df2$Spike_Level <- as.character(df2$Spike_Level)
    
      ggplot(df2, aes(fill=Sample, y=Intensity, x=Spike_Level)) + 
        geom_bar(position="dodge", stat="identity")+
        theme(legend.position = "none")
      ggsave(file_name, width=8, height=4)
  }
}

#-----------------------------------------------------------------------------------------
barplot_adh <- function(data_cv, data_names, plot_title, file_name) {
  df2 <- data.frame(data_names)
  df2$CV <- data_cv
  data_av <- mean(data_cv)
  colnames(df2) <- c("Peptide", "CV")
  df2$Peptide <- factor(df2$Peptide, levels = df2$Peptide)
  
  ggplot(data=df2, aes(x=Peptide, y=CV)) +
    geom_bar(stat="identity", fill=rainbow(n=12))+ theme_classic() + 
    ggtitle(plot_title) + 
    xlab(NULL)+
    ylab(NULL)+
    #scale_y_discrete(labels = NULL) +
    coord_cartesian(ylim=NULL, expand = TRUE) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(size=5, angle = 90,  color="black"),
          axis.text.y = element_text(size=5,  color="black"),
    )+
    geom_hline(yintercept = data_av, color="blue")+
    geom_text(aes(2, data_av+5, label = str_c("Average = ", trunc(round(data_av,2)), "%"), vjust = -1))
  ggsave(file_name, width=8, height=4)
  return("done")
}



#-----------------------------------------------------------------------------------------

protein_qc_plots<- function(data_in, plot_title, plot_dir) {
 
  dpmsr_set$y$protein_list <<- list()
  
  if(!is.null(dpmsr_set$x$adh_list)){
    adh_plot <-subset(data_in, Accession %in% dpmsr_set$x$adh_list)  
    adh_plot <- adh_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(adh_plot, str_c("ADH_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "ADH")
  }
  
  if(!is.null(dpmsr_set$x$bira_list)){
    bira_plot <-subset(data_in, Accession %in% dpmsr_set$x$bira_list)  
    bira_plot <- bira_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(bira_plot, str_c("BirA_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "BirA")
  }
  
  if(!is.null(dpmsr_set$x$bait_list)){
    bait_plot <-subset(data_in, Accession %in% dpmsr_set$x$bait_list)  
    bait_plot <- bait_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(bait_plot, str_c("Bait_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "Bait")
  }
  
  if(!is.null(dpmsr_set$x$avidin_list)){
    avidin_plot <-subset(data_in, Accession %in% dpmsr_set$x$avidin_list)  
    avidin_plot <- avidin_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(avidin_plot, str_c("Avidin_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "Avidin")
  }
  
  if(!is.null(dpmsr_set$x$carbox_list)){
    carbox_plot <-subset(data_in, Accession %in% dpmsr_set$x$carbox_list)  
    carbox_plot <- carbox_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(carbox_plot, str_c("Carbox_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "Carbox")
  }
  
  if(!is.null(dpmsr_set$x$protein1_list)){
    protein1_plot <-subset(data_in, Accession %in% dpmsr_set$x$protein1_list)  
    protein1_plot <- protein1_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(protein1_plot, str_c("Protein1_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "Protein1")
  }
  
  if(!is.null(dpmsr_set$x$protein2_list)){
    protein2_plot <-subset(data_in, Accession %in% dpmsr_set$x$protein2_list)  
    protein2_plot <- protein2_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(protein2_plot, str_c("Protein2_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "Protein2")
  }
  
  if(!is.null(dpmsr_set$x$protein3_list)){
    protein3_plot <-subset(data_in, Accession %in% dpmsr_set$x$protein3_list)  
    protein3_plot <- protein3_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(protein3_plot, str_c("Protein3_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "Protein3")
  }
  
  if(!is.null(dpmsr_set$x$protein4_list)){
    protein4_plot <-subset(data_in, Accession %in% dpmsr_set$x$protein4_list)  
    protein4_plot <- protein4_plot[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    bar_plot(protein4_plot, str_c("Protein4_", plot_title), plot_dir)
    dpmsr_set$y$protein_list <<- c(dpmsr_set$y$protein_list, "Protein4")
  }
}
  
 #------------------------------------------------------------------------------------------

cv_stats_plot <- function(summary_cv, total_cv) {
  # create summary table for CV's for groups under different normalize conditions
  Simple_Excel(summary_cv, str_c(file_prefix1, "_Average_CV.xlsx", collapse = " "))
  
  total_cv <- total_cv[2:ncol(total_cv)]
  png(filename=str_c(output_dir, "CV_Comparison", "_boxplot.png"), width = 800, height = 600)
  data_box <- log2(total_cv)
  data_box[data_box ==-Inf ] <- NA
  boxplot(data_box, 
          col = color_choices[1:group_number], 
          notch = TRUE, 
          boxwex = 0.8,
          main = c("CV Comparison", group_title),
          axes=TRUE,
          horizontal = TRUE,
          las=1,
          par(mar=c(8,10,4,2)))
  dev.off()
}

