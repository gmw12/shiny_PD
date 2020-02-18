project_overview <- function(){
  
  df_peptide <- dpmsr_set$data$data_raw_peptide  
  
  count_peptides <- nrow(df_peptide)
  
  df_peptide$Oxidation <- grepl("Oxidation", df_peptide$Modifications)
  df_peptide$Deamidated <- grepl("Deamidated", df_peptide$Modifications)
  df_peptide$Carbamidomethyl <- grepl("Carbamidomethyl", df_peptide$Modifications)
  df_peptide$Phospho <- grepl("Phospho", df_peptide$Modifications)
  df_peptide$GG <- grepl("GlyGly", df_peptide$Modifications)
  df_peptide$semi <- grepl("^\\[.*(R|K).*\\]", df_peptide$Annotated.Sequence) * ( grepl("(R|K).\\[.*\\]$", df_peptide$Annotated.Sequence) | 
                                                                                     grepl("\\[-\\]$", df_peptide$Annotated.Sequence) )

  stat_list <-list()
  stat_list[["MSMS"]] <- nrow(dpmsr_set$data$data_raw_msms)
  stat_list[["PSM"]] <- nrow(dpmsr_set$data$data_raw_psm)
  stat_list[["PSM/MSMS%"]] <-  nrow(dpmsr_set$data$data_raw_psm)/nrow(dpmsr_set$data$data_raw_msms)*100
  stat_list[["Peptides"]] <- nrow(dpmsr_set$data$data_raw_peptide)
  stat_list[["Proteins"]] <- nrow(dpmsr_set$data$data_raw_protein)
  stat_list[["Oxidation%"]] <- round(length(df_peptide$Oxidation[df_peptide$Oxidation==TRUE])/count_peptides*100,1)
  stat_list[["Deamidated%"]] <- round(length(df_peptide$Deamidated [df_peptide$Deamidated ==TRUE])/count_peptides*100,1)
  stat_list[["Carbamidomethyl%"]] <- round(length(df_peptide$Carbamidomethyl[df_peptide$Carbamidomethyl==TRUE])/count_peptides*100,1)
  stat_list[["Phospho"]] <- length(df_peptide$Phospho[df_peptide$Phospho==TRUE])
  stat_list[["GlyGly"]] <- length(df_peptide$GG[df_peptide$GG==TRUE])
  stat_list[["Missed%"]] <- round(length(dpmsr_set$data$data_raw_peptide$Number.of.Missed.Cleavages[dpmsr_set$data$data_raw_peptide$Number.of.Missed.Cleavages>0])/count_peptides*100,1)
  stat_list[["Semi%"]] <- round(length(df_peptide$semi[df_peptide$semi==0]) /count_peptides *100,1) 
  
  dpmsr_set$overview <<- stat_list

  #---------------------------------------------------------------------
  
  df_ms <- dpmsr_set$data$data_raw_msms
  df_psm <- dpmsr_set$data$data_raw_psm
  df_peptide <- dpmsr_set$data$data_raw_peptide
  
  
  ms_bins_start <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,200)
  ms_bins_stop <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,200,300)
  
  psm_inj_summary <- data.frame(ms_bins_start)
  psm_inj_summary$ms_bins_stop <- ms_bins_stop
  psm_inj_summary$name <- str_c(psm_inj_summary$ms_bins_start, "_", psm_inj_summary$ms_bins_stop)
  psm_inj_summary$total <- 1
  
  for(i in 1:length(ms_bins_start)){
    psm_inj_summary$total[i] <- length(df_psm$Ion.Inject.Time.in.ms[df_psm$Ion.Inject.Time.in.ms > psm_inj_summary$ms_bins_start[i] & 
                                                                      df_psm$Ion.Inject.Time.in.ms <= psm_inj_summary$ms_bins_stop[i] ])
  }
  
  ms_inj_summary <- data.frame(ms_bins_start)
  ms_inj_summary$ms_bins_stop <- ms_bins_stop
  ms_inj_summary$name <- str_c(ms_inj_summary$ms_bins_start, "_", ms_inj_summary$ms_bins_stop)
  ms_inj_summary$total <- 1
  
  for(i in 1:length(ms_bins_start)){
    ms_inj_summary$total[i] <- length(df_ms$Ion.Inject.Time.in.ms[df_ms$Ion.Inject.Time.in.ms > ms_inj_summary$ms_bins_start[i] & 
                                                                    df_ms$Ion.Inject.Time.in.ms <= ms_inj_summary$ms_bins_stop[i] ])
  }
  
  df<- psm_inj_summary
  ggplot(df, aes(x=name, y=total)) + 
    geom_bar(stat="identity", width=.5, fill="blue") + 
    scale_x_discrete(limits = df$name) +
    labs(title="Project Ion Inject Summary", x="Inject Time(ms)", y="Count") + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
    file_name <- str_c(dpmsr_set$file$qc_dir, "PSM_inj_summary.png")
    ggsave(file_name, width=5, height=3)
    
#----------------------------------------    
    
    df2 <- data.frame(cbind(df_psm$Rank,  df_psm$Delta.M.in.ppm, df_psm$Ion.Inject.Time.in.ms ))
    df2$psm <- "All"
    df2$psm <- as.factor(df2$psm)

    
    g <- ggplot(df2, aes(psm, X2))
    g + geom_violin(fill="blue") + 
      labs(title="PSM Mass Accuracy", 
           x="PSM's",
           y="ppm")    
    file_name <- str_c(dpmsr_set$file$qc_dir, "Mass_Accuracy.png")
    ggsave(file_name, width=5, height=3)

}
