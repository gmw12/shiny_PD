project_overview <- function(){
  require(Peptides)
  
  df_peptide <- dpmsr_set$data$data_raw_peptide  
  
  count_peptides <- nrow(df_peptide)
  
  df_peptide$Oxidation <- grepl("Oxidation", df_peptide$Modifications)
  df_peptide$Deamidated <- grepl("Deamidated", df_peptide$Modifications)
  df_peptide$Carbamidomethyl <- grepl("Carbamidomethyl", df_peptide$Modifications)
  df_peptide$Phospho <- grepl("Phospho", df_peptide$Modifications)
  df_peptide$GG <- grepl("GlyGly", df_peptide$Modifications)
  if("Annotated Sequence" %in% colnames(df_peptide)){
    df_peptide$semi <- grepl("^\\[.*(R|K).*\\]", df_peptide$Annotated.Sequence) * 
                            ( grepl("(R|K).\\[.*\\]$", df_peptide$Annotated.Sequence) | 
                            grepl("\\[-\\]$", df_peptide$Annotated.Sequence) )
  } else {
    df_peptide$semi <- 1
  }
  
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
    geom_bar(stat="identity", width=.5, fill="darkorange") + 
    scale_x_discrete(limits = df$name) +
    labs(title="Project Ion Inject Summary", x="Inject Time(ms)", y="Count") + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6))
    file_name <- str_c(dpmsr_set$file$qc_dir, "PSM_inj_summary.png")
    ggsave(file_name, width=5, height=3)
    
    
#----------------------------------------    
    
    df2 <- data.frame(cbind(df_psm$Rank,  df_psm$Delta.M.in.ppm, df_psm$Ion.Inject.Time.in.ms))
    df2$psm <- "All"
    df2$psm <- as.factor(df2$psm)

    
    g <- ggplot(df2, aes(psm, X2))
    g + geom_violin(fill="blue") + 
      labs(title="PSM Mass Accuracy", 
           x="PSM's",
           y="ppm")    
    file_name <- str_c(dpmsr_set$file$qc_dir, "Mass_Accuracy.png")
    ggsave(file_name, width=5, height=3)


    #--RT--------------------------------------    
    df2 <- data.frame(df_peptide$RT.in.min.by.Search.Engine.Mascot)
    colnames(df2) <- c("RT")
    
    g <- ggplot(df2, aes(x=RT))
    g + geom_density(fill="darkgreen") + theme_classic() +
      labs(title="Peptide RT", 
           x="RT, min",
           y="Density")   +
      theme(plot.title = element_text(hjust = 0.5))  
    file_name <- str_c(dpmsr_set$file$qc_dir, "Peptide_RT.png")
    ggsave(file_name, width=5, height=3)
    
    
    
    #--RT width--------------------------------------    
    df2 <- data.frame(dpmsr_set$data$data_features)
    df2$RT_width <- (df2$Right.RT.in.min-df2$Left.RT.in.min)*60
    df2 <- df2[(df2$RT_width < 100),]
    
    
    g <- ggplot(df2, aes(x=RT_width))
    g + geom_density(fill="purple") + theme_classic() +
      labs(title="Feature Peak Width", 
           x="Peak Width, sec",
           y="Density")   +
      theme(plot.title = element_text(hjust = 0.5))  
    file_name <- str_c(dpmsr_set$file$qc_dir, "Feature_Peak_Width.png")
    ggsave(file_name, width=5, height=3)
    
    #--MZ--------------------------------------     
    df2 <- data.frame(df_peptide$mz.in.Da.by.Search.Engine.Mascot)
    colnames(df2) <- "Mass"
    g <- ggplot(df2, aes(x=Mass))
    g + geom_density(fill="darkred") + theme_classic() +
      labs(title="Peptide mz", 
           x="m/z",
           y="Density")   +
      theme(plot.title = element_text(hjust = 0.5))  
    file_name <- str_c(dpmsr_set$file$qc_dir, "Peptide_MZ.png")
    ggsave(file_name, width=5, height=3)
    
 #--Charge----------------------------------------------------------   
    df3 <- data.frame(df_peptide$Charge.by.Search.Engine.Mascot)
    df3$count <- 1
    colnames(df3) <- c("Charge", "Count")
    df3 <- df3 %>% group_by(Charge) %>% summarise_all(list(sum))
    df3$Charge <- factor(df3$Charge, levels=df3$Charge)
    
    g <- ggplot(data=df3, aes(x=Charge, y=Count, fill=Charge) )
    g +  geom_bar(stat = "identity")  + theme_classic() +
      scale_fill_brewer(palette = "Set1") + 
      ggtitle("Peptide ID by Charge") +
      theme(plot.title = element_text(hjust = 0.5)) 
    file_name <- str_c(dpmsr_set$file$qc_dir, "Peptide_Charge.png")
    ggsave(file_name, width=5, height=3)



#--FAIMS----------------------------------------------------------   
  FAIMS_overview <- function(){  
    df3 <- data.frame(dpmsr_set$data$data_raw_msms$Comp.Voltage.in.V)
    df3$count <- 1
    colnames(df3) <- c("CV", "Count")
    df3 <- df3 %>% group_by(CV) %>% summarise_all(list(sum))
    df3$CV <- factor(df3$CV, levels=df3$CV)
    
    g <- ggplot(data=df3, aes(x=CV, y=Count, fill=CV) )
    g +  geom_bar(stat = "identity")  + theme_classic() +
      scale_fill_brewer(palette = "Set1") + 
      ggtitle("PSM ID by CV") +
      theme(plot.title = element_text(hjust = 0.5)) 
    file_name <- str_c(dpmsr_set$file$qc_dir, "FAIMS_PSM_CV.png")
    ggsave(file_name, width=5, height=3)
}

 try(FAIMS_overview(), silent=TRUE)
    
    #--PI--------------------------------------     
    df2 <- data.frame(pI(dpmsr_set$data$data_raw_peptide$Sequence, pKscale="EMBOSS"))
    colnames(df2) <- "PI"
    g <- ggplot(df2, aes(x=PI))
    g + geom_density(fill="darkgray") + theme_classic() +
      labs(title="Peptide PI", 
           x="PI",
           y="Density")   +
      theme(plot.title = element_text(hjust = 0.5))  
    file_name <- str_c(dpmsr_set$file$qc_dir, "Peptide_PI.png")
    ggsave(file_name, width=5, height=3)
    
    #--PI--------------------------------------     
    df2 <- data.frame(t(data.frame(crucianiProperties(dpmsr_set$data$data_raw_peptide$Sequence))))
    rownames(df2) <- NULL
    colnames(df2) <- c("Polarity", "Hydrophobicity", "H-bonding")
    df3 <- gather(df2)
    g <- ggplot(df3, aes(x=value, fill=key))
    g + geom_density(alpha=0.4) + theme_classic() +
      labs(title="Peptide Properities", 
           x="Scale",
           y="Density")   +
      theme(plot.title = element_text(hjust = 0.5))  
    file_name <- str_c(dpmsr_set$file$qc_dir, "Peptide_Cruc.png")
    ggsave(file_name, width=5, height=3)  
    
    #--AA Summary----------------------------------------------------------   
    df2 <- data.frame(aaComp(dpmsr_set$data$data_raw_peptide$Sequence))
    df2 <- t(data.frame(df2))
    toDelete <- seq(1, nrow(df2), 2)
    df2 <- df2[toDelete,]
    df3 <- data.frame(colnames(df2))
    df3$sum <- colSums(df2)
    colnames(df3) <- c("class", "sum")
    df3$class <- factor(df3$class, levels=df3$class)
    
    g <- ggplot(data=df3, aes(x=class, y=sum, fill=class) )
    g +  geom_bar(stat = "identity")  + theme_classic() +
      scale_fill_brewer(palette = "Set1") + 
      ggtitle("AA Type") +
      theme(plot.title = element_text(hjust = 0.5)) 
    file_name <- str_c(dpmsr_set$file$qc_dir, "Peptide_AAInfo.png")
    ggsave(file_name, width=5, height=3)
    
    #--aliphatic index--------------------------------------     
    df2 <- data.frame(aIndex(dpmsr_set$data$data_raw_peptide$Sequence))
    colnames(df2) <- "AI"
    g <- ggplot(df2, aes(x=AI))
    g + geom_density(fill="lightblue") + theme_classic() +
      labs(title="Peptide Aliphatic Index", 
           x="Aliphatic Index",
           y="Density")   +
      theme(plot.title = element_text(hjust = 0.5))  
    file_name <- str_c(dpmsr_set$file$qc_dir, "Peptide_AI.png")
    ggsave(file_name, width=5, height=3)
    
  
    dpmsr_set$data$data_raw_msms <<- NULL
    dpmsr_set$data$data_raw_psm <<- NULL
    dpmsr_set$data$data_raw_decoypsm <<- NULL
    dpmsr_set$data$data_features <<- NULL
}











