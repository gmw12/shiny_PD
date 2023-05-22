project_overview <- function(session, input, output){
  cat(file=stderr(), "Starting project_overview...", "\n")
  library(Peptides)
  start_time <- Sys.time()
  plot_dir <- dpmsr_set$file$qc_dir
  
  df_peptide <- dpmsr_set$data$data_raw_peptide  
  
  #Gather peptide info for stat_list
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
  
  
  cat(file=stderr(), "Starting project_overview... Stat List", "\n")
  
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
  test_data <- dpmsr_set$data$data_raw_peptide %>% dplyr::select(contains("Abundance.F"))
  total_missing <- table(is.na(test_data))
  stat_list[["Peptide_IDs"]] <- total_missing[1]
  stat_list[["Peptide_Missing"]] <- total_missing[2]
  stat_list[["Missing%"]] <- round(total_missing[2]/total_missing[1]*100,1)
  dpmsr_set$overview <<- stat_list

  #---------------------------------------------------------------------
  cat(file=stderr(), "Starting project_overview... data processing", "\n")
  
  df_ms <- dpmsr_set$data$data_raw_msms
  df_psm <- dpmsr_set$data$data_raw_psm
  df_peptide <- dpmsr_set$data$data_raw_peptide
  
  ms_bins <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,200,300)
  
  psm_inj_summary <- data.frame(ms_bins[1:22])
  psm_inj_summary$ms_bins_stop <- ms_bins[2:23]
  psm_inj_summary$name <- str_c(psm_inj_summary$ms_bins_start, "_", psm_inj_summary$ms_bins_stop)
  psm_inj_summary$total <- tapply(df_psm$Ion.Inject.Time.in.ms, cut(df_psm$Ion.Inject.Time.in.ms, breaks = ms_bins), length)
  
  ms_inj_summary <- data.frame(ms_bins[1:22])
  ms_inj_summary$ms_bins_stop <- ms_bins[2:23]
  ms_inj_summary$name <- str_c(ms_inj_summary$ms_bins_start, "_", ms_inj_summary$ms_bins_stop)
  ms_inj_summary$total <- tapply(df_ms$Ion.Inject.Time.in.ms, cut(df_ms$Ion.Inject.Time.in.ms, breaks = ms_bins), length)
  
  feature_bins <- c(0,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
  feature_width_summary <- data.frame(feature_bins[1:19])
  feature_width_summary$feature_bins_stop <- feature_bins[2:20]
  feature_width_summary$name <- str_c(feature_width_summary$feature_bins_start, "_", feature_width_summary$feature_bins_stop)
  feature_width_summary$total <- 0
  
  if(!is.null(dpmsr_set$data$data_features)){
    df_features <- (dpmsr_set$data$data_features$Right.RT.in.min - dpmsr_set$data$data_features$Left.RT.in.min)*60
    df_features <- df_features[df_features <100]
    feature_width_summary$total <- tapply(df_features, cut(df_features, breaks = feature_bins), length)
  }

  #----------------------------------------    
  cat(file=stderr(), "Starting project_overview... PSM_Inj", "\n")
  
  psm_inj <- function(df, plot_dir){  
    ggplot2::ggplot(df, ggplot2::aes(x=name, y=total)) + 
      ggplot2::geom_bar(stat="identity", width=.5, fill="darkorange") + 
      ggplot2::scale_x_discrete(limits = df$name) +
      ggplot2::labs(title="Project Ion Inject Summary", x="Inject Time(ms)", y="Count") + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=65, vjust=0.6))
      file_name <- stringr::str_c(plot_dir, "PSM_inj_summary.png")
      ggplot2::ggsave(file_name, width=5, height=3)
  }  
  
    test_rt <- callr::r_bg(psm_inj, args=list(psm_inj_summary, plot_dir), supervise = TRUE)

    #---Mass Accuracy-------------------------------------    
    cat(file=stderr(), "Starting project_overview... MassAccuracy", "\n")
    
    mass_accuracy <- function(df2, plot_dir){  
      df2$psm <- "All"
      df2$psm <- as.factor(df2$psm)
      
      g <- ggplot2::ggplot(df2[df2$X2 < 10,], ggplot2::aes(psm, X2))
      g + ggplot2::geom_violin(fill="blue") + 
        ggplot2::labs(title="PSM Mass Accuracy", 
             x="PSM's",
             y="ppm")    
      file_name <- stringr::str_c(plot_dir, "Mass_Accuracy.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }
    
    df2 <- data.frame(cbind(df_psm$Rank, df_psm$Delta.M.in.ppm, df_psm$Ion.Inject.Time.in.ms))
    test_ma <- callr::r_bg(mass_accuracy, args=list(df2, plot_dir), supervise = TRUE)

    
    #--RT--------------------------------------    
    cat(file=stderr(), "Starting project_overview... RT", "\n")
    
    rt_summary <- function(data, plot_dir){  
      df2 <- data |> dplyr::select(dplyr::contains('RT.in.min.by.Search.Engine.'))
      colnames(df2) <- c("RT")
      
      g <- ggplot2::ggplot(df2, ggplot2::aes(x=RT))
      g + ggplot2::geom_density(fill="darkgreen") + ggplot2::theme_classic() +
        ggplot2::labs(title="Peptide RT", 
             x="RT, min",
             y="Density")   +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  
      file_name <- stringr::str_c(plot_dir, "Peptide_RT.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }
  
    test_rt <- callr::r_bg(rt_summary, args=list(df_peptide, plot_dir), supervise = TRUE)
    
    
#--RT width--------------------------------------    
    cat(file=stderr(), "Starting project_overview... RT Width", "\n")
  
    rtw_summary <- function(df, plot_dir){  
      ggplot2::ggplot(df, ggplot2::aes(x=name, y=total)) + 
        ggplot2::geom_bar(stat="identity", width=.5, fill="darkorange") + 
        ggplot2::scale_x_discrete(limits = df$name) +
        ggplot2::labs(title="Feature Peak Width", x="seconds", y="Count") + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=65, vjust=0.6))
      file_name <- stringr::str_c(plot_dir, "Feature_Peak_Width.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }  
    
    test_rt <- callr::r_bg(psm_inj, args=list(psm_inj_summary, plot_dir), supervise = TRUE)  
    
  rtw_summary2 <- function(df2, plot_dir){    
      g <- ggplot2::ggplot(df2, ggplot2::aes(x=RT_width))
      g + ggplot2::geom_density(fill="purple") + ggplot2::theme_classic() +
        ggplot2::labs(title="Feature Peak Width", 
             x="Peak Width, sec",
             y="Density")   +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  
      file_name <- stringr::str_c(plot_dir, "Feature_Peak_Width.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }
  
  test_rtw <- callr::r_bg(rtw_summary, args=list(feature_width_summary, plot_dir), supervise = TRUE)

  
    #--MZ--------------------------------------  
    cat(file=stderr(), "Starting project_overview... MZ", "\n")
    
    mz_summary <- function(data, plot_dir){
      df2 <- data |> dplyr::select(dplyr::starts_with('mz.in.Da.by.Search.Engine.'))
      colnames(df2) <- "Mass"
      g <- ggplot2::ggplot(df2, ggplot2::aes(x=Mass))
      g + ggplot2::geom_density(fill="darkred") + ggplot2::theme_classic() +
        ggplot2::labs(title="Peptide mz", 
             x="m/z",
             y="Density")   +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  
      file_name <- stringr::str_c(plot_dir, "Peptide_MZ.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }
    
    test_mz <- callr::r_bg(mz_summary, args=list(df_peptide, plot_dir), supervise = TRUE)
    
 #--Charge----------------------------------------------------------   
    cat(file=stderr(), "Starting project_overview... Charge", "\n")
    
    charge_summary <- function(data, plot_dir){
      df3 <- data |> dplyr::select(dplyr::contains('Charge.by.Search.Engine.'))
      df3$count <- 1
      colnames(df3) <- c("Charge", "Count")
      df3 <- df3 |> dplyr::group_by(Charge) |> dplyr::summarise_all(list(sum))
      df3$Charge <- factor(df3$Charge, levels=df3$Charge)
      
      g <- ggplot2::ggplot(data=df3, ggplot2::aes(x=Charge, y=Count, fill=Charge))
      g +  ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_brewer(palette = "Set1") + 
        ggplot2::ggtitle("Peptide ID by Charge") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) 
      file_name <- stringr::str_c(plot_dir, "Peptide_Charge.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }
    
    test_charge <- callr::r_bg(charge_summary, args=list(df_peptide, plot_dir), supervise = TRUE)

  #--FAIMS----------------------------------------------------------  
  cat(file=stderr(), "Starting project_overview... FAIMS", "\n")
    
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
    cat(file=stderr(), "Starting project_overview... PI", "\n")
  
    pi_summary <- function(data, plot_dir){
      df2 <- data.frame(Peptides::pI(data, pKscale="EMBOSS"))
      colnames(df2) <- "PI"
      g <- ggplot2::ggplot(df2, ggplot2::aes(x=PI))
      g + ggplot2::geom_density(fill="darkgray") + ggplot2::theme_classic() +
        ggplot2::labs(title="Peptide PI", 
             x="PI",
             y="Density")   +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  
      file_name <- stringr::str_c(plot_dir, "Peptide_PI.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }
    
    test_pi <- callr::r_bg(pi_summary, args=list(dpmsr_set$data$data_raw_peptide$Sequence, plot_dir), supervise = TRUE)
    
    #--Cruc--------------------------------------     
    cat(file=stderr(), "Starting project_overview... Cruc", "\n")
    
    cruc_summary <- function(data, plot_dir){
      df2 <- data.frame(t(data.frame(Peptides::crucianiProperties(data))))
      rownames(df2) <- NULL
      colnames(df2) <- c("Polarity", "Hydrophobicity", "H-bonding")
      df3 <- tidyr::gather(df2)
      g <- ggplot2::ggplot(df3, ggplot2::aes(x=value, fill=key))
      g + ggplot2::geom_density(alpha=0.4) + ggplot2::theme_classic() +
        ggplot2::labs(title="Peptide Properities", 
             x="Scale",
             y="Density")   +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  
      file_name <- stringr::str_c(plot_dir, "Peptide_Cruc.png")
      ggplot2::ggsave(file_name, width=5, height=3) 
    }
    
    test_cruc <- callr::r_bg(cruc_summary, args=list(dpmsr_set$data$data_raw_peptide$Sequence, plot_dir), supervise = TRUE)
    
    #--AA Summary----------------------------------------------------------   
    cat(file=stderr(), "Starting project_overview... AA Summary", "\n")
    
    aa_summary <- function(data, plot_dir){
      df <- data.frame(Peptides::aaComp(data))
      df <- t(data.frame(df))
      toDelete <- seq(1, nrow(df), 2)
      df <- df[toDelete,]
      df2 <- data.frame(colnames(df))
      df2$sum <- colSums(df)
      colnames(df2) <- c("class", "sum")
      df2$class <- factor(df2$class, levels=df2$class)
      
      g <- ggplot2::ggplot(data=df2, ggplot2::aes(x=class, y=sum, fill=class) )
      g +  ggplot2::geom_bar(stat = "identity")  + ggplot2::theme_classic() +
        ggplot2::scale_fill_brewer(palette = "Set1") + 
        ggplot2::ggtitle("AA Type") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) 
      file_name <- stringr::str_c(plot_dir, "Peptide_AAInfo.png")
      ggplot2::ggsave(file_name, width=5, height=3)
  }
    test_aa <- callr::r_bg(aa_summary, args=list(dpmsr_set$data$data_raw_peptide$Sequence, plot_dir), supervise = TRUE)
    
    #--aliphatic index--------------------------------------     
    cat(file=stderr(), "project_overview... aliphatic index", "\n")

    aliphatic_index <- function(data, plot_dir){
      df <- data.frame(Peptides::aIndex(data))
      colnames(df) <- "AI"
      g <- ggplot2::ggplot(df, ggplot2::aes(x=AI))
      g + ggplot2::geom_density(fill="lightblue") + ggplot2::theme_classic() +
        ggplot2::labs(title="Peptide Aliphatic Index", 
             x="Aliphatic Index",
             y="Density")   +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  
      file_name <- stringr::str_c(plot_dir, "Peptide_AI.png")
      ggplot2::ggsave(file_name, width=5, height=3)
    }
  
    test_ai <- callr::r_bg(aliphatic_index, args=list(dpmsr_set$data$data_raw_peptide$Sequence, plot_dir), supervise = TRUE)
    
    cat(file=stderr(), "Starting project_overview... trim dpmsr_set", "\n")
    dpmsr_set$data$data_raw_msms <<- NULL
    dpmsr_set$data$data_raw_psm <<- NULL
    dpmsr_set$data$data_raw_decoypsm <<- NULL
    dpmsr_set$data$data_raw_decoypeptide <<- NULL
    dpmsr_set$data$data_raw_decoyprotein <<- NULL
    dpmsr_set$data$data_features <<- NULL
    
    detach(package:Peptides, unload=TRUE)
    gc()
    cat(file=stderr(), "Project_overview... complete", "\n")
    cat(file=stderr(), str_c("Overview load time ---> ", Sys.time()-start_time), "\n")
}











