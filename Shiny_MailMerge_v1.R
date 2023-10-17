report_templates <- function(session, input){
  
    df <- as.character(dpmsr_set$x$Project)
    df <- data.frame(df)
    colnames(df) <- "Project"
    df$Report_date <- as.character(dpmsr_set$x$Report_date)
    df$Client <- as.character(dpmsr_set$x$Client)
    df$PI <- as.character(dpmsr_set$x$PI )
    df$Database <- as.character(dpmsr_set$x$Database)
    df$Species <- as.character(dpmsr_set$x$Species) 
    df$Database_date <- as.character(dpmsr_set$x$Database_date)
    df$Sample_number <- as.character(dpmsr_set$x$Sample_number)
    df$Replicate_number <- as.character(dpmsr_set$x$Replicate_number)
    df$Injection_volume <- as.character(dpmsr_set$x$Injection_volume)
    df$Injection_percent <- as.character(as.numeric(dpmsr_set$x$Injection_volume)/as.numeric(dpmsr_set$x$Final_volume)*100)
    df$SPQC_aliquot <- as.character(dpmsr_set$x$SPQC_aliquot)
    df$ADH_recon_conc <- as.character(dpmsr_set$x$ADH_recon_conc)
    df$Gradient_start <- as.character(dpmsr_set$x$Gradient_start) 
    df$Gradient_end <- as.character(dpmsr_set$x$Gradient_end)
    df$Gradient_length <- as.character(dpmsr_set$x$Project)
    df$Final_volume <- as.character(dpmsr_set$x$Final_volume)
    df$Strap_binding_buffer <- as.character(trunc(as.numeric(dpmsr_set$x$Strap_binding_buffer),0))
    df$Spike_1 <- as.character(dpmsr_set$x$Spike_1)
    df$Spike_2 <- as.character(dpmsr_set$x$Spike_2)
    df$Spike_unit <- as.character(dpmsr_set$x$Spike_unit) 
   
    sample_groups <- dpmsr_set$y$sample_groups[dpmsr_set$y$sample_groups$Group != "SPQC",]
    df$Sample_groups <- as.character(str_c(sample_groups$Group, collapse = ", "))
    comma_pos <- str_locate_all(df$Sample_groups, ",")[[1]]
    comma_pos_last <- comma_pos[nrow(comma_pos), ] 
    str_sub(df$Sample_groups, comma_pos_last[1], comma_pos_last[2]) <- " and" 

    df$Total_samples <- as.character(dpmsr_set$y$sample_number)
    df$SPQC_injections <- dpmsr_set$y$sample_groups[dpmsr_set$y$sample_groups$Group=="SPQC",]$Count
    df$MSMS_count <- dpmsr_set$overview$MSMS
    df$PSM_count <- dpmsr_set$overview$PSM
    df$Peptide_count <- nrow(dpmsr_set$data$data_to_norm)
    if (dpmsr_set$x$final_data_output == "Protein"){
        df$Protein_count <- nrow(dpmsr_set$data$final$impute)
    }
    
    #test <- input$select_final_data_report
    test <- "sltmm"
    df$SPQC_CV <- as.character(round(as.numeric(dpmsr_set$data$summary_cv[dpmsr_set$data$summary_cv$Group == "SPQC",][test]),1))
    
    sample_cv <- dpmsr_set$data$summary_cv[dpmsr_set$data$summary_cv$Group != "SPQC",]
    sample_cv[test] <- round(sample_cv[test],1)
    sample_cv[test] <- lapply(sample_cv[test], as.character)
    df$Sample_cv <- str_c(unlist(sample_cv[test]), collapse = ", ")
    comma_pos <- str_locate_all(df$Sample_cv, ",")[[1]]
    comma_pos_last <- comma_pos[nrow(comma_pos), ] 
    str_sub(df$Sample_cv, comma_pos_last[1], comma_pos_last[2]) <- " and" 
    
    Simple_Excel(df, "Mail_Merge", str_c(dpmsr_set$file$extra_dir, "Mail_Merge.xlsx", collapse = " "))

}
