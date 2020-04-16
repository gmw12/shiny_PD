
#--------------------------------------------------------------------------------------
stat_prep <- function(){
  ncores <- detectCores()
  data_list <- dpmsr_set$data$impute
  dpmsr_set$data$final <<- mclapply(data_list, stat_prep_parallel, mc.cores = ncores)
  
  #check that the parallel processing went through, if not do it again one at a time
  check_stats_prep_parallel(data_list)
  
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    dpmsr_set$data$finalraw <<- collapse_peptide(dpmsr_set$data$normalized$impute)
    Simple_Excel(dpmsr_set$data$finalraw, str_c(dpmsr_set$file$extra_prefix, "_final_raw_protein.xlsx", collapse = " "))
  } else{
    dpmsr_set$data$finalraw <<- dpmsr_set$data$normalized$impute
    Simple_Excel(dpmsr_set$data$finalraw, str_c(dpmsr_set$file$extra_prefix, "_final_raw_peptide.xlsx", collapse = " "))
  }
}
#--------------------------------------------------------------------------------------
stat_prep_parallel <- function(data_list){
  stat_prep_collapse(data_list)
}

#--------------------------------------------------------------------------------------
stat_prep_collapse <- function(data_in) {
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    data_in <- collapse_peptide(data_in)
  }
  
  if (as.logical(dpmsr_set$x$peptide_ptm_out)){
    data_in <- data_in [grep(dpmsr_set$x$peptide_report_grep, data_in$Modifications),]
  }
  
  info_columns_final <- ncol(data_in)-dpmsr_set$y$sample_number
  annotate_in <- data_in[1:info_columns_final]
  data_in <- data_in[(info_columns_final+1):ncol(data_in)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  #generate %CV's for each group
  for(i in 1:dpmsr_set$y$group_number) 
  {
    stat_df[ , str_c(dpmsr_set$y$sample_groups$Group[i], "_CV")] <- percentCV_gw(data_in[dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i]])
  } 
  
  stat_df <- stat_df[2:ncol(stat_df)]
  colnames(data_in) <- dpmsr_set$design$Header2
  data_table <- cbind(annotate_in, data_in, stat_df)
  
  
  return(data_table)
}

#--------------------------------------------------------------------------------------
check_stats_prep_parallel <- function(data_list){
  for(data_name in names(data_list)){
    if(is.null(dpmsr_set$data$final[[data_name]]    )){
      cat(file=stderr(), str_c("Stats Prep Parallel function...", data_name), "\n")
      dpmsr_set$data$final[[data_name]] <<- stat_prep_collapse(dpmsr_set$data$impute[[data_name]])
    }
  }
}

#--------------------------------------------------------------------------------------

 #create data frame for comparisons
set_stat_groups <- function(session, input, output){
  cat(file=stderr(), "set_stat_groups....1", "\n")
  stats_data <- dpmsr_set$data$impute[[input$select_final_data_stats]]
  stats_data <- stats_data[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
  
  if(input$select_final_data_stats == "impute"){
    colnames(stats_data) <- dpmsr_set$design$Header3
  }else{
    colnames(stats_data) <- dpmsr_set$design$Header2
  }
  comp_groups <- data.frame(seq(from=1, to=as.numeric(input$comp_number)))
  colnames(comp_groups) <- "CompNumber"
  color_choices <- distinctColorPalette(as.numeric(input$comp_number)*2)
  color_choices_N <- color_choices[c(TRUE, FALSE)]
  color_choices_D <- color_choices[c(FALSE, TRUE)]
  cat(file=stderr(), "set_stats_groups....2", "\n")
  
  for (i in 1:input$comp_number){
    cat(file=stderr(), str_c("stats_comp ", i, " of ", input$comp_number), "\n")
    stats_data_N <- stats_data
    stats_data_D <- stats_data
    
    for(stats_group in input[[str_c('comp_',i,'N')]]){stats_data_N <- stats_data_N %>% dplyr::select(contains(stats_group))   }
    for(stats_group in input[[str_c('comp_',i,'D')]]){stats_data_D <- stats_data_D %>% dplyr::select(contains(stats_group))   }
    
    comp_groups$comp_N[i] <- paste(input[[str_c('comp_',i,'N')]], collapse = "_")
    comp_groups$comp_D[i] <- paste(input[[str_c('comp_',i,'D')]], collapse = "_")
    comp_groups$comp_name[i] <- str_c(comp_groups$comp_N[i], "_v_", comp_groups$comp_D[i])
    comp_groups$com_N_headers[i] <- list(colnames(stats_data_N))
    comp_groups$com_D_headers[i] <- list(colnames(stats_data_D))
    comp_groups$fc[i] <- str_c(comp_groups$comp_name[i], "_FC")
    comp_groups$fc2[i] <- str_c(comp_groups$comp_name[i], "_FC2")
    comp_groups$pval[i] <- str_c(comp_groups$comp_name[i], "_pval")
    comp_groups$limma_pval[i] <- str_c(comp_groups$comp_name[i], "_limma_pval")
    comp_groups$exactTest[i] <- str_c(comp_groups$comp_name[i], "_exactTest")
    comp_groups$mf[i] <- str_c(comp_groups$comp_name[i], "_MF")
    comp_groups$N_count[i] <- ncol(stats_data_N)
    comp_groups$D_count[i] <- ncol(stats_data_D)
    comp_groups$N_color[1] <- color_choices_N[i]
    comp_groups$D_color[1] <- color_choices_D[i]
    updatePickerInput(session, inputId = str_c('comp_',i,'N'), label=str_c("Comp",i,"_N  selected -> ", ncol(stats_data_N)) )
    updatePickerInput(session, inputId = str_c('comp_',i,'D'), label=str_c("Comp",i,"_D  selected -> ", ncol(stats_data_D)) )
    updateTextInput(session, inputId = str_c("volcano",i,"_stats_plot_title"),  value = str_c("Volcano ", comp_groups$comp_name[i]))
    dpmsr_set$y$stats[[str_c("comp", i, "_N")]] <<- input[[str_c('comp_',i,'N')]]
    dpmsr_set$y$stats[[str_c("comp", i, "_D")]] <<- input[[str_c('comp_',i,'D')]]
    comp_groups$sample_numbers_N[i] <- list(match(colnames(stats_data_N), names(stats_data)))
    comp_groups$sample_numbers_D[i] <- list(match(colnames(stats_data_D), names(stats_data)))
  }
  
  cat(file=stderr(), "set_stats_groups....4", "\n")
  dpmsr_set$y$stats$groups <<- comp_groups
  dpmsr_set$y$stats$comp_number <<- input$comp_number
  
  updatePickerInput(session, "stats_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name)
  updateTextInput(session, "stats_barplot_title", value = str_c("Barplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  updateTextInput(session, "stats_boxplot_title", value = str_c("Boxplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  updateSelectInput(session, "stats_select_data_comp", choices = dpmsr_set$y$stats$groups$comp_name, 
                    selected = dpmsr_set$y$stats$groups$comp_name[1])
  
}


#--------------------------------------------------------------------------------------------------------------------------------
#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_calc <- function(session, input, output) {
  data_in <- dpmsr_set$data$impute[[input$select_final_data_stats]]
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    data_in <- collapse_peptide(data_in)
  }
  if (as.logical(dpmsr_set$x$peptide_ptm_out)){
    data_in <- data_in [grep(dpmsr_set$x$peptide_report_grep, data_in$Modifications),]
  }
  
  annotate_in <- data_in[1:dpmsr_set$y$info_columns_final]
  data_in <- data_in[(dpmsr_set$y$info_columns_final+1):ncol(data_in)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  imputed_df <- dpmsr_set$data$Protein_imputed_df[3:ncol(dpmsr_set$data$Protein_imputed_df)]  
  
  for(i in 1:nrow(dpmsr_set$y$stats$groups)) 
  {
    stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_N[i], "_CV")] <- percentCV_gw(data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ])
    stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_D[i], "_CV")] <- percentCV_gw(data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ])
  } 

  #generate pvalue and FC for each comparison
  for(i in 1:nrow(dpmsr_set$y$stats$groups))
  {
    comp_N_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    stat_df[ ,dpmsr_set$y$stats$groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
    #stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
    #stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
    comp_N_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    stat_df[ ,dpmsr_set$y$stats$groups$mf[i]] <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
  } 
  
  # Create tables for excel--------------------------------------------------
  if(input$select_final_data_stats == "impute"){
    colnames(data_in) <- dpmsr_set$design$Header3
  }else{
    colnames(data_in) <- dpmsr_set$design$Header2
  }
  data_table <- cbind(annotate_in, data_in, stat_df[2:ncol(stat_df)])
  dpmsr_set$data$stats$final <<- data_table
  dpmsr_set$data$stats$final_comp <<- input$select_final_data_stats
  dpmsr_set$y$stats$pvalue_cutoff <<- input$stats_pvalue_cutoff
  dpmsr_set$y$stats$foldchange_cutoff <<- input$stats_foldchange_cutoff
  
  return()
}

#----------------------------------------------------------------------------------------
# create final excel documents
stats_Final_Excel <- function() {
  require(openxlsx)
    
    filename <- str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp, "//Final_", 
                      dpmsr_set$data$stats$final_comp,  "_stats.xlsx")
    df <- dpmsr_set$data$stats$final
    #remove FC2 from df for excel
    df <- df[,-grep(pattern="_FC2", colnames(df))]
    
    if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
      df_raw <- dpmsr_set$data$finalraw
    }else{
      df_raw <- dpmsr_set$data$finalraw
      #if PTM out need to reduce raw data frame for excel
      if (as.logical(dpmsr_set$x$peptide_ptm_out)){
        df_raw <- df_raw[grep(dpmsr_set$x$peptide_report_grep, df_raw$Modifications),]
      }
    }
    
    # save the option to save normalized table with raw
    #df2 <- cbind(df_raw, df[(dpmsr_set$y$info_columns_final+1):ncol(df)])
    
    df2 <- df
    
    if (as.logical(dpmsr_set$x$accession_report_out)){
      df <-subset(df, Accession %in% dpmsr_set$x$accession_report_list )
      df2 <-subset(df2, Accession %in% dpmsr_set$x$accession_report_list )
    }
    
    nextsheet <- 1
    
    wb <- createWorkbook()
    
    if(dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Protein"){
      
      raw_protein <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Protein_to_Protein_Raw.xlsx"))
      raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_ProteinPeptide_to_Peptide_Raw.xlsx"))
      
      if(as.logical(dpmsr_set$x$peptide_isoform) && dpmsr_set$x$raw_data_input=="Peptide"){
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Isoform_to_Isoform_Raw.xlsx"))
      }
      
      if(!as.logical(dpmsr_set$x$peptide_isoform) && dpmsr_set$x$raw_data_input=="Peptide"){
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Peptide_to_Peptide_Raw.xlsx"))
      }
      
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet +1
      addWorksheet(wb, "Raw Peptide Data")
      writeData(wb, sheet = nextsheet, raw_peptide)
      nextsheet <- nextsheet +1
      addWorksheet(wb, "Raw Protein Data")
      writeData(wb, sheet = nextsheet, raw_protein) 
      nextsheet <- nextsheet +1
    }else if (dpmsr_set$x$raw_data_input=="Peptide"){
      if(as.logical(dpmsr_set$x$peptide_isoform)){
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Isoform_to_Isoform_Raw.xlsx"))
      }else{
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Peptide_to_Peptide_Raw.xlsx"))
      }
      addWorksheet(wb, "Raw Peptide Data")
      writeData(wb, sheet = nextsheet, raw_peptide) 
      nextsheet <- nextsheet +1
    }
    
    #addWorksheet(wb, deparse(substitute(df2)))
    addWorksheet(wb, "Normalized Data")
    writeData(wb, sheet = nextsheet, df2)  
    nextsheet <- nextsheet +1
    z=1
    for(i in 1:as.numeric(dpmsr_set$y$stats$comp_number))  {
      comp_string <- dpmsr_set$y$stats$groups$comp_name[i]
      df_N <- df %>% dplyr::select(unlist(dpmsr_set$y$stats$groups$com_N_headers[i]))  
      df_D <- df %>% dplyr::select(unlist(dpmsr_set$y$stats$groups$com_D_headers[i]) ) 
      df_Comp <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$comp_name [i]) ) 
      df_final <- cbind(df[1:(dpmsr_set$y$info_columns_final)], df_N, df_D, df_Comp)
      
      filtered_df <- subset(df_final, df_final[ , dpmsr_set$y$stats$groups$pval[i]] <= as.numeric(dpmsr_set$x$pvalue_cutoff) &  
                              (df_final[ , dpmsr_set$y$stats$groups$fc[i]] >= as.numeric(dpmsr_set$x$foldchange_cutoff) | 
                                 df_final[ , dpmsr_set$y$stats$groups$fc[i]] <= -as.numeric(dpmsr_set$x$foldchange_cutoff))) 
      filtered_df <- subset(filtered_df, filtered_df[ , dpmsr_set$y$stats$groups$mf[i]] >= as.numeric(dpmsr_set$x$missing_factor))
      addWorksheet(wb, comp_string)
      writeData(wb, sheet = nextsheet, filtered_df)
      nextsheet <- nextsheet +1
      z <- z+2
    }
    cat(file=stderr(), "writting excel to disk...", "\n")
    saveWorkbook(wb, filename, overwrite = TRUE)
    
}






