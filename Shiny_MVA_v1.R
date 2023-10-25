
#--------------------------------------------------------------------------------------
stat_prep <- function(){
  cat(file = stderr(), ("stat_prep function started..."), "\n")
  ncores <- detectCores()
  if (is.na(ncores) | ncores < 1) {ncores <- 1}

  data_list <- dpmsr_set$data$impute
  #clear previous data
  dpmsr_set$data$final <- NULL
  
  #iq:fast_maxlfq function has issue with parallel processing, only works on first pass (c++ ???)
  #directlfq is python code 
  cat(file = stderr(), (str_c("stat_prep function started...1  cores=", ncores)), "\n")
  if (dpmsr_set$x$rollup_method != "IQ_MaxLFQ" && dpmsr_set$x$rollup_method != "DirectLFQ") {
    dpmsr_set$data$final <<- mclapply(data_list, stat_prep_parallel, mc.cores = ncores)
  }
  
  #check that the parallel processing went through, if not do it again one at a time, MaxLFQ will run here
  cat(file = stderr(), ("stat_prep function started...2"), "\n")
  check_stats_prep_parallel(data_list)
  
  cat(file = stderr(), ("stat_prep function started...3"), "\n")
  if (dpmsr_set$y$state == "Peptide" && dpmsr_set$x$final_data_output == "Protein") {
    dpmsr_set$data$finalraw <<- collapse_peptide(dpmsr_set$data$normalized$impute)
    Simple_Excel(dpmsr_set$data$finalraw, "final_raw_protein", str_c(dpmsr_set$file$extra_prefix, "_final_raw_protein.xlsx", collapse = " "))
  } else{
    dpmsr_set$data$finalraw <<- dpmsr_set$data$normalized$impute
    Simple_Excel(dpmsr_set$data$finalraw, "final_raw_peptide", str_c(dpmsr_set$file$extra_prefix,  "_final_raw_peptide.xlsx", collapse = " "))
  }
  
  
  cat(file = stderr(), ("stat_prep function started...4"), "\n")
  
  if (is.null(dpmsr_set$data$directlfq$directlfq_protein_impute)) {
    cat(file = stderr(), ("no directlfq data..."), "\n")
  }else{
    dpmsr_set$data$final$directlfq_full <<- stat_prep_collapse(dpmsr_set$data$directlfq$directlfq_protein_impute, directlfq = TRUE) 
  }
  
  cat(file = stderr(), ("stat_prep function complete..."), "\n")
}

#--------------------------------------------------------------------------------------
save_final_nostats <- function(){
  #save to excel final without stats
  for (name in names(dpmsr_set$data$final)) {
    cat(file = stderr(), str_c("Saving final excel w/o stats....", name), "\n")
    filename <- str_c(dpmsr_set$file$output_dir, name, "//Final_", name, "_noStats.xlsx")
    Simple_Excel(dpmsr_set$data$final[[name]], name, filename)
  }
  
  for (name in names(dpmsr_set$data$impute)) {
    cat(file = stderr(), str_c("Saving imputed peptide excel w/o stats....", name), "\n")
    filename <- str_c(dpmsr_set$file$output_dir, name, "//Imputed_Peptide_", name, ".xlsx")
    Simple_Excel(dpmsr_set$data$impute[[name]], name, filename)
  }
}

#--------------------------------------------------------------------------------------
stat_prep_parallel <- function(data_list){
  stat_prep_collapse(data_list)
}

#--------------------------------------------------------------------------------------
stat_prep_collapse <- function(data_in, directlfq=FALSE) {
  
  if (dpmsr_set$y$state != "Protein" && dpmsr_set$x$final_data_output == "Protein" && !directlfq) {
    data_in <- collapse_peptide(data_in)
    }
  
  if (dpmsr_set$y$state != "Protein" && dpmsr_set$x$final_data_output != "Protein" && !directlfq) {
    data_in <- collapse_precursor(data_in)
  }
  
  
  #this sets the final table for graphs to use only those that will be reported
  if (as.logical(dpmsr_set$x$peptide_ptm_norm)) {
    if ("Modifications" %in% colnames(data_in)) {
      data_in <- data_in[grep(dpmsr_set$x$peptide_norm_grep, data_in$Modifications),]
    } else {
      data_in <- data_in[grep(dpmsr_set$x$peptide_norm_grep, data_in$Sequence),]
    } 
  }
  
  info_columns_final <- ncol(data_in) - dpmsr_set$y$sample_number
  annotate_in <- data_in[1:info_columns_final]

  data_in <- data_in[(info_columns_final + 1):ncol(data_in)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  #generate %CV's for each group
  for (i in 1:dpmsr_set$y$group_number) 
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
  for (data_name in names(data_list)) {
    if (is.null(dpmsr_set$data$final[[data_name]] )) {
      cat(file = stderr(), str_c("Stats Prep Parallel function...", data_name), "\n")
      dpmsr_set$data$final[[data_name]] <<- stat_prep_collapse(dpmsr_set$data$impute[[data_name]])
    }
  }
}

#--------------------------------------------------------------------------------------

 #create data frame for comparisons
set_stat_groups <- function(session, input, output){
  cat(file = stderr(), "set_stat_groups....1", "\n")
  
  # test for comp descriptions over 31 characters - limit for excel sheet name
  comp_warning <- FALSE
  comp_warning_list <- list()
  dpmsr_set$y$uniquegroups <<- unique(dpmsr_set$design$Group)
  
  stats_data <- dpmsr_set$data$final[[input$select_final_data_stats]]
  stats_data <- stats_data[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number)]
  
  if (input$select_final_data_stats == "impute") {
    colnames(stats_data) <- dpmsr_set$design$Header3
  }else{
    colnames(stats_data) <- dpmsr_set$design$Header2
  }
  
  #save spqc information
  dpmsr_set$y$stats$comp_spqc <<- input$comp_spqc
  spqc_df <- stats_data %>% dplyr::select(contains(input$comp_spqc))
  dpmsr_set$y$stats$comp_spqc_sample_numbers <<- list(match(colnames(spqc_df), names(stats_data)))
  dpmsr_set$y$stats$comp_spqc_number <<- length(unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers))
  dpmsr_set$y$stats$comp_spqc_sample_headers <<- list(colnames(spqc_df))
  
  comp_groups <- data.frame(seq(from = 1, to = as.numeric(input$comp_number)))
  colnames(comp_groups) <- "CompNumber"
  color_choices <- distinctColorPalette(as.numeric(input$comp_number)*2)
  color_choices_N <- color_choices[c(TRUE, FALSE)]
  color_choices_D <- color_choices[c(FALSE, TRUE)]
  cat(file = stderr(), "set_stat_groups....2", "\n")
  
  for (i in 1:input$comp_number) {
    cat(file = stderr(), str_c("begin stat_comp ", i, " of ", input$comp_number), "\n")
    stats_data_N <- stats_data
    stats_data_D <- stats_data
    
    #REBUILDING COMPARISON SELECTOR
    # using contains instead of matches, trying to find all inclusive sample group with not SPQC type
    #for(stats_group in input[[str_c('comp_',i,'N')]]){stats_data_N <- stats_data_N %>% dplyr::select(contains(stats_group))   }
    #for(stats_group in input[[str_c('comp_',i,'D')]]){stats_data_D <- stats_data_D %>% dplyr::select(contains(stats_group))   }
 
    # adding section to reduce complexity of grouping, Min.Samples will not allow other factors into the comparision
    factorsN <- input[[str_c('comp_',i,'N')]]
    factorsD <- input[[str_c('comp_',i,'D')]]
    
    #trouble shooting
    #dfD <<- stats_data_D
    #dfN <<- stats_data_N
    #fN <<- factorsN
    #fD <<- factorsD
    # stats_data_D <- dfD
    # stats_data_N <- dfN
    # factorsN <- fN
    # factorsD <- fD
    # rm(stats_data_D, stats_data_N, factorsN, factorsD, n_max, n_min, d_max, d_min)
    
    # find sample position numbers
    samples_N <- which(dpmsr_set$design$Group %in% as.list(factorsN))
    samples_D <- which(dpmsr_set$design$Group %in% as.list(factorsD))
    cat(file = stderr(), str_c("samples_N = ", samples_N), "\n")
    cat(file = stderr(), str_c("samples_D = ", samples_D), "\n")
    
    
    # reduce dataframe to samples of interest
    stats_data_N <- stats_data_N[c(samples_N)]
    stats_data_D <- stats_data_D[c(samples_D)]
    cat(file = stderr(), str_c("stats_data_N = ", ncol(stats_data_N)), "\n")
    cat(file = stderr(), str_c("stats_data_D = ", ncol(stats_data_D)), "\n")
    
    # set comp group names
    comp_groups$comp_N[i] <- paste(unique(unlist(str_split(factorsN, "_"))), collapse = "_")
    comp_groups$comp_D[i] <- paste(unique(unlist(str_split(factorsD, "_"))), collapse = "_")
    comp_groups$comp_name[i] <- str_c(comp_groups$comp_N[i], "_v_", comp_groups$comp_D[i])
    if (nchar(comp_groups$comp_name[i]) > 31) {
      comp_groups$comp_name[i] <- substr(comp_groups$comp_name[i], 1, 31)
      comp_warning <- TRUE
      comp_warning_list = c(comp_warning_list,i)
      cat(file = stderr(), str_c("Comparison ", i," is over 31 characters"), "\n")
    }
    
    if (dpmsr_set$x$primary_group) {
      comp_groups$primary_comp_N[i] <- str_split(comp_groups$comp_N[i], "_")[[1]][1]
      comp_groups$primary_comp_D[i] <- str_split(comp_groups$comp_D[i], "_")[[1]][1]
    }
    

    # set column headers for comparison stat columns
    # comp_groups$fc[i] <- str_c(comp_groups$comp_name[i], "_FC")
    # comp_groups$fc2[i] <- str_c(comp_groups$comp_name[i], "_FC2")
    # comp_groups$pval[i] <- str_c(comp_groups$comp_name[i], "_pval")
    # comp_groups$limma_pval[i] <- str_c(comp_groups$comp_name[i], "_limma_pval")
    # comp_groups$exactTest[i] <- str_c(comp_groups$comp_name[i], "_exactTest")
    # comp_groups$adjpval[i] <- str_c(comp_groups$comp_name[i], "_adjpval")
    # comp_groups$cohensd[i] <- str_c(comp_groups$comp_name[i], "_cohensd")
    # comp_groups$mf[i] <- str_c(comp_groups$comp_name[i], "_MF")
    # comp_groups$cv_N[i] <- str_c(comp_groups$comp_N[i], "_CV")  
    # comp_groups$cv_D[i] <- str_c(comp_groups$comp_D[i], "_CV")  
    
    comp_groups$N_count[i] <- ncol(stats_data_N)
    comp_groups$D_count[i] <- ncol(stats_data_D)
    comp_groups$N_color[i] <- color_choices_N[i] #first i was 1
    comp_groups$D_color[i] <- color_choices_D[i] #first i was 1
    
    #update screen to display the number of samples in the comparison
    updatePickerInput(session, inputId = str_c('comp_',i,'N'), label = str_c("Numerator selected -> ", ncol(stats_data_N)) )
    updatePickerInput(session, inputId = str_c('comp_',i,'D'), label = str_c("Denominator selected -> ", ncol(stats_data_D)) )
    
    updateTextInput(session, inputId = str_c("comp",i,"_name"),  value = comp_groups$comp_name[i])
    
    #save original sample selection
    dpmsr_set$y$stats[[str_c("comp", i, "_N")]] <<- input[[str_c('comp_',i,'N')]]
    dpmsr_set$y$stats[[str_c("comp", i, "_D")]] <<- input[[str_c('comp_',i,'D')]]
    
    comp_groups$sample_numbers_N[i] <- list(samples_N)
    comp_groups$sample_numbers_D[i] <- list(samples_D)
    
    cat(file = stderr(), str_c("end stat_comp ", i, " of ", input$comp_number), "\n")
  }
  
  cat(file = stderr(), "set_stats_groups....3", "\n")
  dpmsr_set$y$stats$groups <<- comp_groups
  dpmsr_set$y$stats$comp_number <<- input$comp_number
  
# update UI widgets after stat comparisons are selected
  updateTextInput(session, "stats_barplot_title", value = str_c("Barplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  updateTextInput(session, "stats_boxplot_title", value = str_c("Boxplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  update_comparisons(session, input, output)
  
  if (comp_warning) {
    comp_warning_list <- paste(unlist(comp_warning_list), collapse = ",")
    shinyalert("Oops!", str_c("Comparison(", unlist(comp_warning_list), ") too long, max=31"), type = "error")
    }
  
  cat(file = stderr(), "set_stat_groups....end", "\n")
}

#--------------------------------------------------------------------------------------------------------------------------------
#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_calc2 <- function(session, input, output) {
  cat(file = stderr(), "started stat_calc2", "\n")
  cat(file = stderr(), str_c("Final data input is ", dpmsr_set$x$raw_data_input), "\n")
  cat(file = stderr(), str_c("Final data output is ", dpmsr_set$x$final_data_output), "\n")
  
  input_select_final_data_stats <- input$select_final_data_stats
  input_peptide_missing_filter <- input$peptide_missing_filter
  input_peptide_cv_filter <- input$peptide_cv_filter
  input_peptide_missing_factor <- input$peptide_missing_factor
  input_peptide_cv_factor <- input$peptide_cv_factor
  input_checkbox_add_gene_column <- input$checkbox_add_gene_column
  
  #get stat parameters from inputs
  set_stat_parameters(session, input, output)
  
  # Primary Loop
  # Testing - data_in <- dpmsr_set$data$impute$impute;  i=1; input_select_final_data_stats="impute"; input_peptide_missing_filter <- TRUE; input_peptide_cv_filter <- FALSE; input_peptide_missing_factor <- 0.6; input_peptide_cv_factor <- 100; input_checkbox_add_gene_column <- FALSE
  
  for (i in 1:nrow(dpmsr_set$y$stats$groups)) 
  {
    cat(file = stderr(), str_c("Calculating stats...", i, " of ", nrow(dpmsr_set$y$stats$groups)), "\n")
    
    #transfer new stat group name incase customized
    dpmsr_set$y$stats$groups$comp_name[i] <<- dpmsr_set$y$stats[[str_c('comp',i,'_name')]]
    updateTextInput(session, inputId = str_c("volcano",i,"_stats_plot_title"),  value = str_c("Volcano ", dpmsr_set$y$stats$groups$comp_name[i]))
    
    #assign headers, group name may have changed since setting stats
    dpmsr_set$y$stats$groups <<- set_stat_headers(dpmsr_set$y$stats$groups, i)
  
    #get data for stats
    cat(file = stderr(), "stat_calc2...1", "\n")
    if (input_select_final_data_stats == "directlfq_full") {
      cat(file = stderr(), "stat_calc2...1.1", "\n")
      data_in <- dpmsr_set$data$directlfq$directlfq_protein_impute
    # }else if ((dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") && dpmsr_set$x$final_data_output == "Peptide") {
    #   data_in <- dpmsr_set$data$final[[input$select_final_data_stats]]
    #   data_in <- data_in %>% dplyr::select(!ends_with("_CV"))
    } else {
      cat(file = stderr(), "stat_calc2...1.2", "\n")
      data_in <- dpmsr_set$data$impute[[input_select_final_data_stats]]
      #data_in <- dpmsr_set$data$impute$sltmm
    }
    
    info_columns <- ncol(data_in) - dpmsr_set$y$sample_number
    annotate_in <- data_in[1:info_columns]
    data_in <- data_in[(info_columns + 1):ncol(data_in)]
    
    # get impute data for stats
    cat(file = stderr(), "stat_calc2...2", "\n")
    if (dpmsr_set$x$raw_data_input == "Protein" || input_select_final_data_stats == "directlfq_full") {
      cat(file = stderr(), "stat_calc2...2.1", "\n")
      imputed_info_columns <- ncol(dpmsr_set$data$protein_imputed_df) - dpmsr_set$y$sample_number
      imputed_df <- dpmsr_set$data$protein_imputed_df[(imputed_info_columns + 1):ncol(dpmsr_set$data$protein_imputed_df)]  
    }else{
      cat(file = stderr(), "stat_calc2...2.2", "\n")
      imputed_info_columns <- ncol(dpmsr_set$data$peptide_imputed_df) - dpmsr_set$y$sample_number
      imputed_df <- dpmsr_set$data$peptide_imputed_df[(imputed_info_columns + 1):ncol(dpmsr_set$data$peptide_imputed_df)]  
    }
    
    #isolate data for comparison (data and imputed)  
    cat(file = stderr(), "stat_calc2...3", "\n")
    comp_N_data <- data_in[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_data <- data_in[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    spqc_data <- data_in[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    spqc_impute_data <- imputed_df[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    comp_N_imputed <- imputed_df[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_imputed <- imputed_df[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    comp_N_imputed <- mutate_all(comp_N_imputed, function(x) as.numeric(x))
    comp_D_imputed <- mutate_all(comp_D_imputed, function(x) as.numeric(x))
    spqc_impute_data <- mutate_all(spqc_impute_data, function(x) as.numeric(x))
    sample_count <- ncol(comp_N_data) + ncol(comp_D_data) + ncol(spqc_data)
    sample_N_count <- ncol(comp_N_data)
    sample_D_count <- ncol(comp_D_data)
    sample_SPQC_count <- ncol(spqc_data)
    
    #---------------------------------------------------------------------------------------------------------------    
    #if data is peptide/peptide or protein/protein will need to created imputed column now
    
    #reduce peptide PD imputed column to only samples of comparison
    cat(file = stderr(), "stat_calc2...4", "\n")
    if (dpmsr_set$x$final_data_output == "Peptide" && (dpmsr_set$x$raw_data_input != "Precursor" && dpmsr_set$x$raw_data_input != "Precursor_PTM" )) {
      cat(file = stderr(), "stat_calc2...4.1", "\n")
      df_impute_peptide <- cbind(comp_N_imputed, comp_D_imputed, spqc_impute_data)
      df_impute_peptide[df_impute_peptide == 0] <- "-"
      while (ncol(df_impute_peptide) > 1) {
        df_impute_peptide[,1] <- str_c(df_impute_peptide[,1], ".", df_impute_peptide[,2])
        df_impute_peptide[,2] <- NULL
      }
      colnames(df_impute_peptide) <- "Detected_Imputed"
      annotate_in$Detected_Imputed <- df_impute_peptide$Detected_Imputed
      df <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data)
    }
    
    #reduce protein PD imputed column to only samples of comparison
    cat(file = stderr(), "stat_calc2...5", "\n")
    if (dpmsr_set$x$raw_data_input == "Protein" || input_select_final_data_stats == "directlfq_full") {
      cat(file = stderr(), "stat_calc2...5.1", "\n")
      df_impute_protein <- cbind(comp_N_imputed, comp_D_imputed, spqc_impute_data)
      df_impute_protein[df_impute_protein == 0] <- "-"
      while (ncol(df_impute_protein) > 1) {
        df_impute_protein[,1] <- str_c(df_impute_protein[,1], ".", df_impute_protein[,2])
        df_impute_protein[,2] <- NULL
      }
      colnames(df_impute_protein) <- "Detected_Proteins"
      annotate_in$Detected_Proteins <- df_impute_protein$Detected_Proteins
      df <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data)
      df$Protein_Imputed <- NULL
      annotate_in$Protein_Imputed <- NULL
      dpmsr_set$y$info_columns <<- ncol(annotate_in)
      dpmsr_set$y$info_columns_final <<- ncol(annotate_in)
    }
    
    # if(input$select_final_data_stats=="directlfq_full"){
    #   df_impute_protein <- cbind(comp_N_imputed, comp_D_imputed, spqc_impute_data)
    #   df_impute_protein[df_impute_protein==0] <- "-"
    #   while (ncol(df_impute_protein)>1) {
    #     df_impute_protein[,1] <- str_c(df_impute_protein[,1], ".", df_impute_protein[,2])
    #     df_impute_protein[,2] <- NULL
    #   }
    #   colnames(df_impute_protein) <- "Detected_Proteins"
    
    #   annotate_in$Detected_Proteins <- df_impute_protein$Detected_Proteins
    #   df <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data)
    #   df$Protein_Imputed <- NULL
    #   annotate_in$Protein_Imputed <- NULL
    #   dpmsr_set$y$info_columns <<- ncol(annotate_in)
    #   dpmsr_set$y$info_columns_final <<- ncol(annotate_in)
    # }
    
    #------------------------------------------------------------------------------------------------------------------
    
    # if TMT SPQC norm will skip the peptide to protein section, set df now
    if (as.logical(dpmsr_set$x$tmt_spqc_norm)) {
      cat(file = stderr(), "tmt_spqc_norm is set", "\n")
      df <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data)
    }
    
    
    #--------------------------------------------------------------------------------------------------------------
    #section for peptide to protein final data, and precursor to peptide
    cat(file = stderr(), "stat_calc2...6 (peptide to protein - or - precursor to peptide)", "\n")
    if ((dpmsr_set$x$raw_data_input != "Protein" &  dpmsr_set$x$final_data_output == "Protein" & (!as.logical(dpmsr_set$x$tmt_spqc_norm)) & 
      input_select_final_data_stats != "directlfq_full"  ) || ((dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") &&
                                                               dpmsr_set$x$final_data_output == "Peptide")) {
      cat(file = stderr(), "stat_calc2...6.1", "\n")
      mf <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
      N_CV <- percentCV_gw(comp_N_data)
      D_CV <- percentCV_gw(comp_D_data)
      minCV <- pmin(N_CV, D_CV)
      
      df <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data, mf, N_CV, D_CV, minCV)
      df_impute <- cbind(annotate_in, comp_N_imputed, comp_D_imputed, mf, N_CV, D_CV, minCV)
      df_impute_summary <- cbind(annotate_in, comp_N_imputed, comp_D_imputed, spqc_impute_data, mf, N_CV, D_CV, minCV)
    
      if (input_peptide_missing_filter) {
        df <- df[(df$mf >= input_peptide_missing_factor),]
        df_impute <- df_impute[(df_impute$mf >= input_peptide_missing_factor),]  
        df_impute_summary <- df_impute_summary[(df_impute_summary$mf >= input_peptide_missing_factor),] 
        #df <- df[(df$mf >= 0.6),]; df_impute <- df_impute[(df_impute$mf >= 0.6),]; df_impute_summary <- df_impute_summary[(df_impute_summary$mf >= 0.6),]  
      }
  
      if (input_peptide_cv_filter) {
        df <- df[(df$minCV <= input$peptide_cv_factor),]
        df_impute <- df_impute[(df_impute$minCV <= input$peptide_cv_factor),]  
        df_impute_summary <- df_impute_summary[(df_impute_summary$minCV <= input$peptide_cv_factor),]  
      }
    
      cat(file = stderr(), "stat_calc2...7", "\n")
      df$mf <- NULL
      df_impute$mf <- NULL
      df_impute_summary$mf <- NULL
      df$N_CV <- NULL
      df_impute$N_CV <- NULL
      df_impute_summary$N_CV <- NULL
      df$D_CV <- NULL
      df_impute$D_CV <- NULL
      df$D_CV <- NULL
      df_impute_summary$D_CV <- NULL
      df$minCV <- NULL
      df_impute$minCV <- NULL
      df_impute_summary$minCV <- NULL
      
      cat(file = stderr(), "stat_calc2...8", "\n")
      df_impute_summary <- df_impute_summary[(info_columns + 1):ncol(df_impute_summary)] %>% mutate_all(as.character)
      while (ncol(df_impute_summary) > 1) {
        df_impute_summary[,1] <- str_c(df_impute_summary[,1], ".", df_impute_summary[,2])
        df_impute_summary[,2] <- NULL
      }
      
      colnames(df_impute_summary) <- "Detected_Imputed"
      df$Detected_Imputed <- df_impute_summary$Detected_Imputed
      dpmsr_set$data$stats$peptide[[dpmsr_set$y$stats$groups$comp_name[i]]] <<- df


      #testing block
      #test_df <<- df; test_df_impute <<- df_impute; test_info_columns <<- info_columns 
      # df <- test_df; df_impute <- test_df_impute ;info_columns <- test_info_columns
      
      #collapse peptide data and imputed peptide info, this is done with only the stat groups
      cat(file = stderr(), "stat_calc2...9", "\n")
      if (dpmsr_set$x$final_data_output == "Peptide" && (dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM")) {
        cat(file = stderr(), "stat_calc2...9.1", "\n")
        df <- collapse_precursor(df, info_columns, stats = TRUE)
        df_impute <- collapse_precursor(df_impute, info_columns, stats = TRUE)
      }else{
        cat(file = stderr(), "stat_calc2...9.2", "\n")
        #test_df <<- df; test_df_impute <<- df_impute; test_info_colums <<- info_columns
        df <- collapse_peptide(df, info_columns, stats = TRUE)
        df_impute <- collapse_peptide(df_impute, info_columns, stats = TRUE, impute = TRUE)
      }
      
      
      cat(file = stderr(), "stat_calc2...10", "\n")
      #create string column with imputed peptide info
      info_columns <- ncol(df) - sample_count
      test2 <- df_impute[(info_columns + 1):ncol(df_impute)]
      test2[test2 == 0] <- "-"
      test2 <- test2 %>% mutate_all(as.character)
      
      while (ncol(test2) > 1) {
        test2[,1] <- str_c(test2[,1], ".", test2[,2])
        test2[,2] <- NULL
      }
      
      cat(file = stderr(), "stat_calc2...11", "\n")
      if (dpmsr_set$x$raw_data_input == "Precursor_PTM") {
        cat(file = stderr(), "stat_calc2...11.1", "\n")
        df <- add_column(df, test2[,1] , .after = "PTMLocations")
        names(df)[grep("PTMLocations", colnames(df)) + 1] <- "Detected_Imputed"
        imputed_df <- df_impute[(ncol(df_impute) - sample_N_count - sample_D_count + 1):ncol(df_impute)]  
        imputed_df[imputed_df > 1] <- 1
      }else{
        cat(file = stderr(), "stat_calc2...11.2", "\n")
        df <- add_column(df, test2[,1] , .after = "Peptides")
        names(df)[grep("Peptides", colnames(df)) + 1] <- "Detected_Imputed"
        imputed_df <- df_impute[(ncol(df_impute) - sample_N_count - sample_D_count + 1):ncol(df_impute)]   
        imputed_df[imputed_df > 1] <- 1
      }


      cat(file = stderr(), "stat_calc2...12", "\n")
      comp_N_imputed <- imputed_df[1:dpmsr_set$y$stats$groups$N_count[i]]
      comp_D_imputed <- imputed_df[(dpmsr_set$y$stats$groups$N_count[i] + 1):(dpmsr_set$y$stats$groups$N_count[i] + dpmsr_set$y$stats$groups$D_count[i])]
    
      annotate_in <- df[1:(ncol(df) - sample_count)]
      cat(file = stderr(), "stat_calc2...13", "\n")
      
      if (input_checkbox_add_gene_column) {
        add_gene <- str_extract(annotate_in$Description, "GN=\\w*")
        add_gene <- gsub("GN=", "", add_gene)
        annotate_in <- annotate_in %>% add_column(Gene = add_gene, .before = 3)
        dpmsr_set$data$final[[input_select_final_data_stats]] <<- dpmsr_set$data$final[[input_select_final_data_stats]] %>% 
          add_column(Gene = annotate_in$Gene, .before = 3)
      }
      
      cat(file = stderr(), "stat_calc2...14", "\n")
      data_in <- df[(ncol(df)-sample_count+1):(ncol(df))]
      comp_N_data <- data_in[1:dpmsr_set$y$stats$groups$N_count[i]]
      comp_D_data <- data_in[(dpmsr_set$y$stats$groups$N_count[i] + 1):(dpmsr_set$y$stats$groups$N_count[i] + dpmsr_set$y$stats$groups$D_count[i])]
      spqc_data <- data_in[(dpmsr_set$y$stats$groups$N_count[i] + dpmsr_set$y$stats$groups$D_count[i] + 1):ncol(data_in)]
  }
    #end of peptide to protein final data
    #-------------------------------------------------------------------------------------------------------------------------

    #test block
    test_imputed_df <<- imputed_df; test_df <<- df; test_data_in <<- data_in; test_comp_N_data <<- comp_N_data; test_comp_D_data <<- comp_D_data; test_spqc_data <<- spqc_data; test_comp_N_imputed <<- comp_N_imputed; test_comp_D_imputed <<- comp_D_imputed; test_annotate_in <<- annotate_in
    # imputed_df <- test_imputed_df; comp_N_data <- test_comp_N_data; comp_D_data <- test_comp_D_data; comp_N_imputed <- test_comp_N_imputed; comp_D_imputed <- test_comp_D_imputed; spqc_data <- test_spqc_data; df <- test_df; annotate_in <- test_annotate_in; i=1
    
    
    #start df for stats -----------------------------------------
    cat(file = stderr(), "stat_calc2...15", "\n")
    stat_df <- annotate_in[1:1]
    #stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_N[i], "_CV")] <- percentCV_gw(comp_N_data)
    #stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_D[i], "_CV")] <- percentCV_gw(comp_D_data)
    stat_df[ , dpmsr_set$y$stats$groups$cv_N[i]] <- percentCV_gw(comp_N_data)
    stat_df[ , dpmsr_set$y$stats$groups$cv_D[i]] <- percentCV_gw(comp_D_data)
    
    stat_df[ , str_c(dpmsr_set$y$stats$comp_spqc, "_CV")] <- percentCV_gw(spqc_data)
    stat_df[ ,dpmsr_set$y$stats$groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
    
    if (input$checkbox_adjpval) {
      stat_df[ ,dpmsr_set$y$stats$groups$adjpval[i]] <- p.adjust(stat_df[ ,dpmsr_set$y$stats$groups$pval[i]], method = input$padjust_options) 
    }
    
    if (input$checkbox_cohensd) {
        stat_df[ ,dpmsr_set$y$stats$groups$cohensd[i]] <- cohend_gw(comp_N_data, comp_D_data, as.logical(input$checkbox_cohensd))
    }
    
    cat(file = stderr(), "stat_calc2...16", "\n")
    if (input$checkbox_limmapvalue) {
      stat_df[ ,dpmsr_set$y$stats$groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data)
    }
    #stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
    #stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
    
    #only include if not TMT SPQC norm
    if (!as.logical(dpmsr_set$x$tmt_spqc_norm)) {
      stat_df[ ,dpmsr_set$y$stats$groups$mf[i]] <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
    }
    
    cat(file = stderr(), "stat_calc2...17", "\n")
    # Create final tables--------------------------------------------------
    if (input_select_final_data_stats == "impute") {
      cat(file = stderr(), "stat_calc2...17.1", "\n")
      colnames(comp_N_data) <- dpmsr_set$design$Header3[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i])]
      colnames(comp_D_data) <- dpmsr_set$design$Header3[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i])]
      colnames(spqc_data) <- dpmsr_set$design$Header3[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    }else{
      cat(file = stderr(), "stat_calc2...17.2", "\n")
      colnames(comp_N_data) <- dpmsr_set$design$Header2[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i])]
      colnames(comp_D_data) <- dpmsr_set$design$Header2[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i])]
      colnames(spqc_data) <- dpmsr_set$design$Header2[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    }
    
    cat(file = stderr(), "stat_calc2...18", "\n")
    data_table <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data, stat_df[2:ncol(stat_df)])
    
    if (input$select_final_data_stats == "directlfq_full") {
      cat(file = stderr(), "remove Detected_Imputed column from directlfq_full", "\n")
      data_table <- data_table[, !names(data_table) %in% c("Detected_Imputed")]
    }
    
    data_table$Stats <- ""
    data_table_cols <- ncol(data_table)
    
    if (!input$checkbox_filter_adjpval) {
      data_table$pval <- ifelse(data_table[[dpmsr_set$y$stats$groups$pval[i]]]  <= input$pvalue_cutoff, 0, 1)
    }else{
      data_table$pval <- ifelse(data_table[[dpmsr_set$y$stats$groups$adjpval[i]]]  <= input$pvalue_cutoff, 0, 1)
    }
    
    if (!as.logical(dpmsr_set$x$tmt_spqc_norm)) {
      data_table$mf <- ifelse(data_table[[dpmsr_set$y$stats$groups$mf[i]]]  >= input$missing_factor, 0, 1)
    }
    
    if (input$stats_spqc_cv_filter) {
      data_table$spqc <- ifelse(data_table[[str_c(dpmsr_set$y$stats$comp_spqc, "_CV")]] <= input$stats_spqc_cv_filter_factor, 0, 1)
    }
    
    if (input$stats_comp_cv_filter) {
      data_table$cv <- ifelse(data_table[[dpmsr_set$y$stats$groups$cv_N[i]]] <= input$stats_comp_cv_filter_factor |
                                data_table[[dpmsr_set$y$stats$groups$cv_D[i]]] <= input$stats_comp_cv_filter_factor, 0, 1)
    }
    if (input$stats_peptide_minimum) {
      data_table$pepmin <- ifelse(data_table$Peptides >= input$stats_peptide_minimum_factor, 0, 1)
    }
    
    cat(file = stderr(), "stat_calc2...19", "\n")
    
    data_table$sum <- rowSums(data_table[(data_table_cols + 1):ncol(data_table)])
    
    data_table$Stats <- ifelse(data_table$sum == 0 & data_table[[dpmsr_set$y$stats$groups$fc[i]]] >= input$foldchange_cutoff, "Up", 
                         ifelse(data_table$sum == 0 & data_table[[dpmsr_set$y$stats$groups$fc[i]]] <= -input$foldchange_cutoff, "Down", ""))       
    
    test_data_table <<- data_table
    
    data_table <- data_table[1:data_table_cols]
    
    data_table <- data_table[order(data_table$Stats, -data_table[[dpmsr_set$y$stats$groups$pval[i]]], decreasing = TRUE),]
    
    #filter stats by ptm or accession if needed
    if (dpmsr_set$x$final_data_output == "Peptide" &  input$checkbox_report_ptm) {
      if (!is.null(data_table$Modifcations)) {
        data_table <- data_table[grep(dpmsr_set$x$peptide_report_grep, data_table$Modifications),]
      }else{
        data_table <- data_table[grep(dpmsr_set$x$peptide_report_grep, data_table$Sequence),]
      }
    }
    
    if (dpmsr_set$x$final_data_output == "Peptide" & input$checkbox_report_accession) {
      data_table <- subset(data_table, Accession %in% dpmsr_set$x$accession_report_list )
    }
    
    dpmsr_set$data$stats[[dpmsr_set$y$stats$groups$comp_name[i]]] <<- data_table
    
    # if (dpmsr_set$x$data_source == "SP"){
    #   dpmsr_set$data$stats[[dpmsr_set$y$stats$groups$comp_name[i]]] <<- expand_spectronaut(data_table) 
    # } 
    
    update_comparisons(session, input, output)
  } 
  # End of Primary Loop

  
  cat(file = stderr(), str_c(dpmsr_set$y$stats$groups$comp_name[i], " data has ", nrow(dpmsr_set$data$stats[[dpmsr_set$y$stats$groups$comp_name[i]]]), " rows"), "\n")
  cat(file = stderr(), "Calculating stats complete...", "\n")
  return()
}



#---------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# create final excel documents
stats_Final_Excel <- function(session, input, output) {
  cat(file = stderr(), "Creating Excel Output File...1", "\n")
  require(openxlsx)
    
    file_dir <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats) 
    filename <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats, "//", input$final_stats_name)
    
    if (!is_dir(file_dir)) {
      cat(file = stderr(), str_c("create_dir...", file_dir), "\n")
      dir_create(file_dir)
    }
    
    
    df <- dpmsr_set$data$final[[input$select_final_data_stats]]
    
    #testing
    #input_final_stats_name <<- input$final_stats_name
    #input_select_final_data_stats <<- input$select_final_data_stats
    # filename <- str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp, "//", input_final_stats_name)
    # df <- dpmsr_set$data$final[[input_select_final_data_stats]]
    
    #remove FC2 from df for excel
    #df <- df[,-grep(pattern="_FC2", colnames(df))]
    
    if (dpmsr_set$y$state == "Peptide" && dpmsr_set$x$final_data_output == "Protein") {
      df_raw <- dpmsr_set$data$finalraw
    }else if ((dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") && (dpmsr_set$x$final_data_output == "Peptide") ) {
      df_raw <- dpmsr_set$data$data_raw_precursor
      }else{
      df_raw <- dpmsr_set$data$finalraw
      #this may not be needed 
      #if PTM out need to reduce raw data frame for excel
      #if (as.logical(dpmsr_set$x$peptide_ptm_out)){
        #df_raw <- df_raw[grep(dpmsr_set$x$peptide_report_grep, df_raw$Modifications),]
      #}
    }
    
    # save the option to save normalized table with raw
    #df2 <- cbind(df_raw, df[(dpmsr_set$y$info_columns_final+1):ncol(df)])
    cat(file = stderr(), "Creating Excel Output File...2", "\n")
    df2 <- df
    
    
    # subsets data for accession specific reports
    if (as.logical(dpmsr_set$x$accession_report_out)) {
      df <- subset(df, Accession %in% dpmsr_set$x$accession_report_list )
      df2 <- subset(df2, Accession %in% dpmsr_set$x$accession_report_list )
    }
    
    nextsheet <- 1
    
    wb <- createWorkbook()
    
    cat(file = stderr(), "Creating Excel Output File...3", "\n")
    #----------------------------------------
    
    if (site_user == "dpmsr" && dpmsr_set$x$raw_data_input != "Protein" &&
        dpmsr_set$x$raw_data_input != "Precursor" && dpmsr_set$x$raw_data_input != "Precursor_PTM") {
      cat(file = stderr(), "Creating Excel Output File...3.1", "\n")
      
      if (as.logical(dpmsr_set$x$peptide_isoform) ) {
          raw_peptide <- dpmsr_set$data$data_peptide_isoform_start
        }else{
          raw_peptide <- dpmsr_set$data$data_peptide_start
        }
  
      peptide_impute <- dpmsr_set$data$impute$impute
        
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Raw Peptide Data")
      writeData(wb, sheet = nextsheet, raw_peptide)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Imputed Peptide Data")
      writeData(wb, sheet = nextsheet, peptide_impute) 
      nextsheet <- nextsheet + 1
    }
    
    
    
    if (site_user == "dpmsr" && dpmsr_set$x$raw_data_input == "Precursor" && dpmsr_set$x$final_data_output == "Protein") {
      cat(file = stderr(), "Creating Excel Output File...3.2", "\n")
      
      precursor_impute <- dpmsr_set$data$impute$impute
      
      raw_peptide <- collapse_precursor(dpmsr_set$data$data_precursor_start, info_columns = 0, stats = FALSE)
      
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Raw Precursor Data")
      writeData(wb, sheet = nextsheet, df_raw)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Raw Peptide Data")
      writeData(wb, sheet = nextsheet, raw_peptide)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Imputed Precursor Data")
      writeData(wb, sheet = nextsheet, dpmsr_set$data$impute$impute) 
      nextsheet <- nextsheet + 1
    }
    
    
    
    if (site_user == "dpmsr" && dpmsr_set$x$final_data_output == "Peptide" &&
        (dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM")) {
      cat(file = stderr(), "Creating Excel Output File...3.3", "\n")
      
      raw_precursor <- dpmsr_set$data$data_raw_precursor
      precursor_impute <- dpmsr_set$data$impute$impute
      
      
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Raw Precursor Data")
      writeData(wb, sheet = nextsheet, raw_precursor)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Imputed Precursor Data")
      writeData(wb, sheet = nextsheet, precursor_impute) 
      nextsheet <- nextsheet + 1
    }  
    
    
    if (site_user == "dpmsr" && dpmsr_set$x$raw_data_input == "Protein") {
      cat(file = stderr(), "Creating Excel Output File...3.4", "\n")
      raw_protein <- dpmsr_set$data$data_protein_start
      protein_impute <- dpmsr_set$data$impute$impute
      
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Raw Protein Data")
      writeData(wb, sheet = nextsheet, raw_protein)
      nextsheet <- nextsheet + 1
      addWorksheet(wb, "Imputed Protein Data")
      writeData(wb, sheet = nextsheet, protein_impute) 
      nextsheet <- nextsheet + 1
    }
    #----------------------------------------
    
    cat(file = stderr(), "Creating Excel Output File...4", "\n")
    
    #addWorksheet(wb, deparse(substitute(df2)))
    addWorksheet(wb, "Normalized Data")
    writeData(wb, sheet = nextsheet, df2)  
    nextsheet <- nextsheet + 1
    z = 1
    for (i in 1:as.numeric(dpmsr_set$y$stats$comp_number))  {
      comp_string <- dpmsr_set$y$stats$groups$comp_name[i]
      filtered_df <- dpmsr_set$data$stats[[comp_string]]
      cat(file = stderr(), str_c("Saving excel...", comp_string), "\n")
      addWorksheet(wb, comp_string)
      writeData(wb, sheet = nextsheet, filtered_df)
      nextsheet <- nextsheet + 1
      z <- z + 2
    }
    cat(file = stderr(), "writting excel to disk...", "\n")
    saveWorkbook(wb, filename, overwrite = TRUE)
    cat(file = stderr(), str_c("Creating Excel Output File...", filename), "\n")
}

# data table filter ------------------------------------------------------

stats_data_table_filter <- function(session, input, output) {
  cat(file = stderr(), "stats_data_table_filter... start", "\n")
  
  comp_string <- input$stats_select_data_comp
  comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
  sample_number <- dpmsr_set$y$stats$groups$N_count[comp_number] + dpmsr_set$y$stats$groups$D_count[comp_number]
  
  cat(file = stderr(), "stats_data_table_filter... 1", "\n")
  df_filter <- dpmsr_set$data$stats[[comp_string]]
  
  if (is.grouped_df(df_filter)) {df_filter <- ungroup(df_filter)}
  
  #test_df_filter <<- df_filter
  
  
  
  cat(file = stderr(), "stats_data_table_filter... 2", "\n")
  if (input$stats_add_filters) {
    df_filter <- df_filter %>% filter(df_filter$Stats == "Up" | df_filter$Stats == "Down")
  }
  
  cat(file = stderr(), "stats_data_table_filter... 3", "\n")
  
  if (input$stats_data_topn != 0 ) {
    df_filter$sum <- rowSums(df_filter[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + sample_number)])
    df_filter <- df_filter[order(-df_filter$sum),]                      
    df_filter <- df_filter[1:input$stats_data_topn,]
    df_filter$sum <- NULL
  }
  
  cat(file = stderr(), "stats_data_table_filter... 4", "\n")
  
  if (input$stats_data_accession != "0" ) {
    #df_filter <-subset(df_filter, Accession %in% as.character(input$stats_data_accession)  )
    #df_filter <- filter(df_filter, grepl(as.character(input$stats_data_accession)), Accession)
    df_filter <- df_filter[grep(as.character(input$stats_data_accession), df_filter$Accession), ]
  }
  
  if (input$stats_data_description != "0") {
    df_filter <- df_filter[grep(as.character(input$stats_data_description), df_filter$Description), ]
  }
  cat(file = stderr(), "stats_data_table_filter... end", "\n")
  return(df_filter)
}



# peptide zscore ------------------------------------------------------

peptide_zscore <- function(df_peptide, info_columns) {
  info_df <- df_peptide[1:info_columns]
  df <- df_peptide[(info_columns + 1):ncol(df_peptide)]
  df <- log(df, 2)
  z_mean <- rowMeans(df)
  z_stdev  <- apply(df[1:ncol(df)], 1, sd)
  
  df <- apply(df, MARGIN = 2, function(x) (x - z_mean)/z_stdev  )
  df <- data.frame(cbind(info_df, df), stringsAsFactors = FALSE)
  
  return(df)
  
}

# TMT IRS data ------------------------------------------------------

TMT_IRS_protein_data <- function(session, input, output) {
  
  cat(file = stderr(), "TMT IRS Protein stats and plots..." , "\n")
  
  comp_string <- input$stats_oneprotein_plot_comp
  comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string)
  
  df <- dpmsr_set$data$stats[[comp_string]]
  df <- subset(df, Accession %in% as.character(input$stats_oneprotein_accession)  )
  sample_number <- dpmsr_set$y$stats$groups$N_count[comp_number] + dpmsr_set$y$stats$groups$D_count[comp_number]
  
  #add spqc to plots
  if (input$stats_oneprotein_plot_spqc) {
    df <- df[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
  }else{
    df <- df[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + sample_number)]
  }
  
  comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
  #add spqc to plots
  if (input$stats_oneprotein_plot_spqc) {
    comp_rows <- c(comp_rows, dpmsr_set$y$stats$comp_spqc_sample_numbers)
  }
  
  comp_rows <- unlist(comp_rows)
  namex <- dpmsr_set$design$Label[comp_rows]
  color_list <- dpmsr_set$design$colorlist[comp_rows]
  groupx <- dpmsr_set$design$Group[comp_rows]

  
  return(list("df" = df, "namex" = namex, "color_list" = color_list, "comp_string" = comp_string))
}

# one protein data ------------------------------------------------------

oneprotein_data <- function(session, input, output) {
  
  cat(file = stderr(), "One Protein stats and plots..." , "\n")

  comp_string <- input$stats_oneprotein_plot_comp
  comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string)
  cat(file = stderr(), str_c("comp_string=", comp_string, "   comp_number=", comp_number), "\n")
  
  cat(file = stderr(), "One Protein stats and plots...1" , "\n")
  df <- dpmsr_set$data$stats[[comp_string]]

  cat(file = stderr(), str_c("accession=", input$stats_oneprotein_accession), "\n")
  df <- df[grep(as.character(input$stats_oneprotein_accession), df$Accession), ]
  #df <-subset(df, df$Accession %in% as.character(input$stats_oneprotein_accession)  )
  
  sample_number <- dpmsr_set$y$stats$groups$N_count[comp_number] + dpmsr_set$y$stats$groups$D_count[comp_number]
  
  cat(file = stderr(), "One Protein stats and plots...2" , "\n")
  
  if (dpmsr_set$x$final_data_output == "Peptide") {
    df_peptide <- dpmsr_set$data$stats[[comp_string]]
    df_peptide <- df_peptide[1:(dpmsr_set$y$info_columns_final + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
  }else{
    df_peptide <- dpmsr_set$data$stats$peptide[[comp_string]]
  }
  
  df_peptide <- subset(df_peptide, Accession %in% as.character(input$stats_oneprotein_accession)  )
  peptide_info_columns <- ncol(df_peptide) - sample_number - dpmsr_set$y$stats$comp_spqc_number
  
  #add spqc to plots
  cat(file = stderr(), "One Protein stats and plots...3" , "\n")
  if (input$stats_oneprotein_plot_spqc) {
    df <- df[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
  }else{
    df <- df[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + sample_number)]
  }
  
  comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
  #add spqc to plots
  if (input$stats_oneprotein_plot_spqc) {
    comp_rows <- c(comp_rows, dpmsr_set$y$stats$comp_spqc_sample_numbers)
  }
  
  #test_df_peptide <<- df_peptide
  
  cat(file = stderr(), "One Protein stats and plots...4" , "\n")
  comp_rows <- unlist(comp_rows)
  namex <- dpmsr_set$design$Label[comp_rows]
  color_list <- dpmsr_set$design$colorlist[comp_rows]
  groupx <- dpmsr_set$design$Group[comp_rows]
  
  cat(file = stderr(), "One Protein stats and plots...5" , "\n")
  if (dpmsr_set$x$final_data_output != "Peptide") {
    colnames(df_peptide)[(peptide_info_columns + 1):ncol(df_peptide)] <- 
      c(dpmsr_set$design$Header1[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number])],
        dpmsr_set$design$Header1[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[comp_number])],
        dpmsr_set$design$Header1[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)] )
  } 
  
  #sort peptides by intensity, keep highest abundant peptide
  cat(file = stderr(), "One Protein stats and plots...6" , "\n")
  df_peptide$sum <- rowSums(df_peptide[(peptide_info_columns + 1):(peptide_info_columns + sample_number + 1)])
  df_peptide <- df_peptide[order(df_peptide$sum), ]
  df_peptide$sum <- NULL
  df_peptide <- df_peptide[!duplicated(df_peptide$Sequence),]
  
  #test_df_peptide <<- df_peptide
  
  if (input$stats_use_zscore) {df_peptide <- peptide_zscore(df_peptide, peptide_info_columns)}
  
  cat(file = stderr(), "One Protein stats and plots...end" , "\n")
  
  df_list <- list("df" = df, "df_peptide" = df_peptide, "namex" = namex, "color_list" = color_list, "comp_string" = comp_string, "peptide_info_columns" = peptide_info_columns)
  #test_df_list_one_protein <<- df_list
  return(df_list)
}

# one peptide data ------------------------------------------------------

onepeptide_data <- function(session, input, output) {
    
    cat(file = stderr(), "One Peptide stats and plots..." , "\n")
    
    comp_string <- input$stats_onepeptide_plot_comp
    comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string)
    
    cat(file = stderr(), "One Peptide stats and plots...1" , "\n")
    df <- dpmsr_set$data$stats[[comp_string]]
    df <- df[grep(as.character(input$stats_onepeptide_accession), df$Accession), ]
    sample_number <- dpmsr_set$y$stats$groups$N_count[comp_number] + dpmsr_set$y$stats$groups$D_count[comp_number]
    df_peptide <- df
    info_columns <- dpmsr_set$y$info_columns_final
    df_peptide_stats <- df_peptide[(info_columns + sample_number + dpmsr_set$y$stats$comp_spqc_number + 1):ncol(df_peptide)]
    
    #continue filtering of df, used for first plot of a single peptide
    cat(file = stderr(), "One Peptide stats and plots...2" , "\n")
    df <- subset(df, Sequence %in% as.character(input$stats_onepeptide_sequence)  )
    grep_mod <- stringr::str_replace_all(input$stats_onepeptide_modification, "\\[", "\\\\[")
    grep_mod <- stringr::str_replace_all(grep_mod, "\\]", "\\\\]")
    grep_mod <- stringr::str_replace_all(grep_mod,"\\(", "\\\\(")
    grep_mod <- stringr::str_replace_all(grep_mod,"\\)", "\\\\)")
    df <- df %>% filter(stringr::str_detect(Modifications, grep_mod) )
    
    #add spqc to plots
    cat(file = stderr(), "One Peptide stats and plots...3" , "\n")
    if (input$stats_onepeptide_plot_spqc) {
      df_peptide <- df_peptide[1:(info_columns + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
      df <- df[(info_columns + 1):(info_columns + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
      comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number], 
                     dpmsr_set$y$stats$comp_spqc_sample_numbers )
    }else{
      df_peptide <- df_peptide[1:(info_columns + sample_number)]
      df <- df[(info_columns + 1):(info_columns + sample_number)]
      comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
    }

    cat(file = stderr(), "One Peptide stats and plots...4" , "\n")
    comp_rows <- unlist(comp_rows)
    namex <- dpmsr_set$design$Label[comp_rows]
    color_list <- dpmsr_set$design$colorlist[comp_rows]
    groupx <- dpmsr_set$design$Group[comp_rows]


    if (input$stats_onepeptide_use_zscore) {df_peptide <- peptide_zscore(df_peptide, info_columns)}
    
    df_list <- list("df" = df, "df_peptide" = df_peptide, "df_peptide_stats" = df_peptide_stats, "info_columns" = info_columns, "namex" = namex, "color_list" = color_list, "comp_string" = comp_string)
    #test_df_list <<- df_list
    
    cat(file = stderr(), "One Peptide stats and plots...end" , "\n")
    return(df_list)
}
    



# stat filter ------------------------------------------------------  
    stats_filter <- function(df, comp_number) { 
      
      filtered_df <- filter(df, df[[dpmsr_set$y$stats$groups$pval[comp_number]]] <= as.numeric(dpmsr_set$x$pvalue_cutoff) &  
                              (df[[dpmsr_set$y$stats$groups$fc[comp_number]]] >= as.numeric(dpmsr_set$x$foldchange_cutoff) | 
                                 df[[dpmsr_set$y$stats$groups$fc[comp_number]]] <= -as.numeric(dpmsr_set$x$foldchange_cutoff))) 
      
      filtered_df <- filter(filtered_df, filtered_df[[dpmsr_set$y$stats$groups$mf[comp_number] ]] >= as.numeric(dpmsr_set$x$missing_factor))
      
      
      if (dpmsr_set$y$stats$stats_spqc_cv_filter) { 
        filtered_df <- filter(filtered_df, filtered_df[[str_c(dpmsr_set$y$stats$comp_spqc, "_CV")]] <= as.numeric(dpmsr_set$y$stats$stats_spqc_cv_filter_factor ) )
      }
      
      if (dpmsr_set$y$stats$stats_comp_cv_filter) {
        filtered_df <- filter(filtered_df, filtered_df[[dpmsr_set$y$stats$groups$cv_N[comp_number]]] <= as.numeric(dpmsr_set$y$stats$stats_comp_cv_filter_factor) | 
                                filtered_df[[dpmsr_set$y$stats$groups$cv_D[comp_number]]] <= as.numeric(dpmsr_set$y$stats$stats_comp_cv_filter_factor) )
      }
     
      return(filtered_df)
       
    }      
    
    
  peptide_position_lookup <- function(session, input, output, protein_accession)  {
    cat(file = stderr(), "peptide_position_lookup...", "\n")
    
    #----------------------------------------------------------------------------------------------------
    if (dpmsr_set$x$data_source == "PD") {
      # create peptide lookup table
      if (dpmsr_set$x$final_data_output == "Protein") {
        if (class(dpmsr_set$data$data_raw_peptide$Positions.in.Master.Proteins) == "NULL") {
          peptide_pos_lookup <- try(dpmsr_set$data$data_raw_peptide %>% dplyr::select(Sequence, Positions.in.Proteins))
        }else{
          peptide_pos_lookup <- try(dpmsr_set$data$data_raw_peptide %>% dplyr::select(Sequence, Positions.in.Master.Proteins))
        }
      }else{
        if (dpmsr_set$x$peptide_isoform) {
          peptide_pos_lookup <- try(dpmsr_set$data$data_raw_isoform %>% dplyr::select(Sequence, Positions.in.Proteins))
        }else{
          peptide_pos_lookup <- try(dpmsr_set$data$data_raw_peptide %>% dplyr::select(Sequence, Positions.in.Proteins))
        }
      }
      
      cat(file = stderr(), "peptide_position_lookup...1", "\n")
      colnames(peptide_pos_lookup) <- c("Sequence", "Position")
      peptide_pos_lookup$Position <- gsub("\\]; \\[", "xxx",  peptide_pos_lookup$Position)
      s <- strsplit(peptide_pos_lookup$Position, split = "; ")
      peptide_pos_lookup <- data.frame(Sequence = rep(peptide_pos_lookup$Sequence, sapply(s, length)), Position = unlist(s))
      peptide_pos_lookup$Count <- str_count(peptide_pos_lookup$Position, "xxx")
      peptide_pos_lookup <- peptide_pos_lookup %>% separate(Position, c("Accession", "Start_Stop"), sep = " ")
      peptide_pos_lookup$Start_Stop <- gsub("\\[", "",  peptide_pos_lookup$Start_Stop)
      peptide_pos_lookup$Start_Stop <- gsub("\\]", "",  peptide_pos_lookup$Start_Stop)
      peptide_pos_lookup$Start_Stop <- gsub("xxx", ", ",  peptide_pos_lookup$Start_Stop)
      peptide_pos_lookup$AS <- str_c(peptide_pos_lookup$Accession, " ", peptide_pos_lookup$Sequence)
      cat(file = stderr(), "peptide_position_lookup...2", "\n")
      s <- strsplit(peptide_pos_lookup$Start_Stop, split = ", ")
      peptide_pos_lookup <- data.frame(AS = rep(peptide_pos_lookup$AS, sapply(s, length)), Position = unlist(s))
      peptide_pos_lookup <- peptide_pos_lookup %>% separate(AS, c("Accession", "Sequence"), sep = " ")
      peptide_pos_lookup <- peptide_pos_lookup %>% separate(Position, c("Start", "Stop"), sep = "-")
    }
    
    #------------------------------------------------------------------------------------------------------------------
    if (dpmsr_set$x$data_source == "SP") {
      # create peptide lookup table
      if (dpmsr_set$x$final_data_output == "Protein") {
          peptide_pos_lookup <- try(dpmsr_set$data$data_raw_peptide %>% dplyr::select(PG.ProteinAccessions, EG.ModifiedSequence, PEP.PeptidePosition))
      }
      
      cat(file = stderr(), "peptide_position_lookup...3", "\n")
      colnames(peptide_pos_lookup) <- c("Accession", "Sequence", "Start")
      peptide_pos_lookup <- peptide_pos_lookup %>% distinct()
      
      cat(file = stderr(), "peptide_position_lookup...4", "\n")
      peptide_pos_lookup <- expand_spectronaut_protein_lookup(peptide_pos_lookup)
      
      peptide_pos_lookup$Stop <- nchar(peptide_pos_lookup$Sequence) + as.numeric(peptide_pos_lookup$Start)
      
    }    

    peptide_pos_lookup$dup <- str_c(peptide_pos_lookup$Accession, peptide_pos_lookup$Sequence, peptide_pos_lookup$Start, peptide_pos_lookup$Stop)
    peptide_pos_lookup <- peptide_pos_lookup[!duplicated(peptide_pos_lookup$dup),]
    peptide_pos_lookup$dup <- NULL
  
    peptide_pos_lookup <- subset(peptide_pos_lookup, Accession %in% protein_accession )
    #peptide_pos_lookup <- subset(peptide_pos_lookup, Accession %in% "P60202"  )  
    cat(file = stderr(), "peptide_position_lookup...end", "\n")
    return(peptide_pos_lookup)
    
  }

  
  
    
  
#----------------------------------------------------------------------------------------
set_stat_parameters <- function(session, input, output)  {  
    #single shot observers
    dpmsr_set$data$stats$final_comp <<- input$select_final_data_stats
    dpmsr_set$y$stats$pvalue_cutoff <<- input$stats_pvalue_cutoff
    dpmsr_set$y$stats$foldchange_cutoff <<- input$stats_foldchange_cutoff
    dpmsr_set$y$stats$stats_spqc_cv_filter <<- input$stats_spqc_cv_filter
    dpmsr_set$y$stats$stats_spqc_cv_filter_factor <<- input$stats_spqc_cv_filter_factor
    dpmsr_set$y$stats$stats_comp_cv_filter <<- input$stats_spqc_cv_filter
    dpmsr_set$y$stats$stats_comp_cv_filter_factor <<- input$stats_spqc_cv_filter_factor
    dpmsr_set$y$peptide_missing_filter <<- input$peptide_missing_filter
    dpmsr_set$y$peptide_missing_factor <<- input$peptide_missing_factor
    dpmsr_set$y$peptide_cv_filter <<- input$peptide_cv_filter
    dpmsr_set$y$peptide_cv_factor <<- input$peptide_cv_factor
    dpmsr_set$y$stats$stats_peptide_minimum <<- input$stats_peptide_minimum
    dpmsr_set$y$stats$stats_peptide_minimum_factor <<- input$stats_peptide_minimum_factor
}

  #----------------------------------------------------------------------------------------
set_stat_headers <- function(comp_groups, i)  {  
  cat(file = stderr(), "set stat headers...", "\n")
  #set comp groups names
  comp_groups$fc[i] <- str_c(comp_groups$comp_name[i], "_FC")
  comp_groups$fc2[i] <- str_c(comp_groups$comp_name[i], "_FC2")
  comp_groups$pval[i] <- str_c(comp_groups$comp_name[i], "_pval")
  comp_groups$limma_pval[i] <- str_c(comp_groups$comp_name[i], "_limma_pval")
  comp_groups$exactTest[i] <- str_c(comp_groups$comp_name[i], "_exactTest")
  comp_groups$adjpval[i] <- str_c(comp_groups$comp_name[i], "_adjpval")
  comp_groups$cohensd[i] <- str_c(comp_groups$comp_name[i], "_cohensd")
  comp_groups$mf[i] <- str_c(comp_groups$comp_name[i], "_MF")
  
  if (dpmsr_set$x$primary_group) {
    comp_groups$cv_N[i] <- str_c(comp_groups$primary_comp_N[i], "_CV")  
    comp_groups$cv_D[i] <- str_c(comp_groups$primary_comp_D[i], "_CV")
  }else{
    comp_groups$cv_N[i] <- str_c(comp_groups$comp_N[i], "_CV")  
    comp_groups$cv_D[i] <- str_c(comp_groups$comp_D[i], "_CV")
  }
  return(comp_groups)
}
  
  
  
#----------------------------------------------------------------------------------------  
expand_spectronaut <- function(df)  {    
    
  df <- df %>% add_column(ProteinGroup = "", .after = 'ProteinName')
  dpmsr_set$y$info_columns <<- dpmsr_set$y$info_columns + 1
  dpmsr_set$y$info_columns_final <<- dpmsr_set$y$info_columns_final + 1
  df$ProteinGroup <- df$Accession
  
  #expand protein groups
  df$group <- ""
  #find groups
  for (i in 1:nrow(df)) {
    df$group[i] <- grepl(";", df$Accession[i])
  }
  #isolate groups from singles
  test1 <- df[df$group == FALSE,]
  test2 <- df[df$group == TRUE,]
  
  #loop through groups and create single entry for each accession
  for (i in 1:nrow(test2)) {
    accessions <- unlist(str_split(test2$Accession[i], ";"))
    genes <- unlist(str_split(test2$Gene[i], ";"))
    descriptions <- unlist(str_split(test2$Description[i], ";"))
    names <- unlist(str_split(test2$ProteinName[i], ";"))

    for (a in 1:length(accessions)) {
      temp_df <- test2[i,]
      temp_df$Gene <- genes[a]
      temp_df$Accession <- accessions[a]
      temp_df$Description <- descriptions[a]
      temp_df$ProteinName <- names[a]
      test1 <- rbind(test1, temp_df)
    }
  }

  test1$group <- NULL
  return(test1)
}
  
  #----------------------------------------------------------------------------------------  
  expand_spectronaut_protein_lookup <- function(df)  {    

    #expand protein groups
    df$group <- ""
    
    #find groups
    for (i in 1:nrow(df)) {
      df$group[i] <- grepl(";", df$Accession[i])
    }
    
    #isolate groups from singles
    test1 <- df[df$group == FALSE,]
    test2 <- df[df$group == TRUE,]
    
    #loop through groups and create single entry for each accession
    for (i in 1:nrow(test2)) {
      accessions <- unlist(str_split(test2$Accession[i], ";"))
      sequence <- test2$Sequence[i]
      start <- unlist(str_split(test2$Start[i], ";"))
      
      for (a in 1:length(accessions)) {
        temp_df <- test2[i,]
        temp_df$Accession <- accessions[a]
        temp_df$Sequence <- sequence
        temp_df$Start <- start[a]
        test1 <- rbind(test1, temp_df)
      }
    }
    
    test1$Start <- gsub(",.*", "",  test1$Start)
    
    test1$group <- NULL
    return(test1)
  }
  

