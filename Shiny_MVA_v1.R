
#--------------------------------------------------------------------------------------
stat_prep <- function(){
  ncores <- detectCores()
  if (is.na(ncores) | ncores < 1) {ncores <- 1}
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
save_final_nostats <- function(){
  #save to excel final without stats
  for (name in names(dpmsr_set$data$final)){
    cat(file=stderr(), str_c("Saving final excel w/o stats....", name), "\n")
    filename <- str_c(dpmsr_set$file$output_dir, name, "//Final_", name, "_noStats.xlsx")
    Simple_Excel_name(dpmsr_set$data$final[[name]], filename, name)
  }
  
  for (name in names(dpmsr_set$data$impute)){
    cat(file=stderr(), str_c("Saving imputed peptide excel w/o stats....", name), "\n")
    filename <- str_c(dpmsr_set$file$output_dir, name, "//Imputed_Peptide_", name, ".xlsx")
    Simple_Excel_name(dpmsr_set$data$impute[[name]], filename, name)
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
  
  #this sets the final table for graphs to use only those that will be reported
  if (as.logical(dpmsr_set$x$peptide_ptm_norm)){
    data_in <- data_in [grep(dpmsr_set$x$peptide_norm_grep, data_in$Modifications),]
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
  stats_data <- dpmsr_set$data$final[[input$select_final_data_stats]]
  stats_data <- stats_data[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
  
  if(input$select_final_data_stats == "impute"){
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
    
    #REBUILDING COMPARISON SELECTOR
    # using contains instead of matches, trying to find all inclusive sample group with not SPQC type
    #for(stats_group in input[[str_c('comp_',i,'N')]]){stats_data_N <- stats_data_N %>% dplyr::select(contains(stats_group))   }
    #for(stats_group in input[[str_c('comp_',i,'D')]]){stats_data_D <- stats_data_D %>% dplyr::select(contains(stats_group))   }
 
    for(stats_group in input[[str_c('comp_',i,'N')]]){
      grep_comp <- paste("[ _]", stats_group, "[ _$]", sep="")
      stats_data_N <- stats_data_N[,grepl(grep_comp, colnames(stats_data_N))]
      }
    
    for(stats_group in input[[str_c('comp_',i,'D')]]){
      grep_comp <- paste("[ _]", stats_group, "[ _$]", sep="")
      stats_data_D <- stats_data_D[,grepl(grep_comp, colnames(stats_data_D))]
     }

      
    #work around to get all samples into a comparison for graphs
    if( input[[str_c('comp_',i,'N')]] == "All.Samples"){
      stats_data_N <- stats_data %>% dplyr::select(!contains(input$comp_spqc))
    }
    if( input[[str_c('comp_',i,'D')]] == "All.Samples"){
      stats_data_D <- stats_data %>% dplyr::select(!contains(input$comp_spqc))
    }
    
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
    comp_groups$adjpval[i] <- str_c(comp_groups$comp_name[i], "_adjpval")
    comp_groups$cohensd[i] <- str_c(comp_groups$comp_name[i], "_cohensd")
    comp_groups$mf[i] <- str_c(comp_groups$comp_name[i], "_MF")
    comp_groups$cv_N[i] <- str_c(comp_groups$comp_N[i], "_CV")  
    comp_groups$cv_D[i] <- str_c(comp_groups$comp_D[i], "_CV")  
    comp_groups$N_count[i] <- ncol(stats_data_N)
    comp_groups$D_count[i] <- ncol(stats_data_D)
    comp_groups$N_color[i] <- color_choices_N[i] #first i was 1
    comp_groups$D_color[i] <- color_choices_D[i] #first i was 1
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
  
# update UI widgets after stat comparisons are selected
  updateTextInput(session, "stats_barplot_title", value = str_c("Barplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  updateTextInput(session, "stats_boxplot_title", value = str_c("Boxplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  update_comparisons(session, input, output)
}


#----no longer used----------------------------------------------------------------------------------------------------------------------------
#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
# stat_calc <- function(session, input, output) {
#   data_in <- dpmsr_set$data$final[[input$select_final_data_stats]]
# 
#   annotate_in <- data_in[1:dpmsr_set$y$info_columns_final]
#   data_in <- data_in[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
#   #start df for stats
#   stat_df <- annotate_in[1:1]
#   imputed_df <- dpmsr_set$data$Protein_imputed_df[3:ncol(dpmsr_set$data$Protein_imputed_df)]  
#   
#   for(i in 1:nrow(dpmsr_set$y$stats$groups)) 
#   {
#     stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_N[i], "_CV")] <- percentCV_gw(data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ])
#     stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_D[i], "_CV")] <- percentCV_gw(data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ])
#   } 
# 
#   stat_df[ , str_c(dpmsr_set$y$stats$comp_spqc, "_CV")] <- percentCV_gw(data_in[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)])
#   
#   #generate pvalue and FC for each comparison
#   for(i in 1:nrow(dpmsr_set$y$stats$groups))
#   {
#     comp_N_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
#     comp_D_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
#     stat_df[ ,dpmsr_set$y$stats$groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
#     stat_df[ ,dpmsr_set$y$stats$groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
#     stat_df[ ,dpmsr_set$y$stats$groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
#     #stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
#     #stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
#     comp_N_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
#     comp_D_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
#     stat_df[ ,dpmsr_set$y$stats$groups$mf[i]] <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
#   } 
#   
#   # Create final tables--------------------------------------------------
#   if(input$select_final_data_stats == "impute"){
#     colnames(data_in) <- dpmsr_set$design$Header3
#   }else{
#     colnames(data_in) <- dpmsr_set$design$Header2
#   }
#   data_table <- cbind(annotate_in, data_in, stat_df[2:ncol(stat_df)])
#   dpmsr_set$data$stats$final <<- data_table
#   
#   #single shot observers
#   dpmsr_set$data$stats$final_comp <<- input$select_final_data_stats
#   dpmsr_set$y$stats$pvalue_cutoff <<- input$stats_pvalue_cutoff
#   dpmsr_set$y$stats$foldchange_cutoff <<- input$stats_foldchange_cutoff
#   dpmsr_set$y$stats$stats_spqc_cv_filter <<- input$stats_spqc_cv_filter
#   dpmsr_set$y$stats$stats_spqc_cv_filter_factor <<- input$stats_spqc_cv_filter_factor
#   dpmsr_set$y$stats$stats_comp_cv_filter <<- input$stats_spqc_cv_filter
#   dpmsr_set$y$stats$stats_comp_cv_filter_factor <<- input$stats_spqc_cv_filter_factor
#   return()
# }



#--------------------------------------------------------------------------------------------------------------------------------
#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_calc2 <- function(session, input, output) {
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

  
  for(i in 1:nrow(dpmsr_set$y$stats$groups)) 
  {
    cat(file=stderr(), str_c("Calculating stats...", i, " of ", nrow(dpmsr_set$y$stats$groups)), "\n")
    data_in <- dpmsr_set$data$impute[[input$select_final_data_stats]]
    info_columns <- ncol(data_in) - dpmsr_set$y$sample_number
    annotate_in <- data_in[1:info_columns]
    data_in <- data_in[(info_columns+1):ncol(data_in)]
    
    imputed_info_columns <- ncol(dpmsr_set$data$peptide_imputed_df) - dpmsr_set$y$sample_number
    imputed_df <- dpmsr_set$data$peptide_imputed_df[(imputed_info_columns+1):ncol(dpmsr_set$data$peptide_imputed_df)]  
    
    
    comp_N_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    spqc_data <- data_in[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    spqc_impute_data <- imputed_df[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    comp_N_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    comp_N_imputed <- mutate_all(comp_N_imputed, function(x) as.numeric(x))
    comp_D_imputed <- mutate_all(comp_D_imputed, function(x) as.numeric(x))
    spqc_impute_data <- mutate_all(spqc_impute_data, function(x) as.numeric(x))
    
    sample_count <- ncol(comp_N_data) + ncol(comp_D_data) + ncol(spqc_data)
    
    #--reduce peptide PD imputed column to only samples of comparison----------------------------------------------------------
    if(dpmsr_set$x$final_data_output == "Peptide"){
      df_impute_peptide <- cbind(comp_N_imputed, comp_D_imputed, spqc_impute_data)
      df_impute_peptide[df_impute_peptide==0] <- "-"
      while (ncol(df_impute_peptide)>1) {
        df_impute_peptide[,1] <- str_c(df_impute_peptide[,1], ".", df_impute_peptide[,2])
        df_impute_peptide[,2] <- NULL
      }
      colnames(df_impute_peptide) <- "PD_Detected_Peptides"
      annotate_in$PD_Detected_Peptides <- df_impute_peptide$PD_Detected_Peptides
    }
    
    #--section for peptide to protein final data----------------------------------------------------------
    if(dpmsr_set$x$final_data_output == "Protein"){

      cat(file=stderr(), str_c("Final data output is ", dpmsr_set$x$final_data_output), "\n")
      
      mf <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
      N_CV <- percentCV_gw(comp_N_data)
      D_CV <- percentCV_gw(comp_D_data)
      minCV <- pmin(N_CV, D_CV)
      
      
      # exclude TMT SPQC norm - can not have imputed peptide level information
      if (as.logical(dpmsr_set$x$tmt_spqc_norm)){
        df <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data)
      }else{
        df <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data, mf, N_CV, D_CV, minCV)
        df_impute <- cbind(annotate_in, comp_N_imputed, comp_D_imputed, mf, N_CV, D_CV, minCV)
        df_impute_summary <- cbind(annotate_in, comp_N_imputed, comp_D_imputed, spqc_impute_data, mf, N_CV, D_CV, minCV)
      
      
        if(input$peptide_missing_filter){
          df <- df[(df$mf >= input$peptide_missing_factor),]
          df_impute <- df_impute[(df_impute$mf >= input$peptide_missing_factor),]  
          df_impute_summary <- df_impute_summary[(df_impute_summary$mf >= input$peptide_missing_factor),]  
        }
    
        if(input$peptide_cv_filter){
          df <- df[(df$minCV <= input$peptide_cv_factor),]
          df_impute <- df_impute[(df_impute$minCV <= input$peptide_cv_factor),]  
          df_impute_summary <- df_impute_summary[(df_impute_summary$minCV <= input$peptide_cv_factor),]  
        }
      
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
        
        df_impute_summary <- df_impute_summary[(info_columns+1):ncol(df_impute_summary)] %>% mutate_all(as.character)
        while (ncol(df_impute_summary)>1) {
          df_impute_summary[,1] <- str_c(df_impute_summary[,1], ".", df_impute_summary[,2])
          df_impute_summary[,2] <- NULL
        }
        colnames(df_impute_summary) <- "PD_Detected_Peptides"
        
        df$PD_Detected_Peptides <- df_impute_summary
        dpmsr_set$data$stats$peptide[[dpmsr_set$y$stats$groups$comp_name[i]]] <<- df
        
        #collapse peptide data and imputed peptide info, this is done with only the stat groups
        df <- collapse_peptide_stats(df, info_columns)
        df_impute <- collapse_peptide_stats(df_impute, info_columns)
      
        #xdf <<- df
        #xdf_impute <<- df_impute
        
        #create string column with imputed peptide info
        info_columns <- ncol(df) - sample_count
        test2 <- df_impute[(info_columns+1):ncol(df_impute)]
        
        test2[test2==0] <- "-"
        test2 <- test2 %>% mutate_all(as.character)
        while (ncol(test2)>1) {
          test2[,1] <- str_c(test2[,1], ".", test2[,2])
          test2[,2] <- NULL
        }
        
        df <- add_column(df, test2[,1] , .after = "Peptides")
        names(df)[4] <- "PD_Detected_Peptides"
        #xdf2 <<- df
        
        #xdfimpute2 <<- df_impute
        imputed_df <- df_impute[4:ncol(df_impute)]  
        imputed_df[imputed_df>1] <- 1
        comp_N_imputed <- imputed_df[1:dpmsr_set$y$stats$groups$N_count[i]]
        comp_D_imputed <- imputed_df[(dpmsr_set$y$stats$groups$N_count[i]+1):(dpmsr_set$y$stats$groups$N_count[i]+dpmsr_set$y$stats$groups$D_count[i])]
        
      }
      
      #xdf2 <<- df
      
      annotate_in <- df[1:dpmsr_set$y$info_columns_final]
      
      if(input$checkbox_add_gene_column) {
        add_gene <- str_extract(annotate_in$Description, "GN=\\w*")
        add_gene <- gsub("GN=", "", add_gene)
        annotate_in <- annotate_in %>% add_column(Gene = add_gene, .before = 3)
        dpmsr_set$data$final[[input$select_final_data_stats]] <<- dpmsr_set$data$final[[input$select_final_data_stats]] %>% 
          add_column(Gene = annotate_in$Gene, .before = 3)
        
      }
      
      data_in <- df[(dpmsr_set$y$info_columns_final+1):(ncol(df))]
      
      comp_N_data <- data_in[1:dpmsr_set$y$stats$groups$N_count[i]]
      comp_D_data <- data_in[(dpmsr_set$y$stats$groups$N_count[i]+1):(dpmsr_set$y$stats$groups$N_count[i]+dpmsr_set$y$stats$groups$D_count[i])]
      spqc_data <- data_in[(dpmsr_set$y$stats$groups$N_count[i]+dpmsr_set$y$stats$groups$D_count[i]+1):ncol(data_in)]
    }
    
    #--end of peptide to protein final data-----------------------------------------------
    
    
    #start df for stats -----------------------------------------
    stat_df <- annotate_in[1:1]
    stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_N[i], "_CV")] <- percentCV_gw(comp_N_data)
    stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_D[i], "_CV")] <- percentCV_gw(comp_D_data)
    stat_df[ , str_c(dpmsr_set$y$stats$comp_spqc, "_CV")] <- percentCV_gw(spqc_data)
    stat_df[ ,dpmsr_set$y$stats$groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
    
    if(input$checkbox_adjpval){
      stat_df[ ,dpmsr_set$y$stats$groups$adjpval[i]] <- p.adjust(stat_df[ ,dpmsr_set$y$stats$groups$pval[i]], method = input$padjust_options) 
    }
    
    if(input$checkbox_cohensd){
        stat_df[ ,dpmsr_set$y$stats$groups$cohensd[i]] <- cohend_gw(comp_N_data, comp_D_data, as.logical(input$checkbox_cohensd))
    }
    
    if(input$checkbox_limmapvalue){
      stat_df[ ,dpmsr_set$y$stats$groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data)
    }
    #stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
    #stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
    
    #only include if not TMT SPQC norm
    if (!as.logical(dpmsr_set$x$tmt_spqc_norm)){
      stat_df[ ,dpmsr_set$y$stats$groups$mf[i]] <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
    }
    
    # Create final tables--------------------------------------------------
    if(input$select_final_data_stats == "impute"){
      colnames(comp_N_data) <- dpmsr_set$design$Header3[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i])]
      colnames(comp_D_data) <- dpmsr_set$design$Header3[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i])]
      colnames(spqc_data) <- dpmsr_set$design$Header3[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    }else{
      colnames(comp_N_data) <- dpmsr_set$design$Header2[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i])]
      colnames(comp_D_data) <- dpmsr_set$design$Header2[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i])]
      colnames(spqc_data) <- dpmsr_set$design$Header2[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
    }
    
    
    data_table <- cbind(annotate_in, comp_N_data, comp_D_data, spqc_data, stat_df[2:ncol(stat_df)])
    
    data_table$Stats <- ""
    data_table_cols <- ncol(data_table)
    
    if (!input$checkbox_filter_adjpval) {
      data_table$pval <- ifelse(data_table[[dpmsr_set$y$stats$groups$pval[i]]]  <= input$pvalue_cutoff, 0, 1)
    }else{
      data_table$pval <- ifelse(data_table[[dpmsr_set$y$stats$groups$adjpval[i]]]  <= input$pvalue_cutoff, 0, 1)
    }
    
    if (!as.logical(dpmsr_set$x$tmt_spqc_norm)){
      data_table$mf <- ifelse(data_table[[dpmsr_set$y$stats$groups$mf[i]]]  >= input$missing_factor, 0, 1)
    }
    
    if(input$stats_spqc_cv_filter){
      data_table$spqc <- ifelse(data_table[[str_c(dpmsr_set$y$stats$comp_spqc, "_CV")]] <= input$stats_spqc_cv_filter_factor, 0, 1)
    }
    if(input$stats_comp_cv_filter){
      data_table$cv <- ifelse(data_table[[str_c(dpmsr_set$y$stats$groups$comp_N[i], "_CV")]] <= input$stats_comp_cv_filter_factor |
                                data_table[[str_c(dpmsr_set$y$stats$groups$comp_D[i], "_CV")]] <= input$stats_comp_cv_filter_factor, 0, 1)
    }
    if(input$stats_peptide_minimum){
      data_table$pepmin <- ifelse(data_table$Peptides >= input$stats_peptide_minimum_factor, 0, 1)
    }
    
    data_table$sum <- rowSums(data_table[(data_table_cols+1):ncol(data_table)])
    
    data_table$Stats <- ifelse(data_table$sum == 0 & data_table[[dpmsr_set$y$stats$groups$fc[i]]] >= input$foldchange_cutoff, "Up", 
                         ifelse(data_table$sum == 0 & data_table[[dpmsr_set$y$stats$groups$fc[i]]] <= -input$foldchange_cutoff, "Down", ""))       
    
    data_table <- data_table[1:data_table_cols]
    
    data_table <- data_table[order(data_table$Stats, -data_table[[dpmsr_set$y$stats$groups$pval[i]]], decreasing = TRUE),]
    
    #filter stats by ptm or accession if needed
    if(dpmsr_set$x$final_data_output == "Peptide" &  input$checkbox_report_ptm){
      data_table <- data_table[grep(dpmsr_set$x$peptide_report_grep, data_table$Modifications),]
    }
    
    if(dpmsr_set$x$final_data_output == "Peptide" & input$checkbox_report_accession){
      data_table <-subset(data_table, Accession %in% dpmsr_set$x$accession_report_list )
    }
    
    dpmsr_set$data$stats[[dpmsr_set$y$stats$groups$comp_name[i]]] <<- data_table
  } 
    
  cat(file=stderr(), str_c(dpmsr_set$y$stats$groups$comp_name[i], " data has ", nrow(dpmsr_set$data$stats[[dpmsr_set$y$stats$groups$comp_name[i]]]), " rows"), "\n")
  cat(file=stderr(), "Calculating stats complete...", "\n")
  return()
}




#----------------------------------------------------------------------------------------
# create final excel documents
stats_Final_Excel <- function(session, input, output) {
  cat(file=stderr(), "Creating Excel Output File...1", "\n")
  require(openxlsx)
    
    filename <- str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp, "//", input$final_stats_name)
    df <- dpmsr_set$data$final[[input$select_final_data_stats]]
    
    #testing
    input_final_stats_name <<- input$final_stats_name
    input_select_final_data_stats <<- input$select_final_data_stats
    # filename <- str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp, "//", input_final_stats_name)
    # df <- dpmsr_set$data$final[[input_select_final_data_stats]]
    
    #remove FC2 from df for excel
    #df <- df[,-grep(pattern="_FC2", colnames(df))]
    
    if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
      df_raw <- dpmsr_set$data$finalraw
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
    cat(file=stderr(), "Creating Excel Output File...2", "\n")
    df2 <- df
    
    
    # subsets data for accession specific reports
    if (as.logical(dpmsr_set$x$accession_report_out)){
      df <-subset(df, Accession %in% dpmsr_set$x$accession_report_list )
      df2 <-subset(df2, Accession %in% dpmsr_set$x$accession_report_list )
    }
    
    nextsheet <- 1
    
    wb <- createWorkbook()
    
    cat(file=stderr(), "Creating Excel Output File...3", "\n")
    #----------------------------------------
    
    if(site_user == "dpmsr" && dpmsr_set$x$raw_data_input != "Protein"){
  
     #raw_protein <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Protein_to_Protein_Raw.xlsx"))
     #raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_ProteinPeptide_to_Peptide_Raw.xlsx"))
     #raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Isoform_to_Isoform_Raw.xlsx"))
     #raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Peptide_to_Peptide_Raw.xlsx"))
  
    if(as.logical(dpmsr_set$x$peptide_isoform) ){
        raw_peptide <- dpmsr_set$data$data_peptide_isoform_start
      }else{
        raw_peptide <- dpmsr_set$data$data_peptide_start
      }

      peptide_impute <- dpmsr_set$data$impute$impute
      
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet +1
      addWorksheet(wb, "Raw Peptide Data")
      writeData(wb, sheet = nextsheet, raw_peptide)
      nextsheet <- nextsheet +1
      addWorksheet(wb, "Imputed Peptide Data")
      writeData(wb, sheet = nextsheet, peptide_impute) 
      nextsheet <- nextsheet +1
    }
    
    if(site_user == "dpmsr" && dpmsr_set$x$raw_data_input == "Protein"){
      raw_protein <- pmsr_set$data$data_protein_start
      protein_impute <- dpmsr_set$data$impute$impute
      
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet +1
      addWorksheet(wb, "Raw Protein Data")
      writeData(wb, sheet = nextsheet, raw_protein)
      nextsheet <- nextsheet +1
      addWorksheet(wb, "Imputed Protein Data")
      writeData(wb, sheet = nextsheet, protein_impute) 
      nextsheet <- nextsheet +1
    }
    #----------------------------------------
    
    
    
    #addWorksheet(wb, deparse(substitute(df2)))
    addWorksheet(wb, "Normalized Data")
    writeData(wb, sheet = nextsheet, df2)  
    nextsheet <- nextsheet +1
    z=1
    for(i in 1:as.numeric(dpmsr_set$y$stats$comp_number))  {
      comp_string <- dpmsr_set$y$stats$groups$comp_name[i]
      filtered_df <- dpmsr_set$data$stats[[comp_string]]
      cat(file=stderr(), str_c("Saving excel...", comp_string), "\n")
      addWorksheet(wb, comp_string)
      writeData(wb, sheet = nextsheet, filtered_df)
      nextsheet <- nextsheet +1
      z <- z+2
    }
    cat(file=stderr(), "writting excel to disk...", "\n")
    saveWorkbook(wb, filename, overwrite = TRUE)
    cat(file=stderr(), "Creating Excel Output File...end", "\n")
}

# data table filter ------------------------------------------------------

stats_data_table_filter <- function(session, input, output) {
  cat(file=stderr(), "stats_data_table_filter... start", "\n")
  
  comp_string <- input$stats_select_data_comp
  comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
  sample_number <- dpmsr_set$y$stats$groups$N_count[comp_number] + dpmsr_set$y$stats$groups$D_count[comp_number]
  
  cat(file=stderr(), "stats_data_table_filter... 1", "\n")
  df_filter <- dpmsr_set$data$stats[[comp_string]]
  test_df_filter <<- df_filter
  
  cat(file=stderr(), "stats_data_table_filter... 2", "\n")
  if(input$stats_add_filters){
    df_filter <- df_filter %>% filter(df_filter$Stats == "Up" | df_filter$Stats == "Down")
  }
  
  cat(file=stderr(), "stats_data_table_filter... 3", "\n")
  
  if(input$stats_data_topn != 0 ){
    df_filter$sum <- rowSums(df_filter[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+sample_number)])
    df_filter <- df_filter[order(-df_filter$sum),]                      
    df_filter <- df_filter[1:input$stats_data_topn,]
    df_filter$sum <- NULL
  }
  
  cat(file=stderr(), "stats_data_table_filter... 4", "\n")
  
  if(input$stats_data_accession != "0" ){
    df_filter <-subset(df_filter, Accession %in% as.character(input$stats_data_accession)  )
  }
  
  if(input$stats_data_description != "0") {
    df_filter <-df_filter[grep(as.character(input$stats_data_description), df_filter$Description), ]
  }
  cat(file=stderr(), "stats_data_table_filter... end", "\n")
  return(df_filter)
}



# peptide zscore ------------------------------------------------------

peptide_zscore <- function(df_peptide, info_columns) {
  info_df <- df_peptide[1:info_columns]
  df <- df_peptide[(info_columns+1):ncol(df_peptide)]
  df<- log(df, 2)
  z_mean <- rowMeans(df)
  z_stdev  <- apply(df[1:ncol(df)], 1, sd)
  
  df <- apply(df, MARGIN = 2, function(x) (x-z_mean)/z_stdev  )
  df <- data.frame(cbind(info_df, df), stringsAsFactors = FALSE)
  
  return(df)
  
}


# one protein data ------------------------------------------------------

oneprotein_data <- function(session, input, output) {
  
  cat(file=stderr(), "One Protein stats and plots..." , "\n")
  
  comp_string <- input$stats_oneprotein_plot_comp
  comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string)
  
  df <- dpmsr_set$data$stats[[comp_string]]
  df <-subset(df, Accession %in% as.character(input$stats_oneprotein_accession)  )
  sample_number <- dpmsr_set$y$stats$groups$N_count[comp_number] + dpmsr_set$y$stats$groups$D_count[comp_number]
  
  df_peptide <- dpmsr_set$data$stats$peptide[[comp_string ]]    
  df_peptide <- subset(df_peptide, Accession %in% as.character(input$stats_oneprotein_accession)  )
  
  peptide_info_columns <- ncol(df_peptide) - sample_number - dpmsr_set$y$stats$comp_spqc_number
  
  #add spqc to plots
  if(input$stats_oneprotein_plot_spqc){
    df <- df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
  }else{
    df <- df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final + sample_number)]
  }
  
  comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
  #add spqc to plots
  if(input$stats_oneprotein_plot_spqc){
    comp_rows <- c(comp_rows, dpmsr_set$y$stats$comp_spqc_sample_numbers)
  }
  
  comp_rows <- unlist(comp_rows)
  namex <- dpmsr_set$design$Label[comp_rows]
  color_list <- dpmsr_set$design$colorlist[comp_rows]
  groupx <- dpmsr_set$design$Group[comp_rows]
  
  
  colnames(df_peptide)[(peptide_info_columns+1):ncol(df_peptide)] <- 
    c(dpmsr_set$design$Header1[unlist(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number])],
      dpmsr_set$design$Header1[unlist(dpmsr_set$y$stats$groups$sample_numbers_D[comp_number])],
      dpmsr_set$design$Header1[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)] )
  
  
  #sort peptides by intensity, keep highest abundant peptide
  df_peptide$sum <- rowSums(df_peptide[(peptide_info_columns+1):(peptide_info_columns+sample_number+1)])
  df_peptide <- df_peptide[order(df_peptide$sum), ]
  df_peptide$sum <- NULL
  df_peptide <- df_peptide[!duplicated(df_peptide$Sequence),]
  
  if(input$stats_use_zscore){df_peptide <- peptide_zscore(df_peptide, peptide_info_columns)}
  
  return(list("df"=df, "df_peptide"=df_peptide, "namex"=namex, "color_list"=color_list, "comp_string"=comp_string, "peptide_info_columns"=peptide_info_columns))
}

# one peptide data ------------------------------------------------------

onepeptide_data <- function(session, input, output) {
    
    cat(file=stderr(), "One Peptide stats and plots..." , "\n")
    
    comp_string <- input$stats_onepeptide_plot_comp
    comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string)
    
    df <- dpmsr_set$data$stats[[comp_string]]
    df <-df[grep(as.character(input$stats_onepeptide_accession), df$Accession), ]
    sample_number <- dpmsr_set$y$stats$groups$N_count[comp_number] + dpmsr_set$y$stats$groups$D_count[comp_number]
    df_peptide <- df
    info_columns <-dpmsr_set$y$info_columns_final
    df_peptide_stats <- df_peptide[(info_columns + sample_number + dpmsr_set$y$stats$comp_spqc_number+1):ncol(df_peptide)]
    
    #continue filtering of df, used for first plot of a single peptide
    df <-subset(df, Sequence %in% as.character(input$stats_onepeptide_sequence)  )
    grep_mod <- stringr::str_replace_all(input$stats_onepeptide_modification, "\\[", "\\\\[")
    grep_mod <- stringr::str_replace_all(grep_mod, "\\]", "\\\\]")
    grep_mod <- stringr::str_replace_all(grep_mod,"\\(", "\\\\(")
    grep_mod <- stringr::str_replace_all(grep_mod,"\\)", "\\\\)")
    df <- df %>% filter(stringr::str_detect(Modifications, grep_mod) )
    
    #add spqc to plots
    if(input$stats_onepeptide_plot_spqc){
      df_peptide <- df_peptide[1:(info_columns + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
      df <- df[(info_columns+1):(info_columns + sample_number + dpmsr_set$y$stats$comp_spqc_number)]
      comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number], 
                     dpmsr_set$y$stats$comp_spqc_sample_numbers )
    }else{
      df_peptide <- df_peptide[1:(info_columns + sample_number)]
      df <- df[(info_columns+1):(info_columns + sample_number)]
      comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
    }

    comp_rows <- unlist(comp_rows)
    namex <- dpmsr_set$design$Label[comp_rows]
    color_list <- dpmsr_set$design$colorlist[comp_rows]
    groupx <- dpmsr_set$design$Group[comp_rows]


    if(input$stats_onepeptide_use_zscore){df_peptide <- peptide_zscore(df_peptide, info_columns)}
    
    return(list("df"=df, "df_peptide"=df_peptide, "df_peptide_stats"=df_peptide_stats, "info_columns"=info_columns, "namex"=namex, "color_list"=color_list, "comp_string"=comp_string))
}
    



# stat filter ------------------------------------------------------  
    stats_filter <- function(df, comp_number) { 
      
      filtered_df <- filter(df, df[[dpmsr_set$y$stats$groups$pval[comp_number]]] <= as.numeric(dpmsr_set$x$pvalue_cutoff) &  
                              (df[[dpmsr_set$y$stats$groups$fc[comp_number]]] >= as.numeric(dpmsr_set$x$foldchange_cutoff) | 
                                 df[[dpmsr_set$y$stats$groups$fc[comp_number]]] <= -as.numeric(dpmsr_set$x$foldchange_cutoff))) 
      
      filtered_df <- filter(filtered_df, filtered_df[[dpmsr_set$y$stats$groups$mf[comp_number] ]] >= as.numeric(dpmsr_set$x$missing_factor))
      
      
      if(dpmsr_set$y$stats$stats_spqc_cv_filter){
        filtered_df <- filter(filtered_df, filtered_df[[str_c(dpmsr_set$y$stats$comp_spqc, "_CV")]] <= as.numeric(dpmsr_set$y$stats$stats_spqc_cv_filter_factor ) )
      }
      
      if(dpmsr_set$y$stats$stats_comp_cv_filter){
        filtered_df <- filter(filtered_df, filtered_df[[dpmsr_set$y$stats$groups$cv_N[comp_number]]] <= as.numeric(dpmsr_set$y$stats$stats_comp_cv_filter_factor) | 
                                filtered_df[[dpmsr_set$y$stats$groups$cv_D[comp_number]]] <= as.numeric(dpmsr_set$y$stats$stats_comp_cv_filter_factor) )
      }
     
      return(filtered_df)
       
    }      
    
    
  peptide_position_lookup <- function(session, input, output, protein_accession)  {
    # create peptide lookup table
    if(dpmsr_set$x$final_data_output == "Protein"){
      if(class(dpmsr_set$data$data_raw_peptide$Positions.in.Master.Proteins)=="NULL"){
        peptide_pos_lookup <- try(dpmsr_set$data$data_raw_peptide %>% dplyr::select(Sequence, Positions.in.Proteins))
      }else{
        peptide_pos_lookup <- try(dpmsr_set$data$data_raw_peptide %>% dplyr::select(Sequence, Positions.in.Master.Proteins))
      }
    }else{
      peptide_pos_lookup <- dpmsr_set$data$data_raw_isoform %>% dplyr::select(Sequence, Positions.in.Proteins)
    }
    
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
    s <- strsplit(peptide_pos_lookup$Start_Stop, split = ", ")
    peptide_pos_lookup <- data.frame(AS = rep(peptide_pos_lookup$AS, sapply(s, length)), Position = unlist(s))
    peptide_pos_lookup <- peptide_pos_lookup %>% separate(AS, c("Accession", "Sequence"), sep = " ")
    peptide_pos_lookup <- peptide_pos_lookup %>% separate(Position, c("Start", "Stop"), sep = "-")
    peptide_pos_lookup$dup <- str_c(peptide_pos_lookup$Accession, peptide_pos_lookup$Sequence, peptide_pos_lookup$Start, peptide_pos_lookup$Stop)
    peptide_pos_lookup <- peptide_pos_lookup[!duplicated(peptide_pos_lookup$dup),]
    peptide_pos_lookup$dup <- NULL
    
    peptide_pos_lookup <- subset(peptide_pos_lookup, Accession %in% protein_accession )
    #peptide_pos_lookup <- subset(peptide_pos_lookup, Accession %in% "P60202"  )  
  
    return(peptide_pos_lookup)
    
  }
  
  
  
  
  
  
  