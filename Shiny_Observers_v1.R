
#-----------------------------------------------------------
#-----------------------------------------------------------

update_dpmsr_set_from_widgets <- function(session, input, output){
  #--Select type of data to load, and final format--------------------------------------------------    
  observe({
    if (input$radio_input==1){dpmsr_set$x$raw_data_input <<-"Protein_Peptide"}
    else if (input$radio_input==2){dpmsr_set$x$raw_data_input<<-"Protein"}
    else if (input$radio_input==3){dpmsr_set$x$raw_data_input<<-"Peptide"}
    else if (input$radio_input==4){dpmsr_set$x$raw_data_input<<-"PSM_FDR"}
  })
  
  observe({
    if (input$radio_output=="1"){
      dpmsr_set$x$final_data_output <<- "Protein"
      cat(file=stderr(), str_c("output test  ", input$radio_output), "\n")
      }
    else if (input$radio_output=="2"){
      dpmsr_set$x$final_data_output <<-"Peptide"
      updateCheckboxInput(session, "peptide_missing_filter",  value = FALSE)
      updateCheckboxInput(session, "peptide_cv_filter",  value = FALSE)
      updateCheckboxInput(session, "stats_peptide_minimum",  value = FALSE)
      cat(file=stderr(), str_c("output test  ", input$radio_output), "\n")
      }
  }) 
  

  
  observe({
    dpmsr_set$x$peptides_to_use <<- input$razor
  })
  
  observe({
    if (input$checkbox_isoform){dpmsr_set$x$peptide_isoform <<- TRUE}else{dpmsr_set$x$peptide_isoform <<- FALSE}
  })
  

  observe({
    dpmsr_set$x$file_prefix <<-  input$fileprefix
  })
  
  #-Filters-----------------------------------------------------------------------------------------------------    
  observe({
    if (input$checkbox_require_x){dpmsr_set$x$require_x <<-TRUE}else{dpmsr_set$x$require_x <<-FALSE}
  })
  
  observe({
    dpmsr_set$x$require_x_cutoff <<-  input$require_x_cutoff
  })
  
  observe({
    dpmsr_set$x$minimum_measured_all <<-  input$minimum_measured_all
  })
  
  observe({
    if (input$checkbox_filtercv){dpmsr_set$x$filter_cv <<-TRUE}else{dpmsr_set$x$filter_cv <<-FALSE}
  })
  
  observe({
    dpmsr_set$x$filter_cv_group <<-  input$text_filtercvgroup
  })
  
  observe({
    dpmsr_set$x$filter_cv_cutoff <<-  input$number_filtercvcutoff
  })
  
  observe({
    if (input$checkbox_norm_ptm){dpmsr_set$x$peptide_ptm_norm <<- TRUE
    }else{dpmsr_set$x$peptide_ptm_norm <<- FALSE}
  })
  
  observe({
    dpmsr_set$x$peptide_norm_grep <<-  input$peptide_norm_grep
  })
  
  observe({
    if (input$checkbox_impute_ptm){dpmsr_set$x$peptide_ptm_impute <<- TRUE
    }else{dpmsr_set$x$peptide_ptm_impute <<- FALSE}
  })
  
  observe({
    dpmsr_set$x$peptide_impute_grep <<-  input$peptide_impute_grep
  })
  
  observe({
    if (input$checkbox_report_ptm){dpmsr_set$x$peptide_report_ptm <<- TRUE
    }else{dpmsr_set$x$peptide_report_ptm <<- FALSE}
  })
  
  observe({
    dpmsr_set$x$peptide_report_grep <<-  input$peptide_report_grep
  })
  
  #-Normalize-----------------------------------------------------------------------------------------------------      
  
  observe({
    if (input$checkbox_n1){dpmsr_set$x$sl <<-TRUE}else{dpmsr_set$x$sl <<-FALSE}
  })
  observe({
    if (input$checkbox_n2){dpmsr_set$x$tmm <<-TRUE}else{dpmsr_set$x$tmm <<-FALSE}
  })
  observe({
    if (input$checkbox_n3){dpmsr_set$x$sltmm <<-TRUE}else{dpmsr_set$x$sltmm <<-FALSE}
  })
  observe({
    if (input$checkbox_n4){dpmsr_set$x$quantile <<-TRUE}else{dpmsr_set$x$quantile <<-FALSE}
  })
  observe({
    if (input$checkbox_n5){dpmsr_set$x$lr <<-TRUE}else{dpmsr_set$x$lr <<-FALSE}
  })
  observe({
    if (input$checkbox_n6){dpmsr_set$x$loess <<-TRUE}else{dpmsr_set$x$loess <<-FALSE}
  })
  observe({
    if (input$checkbox_n7){dpmsr_set$x$vsn <<-TRUE}else{dpmsr_set$x$vsn <<-FALSE}
  })
  observe({
    if (input$checkbox_n8){dpmsr_set$x$ti <<-TRUE}else{dpmsr_set$x$ti <<-FALSE}
  })
  observe({
    if (input$checkbox_n9){dpmsr_set$x$mi <<-TRUE}else{dpmsr_set$x$mi <<-FALSE}
  })    
  observe({
    if (input$checkbox_n10){dpmsr_set$x$ai <<-TRUE}else{dpmsr_set$x$ai <<-FALSE}
  })
  observe({
    if (input$checkbox_n11){dpmsr_set$x$protein <<-TRUE}else{dpmsr_set$x$protein <<-FALSE}
  })
  
  observe({
    dpmsr_set$x$protein_norm_list <<-  input$protein_norm_list
  })
  
  #-impute-----------------------------------------------------------------------------------------------------      
  observe({
    if (input$radio_impute ==1){dpmsr_set$x$impute_method <<-"Duke"}
    else if (input$radio_impute ==2){dpmsr_set$x$impute_method <<-"Floor"}
    else if (input$radio_impute ==3){dpmsr_set$x$impute_method <<-"Minimum"}
    else if (input$radio_impute ==4){dpmsr_set$x$impute_method <<-"Average/Group"}
    else if (input$radio_impute ==5){dpmsr_set$x$impute_method <<-"KNN"}
    else if (input$radio_impute ==6){dpmsr_set$x$impute_method <<-"LocalLeastSquares"}
    else if (input$radio_impute ==7){dpmsr_set$x$impute_method <<-"MLE"}    
    else if (input$radio_impute ==8){dpmsr_set$x$impute_method <<-"BottomX"}   
    else if (input$radio_impute ==9){dpmsr_set$x$impute_method <<-"Average/Global"}  
  })
  
  observe({
    dpmsr_set$x$bottom_x <<- input$bottom_x
  })
  
  observe({
    if (input$checkbox_impute_ptm){dpmsr_set$x$peptide_ptm_impute <<- TRUE
    }else{dpmsr_set$x$peptide_ptm_impute <<- FALSE}
  })
  
  observe({
    dpmsr_set$x$area_floor <<-  input$impute_floor
  })
  
  observe({
    if (input$checkbox_misaligned){dpmsr_set$x$duke_misaligned <<-TRUE}else{dpmsr_set$x$duke_misaligned <<- FALSE}
  })
  
  observe({
    dpmsr_set$x$missing_cutoff <<- input$missing_cutoff
  })
  
  observe({
    dpmsr_set$x$misaligned_cutoff <<- input$misaligned_cutoff
  })
  
  observe({
    dpmsr_set$x$int_cutoff_sd <<- input$intensity_cutoff_mean_sd
  })
  
  #-stats-----------------------------------------------------------------------------------------------------      
  
  observe({
    dpmsr_set$x$comp_number <<- input$comp_number
  })
  
  observe({
    if (input$pair_comp){dpmsr_set$x$pair_comp <<-TRUE}else{dpmsr_set$x$pair_comp <<- FALSE}
  })
  
  observe({
    dpmsr_set$x$pvalue_cutoff <<- input$pvalue_cutoff
  })
  
  observe({
    dpmsr_set$x$foldchange_cutoff <<- input$foldchange_cutoff
  })
  
  observe({
    dpmsr_set$x$missing_factor <<- input$missing_factor
  })
  
  observe({
    dpmsr_set$y$peptide_missing_filter <<- input$peptide_missing_filter
  })
  
  observe({
    dpmsr_set$y$peptide_missing_factor <<- input$peptide_missing_factor
  })
  
  observe({
    dpmsr_set$y$peptide_cv_filter <<- input$peptide_cv_filter
  })
  
  observe({
    dpmsr_set$y$peptide_cv_factor <<- input$peptide_cv_factor
  })
  
  observe({
    if (input$checkbox_report_accession){dpmsr_set$x$accession_report_out <<-TRUE
    }else{dpmsr_set$x$accession_report_out <<-FALSE}
  })
  
  observe({
    dpmsr_set$x$accession_report_list <<-  input$report_accession
  })
  
  
  #-Plots-----------------------------------------------------------------------------------------------------    
  observe({
    dpmsr_set$x$adh_list <<-  input$adh_list
  })
  observe({
    dpmsr_set$x$bait_list <<-  input$bait_list
  })
  observe({
    dpmsr_set$x$avidin_list <<-  input$avidin_list
  })
  observe({
    dpmsr_set$x$carbox_list <<-  input$carbox_list
  })
  observe({
    dpmsr_set$x$bira_list <<-  input$bira_list
  })
  observe({
    dpmsr_set$x$protein1_list <<-  input$protein1_list
  })
  observe({
    dpmsr_set$x$protein2_list <<-  input$protein2_list
  })
  observe({
    dpmsr_set$x$protein3_list <<-  input$protein3_list
  })
  observe({
    dpmsr_set$x$protein4_list <<-  input$protein4_list
  })
  observe({
    dpmsr_set$x$qc_spike_id <<-  input$protein_spike_list
  })
  observe({
    dpmsr_set$y$organism <<-  input$select_organism
  })

 
  

  
  
   
}


