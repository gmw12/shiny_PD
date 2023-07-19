
update_widget_startup <- function(session, input, output){
  cat(file=stderr(), "update_widget_startup...1", "\n")
  updateTextInput(session, "fileprefix", value = as.character(dpmsr_set$x$file_prefix))
  
  if(dpmsr_set$x$final_data_output == "Protein"){
    fdo <- 1
  }else if(dpmsr_set$x$final_data_output == "Peptide"){
    fdo <- 2
    }
  updateRadioButtons(session, "radio_output", selected = fdo )
  
  cat(file=stderr(), "update_widget_startup...2", "\n")
  
  if(dpmsr_set$x$raw_data_input == "Protein_Peptide"){rdi <- 1}
    else if(dpmsr_set$x$raw_data_input == "Protein"){rdi <- 2}
    else if(dpmsr_set$x$raw_data_input == "Peptide"){rdi <- 3}
    else if(dpmsr_set$x$raw_data_input == "PSM_FDR"){rdi <- 4}
  updateRadioButtons(session, "radio_input", selected = rdi )
  
  cat(file=stderr(), "update_widget_startup...3", "\n")
  if(is.null(dpmsr_set$x$data_source)){
    dpmsr_set$x$data_source <- "PD"
    updateRadioButtons(session, "data_source", selected = 1) 
  }else{
    try(if(dpmsr_set$x$data_source == "SP"){updateRadioButtons(session, "data_source", selected = 2) })
    try(if(dpmsr_set$x$data_source == "PD"){updateRadioButtons(session, "data_source", selected = 1) })
  }
  
  cat(file=stderr(), "update_widget_startup...4", "\n")
  updateCheckboxInput(session, "primary_group", value = as.logical(dpmsr_set$x$primary_group)) 
  
  updateSelectInput(session, "razor", selected = dpmsr_set$x$peptides_to_use )

  if (as.logical(dpmsr_set$x$peptide_isoform)) {test <- 1}else{test<-0}
  
  updateCheckboxInput(session, "checkbox_isoform", value = test)

  updateCheckboxInput(session, "checkbox_norm_ptm", value = as.logical(dpmsr_set$x$peptide_ptm_norm ))
  updateCheckboxInput(session, "checkbox_impute_ptm", value = as.logical(dpmsr_set$x$peptide_ptm_impute ))
  
  updateTextInput(session, "peptide_norm_grep", value = as.character(dpmsr_set$x$peptide_norm_grep ))
  updateTextInput(session, "peptide_impute_grep", value = as.character(dpmsr_set$x$peptide_impute_grep ))
  
  updateCheckboxInput(session, "checkbox_norm_include", value = as.logical(dpmsr_set$x$checkbox_norm_include))
  updateCheckboxInput(session, "checkbox_norm_exclude", value = as.logical(dpmsr_set$x$checkbox_norm_exclude ))
  updateTextInput(session, "include_norm_grep", value = as.character(dpmsr_set$x$include_norm_grep ))
  updateTextInput(session, "exclude_norm_grep", value = as.character(dpmsr_set$x$exclude_norm_grep ))
  
  updateCheckboxInput(session, "checkbox_tmt", value = as.logical(dpmsr_set$x$tmt_spqc_norm)) 
  updateCheckboxInput(session, "checkbox_report_accession", value = as.logical(dpmsr_set$x$accession_report_out)) 
  
  cat(file=stderr(), "update_widget_startup...5", "\n")
  
  updateTextInput(session, "report_accession", value = as.character(dpmsr_set$x$accession_report_list))
  
  updateTextInput(session, "adh_list", value = as.character(dpmsr_set$x$adh_list ))
  updateTextInput(session, "bait_list", value = as.character(dpmsr_set$x$bait_list ))
  updateTextInput(session, "avidin_list", value = as.character(dpmsr_set$x$avidin_list ))
  updateTextInput(session, "carbox_list", value = as.character(dpmsr_set$x$carbox_list ))
  updateTextInput(session, "bira_list", value = as.character(dpmsr_set$x$bira_list ))
  updateTextInput(session, "protein1_list", value = as.character(dpmsr_set$x$protein1_list ))
  updateTextInput(session, "protein2_list", value = as.character(dpmsr_set$x$protein2_list ))
  updateTextInput(session, "protein3_list", value = as.character(dpmsr_set$x$protein3_list ))
  updateTextInput(session, "protein4_list", value = as.character(dpmsr_set$x$protein4_list ))
  updateTextInput(session, "protein_spike_list", value = as.character(dpmsr_set$x$qc_spike_id))
  
  cat(file=stderr(), "update_widget_startup...end", "\n")
}

#-Filter--------------------------------------------------------------------- 
update_widget_filter <- function(session, input, output){
  cat(file=stderr(), "update_widget_filter...", "\n")
  
  updateNumericInput(session, "minimum_measured_all", value = as.numeric(dpmsr_set$x$minimum_measured_all))
  
  updateCheckboxInput(session, "checkbox_require_x", value = as.logical(dpmsr_set$x$require_x))
  
  updateNumericInput(session, "require_x_cutoff", value = as.numeric(dpmsr_set$x$require_x_cutoff))
  updateCheckboxInput(session, "checkbox_filtercv", value = as.logical(dpmsr_set$x$filter_cv))
  
  
  updateSelectInput(session, "text_filtercvgroup", choices = unique(dpmsr_set$design$Group)  )
  updateSelectInput(session, "text_filtercvgroup", selected = as.character(dpmsr_set$x$filter_cv_group))
  
  updateNumericInput(session, "number_filtercvcutoff", value = as.numeric(dpmsr_set$x$filter_cv_cutoff))


  updateCheckboxInput(session, "checkbox_report_ptm", value = as.logical(dpmsr_set$x$peptide_report_ptm ))
  
  updateTextInput(session, "peptide_report_grep", value = as.character(dpmsr_set$x$peptide_report_grep ))
  
}
  
  
  #-Norm---------------------------------------------------------------------
update_widget_norm <- function(session, input, output){
  cat(file=stderr(), "update_widget_norm...", "\n")
  
  if (as.logical(dpmsr_set$x$sl)) {sl_norm <- 1}else{sl_norm<-0}
  updateCheckboxInput(session, "checkbox_n1", value = sl_norm)
  
  if (as.logical(dpmsr_set$x$tmm)) {tmm_norm <- 1}else{tmm_norm<-0}
  updateCheckboxInput(session, "checkbox_n2", value = tmm_norm)
  
  if (as.logical(dpmsr_set$x$sltmm)) {sltmm_norm <- 1}else{sltmm_norm<-0}
  updateCheckboxInput(session, "checkbox_n3", value = sltmm_norm)
  
  if (as.logical(dpmsr_set$x$directlfq)) {directlfq_norm <- 1}else{directlfq_norm<-0}
  updateCheckboxInput(session, "checkbox_n13", value = directlfq_norm)
  
  if (as.logical(dpmsr_set$x$quantile)) {quantile_norm <- 1}else{quantile_norm<-0}
  updateCheckboxInput(session, "checkbox_n4", value = quantile_norm)
  
  if (as.logical(dpmsr_set$x$lr)) {lr_norm <- 1}else{lr_norm<-0}
  updateCheckboxInput(session, "checkbox_n5", value = lr_norm)
  
  if (as.logical(dpmsr_set$x$loess)) {loess_norm <- 1}else{loess_norm<-0}
  updateCheckboxInput(session, "checkbox_n6", value = loess_norm)
  
  if (as.logical(dpmsr_set$x$vsn)) {vsn_norm <- 1}else{vsn_norm<-0}
  updateCheckboxInput(session, "checkbox_n7", value = vsn_norm)
  
  if (as.logical(dpmsr_set$x$ti)) {ti_norm <- 1}else{ti_norm<-0}
  updateCheckboxInput(session, "checkbox_n8", value = ti_norm)
  
  if (as.logical(dpmsr_set$x$mi)) {mi_norm <- 1}else{mi_norm<-0}
  updateCheckboxInput(session, "checkbox_n9", value = mi_norm)
  
  if (as.logical(dpmsr_set$x$ai)) {ai_norm <- 1}else{ai_norm<-0}
  updateCheckboxInput(session, "checkbox_n10", value = ai_norm)
  
  if (as.logical(dpmsr_set$x$protein)) {protein_norm <- 1}else{protein_norm<-0}
  updateCheckboxInput(session, "checkbox_n11", value = protein_norm)
  
  updateTextInput(session, "protein_norm_list", value = as.character(dpmsr_set$x$protein_norm_list ))
  
  #TMT SPQC----------------------------------------------------------------------------------
  
  updateNumericInput(session, "tmt_sets", value = as.numeric(dpmsr_set$x$tmt_spqc_sets))
  updateNumericInput(session, "tmt_channels", value = as.numeric(dpmsr_set$x$tmt_spqc_channels))
  updateNumericInput(session, "tmt_filter_sd", value = as.numeric(dpmsr_set$x$tmt_spqc_sets))
  
  if (as.logical(dpmsr_set$x$tmt_filter_peptide)) {tmt_filter <- 1}else{tmt_filter<-0}
  updateCheckboxInput(session, "checkbox_tmt_filter", value = tmt_filter)
}

#-Impute---------------------------------------------------------------------
  update_widget_impute <- function(session, input, output){ 
    cat(file=stderr(), "update_widget_impute...", "\n")
      
    if(dpmsr_set$x$impute_method == "Duke"){impute_method <- 1}
      else if(dpmsr_set$x$impute_method == "Floor"){impute_method <- 2}
      else if(dpmsr_set$x$impute_method == "Minimum"){impute_method <- 3}
      else if(dpmsr_set$x$impute_method == "Average/Group"){impute_method <- 4}
      else if(dpmsr_set$x$impute_method == "KNN"){impute_method <- 5}
      else if(dpmsr_set$x$impute_method == "LocalLeastSquares"){impute_method <- 6}
      else if(dpmsr_set$x$impute_method == "MLE"){impute_method <- 7}
      else if(dpmsr_set$x$impute_method == "BottomX"){impute_method <- 8}
      else if(dpmsr_set$x$impute_method == "Average/Global"){impute_method <- 9}
    updateRadioButtons(session, "radio_impute", selected = impute_method )
    
    if (as.logical(dpmsr_set$x$peptide_ptm_impute) ) {ptmreport <- 1}else{ptmreport<-0}
    updateCheckboxInput(session, "checkbox_impute_ptm", value = ptmreport) 
    
    updateNumericInput(session, "impute_floor", value = as.numeric(dpmsr_set$x$area_floor))
    
    updateNumericInput(session, "bottom_x", value = as.numeric(dpmsr_set$x$bottom_x))
    if(class(dpmsr_set$x$TMT_SPQC_bottom_x)=="NULL") {dpmsr_set$x$TMT_SPQC_bottom_x <<- dpmsr_set$x$bottom_x}
    updateNumericInput(session, "TMT_SPQC_bottom_x", value = as.numeric(dpmsr_set$x$TMT_SPQC_bottom_x))
    
    #if (as.logical(dpmsr_set$x$duke_misaligned)) {misaligned <- 1}else{misaligned<-0}
    updateCheckboxInput(session, "checkbox_misaligned", value = as.logical(dpmsr_set$x$duke_misaligned))
    
    updateNumericInput(session, "missing_cutoff", value = as.numeric(dpmsr_set$x$missing_cutoff))
    updateNumericInput(session, "misaligned_cutoff", value = as.numeric(dpmsr_set$x$misaligned_cutoff))
    updateNumericInput(session, "intensity_cutoff_mean_sd", value = as.numeric(dpmsr_set$x$int_cutoff_sd))
    
  }
  
  #-Rollup---------------------------------------------------------------------
  update_widget_rollup <- function(session, input, output){ 
    cat(file=stderr(), "update_widget_rollup...", "\n")
    
    if(dpmsr_set$x$rollup_method == "Sum"){rollup_method <- 1}
    else if(dpmsr_set$x$rollup_method == "Median"){rollup_method <- 2}
    else if(dpmsr_set$x$rollup_method == "Median_Polish"){rollup_method <- 3}
    else if(dpmsr_set$x$rollup_method == "Mean"){rollup_method <- 4}
    else if(dpmsr_set$x$rollup_method == "IQ_MaxLFQ"){rollup_method <- 5}
    else if(dpmsr_set$x$rollup_method == "TopN"){rollup_method <- 6}
    updateRadioButtons(session, "radio_rollup", selected = rollup_method )
    
    updateSelectInput(session, "rollup_topN_count", selected = as.numeric(dpmsr_set$y$rollup_topN_count))
    
  }
  
  
  
    
  #---Plots-------------------------------------------------------------------------------
update_widget_post_processing <- function(session, input, output){ 
  cat(file=stderr(), "update_widget_post_processing...", "\n")
  
  if (as.logical(dpmsr_set$x$sl)) {sl_norm <- 1}else{sl_norm<-0}
  updateCheckboxInput(session, "checkbox_nc1", value = sl_norm)
  
  if (as.logical(dpmsr_set$x$tmm)) {tmm_norm <- 1}else{tmm_norm<-0}
  updateCheckboxInput(session, "checkbox_nc2", value = tmm_norm)
  
  if (as.logical(dpmsr_set$x$sltmm)) {sltmm_norm <- 1}else{sltmm_norm<-0}
  updateCheckboxInput(session, "checkbox_nc3", value = sltmm_norm)
  
  if (as.logical(dpmsr_set$x$directlfq)) {directlfq_norm <- 1}else{directlfq_norm<-0}
  updateCheckboxInput(session, "checkbox_nc13", value = directlfq_norm)
  
  if (as.logical(dpmsr_set$x$quantile)) {quantile_norm <- 1}else{quantile_norm<-0}
  updateCheckboxInput(session, "checkbox_nc4", value = quantile_norm)
  
  if (as.logical(dpmsr_set$x$lr)) {lr_norm <- 1}else{lr_norm<-0}
  updateCheckboxInput(session, "checkbox_nc5", value = lr_norm)
  
  if (as.logical(dpmsr_set$x$loess)) {loess_norm <- 1}else{loess_norm<-0}
  updateCheckboxInput(session, "checkbox_nc6", value = loess_norm)
  
  if (as.logical(dpmsr_set$x$vsn)) {vsn_norm <- 1}else{vsn_norm<-0}
  updateCheckboxInput(session, "checkbox_nc7", value = vsn_norm)
  
  if (as.logical(dpmsr_set$x$ti)) {ti_norm <- 1}else{ti_norm<-0}
  updateCheckboxInput(session, "checkbox_nc8", value = ti_norm)
  
  if (as.logical(dpmsr_set$x$mi)) {mi_norm <- 1}else{mi_norm<-0}
  updateCheckboxInput(session, "checkbox_nc9", value = mi_norm)
  
  if (as.logical(dpmsr_set$x$ai)) {ai_norm <- 1}else{ai_norm<-0}
  updateCheckboxInput(session, "checkbox_nc10", value = ai_norm)
  
  if (as.logical(dpmsr_set$x$protein)) {protein_norm <- 1}else{protein_norm<-0}
  updateCheckboxInput(session, "checkbox_nc11", value = protein_norm)
  
  #----------------------------------------------------------------------
  updateSelectInput(session, "norm_type", choices = names(dpmsr_set$data$final) , selected= "impute")
  updateSelectInput(session, "volcano_select", choices = names(dpmsr_set$data$final), selected = "impute")
  updateSelectInput(session, "protein_select", choices = dpmsr_set$y$protein_list , selected= "ADH")
  updateSelectInput(session, "select_oneprotein_norm", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_onepeptide_norm", choices = names(dpmsr_set$data$final), selected= "impute")
  
#----------------------------------------------------------------------

 # updateSelectInput(session, "select_data_comp", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_final_data_stats", choices = names(dpmsr_set$data$final), selected= "impute")
  
  updateTextInput(session, "final_stats_name", value = str_c("Final_", dpmsr_set$data$stats$final_comp,  "_stats.xlsx"))
  
  # update what was previously selected
  try((
    if(as.numeric(dpmsr_set$y$stats$comp_number) > 0){
      cat(file=stderr(), "Update mva...", "\n")
      updateSelectInput(session, "comp_number", selected = as.numeric(dpmsr_set$y$stats$comp_number))
      updateSelectInput(session, "select_final_data_stats", selected = dpmsr_set$data$stats$final_comp)
      updatePickerInput(session, "comp_spqc", selected= dpmsr_set$y$stats$comp_spqc)
      updatePickerInput(session, "stats_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name)
      updateSelectInput(session, "stats_oneprotein_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name, selected = dpmsr_set$y$stats$groups$comp_name[1])
      updateSelectInput(session, "stats_onepeptide_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name, selected = dpmsr_set$y$stats$groups$comp_name[1])
      updateSelectInput(session, "stats_select_data_comp", choices = dpmsr_set$y$stats$groups$comp_name)
      updateSelectInput(session, "select_data_comp_motif", choices = dpmsr_set$y$stats$groups$comp_name)
      updateSelectInput(session, "select_data_comp_wiki", choices = dpmsr_set$y$stats$groups$comp_name)
      updateSelectInput(session, "select_data_comp_profile", choices = dpmsr_set$y$stats$groups$comp_name)
      updateSelectInput(session, "select_data_comp_go", choices = dpmsr_set$y$stats$groups$comp_name)
      updateSelectInput(session, "select_data_comp_string", choices = dpmsr_set$y$stats$groups$comp_name)
      updateSelectInput(session, "select_data_comp_string_enrich", choices = dpmsr_set$y$stats$groups$comp_name)
      updateCheckboxInput(session, "stats_spqc_cv_filter", value = dpmsr_set$y$stats$stats_spqc_cv_filter)
      updateNumericInput(session, "stats_spqc_cv_filter_factor", value = dpmsr_set$y$stats$stats_spqc_cv_filter_factor)
      updateCheckboxInput(session, "stats_comp_cv_filter", value = dpmsr_set$y$stats$stats_spqc_cv_filter)
      updateNumericInput(session, "stats_comp_cv_filter_factor", value = dpmsr_set$y$stats$stats_spqc_cv_filter_factor)
      updateCheckboxInput(session, "stats_peptide_minimum", value = dpmsr_set$y$stats$stats_peptide_minimum)
      updateNumericInput(session, "stats_peptide_minimum_factor", value = dpmsr_set$y$stats$stats_peptide_minimum_factor)
      updatePickerInput(session, "comp_spqc", selected = dpmsr_set$y$stats$comp_spqc )
      for (i in 1:nrow(dpmsr_set$y$stats$groups)){
        updatePickerInput(session, str_c("comp_", i, "N"), selected= dpmsr_set$y$stats[[str_c("comp",i,"_N")]] )
        updatePickerInput(session, str_c("comp_", i, "D"), selected= dpmsr_set$y$stats[[str_c("comp",i,"_D")]] )
      }
    }
      ), silent=TRUE)
  
}


#-Stats---------------------------------------------------------------------

update_widget_stats <- function(session, input, output){
  cat(file=stderr(), "update_widget_stats...", "\n")
  
  updateCheckboxInput(session, "pair_comp", value = as.logical(dpmsr_set$x$pair_comp))
  updateNumericInput(session, "pvalue_cutoff", value = as.numeric(dpmsr_set$x$pvalue_cutoff  ))
  updateNumericInput(session, "foldchange_cutoff", value = as.numeric(dpmsr_set$x$foldchange_cutoff ))
  updateNumericInput(session, "missing_factor", value = as.numeric(dpmsr_set$x$missing_factor ))
  updateSelectInput(session, "select_final_data_stats", selected = dpmsr_set$data$stats$final_comp)
  updateSelectInput(session, "select_organism", selected = dpmsr_set$y$organism)
  
  updateCheckboxInput(session, "peptide_missing_filter", value = as.logical(dpmsr_set$y$peptide_missing_filter))
  updateNumericInput(session, "peptide_missing_factor", value = as.numeric(dpmsr_set$y$peptide_missing_factor  ))
  
  updateCheckboxInput(session, "peptide_cv_filter", value = as.logical(dpmsr_set$y$peptide_cv_filter))
  updateNumericInput(session, "peptide_cv_factor", value = as.numeric(dpmsr_set$y$peptide_cv_factor  ))
  
  update_stat_choices(session, input, output)  
}

#----------------------------------------------------------------------
update_stat_choices <- function(session, input, output){
  cat(file=stderr(), "update_widget_stat_choices...", "\n")
  #updates choice list only (not what was selected)
  updatePickerInput(session, "comp_1N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_2N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_3N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_4N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_5N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_6N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_7N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_8N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_9N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_10N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_11N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_12N", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_1D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_2D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_3D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_4D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_5D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_6D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_7D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_8D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_9D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_10D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_11D", choices = dpmsr_set$y$uniquegroups, selected= "-")
  updatePickerInput(session, "comp_12D", choices = dpmsr_set$y$uniquegroups, selected= "-")  
  updatePickerInput(session, "comp_spqc", choices = c(dpmsr_set$y$uniquegroups, "None"), selected= "-")  
  
  updateTextInput(session, "comp1_name", value = dpmsr_set$y$stats$comp1_name)
  updateTextInput(session, "comp2_name", value = dpmsr_set$y$stats$comp2_name)
  updateTextInput(session, "comp3_name", value = dpmsr_set$y$stats$comp3_name)
  updateTextInput(session, "comp4_name", value = dpmsr_set$y$stats$comp4_name)
  updateTextInput(session, "comp5_name", value = dpmsr_set$y$stats$comp5_name)
  updateTextInput(session, "comp6_name", value = dpmsr_set$y$stats$comp6_name)
  updateTextInput(session, "comp7_name", value = dpmsr_set$y$stats$comp7_name)
  updateTextInput(session, "comp8_name", value = dpmsr_set$y$stats$comp8_name)
  updateTextInput(session, "comp9_name", value = dpmsr_set$y$stats$comp9_name)
  updateTextInput(session, "comp10_name", value = dpmsr_set$y$stats$comp10_name)  
  updateTextInput(session, "comp11_name", value = dpmsr_set$y$stats$comp11_name)
  updateTextInput(session, "comp12_name", value = dpmsr_set$y$stats$comp12_name)
  
  }

#--All-----------------------------------------------------------------
update_widget_all <- function(session, input, output){ 
  cat(file=stderr(), "update_widget_all...", "\n")
  update_widget_startup(session, input, output)
  update_widget_filter(session, input, output)
  update_widget_norm(session, input, output)
  update_widget_impute(session, input, output)
  update_widget_rollup(session, input, output)
  update_widget_stats(session, input, output)
  update_widget_post_processing(session, input, output)
  create_design_table(session, input, output)
}


#--All-----------------------------------------------------------------
update_comparisons <- function(session, input, output){ 
  cat(file=stderr(), "update_widget_comparisons...", "\n")
  updateSelectInput(session, "select_data_comp_motif", choices = dpmsr_set$y$stats$groups$comp_name)
  updateSelectInput(session, "select_data_comp_wiki", choices = dpmsr_set$y$stats$groups$comp_name)
  updateSelectInput(session, "select_data_comp_profile", choices = dpmsr_set$y$stats$groups$comp_name)
  updateSelectInput(session, "select_data_comp_go", choices = dpmsr_set$y$stats$groups$comp_name)
  updateSelectInput(session, "select_data_comp_string", choices = dpmsr_set$y$stats$groups$comp_name)
  updateSelectInput(session, "select_data_comp_string_enrich", choices = dpmsr_set$y$stats$groups$comp_name)
  updatePickerInput(session, "stats_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name)
  updateSelectInput(session, "stats_oneprotein_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name, selected = dpmsr_set$y$stats$groups$comp_name[1])
  updateSelectInput(session, "stats_onepeptide_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name, selected = dpmsr_set$y$stats$groups$comp_name[1])
  updateSelectInput(session, "stats_select_data_comp", choices = dpmsr_set$y$stats$groups$comp_name)
}





