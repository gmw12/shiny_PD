
update_shiny_defaults <- function(session){
  updateTextInput(session, "fileprefix", value = as.character(dpmsr_set$x$file_prefix))
  
  if(dpmsr_set$x$raw_data_input == "Protein_Peptide"){rdi <- 1}
    else if(dpmsr_set$x$raw_data_input == "Protein"){rdi <- 2}
    else if(dpmsr_set$x$raw_data_input == "Peptide"){rdi <- 3}
    else if(dpmsr_set$x$raw_data_input == "PSM_FDR"){rdi <- 4}
  updateRadioButtons(session, "radio_input", selected = rdi )
  

  updateSelectInput(session, "razor", selected = dpmsr_set$x$peptides_to_use )

  
  if(dpmsr_set$x$final_data_output == "Protein"){fdo <- 1}
    else if(dpmsr_set$x$final_data_output == "Peptide"){fdo <- 2}
  updateRadioButtons(session, "radio_output", selected = fdo )
  
  if (as.logical(dpmsr_set$x$peptide_isoform)) {isoform <- 1}else{isoform<-0}
  updateCheckboxInput(session, "checkbox_isoform", value = isoform) 
  
  if (as.logical(dpmsr_set$x$tmt_spqc_norm)) {tmtnorm <- 1}else{tmtnorm<-0}
  updateCheckboxInput(session, "check_tmt", value = tmtnorm) 
  
  if (as.logical(dpmsr_set$x$peptide_ptm_out) ) {ptmreport <- 1}else{ptmreport<-0}
  updateCheckboxInput(session, "checkbox_out_ptm", value = ptmreport) 
  
  updateTextInput(session, "report_grep", value = as.character(dpmsr_set$x$peptide_report_grep))
  
  if (as.logical(dpmsr_set$x$accession_report_out) ) {accession_report <- 1}else{accession_report<-0}
  updateCheckboxInput(session, "checkbox_report_accession", value = accession_report) 
  
  updateTextInput(session, "report_accession", value = as.character(dpmsr_set$x$accession_report_list))
  #----------------------------------------------------------------------
  
  if (as.logical(dpmsr_set$x$require_2)) {require2 <- 1}else{require2<-0}
  updateCheckboxInput(session, "checkbox_require2", value = require2)

  if (as.logical(dpmsr_set$x$filter_cv)) {require_cv <- 1}else{require_cv<-0}
  updateCheckboxInput(session, "checkbox_filtercv", value = require_cv)

  updateTextInput(session, "text_filtercvgroup", value = as.character(dpmsr_set$x$filter_cv_group))
  updateNumericInput(session, "number_filtercvcutoff", value = as.numeric(dpmsr_set$x$filter_cv_cutoff))
  
  if (as.logical(dpmsr_set$x$peptide_ptm_norm) ) {ptmreport <- 1}else{ptmreport<-0}
  updateCheckboxInput(session, "checkbox_norm_ptm", value = ptmreport) 
  
  updateTextInput(session, "filter_grep", value = as.character(dpmsr_set$x$peptide_grep ))
  
  #----------------------------------------------------------------------
  
  if (as.logical(dpmsr_set$x$sl)) {sl_norm <- 1}else{sl_norm<-0}
  updateCheckboxInput(session, "checkbox_n1", value = sl_norm)
  
  if (as.logical(dpmsr_set$x$tmm)) {tmm_norm <- 1}else{tmm_norm<-0}
  updateCheckboxInput(session, "checkbox_n2", value = tmm_norm)
  
  if (as.logical(dpmsr_set$x$sltmm)) {sltmm_norm <- 1}else{sltmm_norm<-0}
  updateCheckboxInput(session, "checkbox_n3", value = sltmm_norm)
  
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
  
  #----------------------------------------------------------------------
  
  if(dpmsr_set$x$impute_method == "Duke"){impute_method <- 1}
    else if(dpmsr_set$x$impute_method == "Floor"){impute_method <- 2}
    else if(dpmsr_set$x$impute_method == "Minimum"){impute_method <- 3}
    else if(dpmsr_set$x$impute_method == "Aerage"){impute_method <- 4}
    else if(dpmsr_set$x$impute_method == "KNN"){impute_method <- 5}
    else if(dpmsr_set$x$impute_method == "LocalLeastSquares"){impute_method <- 6}
    else if(dpmsr_set$x$impute_method == "MLE"){impute_method <- 7}
    else if(dpmsr_set$x$impute_method == "Bottom5"){impute_method <- 8}
  updateRadioButtons(session, "radio_impute", selected = impute_method )
  
  if (as.logical(dpmsr_set$x$peptide_ptm_impute) ) {ptmreport <- 1}else{ptmreport<-0}
  updateCheckboxInput(session, "checkbox_impute_ptm", value = ptmreport) 
  
  updateNumericInput(session, "impute_floor", value = as.numeric(dpmsr_set$x$area_floor))
  
  if (as.logical(dpmsr_set$x$missing_50)) {missing_50 <- 1}else{missing_50<-0}
  updateCheckboxInput(session, "checkbox_missing_50", value = as.logical(dpmsr_set$x$missing_50))
  
  
  #----------------------------------------------------------------------------------
  
  if (as.logical(dpmsr_set$x$sl)) {sl_norm <- 1}else{sl_norm<-0}
  updateCheckboxInput(session, "checkbox_nc1", value = sl_norm)
  
  if (as.logical(dpmsr_set$x$tmm)) {tmm_norm <- 1}else{tmm_norm<-0}
  updateCheckboxInput(session, "checkbox_nc2", value = tmm_norm)
  
  if (as.logical(dpmsr_set$x$sltmm)) {sltmm_norm <- 1}else{sltmm_norm<-0}
  updateCheckboxInput(session, "checkbox_nc3", value = sltmm_norm)
  
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
  
  updateSliderInput(session, "comp_number", value=as.numeric(dpmsr_set$x$comp_number))
  
  updateSelectInput(session, "comp_1N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp1N)
  updateSelectInput(session, "comp_1D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp1D)
  
  updateSelectInput(session, "comp_2N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp2N)
  updateSelectInput(session, "comp_2D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp2D)
  
  updateSelectInput(session, "comp_3N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp3N)
  updateSelectInput(session, "comp_3D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp3D)
  
  updateSelectInput(session, "comp_4N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp4N)
  updateSelectInput(session, "comp_4D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp4D)
  
  updateSelectInput(session, "comp_5N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp5N)
  updateSelectInput(session, "comp_5D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp5D)
  
  updateSelectInput(session, "comp_6N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp6N)
  updateSelectInput(session, "comp_6D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp6D)
  
  updateCheckboxInput(session, "pair_comp", value = as.logical(dpmsr_set$x$pair_comp))
  updateNumericInput(session, "pvalue_cutoff", value = as.numeric(dpmsr_set$x$pvalue_cutoff  ))
  updateNumericInput(session, "foldchange_cutoff", value = as.numeric(dpmsr_set$x$foldchange_cutoff ))

  #----------------------------------------------------------------------
  updateTextInput(session, "adh_list", value = as.character(dpmsr_set$x$adh_list ))
  updateTextInput(session, "bait_list", value = as.character(dpmsr_set$x$bait_list ))
  updateTextInput(session, "avidin_list", value = as.character(dpmsr_set$x$avidin_list ))
  updateTextInput(session, "carbox_list", value = as.character(dpmsr_set$x$carbox_list ))
  updateTextInput(session, "bira_list", value = as.character(dpmsr_set$x$bira_list ))
  updateTextInput(session, "protein1_list", value = as.character(dpmsr_set$x$protein1_list ))
  updateTextInput(session, "protein2_list", value = as.character(dpmsr_set$x$protein2_list ))
  updateTextInput(session, "protein3_list", value = as.character(dpmsr_set$x$protein3_list ))
  updateTextInput(session, "protein4_list", value = as.character(dpmsr_set$x$protein4_list ))

  updateSelectInput(session, "select_data_comp", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  
  updateSelectInput(session, "select_data_comp_motif", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_motif", choices = names(dpmsr_set$data$final), selected= "impute")
  
  updateSelectInput(session, "select_data_comp_wiki", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_wiki", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_profile", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_profile", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_go", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_go", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_string", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_string", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_string_enrich", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  #updateSelectInput(session, "select_final_data_string_enrich", choices = names(dpmsr_set$data$final), selected= "impute")
}



#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

selection_updates_new_start <- function(session, input, output){
  
  #update pulldown for proteinQC plots
  updateSelectInput(session, "norm_type", choices = names(dpmsr_set$data$final) , selected= "impute")
  updateSelectInput(session, "volcano_select", choices = names(dpmsr_set$data$final), selected = "impute")
  
  updateSelectInput(session, "protein_select", choices = dpmsr_set$y$protein_list , selected= "ADH")
  
  updateSelectInput(session, "select_data_comp", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_motif", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_motif", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_wiki", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_wiki", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_profile", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_profile", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_go", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_go", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_string", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_string", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_string_enrich", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  #updateSelectInput(session, "select_final_data_string_enrich", choices = names(dpmsr_set$data$final), selected= "impute")
  
  updateSelectInput(session, "select_oneprotein_norm", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_onepeptide_norm", choices = names(dpmsr_set$data$final), selected= "impute")
  
  
  updateSliderInput(session, "comp_number", value=as.numeric(dpmsr_set$x$comp_number))
  
  updateSelectInput(session, "comp_1N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp1N)
  updateSelectInput(session, "comp_1D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp1D)
  
  updateSelectInput(session, "comp_2N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp2N)
  updateSelectInput(session, "comp_2D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp2D)
  
  updateSelectInput(session, "comp_3N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp3N)
  updateSelectInput(session, "comp_3D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp3D)
  
  updateSelectInput(session, "comp_4N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp4N)
  updateSelectInput(session, "comp_4D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp4D)
  
  updateSelectInput(session, "comp_5N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp5N)
  updateSelectInput(session, "comp_5D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp5D)
  
  updateSelectInput(session, "comp_6N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp6N)
  updateSelectInput(session, "comp_6D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp6D)
  
  updateCheckboxInput(session, "pair_comp", value = as.logical(dpmsr_set$x$pair_comp))
  updateNumericInput(session, "pvalue_cutoff", value = as.numeric(dpmsr_set$x$pvalue_cutoff  ))
  updateNumericInput(session, "foldchange_cutoff", value = as.numeric(dpmsr_set$x$foldchange_cutoff ))
  
  
  
  
  
  
  
}


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

selection_update_after_stats <- function(session, input, output){
  #update pulldown for proteinQC plots
  updateSelectInput(session, "norm_type", choices = names(dpmsr_set$data$final) , selected= "impute")
  updateSelectInput(session, "volcano_select", choices = names(dpmsr_set$data$final), selected = "impute")
  updateSelectInput(session, "select_data_comp", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_motif", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_motif", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_wiki", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_wiki", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_profile", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_profile", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_go", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_go", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_string", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  updateSelectInput(session, "select_final_data_string", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_data_comp_string_enrich", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
  #updateSelectInput(session, "select_final_data_string_enrich", choices = names(dpmsr_set$data$final), selected= "impute")
  
  updateSelectInput(session, "select_oneprotein_norm", choices = names(dpmsr_set$data$final), selected= "impute")
  updateSelectInput(session, "select_onepeptide_norm", choices = names(dpmsr_set$data$final), selected= "impute")
  
  
  
  updateTabsetPanel(session, "nlp1", selected = "tp8") 
  
  updateSliderInput(session, "comp_number", value=as.numeric(dpmsr_set$x$comp_number))
  
  updateSelectInput(session, "comp_1N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp1N)
  updateSelectInput(session, "comp_1D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp1D)
  
  updateSelectInput(session, "comp_2N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp2N)
  updateSelectInput(session, "comp_2D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp2D)
  
  updateSelectInput(session, "comp_3N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp3N)
  updateSelectInput(session, "comp_3D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp3D)
  
  updateSelectInput(session, "comp_4N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp4N)
  updateSelectInput(session, "comp_4D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp4D)
  
  updateSelectInput(session, "comp_5N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp5N)
  updateSelectInput(session, "comp_5D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp5D)
  
  updateSelectInput(session, "comp_6N", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp6N)
  updateSelectInput(session, "comp_6D", choices = dpmsr_set$y$sample_groups$Group, selected= dpmsr_set$x$comp6D)
  
  updateCheckboxInput(session, "pair_comp", value = as.logical(dpmsr_set$x$pair_comp))
  updateNumericInput(session, "pvalue_cutoff", value = as.numeric(dpmsr_set$x$pvalue_cutoff  ))
  updateNumericInput(session, "foldchange_cutoff", value = as.numeric(dpmsr_set$x$foldchange_cutoff ))
  
}






