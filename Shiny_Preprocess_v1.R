#----------------------------------------------------------------------------------------
preprocess_order <- function(){
  # data stat info
  cat(file=stderr(), "preprocess_order()...1", "\n")
  
  if(dpmsr_set$x$raw_data_input!="Protein"){
    dpmsr_set$y$total_columns <<- ncol(dpmsr_set$data$data_peptide_start)
  }else{
    dpmsr_set$y$total_columns <<- ncol(dpmsr_set$data$data_protein_start)
  }
  
  dpmsr_set$y$info_columns <<- dpmsr_set$y$total_columns - dpmsr_set$y$sample_number
  cat(file=stderr(), str_c("dpmsr_set$y$info_columns=", dpmsr_set$y$info_columns ), "\n")
  
  cat(file=stderr(), "preprocess_order()...2", "\n")
  if ((dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") 
      && dpmsr_set$x$final_data_output == "Protein"){
    dpmsr_set$y$info_columns_final <<- 3  #reference collapse peptide function in transform
  }else{
    dpmsr_set$y$info_columns_final <<- dpmsr_set$y$info_columns
  }
  
  cat(file=stderr(), "preprocess_order()...3", "\n")
  if (dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") {
    df <- order_columns(dpmsr_set$data$data_peptide_start)
    colnames(df)[(dpmsr_set$y$info_columns+1):ncol(df)] <- dpmsr_set$design$Header1
    dpmsr_set$data$data_peptide_start <<- df
  } 
  
  cat(file=stderr(), "preprocess_order()...4", "\n")
  if (dpmsr_set$x$raw_data_input=="Protein") {
    df <- order_columns(dpmsr_set$data$data_protein_start)
    colnames(df)[(dpmsr_set$y$info_columns+1):ncol(df)] <- dpmsr_set$design$Header1
    dpmsr_set$data$data_protein_start <<- df
  } 
  cat(file=stderr(), "preprocess_order()...5", "\n")
  if (as.logical(dpmsr_set$x$peptide_isoform)) {
    df <- order_columns(dpmsr_set$data$data_peptide_isoform_start)
    colnames(df)[(dpmsr_set$y$info_columns+1):ncol(df)] <- dpmsr_set$design$Header1
    dpmsr_set$data$data_peptide_isoform_start <<- df
  } 
  cat(file=stderr(), "preprocess_order()...end", "\n")
}

#----------------------------------------------------------------------------------------
preprocess_filter <- function(session, input, output){
  # data stat info

  if (dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") {
    dpmsr_set$data$data_peptide <<- filter_data(session, input, output, dpmsr_set$data$data_peptide_start)
  } 
  
  if (dpmsr_set$x$raw_data_input=="Protein") {
    dpmsr_set$data$data_protein <<- filter_data(session, input, output, dpmsr_set$data$data_protein_start)
  } 
  
  if (as.logical(dpmsr_set$x$peptide_isoform)) {
    dpmsr_set$data$data_peptide_isoform <<- filter_data(session, input, output, dpmsr_set$data$data_peptide_isoform_start)
  } 
  
}

#----------------------------------------------------------------------------------------
# Rearrange columns if raw data is psm, PD does not organize
order_columns <- function(df){
  annotate_df <- df[1:dpmsr_set$y$info_columns]
  df <- df[(dpmsr_set$y$info_columns+1):dpmsr_set$y$total_columns]
  df <- df[, (dpmsr_set$design$PD_Order)]
  df <- cbind(annotate_df, df)
  return(df)
}


#----------------------------------------------------------------------------------------
check_sample_id <- function() {
  cat(file=stderr(), "check_sample_id...1", "\n")
  sample_ids <- dpmsr_set$design[order(dpmsr_set$design$PD_Order),]
  sample_ids <- sample_ids$ID
  sample_ids <- gsub("_.+", "", sample_ids)
  
  if(dpmsr_set$x$raw_data_input == "Protein"){
    if(dpmsr_set$x$data_source == "SP"){
    test_data <- colnames(dpmsr_set$data$data_raw_protein %>% dplyr::select(contains("Quantity")))
    }else{
      test_data <- colnames(dpmsr_set$data$data_raw_peptide %>% dplyr::select(contains("Abundance.F")))
    }
  }else{
    test_data <- colnames(dpmsr_set$data$data_raw_peptide %>% dplyr::select(contains("Abundance.F")))
    }
  
  cat(file=stderr(), "check_sample_id...2", "\n")
  for (i in 1:length(sample_ids)){
    confirm_id <- grepl(sample_ids[i], test_data[i])
    if(!confirm_id){
      shinyalert("Oops!", "Sample ID order does not match sample list!", type = "error")
    }
  }
  
  cat(file=stderr(), "check_sample_id...end", "\n")
  }

#----------------------------------------------------------------------------------------
garbage_cleanup <- function(){
  # reduce size of objects in dpmsr_set
  cat(file=stderr(), "garbage cleanup...", "\n")
  rm(list=ls())
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
}

