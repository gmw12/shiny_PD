#----------------------------------------------------------------------------------------
preprocess_order <- function(){
  # data stat info
  cat(file = stderr(), "preprocess_order()...start", "\n")
  
  if (dpmsr_set$x$raw_data_input == "Protein") {
    cat(file = stderr(), "preprocess_order()...1.1", "\n")
    dpmsr_set$y$total_columns <<- ncol(dpmsr_set$data$data_protein_start)
  }else if (dpmsr_set$x$raw_data_input == "Peptide") {
    cat(file = stderr(), "preprocess_order()...1.2", "\n")
    dpmsr_set$y$total_columns <<- ncol(dpmsr_set$data$data_peptide_start)
  }else if (dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") {
    cat(file = stderr(), "preprocess_order()...1.3", "\n")
    dpmsr_set$y$total_columns <<- ncol(dpmsr_set$data$data_precursor_start)
  }
  
  dpmsr_set$y$info_columns <<- dpmsr_set$y$total_columns - dpmsr_set$y$sample_number
  cat(file = stderr(), str_c("dpmsr_set$y$info_columns=", dpmsr_set$y$info_columns ), "\n")
  
  if ((dpmsr_set$x$raw_data_input == "Protein_Peptide" || dpmsr_set$x$raw_data_input == "Peptide") 
      && dpmsr_set$x$final_data_output == "Protein") {
    cat(file = stderr(), "preprocess_order()...2", "\n")
    dpmsr_set$y$info_columns_final <<- 3  #reference collapse peptide function in transform
  }else{
    dpmsr_set$y$info_columns_final <<- dpmsr_set$y$info_columns
  }
  
  if (dpmsr_set$x$raw_data_input == "Protein_Peptide" || dpmsr_set$x$raw_data_input == "Peptide") {
    cat(file = stderr(), "preprocess_order()...3", "\n")
    cat(file = stderr(), "preprocess_order()...3", "\n")
    df <- order_columns(dpmsr_set$data$data_peptide_start)
    colnames(df)[(dpmsr_set$y$info_columns + 1):ncol(df)] <- dpmsr_set$design$Header1
    dpmsr_set$data$data_peptide_start <<- df
  } 
  
  if (dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") {
    cat(file = stderr(), "preprocess_order()...4", "\n")
    df <- order_columns(dpmsr_set$data$data_precursor_start)
    colnames(df)[(dpmsr_set$y$info_columns + 1):ncol(df)] <- dpmsr_set$design$Header1
    dpmsr_set$data$data_precursor_start <<- df
  } 
  
  if (dpmsr_set$x$raw_data_input == "Protein") {
    cat(file = stderr(), "preprocess_order()...5", "\n")
    df <- order_columns(dpmsr_set$data$data_protein_start)
    colnames(df)[(dpmsr_set$y$info_columns + 1):ncol(df)] <- dpmsr_set$design$Header1
    dpmsr_set$data$data_protein_start <<- df
  } 
  
  if (as.logical(dpmsr_set$x$peptide_isoform)) {
    cat(file = stderr(), "preprocess_order()...6", "\n")
    df <- order_columns(dpmsr_set$data$data_peptide_isoform_start)
    colnames(df)[(dpmsr_set$y$info_columns + 1):ncol(df)] <- dpmsr_set$design$Header1
    dpmsr_set$data$data_peptide_isoform_start <<- df
  } 
  
  cat(file = stderr(), "preprocess_order()...end", "\n")
}

#----------------------------------------------------------------------------------------
preprocess_filter <- function(session, input, output){
  # data stat info

  if (dpmsr_set$x$raw_data_input == "Protein_Peptide" || dpmsr_set$x$raw_data_input == "Peptide") {
    cat(file = stderr(), "preprocess filter protein_peptide or peptide...", "\n")
    dpmsr_set$data$data_peptide <<- filter_data(session, input, output, dpmsr_set$data$data_peptide_start)
  } 
  
  if (dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") {
    cat(file = stderr(), "preprocess filter precursor...", "\n")
    dpmsr_set$data$data_precursor <<- filter_data(session, input, output, dpmsr_set$data$data_precursor_start)
  } 
  
  if (dpmsr_set$x$raw_data_input == "Protein") {
    cat(file = stderr(), "preprocess filter protein...", "\n")
    dpmsr_set$data$data_protein <<- filter_data(session, input, output, dpmsr_set$data$data_protein_start)
  } 
  
  if (as.logical(dpmsr_set$x$peptide_isoform)) {
    cat(file = stderr(), "preprocess filter peptide_isoform...", "\n")
    dpmsr_set$data$data_peptide_isoform <<- filter_data(session, input, output, dpmsr_set$data$data_peptide_isoform_start)
  } 
  
}

#----------------------------------------------------------------------------------------
# Rearrange columns if raw data is psm, PD does not organize
order_columns <- function(df){
  cat(file = stderr(), "order columns start", "\n")
  info_columns <- ncol(df) - dpmsr_set$y$sample_number
  annotate_df <- df[, 1:info_columns]
  df <- df[, (info_columns + 1):ncol(df)]
  df <- df[, (dpmsr_set$design$PD_Order)]
  #make sure data is numeric
  df <- mutate_all(df, function(x) as.numeric(as.character(x)))
  df <- cbind(annotate_df, df)
  cat(file = stderr(), "order columns end", "\n")
  return(df)
}


#----------------------------------------------------------------------------------------
check_sample_id <- function() {
  cat(file = stderr(), "check_sample_id...1", "\n")
  sample_ids <- dpmsr_set$design[order(dpmsr_set$design$PD_Order),]
  sample_ids <- sample_ids$ID
  sample_ids <- gsub("_.+", "", sample_ids)
  
  if (dpmsr_set$x$data_source == "SP") {
    if (dpmsr_set$x$raw_data_input == "Protein") {
      test_data <- colnames(dpmsr_set$data$data_raw_protein %>% dplyr::select(contains("Quantity")))
    }else if (dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") {
      test_data <- colnames(dpmsr_set$data$data_raw_precursor %>% dplyr::select(contains("Quantity")))
    }else {
      test_data <- colnames(dpmsr_set$data$data_raw_peptide %>% dplyr::select(contains("Quantity")))
    }
  }else{
    test_data <- colnames(dpmsr_set$data$data_raw_peptide %>% dplyr::select(contains("Abundance.F")))
  }
  
    
  cat(file = stderr(), "check_sample_id...2", "\n")
  for (i in 1:length(sample_ids)) {
    confirm_id <- grepl(sample_ids[i], test_data[i])
    if (!confirm_id) {
      cat(file = stderr(), str_c("count=", i, "  ", sample_ids[i], " vs ", test_data[i]), "\n")
      shinyalert("Oops!", "Sample ID order does not match sample list!", type = "error")
    }
  }
  
  cat(file = stderr(), "check_sample_id...end", "\n")
  }

#----------------------------------------------------------------------------------------
garbage_cleanup <- function(){
  # reduce size of objects in dpmsr_set
  cat(file = stderr(), "garbage cleanup...", "\n")
  rm(list=ls())
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
}

