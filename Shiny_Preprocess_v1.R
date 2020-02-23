#----------------------------------------------------------------------------------------
preprocess_data <- function(){
  # data stat info
  
  dpmsr_set$y$total_columns <<- ncol(dpmsr_set$data$data_peptide)
  dpmsr_set$y$info_columns <<- dpmsr_set$y$total_columns - dpmsr_set$y$sample_number
  
  if ((dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") 
          && dpmsr_set$x$final_data_output == "Protein"){
    dpmsr_set$y$info_columns_final <<- 3  #reference collapse peptide function in transform
  }else{
    dpmsr_set$y$info_columns_final <<- dpmsr_set$y$info_columns
    }

  if (dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") {
    df <- order_columns(dpmsr_set$data$data_peptide)
    colnames(df)[(dpmsr_set$y$info_columns+1):ncol(df)] <- dpmsr_set$design$Header1
    df <- filter_data(df)
    dpmsr_set$data$data_peptide <<- df
  } 
  
  if (dpmsr_set$x$raw_data_input=="Protein") {
    df <- order_columns(dpmsr_set$data$data_protein)
    colnames(df)[(dpmsr_set$y$info_columns+1):ncol(df)] <- dpmsr_set$design$Header1
    df <- filter_data(df)
    dpmsr_set$data$data_protein <<- df
  } 
  
  if (as.logical(dpmsr_set$x$peptide_isoform)) {
    df <- order_columns(dpmsr_set$data$data_peptide_isoform)
    colnames(df)[(dpmsr_set$y$info_columns+1):ncol(df)] <- dpmsr_set$design$Header1
    df <- filter_data(df)
    dpmsr_set$data$data_peptide_isoform <<- df
  } 
  
}


#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
# create column to flag and delete rows with no data or only one data value
#require row to contain at least one group that has at least 2 values
filter_data <- function(df){
  df[df==0] <- NA
  df$na_count <- apply(df[(dpmsr_set$y$info_columns+1):ncol(df)], 1, function(x) sum(!is.na(x)))
  df <- subset(df, df$na_count > 1)
  df <- df[1:dpmsr_set$y$total_columns]
  Simple_Excel(df, str_c(dpmsr_set$file$extra_prefix, "_raw_filter.xlsx", collapse = " "))
  if(as.logical(dpmsr_set$x$require_2)){
    for(i in 1:dpmsr_set$y$group_number) 
    {
      df$na_count <- apply(df[c((dpmsr_set$y$sample_groups$start[i]+dpmsr_set$y$info_columns):
                                  (dpmsr_set$y$sample_groups$end[i]+dpmsr_set$y$info_columns))], 1, function(x) sum(!is.na(x)))
      colnames(df)[colnames(df)=="na_count"] <- dpmsr_set$y$sample_groups$Group[i]
    }
    df$max <- apply(df[c(dpmsr_set$y$total_columns+1):(ncol(df))], 1, function(x) max(x))
    df <- subset(df, df$max > 1)
    df <- df[1:dpmsr_set$y$total_columns]
    Simple_Excel(df, str_c(dpmsr_set$file$extra_prefix, "_raw_filter_require2.xlsx", collapse = " "))
  }
  #optional filter by cv of specific sample group
  if(as.logical(dpmsr_set$x$filter_cv)){df<-filter_by_cv(df)}
  return(df)
}

#----------------------------------------------------------------------------------------

#filter by CV of a specific group (QCPool)
filter_by_cv <- function(df){
  start <- dpmsr_set$y$info_columns + dpmsr_set$y$sample_groups$start[dpmsr_set$y$sample_groups$Group == dpmsr_set$x$filter_cv_group]
  end <- dpmsr_set$y$info_columns + dpmsr_set$y$sample_groups$end[dpmsr_set$y$sample_groups$Group == dpmsr_set$x$filter_cv_group]
  df$filterCV <- percentCV_gw(df[start:end])
  df <- subset(df, df$filterCV < dpmsr_set$x$filter_cv_cutoff)
  df <- df[-ncol(df)]
  Simple_Excel(df, str_c(dpmsr_set$file$extra_prefix, "_raw_filter_byCV.xlsx", collapse = " "))
  return(df)
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

#these  need a place 
#annotate_df <- data_raw[1:info_columns]
#data_to_normalize <- data_raw[(info_columns+1):ncol(data_raw)]
#----------------------------------------------------------------------------------------
clean_headers <- function(col_headers){
  #colnames(forward_peptide)[(names(forward_peptide) == "Annotated Sequence")] <- "Sequence"
  #----- edit column headers
  col_headers <- str_replace(col_headers, "Protein FDR Confidence: Mascot", "Confidence")
  col_headers <- str_replace(col_headers, "Protein FDR Confidence: Combined", "Confidence")
  col_headers <- str_replace(col_headers, "Annotated Sequence", "Sequence")
  col_headers <- str_replace(col_headers, "Master Protein Accessions", "Accession")
  col_headers <- str_replace(col_headers," \\(by Search Engine\\): Mascot", "")
  col_headers <- str_replace(col_headers,"\\[", "")
  col_headers <- str_replace(col_headers,"\\]", "")
  col_headers <- c(col_headers, dpmsr_set$design$Header1)
  return(col_headers)
}








