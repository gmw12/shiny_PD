#---------------------------------------------------------------------------------------------------------

stat_prep_rollup <- function(session, input, output){ 
  
  cat(file = stderr(), "stat_prep_rollup function...", "\n")
  showModal(modalDialog("Stat prep...", footer = NULL)) 
  stat_prep()
  removeModal()
  cat(file = stderr(), "define final info columns", "\n")
  set_final_info_columns()
  cat(file = stderr(), "qc_apply", "\n")
  showModal(modalDialog("QC stats...", footer = NULL)) 
  qc_apply()
  removeModal()
  cat(file = stderr(), "qc - create_plots", "\n")
  showModal(modalDialog("QC plots...", footer = NULL)) 
  create_qc_plots()
  cat(file = stderr(), "qc - render_plots", "\n")
  qc_render(session, input, output)
  update_widget_post_processing(session, input, output)
  cat(file = stderr(), "qc - save_final_nostats", "\n")
  save_final_nostats()
  removeModal()
}

#--- collapse precursor to peptide-------------------------------------------------------------
collapse_precursor <- function(precursor_data, info_columns = 0, stats = FALSE) {
  cat(file = stderr(), "precursor rollup_sum triggered...", "\n")
  
  #precursor_data <- dpmsr_set$data$impute$sltmm
  #precursor_data <- df_raw
  precursor_data <- dpmsr_set$data$data_precursor_start
  
  #drop columns 
  precursor_data$PrecursorId <- NULL
  precursor_data$Detected_Imputed <- NULL
  
  if (info_columns == 0) {
    info_columns <- ncol(precursor_data) - dpmsr_set$y$sample_number
  }
  
  if (stats == TRUE) {
    info_columns = info_columns - 2
  }
  
  
  #save sample column number to reset info columns below
  sample_columns <- ncol(precursor_data) - info_columns
  
  columns = names(precursor_data)[1:info_columns]
  
  peptide_data <- precursor_data %>%
    group_by_at(vars(one_of(columns))) %>%
    summarise_all(list(sum))
  
  peptide_data <- data.frame(ungroup(peptide_data))
  
  if (nrow(peptide_data) == length(dpmsr_set$data$peptide_missing)) {
    peptide_data <- add_column(peptide_data, dpmsr_set$data$peptide_missing, .after = "Sequence")
    #col_pos <- grep("Sequence", colnames(peptide_data)) + 1
    names(peptide_data)[grep("Sequence", colnames(peptide_data)) + 1] <- "Detected_Imputed"
  } 
  
  cat(file = stderr(), "precursor rollup_sum triggered...end", "\n")
  return(peptide_data)
   
}



#--- collapse peptide to protein-------------------------------------------------------------
collapse_peptide <- function(peptide_data, info_columns=0, stats=FALSE, impute=FALSE){

  cat(file = stderr(), "starting collapse_peptide...", "\n")
  
  if (info_columns == 0) {
    info_columns <- ncol(peptide_data) - dpmsr_set$y$sample_number
  }
  
  #save sample column number to reset info columns below
  sample_columns <- ncol(peptide_data) - info_columns
  
  peptide_data <- collapse_peptide_setup(peptide_data, info_columns)
  #reset info_columns after function call (will reduce info columns to 3 or 4)
  info_columns <- ncol(peptide_data) - sample_columns
  
  if (!impute) {
    if (dpmsr_set$x$rollup_method == "Sum") {
      protein_data <- rollup_sum(peptide_data, info_columns)
    }else if (dpmsr_set$x$rollup_method == "Median") {
      protein_data <- rollup_median(peptide_data, info_columns)
    }else if (dpmsr_set$x$rollup_method == "Median_Polish") {
      protein_data <- rollup_median_polish(peptide_data, info_columns)
    }else if (dpmsr_set$x$rollup_method == "Mean") {
      protein_data <- rollup_mean(peptide_data, info_columns)
    }else if (dpmsr_set$x$rollup_method == "IQ_MaxLFQ") {
      protein_data <- rollup_maxlfq(peptide_data, info_columns)
    }else if (dpmsr_set$x$rollup_method == "TopN") {
      protein_data <- rollup_topn(peptide_data, info_columns)
    }else if (dpmsr_set$x$rollup_method == "DirectLFQ") {
      protein_data <- rollup_directlfq(peptide_data, info_columns)
    }
  }else{
    protein_data <- rollup_sum(peptide_data, info_columns)
  }
  
  protein_data <- data.frame(ungroup(protein_data))
  
  #add imputed column info
  if (!stats) {
    cat(file = stderr(), "collapse peptide to protein... 4", "\n")
    if ((dpmsr_set$x$raw_data_input == "Protein_Peptide" || dpmsr_set$x$raw_data_input == "Peptide") 
        && dpmsr_set$x$final_data_output == "Protein" && !as.logical(dpmsr_set$x$tmt_spqc_norm)   )
      {
        cat(file = stderr(), "collapse peptide to protein... 4.1", "\n")
        protein_data <- add_column(protein_data, dpmsr_set$data$protein_missing, .after = "Peptides")
        #dpmsr_set$y$info_columns_final <<- ncol(protein_data)-dpmsr_set$y$sample_number
        dpmsr_set$y$info_columns_final <<- ncol(protein_data) - sample_columns
        names(protein_data)[grep("Peptides", colnames(protein_data)) + 1] <- "Detected_Imputed"
      }
    
    if (dpmsr_set$x$raw_data_input == "Precursor" && dpmsr_set$x$final_data_output == "Protein" && !as.logical(dpmsr_set$x$tmt_spqc_norm))
    {
      cat(file = stderr(), "collapse to protein... 4.2", "\n")
      protein_data <- add_column(protein_data, dpmsr_set$data$protein_missing, .after = "Peptides")
      #dpmsr_set$y$info_columns_final <<- ncol(protein_data)-dpmsr_set$y$sample_number
      dpmsr_set$y$info_columns_final <<- ncol(protein_data) - sample_columns
      names(protein_data)[grep("Peptides", colnames(protein_data)) + 1] <- "Detected_Imputed"
    }
    
  }

  
  cat(file = stderr(), "finished collapse_peptide...", "\n")
  
  return(protein_data) 
  
  
}


#-------------------------------------------------------------------------------------------

collapse_peptide_setup <- function(peptide_data, info_columns, directlfq=FALSE){
  
  #test_peptide_data <<- peptide_data
  #test_info_columns <<- info_columns
  #peptide_data <- test
  #info_columns <- test_info
  
  cat(file = stderr(), "starting collapse_peptide_setup...", "\n")
  
  peptide_annotate <- peptide_data[1:info_columns]
  peptide_data <- peptide_data[(info_columns + 1):ncol(peptide_data)]
  peptide_data[is.na(peptide_data)] <- 0
  
  # Issue with previous version of dpmsr_set file not have unique peptide information
  if ("Unique" %in% colnames(peptide_annotate))
  {
    cat(file = stderr(), "Unique column found in peptide data...", "\n")
    peptide_annotate <- peptide_annotate[, c("Accession", "Description", "Genes", "Unique")]
  }else{
    cat(file = stderr(), "Unique column NOT found in peptide data...", "\n")
    peptide_annotate <- peptide_annotate[, c("Accession", "Description", "Genes")]
  }
  
  
  #count number of peptides for each protein
  peptide_annotate$Peptides <- 1
  peptide_annotate$Peptides <- as.numeric(peptide_annotate$Peptides)
  
  
  # Issue with previous version of dpmsr_set file not have unique peptide information
  if ("Unique" %in% colnames(peptide_annotate))
  {
    #count number of unique peptides for each protein
    peptide_annotate$Unique[peptide_annotate$Unique == "Unique"] <- 1
    peptide_annotate$Unique[peptide_annotate$Unique != 1] <- 0
    peptide_annotate$Unique <- as.numeric(peptide_annotate$Unique)
    peptide_annotate <- peptide_annotate[, c("Accession", "Description", "Genes", "Peptides", "Unique")]
  }else{
    peptide_annotate <- peptide_annotate[, c("Accession", "Description", "Genes", "Peptides")]
  }
  
  #combine data
  peptide_data <- cbind(peptide_annotate, peptide_data)
  
  cat(file = stderr(), "finished collapse_peptide_setup...", "\n")
  return(peptide_data)
} 




#--------------------------------------------------------------------------------
rollup_sum <- function(peptide_data, info_columns){
  cat(file = stderr(), "rollup_sum triggered...", "\n")
  
  protein_data <- peptide_data %>% group_by(Accession, Description, Genes) %>% summarise_all(list(sum))
  protein_data <- data.frame(ungroup(protein_data))
  
  return(protein_data)
}

#--------------------------------------------------------------------------------
rollup_median <- function(peptide_data, info_columns){
  cat(file = stderr(), "rollup_median triggered...", "\n")
  
  df <- peptide_data[,1:info_columns]
  df <- df %>% group_by(Accession, Description, Genes) %>% summarise_all(list(sum))
  df <- data.frame(ungroup(df))
  
  #select data and log
  df_data <- peptide_data[,(info_columns + 1):ncol(peptide_data)]
  df_data[df_data == 0] <- NA
  df_data <- log2(df_data)
  df_data$Accession <- peptide_data$Accession
  df_data$Description <- peptide_data$Description
  df_data$Genes <- peptide_data$Genes
  
  df_data <- df_data %>% group_by(Accession, Description, Genes) %>% summarise_all(list(median))
  df_data <- data.frame(ungroup(df_data))
  df_data[,3:ncol(df_data)] <- 2^df_data[,3:ncol(df_data)] 
  df_data[is.na(df_data)] <- 0
  
  df_data <- cbind(df, df_data[,3:(ncol(df_data)) ])
  
  return(df_data)
}

#--------------------------------------------------------------------------------
rollup_mean <- function(peptide_data, info_columns){
  cat(file = stderr(), "rollup_mean triggered...", "\n")
  
  df <- peptide_data[,1:info_columns]
  df <- df %>% group_by(Accession, Description) %>% summarise_all(list(sum))
  
  #select data and log
  df_data <- peptide_data[,(info_columns + 1):ncol(peptide_data)]
  df_data[df_data == 0] <- NA
  df_data <- log2(df_data)
  
  df_data$Accession <- peptide_data$Accession
  df_data$Description <- peptide_data$Description
  df_data$Genes <- peptide_data$Genes
  
  df_data <- df_data %>% group_by(Accession, Description, Genes) %>% summarise_all(list(mean))
  df_data <- data.frame(ungroup(df_data))
  df_data[,3:ncol(df_data)] <- 2^df_data[,3:ncol(df_data)]
  df_data[is.na(df_data)] <- 0
  
  df_data <- cbind(df, df_data[,3:(ncol(df_data)) ])
  
  return(df_data)
}

#--------------------------------------------------------------------------------
rollup_median_polish <- function(peptide_data, info_columns){
  cat(file = stderr(), "rollup_median_polish triggered...", "\n")

  #group and sum data (summed samples will be over written below)
  df <- peptide_data %>% group_by(Accession, Description, Genes) %>% summarise_all(list(sum))
  df <- data.frame(ungroup(df))
  
  #select data and log
   df_data <- peptide_data[,(info_columns + 1):ncol(peptide_data)]
   df_data[df_data == 0] <- NA
   df_data <- log2(df_data)
  
  #loop through each grouped protein and rollup
  for (i in 1:nrow(df)) {
    df_temp <- df_data[which(peptide_data$Accession == df$Accession[i] & peptide_data$Description == df$Description[i] & peptide_data$Genes == df$Genes[i]),] 
    df_temp <- medpolish(df_temp, trace.iter = FALSE, na.rm = TRUE)
    data_out <- df_temp$overall + df_temp$col
    data_out <- replace_na(data_out, 0)
    df[i,(info_columns + 1):ncol(df)] <- 2^data_out
  }
  
   
  return(df)
}


#--------------------------------------------------------------------------------
rollup_topn <- function(peptide_data, info_columns){
  cat(file = stderr(), "rollup_topn triggered...", "\n")
  
  #group and sum data (summed samples will be over written below)
  df <- peptide_data %>% group_by(Accession, Description, Genes) %>% summarise_all(list(sum))
  df <- data.frame(ungroup(df))
  
  #select data and log
  df_data <- peptide_data[,(info_columns + 1):ncol(peptide_data)]
  df_data[df_data == 0] <- NA
  df_data <- log2(df_data)

  #loop through each grouped protein and rollup
  for (i in 1:nrow(df)) {
   df_temp <- df_data[which(peptide_data$Accession == df$Accession[i] & peptide_data$Description == df$Description[i] & peptide_data$Genes == df$Genes[i]),] 
   data_mean_row <- rowMeans(df_temp, na.rm = TRUE)
   data_mean_row_sorted <- sort(data_mean_row, decreasing = TRUE, na.last = TRUE, index.return = TRUE)
   data_out <- colMeans(df_temp[data_mean_row_sorted$ix[1:min(dpmsr_set$y$rollup_topN_count, length(data_mean_row))],], na.rm = TRUE)
   data_out <- replace_na(data_out, 0)
   df[i,(info_columns + 1):ncol(df)] <- 2^data_out
  }
  
  return(df)
}
  

#--------------------------------------------------------------------------------
rollup_maxlfq <- function(peptide_data, info_columns){
  cat(file = stderr(), "rollup_maxlfq triggered...", "\n")
  require(iq)
  
  #peptide_data <<- peptide_data
  #info_columns <<- info_columns
  #peptide_data <- peptide_data[-20]

  samples_columns <- ncol(peptide_data) - info_columns
  
  #group and sum data (summed samples will be over written below)
  df <- peptide_data %>% group_by(Accession, Description, Genes) %>% summarise_all(list(sum))

  #select data and log
  df_data <- peptide_data[,(info_columns + 1):ncol(peptide_data)]
  df_data[df_data == 0] <- NA
  df_data <- log2(df_data)
  
  #build iq:maxlfq data input, build first sample manually then loop
  peptide_data$id <- 1:nrow(peptide_data)
  df_maxlfq <- data.frame(peptide_data$Accession)
  sample_list <- colnames(peptide_data[info_columns + 1])
  df_maxlfq$sample_list <- sample_list
  quant <- df_data[,1]
  df_maxlfq$quant <- quant
  df_maxlfq$id <- peptide_data$id
  colnames(df_maxlfq) <- c("protein_list", "sample_list", "quant", "id")
  
  for (i in 2:ncol(df_data)) {
    df_temp <- data.frame(peptide_data$Accession)
    sample_list <- colnames(peptide_data[info_columns + i])
    df_temp$sample_list <- sample_list
    quant <- df_data[,i]
    df_temp$quant <- quant
    df_temp$id <- peptide_data$id
    colnames(df_temp) <- c("protein_list", "sample_list", "quant", "id") 
    df_maxlfq <- rbind(df_maxlfq, df_temp)
  }
  
  #remove missing values, iq functions below do not like them!
  df_maxlfq <- df_maxlfq[!is.na(df_maxlfq$quant),]
    
  protein_table <- iq::fast_MaxLFQ(df_maxlfq)
  result <- data.frame(2^protein_table$estimate)
  #result[is.na(result)] <- 0
    
  #combine annotation data with protein rollup, (replacing summed data)
  result <- result[order(row.names(result)),]
  df <- df[order(df$Accession),]
  
  df[(info_columns + 1):ncol(df)] <- result
  
    
return(df)  
}


#--------------------------------------------------------------------------------------
# DirectLFQ from Mann ()
rollup_directlfq <- function(peptide_data, info_columns){
  cat(file = stderr(), str_c("directlfq_rollup..."), "\n")
  
  peptide_data <<- peptide_data
  info_columns <<- info_columns

  #group and sum data (summed samples will be over written below)
  df <- peptide_data %>% group_by(Accession, Description, Genes) %>% summarise_all(list(sum))
  sample_number <- ncol(peptide_data) - info_columns
  info_column_names <- colnames(peptide_data)[1:info_columns]

  require(reticulate)
  use_python("/home/dpmsr/miniconda3/envs/directlfq/bin/python3")
  directlfq <- import("directlfq.lfq_manager")
  
  protein <- peptide_data$Accession
  df_data <- add_column(peptide_data, protein, .before = 1)
  
  # if (dpmsr_set$x$data_source == "PD"){
  #   ion <- str_c(df_data$Sequence, "_", df_data$Modifications, "_", 1:nrow(df_data))
  # }else if (dpmsr_set$x$data_source == "SP"){
  #   ion <- df_data$PrecursorId
  # }
  # 
  
  ion <- 1:nrow(peptide_data)
  
  df_data <- add_column(df_data, ion, .after = 1)
  
  data_in <- df_data[, !(names(df_data) %in% info_column_names)]
  data_in[data_in == 0] <- NA
  
  directlfq_file <- str_c(dpmsr_set$file$extra_dir, "directlfq.aq_reformat.tsv")
  directlfq_protein_file <- str_c(dpmsr_set$file$extra_dir, "directlfq.aq_reformat.tsv.protein_intensities.tsv")
  directlfq_ion_file <- str_c(dpmsr_set$file$extra_dir, "directlfq.aq_reformat.tsv.ion_intensities.tsv")
  write.table(data_in, file = directlfq_file, sep = "\t", row.names = FALSE)
  
  directlfq$run_lfq(directlfq_file)
  data_out <- Simple_fread((directlfq_protein_file))
  colnames(data_out)[which(names(data_out) == "protein")] <- "Accession"
  
  data_out <- df[1:(ncol(df) - sample_number)] %>% left_join(data_out, by = "Accession")
  data_out <- data_out[, !(names(data_out) %in% c("V1"))]
  data_out[data_out == 0] <- NA
  
  file.remove(directlfq_file)
  file.remove(directlfq_ion_file)
  file.remove(directlfq_protein_file)
  
  return(ungroup(data_out))
}

















