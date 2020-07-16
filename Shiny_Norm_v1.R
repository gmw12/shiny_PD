norm_prep <- function(){
  
  if ((dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") 
      && !as.logical(dpmsr_set$x$peptide_isoform)) {
    dpmsr_set$data$norm_data <<- remove_duplicates(dpmsr_set$data$data_peptide)  
    if (as.logical(dpmsr_set$x$peptide_ptm_norm)){
      dpmsr_set$data$norm_data <<- dpmsr_set$data$norm_data[grep(dpmsr_set$x$peptide_norm_grep, 
                                      dpmsr_set$data$norm_data$Modifications),]
    }
    Simple_Excel(dpmsr_set$data$norm_data,  str_c(dpmsr_set$file$extra_prefix, "_norm_data.xlsx", collapse = " "))
    dpmsr_set$data$data_to_norm <<- dpmsr_set$data$data_peptide 
  }
  
  
  if ((dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") 
      && as.logical(dpmsr_set$x$peptide_isoform)) {
    dpmsr_set$data$norm_data <<- remove_duplicates(dpmsr_set$data$data_peptide)
    if (as.logical(dpmsr_set$x$peptide_ptm_norm)){
      dpmsr_set$data$norm_data <<- dpmsr_set$data$norm_data[grep(dpmsr_set$x$peptide_norm_grep, 
                                                                 dpmsr_set$data$norm_data$Modifications),]
    }
    Simple_Excel(dpmsr_set$data$norm_data,  str_c(dpmsr_set$file$extra_prefix, "_norm_data.xlsx", collapse = " "))
    dpmsr_set$data$data_to_norm <<- dpmsr_set$data$data_peptide_isoform
  }
  
  
  if (dpmsr_set$x$raw_data_input=="Protein") {
    dpmsr_set$data$norm_data <<- dpmsr_set$data$data_protein
    Simple_Excel(dpmsr_set$data$norm_data,  str_c(dpmsr_set$file$extra_prefix, "_norm_data.xlsx", collapse = " "))
    dpmsr_set$data$data_to_norm <<- dpmsr_set$data$data_protein
  }

  #add column for missing value (impute) statistics before datasets expand
  dpmsr_set$data$norm_data <<- add_imputed_column(dpmsr_set$data$norm_data)
  dpmsr_set$data$data_to_norm <<- add_imputed_column(dpmsr_set$data$data_to_norm)
  dpmsr_set$y$info_columns <<- ncol(dpmsr_set$data$norm_data)-dpmsr_set$y$sample_number
  dpmsr_set$y$info_columns_final <<- dpmsr_set$y$info_columns
}


#--------------------------------------------------------------------------------------
apply_norm <- function(){
  ncores <- detectCores()
  if (is.na(ncores)) {ncores <- 1}
  norm_list <- dpmsr_set$y$norm_list
  dpmsr_set$data$normalized <<- mclapply(norm_list, norm_parallel, mc.cores = ncores)
  }


#--------------------------------------------------------------------------------------
norm_parallel <- function(norm_type){
  norm_data <- dpmsr_set$data$norm_data
  data_to_norm <- dpmsr_set$data$data_to_norm
  info_columns <- ncol(norm_data) - dpmsr_set$y$sample_number
 if(norm_type==1){data_sl <- sl_normalize(norm_data, data_to_norm, "SL_Norm", info_columns)}
  else if(norm_type==4){data_quantile <- quantile_normalize(data_to_norm, "Quantile_Norm", info_columns)}
  else if(norm_type==5){data_lr <- lr_normalize(data_to_norm, "LinReg_Norm", info_columns)}
  else if(norm_type==6){data_loess <- loess_normalize(data_to_norm, "LOESS_Norm", info_columns)} 
  else if(norm_type==7){data_vsn <- vsn_normalize(data_to_norm, "VSN_Norm", info_columns)}
  else if(norm_type==8){data_ti <- ti_normalize(norm_data, data_to_norm, "TI_Norm", info_columns)} 
  else if(norm_type==9){data_mi <- mi_normalize(norm_data, data_to_norm, "MI_Norm", info_columns)} 
  else if(norm_type==10){data_ai <- ai_normalize(norm_data, data_to_norm, "AI_Norm", info_columns)} 
  else if(norm_type==99){data_impute <- data_to_norm} 
}



#--------------------------------------------------------------------------------------
# global scaling value, sample loading normalization
sl_normalize <- function(norm_data, data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  norm_data <- norm_data[(info_columns+1):ncol(norm_data)]
  excel_name <- "_Peptide_SL_Norm.xlsx"
  target <- mean(colSums(norm_data, na.rm=TRUE))
  norm_facs <- target / colSums(norm_data,na.rm=TRUE)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_sl_norm.xlsx", collapse = " "))
  return(data_out)
}

# average global scaling value, sample loading normalization
ai_normalize <- function(norm_data, data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  norm_data <- norm_data[(info_columns+1):ncol(norm_data)]
  excel_name <- "_Peptide_AI_Norm.xlsx"
  target <- mean(colMeans(norm_data, na.rm=TRUE))
  norm_facs <- target / colMeans(norm_data, na.rm=TRUE)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_ai_norm.xlsx", collapse = " "))
  return(data_out)
}

# intensity / sum of intensities then * median of sum of intensitites
ti_normalize <- function(norm_data, data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  norm_data <- norm_data[(info_columns+1):ncol(norm_data)]
  excel_name <- "_Peptide_TI_Norm.xlsx"
  target <- median(colSums(norm_data, na.rm=TRUE))
  norm_facs <- target / colSums(norm_data, na.rm=TRUE)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_ti_norm.xlsx", collapse = " "))
  return(data_out)
}

# intensity / sum of intensities then * median of sum of intensitites
mi_normalize <- function(norm_data, data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  norm_data <- norm_data[(info_columns+1):ncol(norm_data)]
  excel_name <- "_Peptide_MI_Norm.xlsx"
  intensity_medians <- colMedians(data.matrix(norm_data), na.rm=TRUE)
  target <- mean(intensity_medians)
  norm_facs <- target / intensity_medians
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_mi_norm.xlsx", collapse = " "))
  return(data_out)
}


# global scaling value, sample loading normalization
vsn_normalize <- function(data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  excel_name <- "_Peptide_VSN_Norm.xlsx"
  data_to_norm[data_to_norm==0] <- NA
  data_out <- normalizeVSN(data.matrix(data_to_norm))
  data_out < data.frame(data_out)
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  #data_out[is.na(data_out)] <- 0.0
  data_out[data_out==0] <- NA
  data_out <- data_out / 10
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_vsn_norm.xlsx", collapse = " "))
  return(data_out)
}


# global scaling value, sample loading normalization
quantile_normalize <- function(data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  excel_name <- "_Peptide_Quantile_Norm.xlsx"
  data_out <- normalize.quantiles(data.matrix(data_to_norm))
  data_out <- data.frame(data_out)
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_quantile_norm.xlsx", collapse = " "))
  return(data_out)
}


# global scaling value, sample loading normalization
loess_normalize <- function(data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  excel_name <- "_Peptide_LOESS_Norm.xlsx"
  data_to_norm <- log2(data_to_norm)
  data_out <- normalizeCyclicLoess(data_to_norm, weights = NULL, span=0.7, iterations = 3, method = "fast")
  data_out <- data.frame(data_out)
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  #data_out[is.na(data_out)] <- 0.0
  data_out[data_out==0] <- NA
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_loess_norm.xlsx", collapse = " "))
  return(data_out)
}


#linear regression normalization
lr_normalize <- function(data_in, data_title, info_columns) {
  annotation_data <- data_in[1:info_columns]
  data_in <- data_in[(info_columns+1):ncol(data_in)]
  excel_name <- "_Peptide_LR_Norm.xlsx"
  #normalize lr on data with no missing values, create new data frame
  data_nomissing <- data_in
  data_nomissing$missingvalues <- rowSums(data_nomissing ==0)
  data_nomissing <- subset(data_nomissing, missingvalues==0)
  data_nomissing <- data_nomissing[1:dpmsr_set$y$sample_number]
  #log2 data
  data_nomissing <- log(data_nomissing,2)
  data_nomissing[data_nomissing==-Inf] = 0  # fix log2 of 0}
  data_in <- log(data_in,2)
  data_in[data_in==-Inf] = 0  # fix log2 of 0}  
  data_in[data_in==0] <- NA
  data_out <- data_in
  #reorders data by intensity independent of indentification
  for(i in 1:dpmsr_set$y$sample_number){
    temp <- data.frame(data_nomissing[,i])
    colnames(temp) <- "test"
    temp <- arrange(temp, test)
    data_nomissing[,i] <- temp
  }
  colnames(data_nomissing) <- seq(from=1, to=dpmsr_set$y$sample_number)
  data_nomissing$avg <- apply(data_nomissing, 1, FUN = function(x) {median(x[x > 0])})
  for(i in 1:dpmsr_set$y$sample_number){
    data_test <- data.frame(cbind(data_nomissing[,i], data_nomissing$avg))
    colnames(data_test) <- c("x", "y")
    LMfit <- rlm(x~y, data_test, na.action=na.exclude)
    Coeffs <- LMfit$coefficients
    m <- Coeffs[2] # y = mX + b
    b <- Coeffs[1] 
    normdata <- (data_in[,i] - b) / m
    data_out[,i] <- normdata
  }
  data_out <- data.frame(2^data_out)
  data_out[data_out==-Inf] = 0  # fix log2 of 0  
  data_out[data_out==0] <- NA
  #data_out[is.na(data_out)] <- 0.0
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_lr_norm.xlsx", collapse = " "))
  return(data_out)
}


# global scaling value, sample loading normalization
protein_normalize <- function(data_to_norm, data_title, info_columns){
  protein_norm_raw <-subset(data_to_norm, Accession %in% dpmsr_set$x$protein_norm_list)
  protein_norm_raw <- protein_norm_raw[(info_columns+1):(info_columns+dpmsr_set$y$sample_number)]
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  #protein_norm_raw$missings <- rowSums(protein_norm_raw == 0.0)
  #protein_norm_raw <- subset(protein_norm_raw, missings==0)
  #protein_norm_raw <- protein_norm_raw[1:sample_number]
  target <- mean(colSums(protein_norm_raw, na.rm=TRUE))
  norm_facs <- target / colSums(protein_norm_raw, na.rm=TRUE)
  data_out <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  #data_out <- finish_norm(data_out, excel_name, data_title)
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_protein_norm.xlsx", collapse = " "))
  return(data_out)
}


#TMM Normalized 
tmm_normalize <- function(norm_data, data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  norm_data <- norm_data[(info_columns+1):ncol(norm_data)]
  #excel_name <- "_Peptide_TMM1_Norm.xlsx"
  #Simple_Excel(cbind(annotate_df, tmm), str_c(file_prefix3, "_Peptide_TMM_Norm_Impute.xlsx", collapse = " ")) 
  tmm <- norm_data
  tmm[is.na(tmm)] <- 0.0
  tmm_factor <- calcNormFactors(tmm, method = "TMM", sumTrim = 0.1)
  data_out <- sweep(tmm, 2, tmm_factor, FUN = "/") # this is data after SL and TMM on original scale
  #Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_TMM2_Norm.xlsx", collapse = " "))
  #Plot_All_gw(data_out, data_title)
  #data_out <- cbind(annotate_df, data_out)
  #if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  #colnames(data_out) <- sample_header_final_norm
  #Simple_Excel(data_out, str_c(file_prefix3, "_", data_title, "_final.xlsx", collapse = " "))
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_tmm_norm.xlsx", collapse = " "))
  return(data_out)
}


# SL TMM Normalized 
sltmm_normalize <- function(norm_data, data_to_norm, data_title, info_columns){
  annotation_data <- data_to_norm[1:info_columns]
  data_to_norm <- data_to_norm[(info_columns+1):ncol(data_to_norm)]
  norm_data <- norm_data[(info_columns+1):ncol(norm_data)]
  #excel_name <- "_Peptide_SLTMM1_Norm.xlsx"
  target <- mean(colSums(norm_data))
  norm_facs <- target / colSums(norm_data)
  sl <- sweep(data_to_norm, 2, norm_facs, FUN = "*")
  Simple_Excel(cbind(annotate_df, sl), str_c(file_prefix3, "_Peptide_SLTMM1_Norm.xlsx", collapse = " "))
  sl <- impute_only(sl)
  Simple_Excel(cbind(annotate_df, sl), str_c(file_prefix3, "_Peptide_SLTMM_Norm_Impute.xlsx", collapse = " ")) 
  #sl_tmm <- calcNormFactors(norm_data, method = "TMM", sumTrim = 0.1) # this is data after SL and TMM on original scale
  sl_tmm <- calcNormFactors(sl, method = "TMM", sumTrim = 0.1) # this is data after SL and TMM on original scaledata_out <- sweep(sl, 2, sl_tmm, FUN = "/")
  data_out <- sweep(sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
  Simple_Excel(cbind(annotate_df, data_out), str_c(file_prefix3, "_Peptide_SLTMM2_Norm.xlsx", collapse = " "))
  Plot_All_gw(data_out, data_title)
  data_out <- cbind(annotate_df, data_out)
  if (peptide_to_protein){data_out <- collapse_peptide(data_out)}
  colnames(data_out) <- sample_header_final_norm
  data_out <- cbind(annotation_data, data_out)
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_sltmm_norm.xlsx", collapse = " "))
  return(data_out)
}



remove_duplicates <- function(data_in){
  data_in$Modifications[is.na(data_in$Modifications)] <- ""
  data_in$dup <- str_c(data_in$Sequence, "_", data_in$Modifications)
  data_out <- distinct(data_in, dup, .keep_all = TRUE)
  data_out$dup <- NULL
  Simple_Excel(data_out,  str_c(dpmsr_set$file$extra_prefix, "_remove_duplicates.xlsx", collapse = " "))
  #data_out <- data_out[(info_columns+1):ncol(data_out)]
  return(data_out)
}




