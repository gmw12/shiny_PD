protein_to_peptide <- function(){
  cat(file = stderr(), "protein_to_peptide", "\n")
  protein <- dpmsr_set$data$data_raw_protein
  peptide_groups <- dpmsr_set$data$data_raw_peptide
  
  #add columns to preserve peptide to protein links
  peptide_groups$Proteins <- peptide_groups$Protein.Accessions
  peptide_groups$Unique <- peptide_groups$Quan.Info
  peptide_groups$Unique[peptide_groups$Unique == ""] <- "Unique"
  
  #protein raw has all confidence proteins - limit to high master
  protein_master <- subset(protein, Master %in% ("IsMasterProtein"))
  protein_high_master <- subset(protein_master, Protein.FDR.Confidence.Combined %in% ("High"))
  master_accessions <- protein_high_master$Accession 
  
  #PD will label which proteins get the razor peptides, 
  protein_razor <- subset(protein, Number.of.Razor.Peptides > 0)
  razor_accessions <- protein_razor$Accession
   
  #gather peptides that are shared
  peptide_shared <- subset(peptide_groups,  Quan.Info %in% ("NotUnique"))
  #gather peptides that have no quant values
  peptide_noquan <- subset(peptide_groups,  Quan.Info %in% ("NoQuanValues"))
  #gather unique peptides
  peptide_unique <- peptide_groups[peptide_groups$Quan.Info == "",]
  
  #expand shared peptides so that each protein has peptides listed separately
  peptide_shared_expand  <- peptide_shared %>% 
    mutate(Master.Protein.Accessions = strsplit(as.character(Master.Protein.Accessions), "; ", fixed = TRUE)) %>% 
    unnest(Master.Protein.Accessions)
  
  if (dpmsr_set$x$peptides_to_use == "Razor") {
    #reduce df to only peptides that have proteins that PD lists as having "razor" peptides
    peptide_shared_expand <- subset(peptide_shared_expand, Master.Protein.Accessions %in% razor_accessions )
    #gather df for razor proteins
    protein_razor_lookup <- protein_razor %>% dplyr::select(Accession, Description, Number.of.Peptides, 
                                                            Coverage.in.Percent, Number.of.Unique.Peptides, Number.of.Razor.Peptides)
    #add columns from protein to df
    peptide_shared_expand <- merge(peptide_shared_expand, protein_razor_lookup, by.x = "Master.Protein.Accessions", by.y = "Accession")
    #create column to check for duplicated peptides
    peptide_shared_expand$duplicated_test <- str_c(peptide_shared_expand$Annotated.Sequence, peptide_shared_expand$Modifications)
    peptide_shared_expand <- peptide_shared_expand[order(peptide_shared_expand$duplicated_test, -peptide_shared_expand$Number.of.Peptides, 
                                                         peptide_shared_expand$Coverage.in.Percent, 
                                                         -peptide_shared_expand$Number.of.Razor.Peptides),]
    #remove duplicated peptides
    peptide_final <- peptide_shared_expand[!duplicated(peptide_shared_expand$duplicated_test),]
    peptide_final$Master.Protein.Descriptions <- peptide_final$Description
    #remove extra columns
    peptide_final <- peptide_final[1:(ncol(peptide_groups))]
    #combine unique and razor/shared peptides
    peptide_final <- rbind(peptide_unique, peptide_final)
  }else if (dpmsr_set$x$peptides_to_use == "Shared") {
    peptide_final <- rbind(peptide_unique, peptide_shared_expand)
  }else{
    peptide_final <- peptide_unique
  }
  
  peptide_final <- peptide_final[order(peptide_final$Master.Protein.Accessions, peptide_final$Sequence),]
  peptide_out <- peptide_final %>% dplyr::select(Confidence, Master.Protein.Accessions, Master.Protein.Descriptions, Proteins, 
                                                 Sequence, Modifications, Unique,
                                                 contains('RT.in.min.by.Search.Engine.'), 
                                                 starts_with('mz.in.Da.by.Search.Engine.'), 
                                                 contains('Charge.by.Search.Engine.'), 
                                                 contains('Percolator.SVM'), 
                                                 contains("Percolator.q.Value"), contains("Abundance.F"))
  
  
  if (ncol(peptide_out) != (12 + dpmsr_set$y$sample_number))
  {
    shinyalert("Oops!", str_c("Number of columns extracted is not as expected ", ncol(peptide_out), "/", (10 + dpmsr_set$y$sample_number)), type = "error")  
  }
  
  colnames(peptide_out)[1:12] <- c("Confidence", "Accession", "Description", "All.Proteins", "Sequence", "Modifications", "Unique", "Retention.Time","Da","mz", "Ion.Score", "q-Value")
  peptide_out <- subset(peptide_out, Accession %in% master_accessions )
  Simple_Excel(peptide_out, "Protein_Peptide_Raw", str_c(dpmsr_set$file$extra_prefix,"_ProteinPeptide_to_Peptide_Raw.xlsx", collapse = " "))
  return(peptide_out)
}


#----------------------------------------------------------------------------------------
protein_to_protein <- function(){
  cat(file = stderr(), "protein_to_protein", "\n")
  protein <- dpmsr_set$data$data_raw_protein
  
  if (dpmsr_set$x$data_source == "PD") {
    cat(file = stderr(), "data type  -> PD", "\n")
    protein <- subset(protein, Master %in% ("IsMasterProtein"))
    protein <- subset(protein, Protein.FDR.Confidence.Combined %in% ("High"))
    protein_out <- protein %>% dplyr::select(Accession, Description, Number.of.Protein.Unique.Peptides, 
                                             contains("Abundance"), -contains("Abundance.Count"))
    colnames(protein_out)[1:3] <- c("Accession", "Description", "Unique.Peptides")
  }
  else if (dpmsr_set$x$data_source == "SP") { 
    cat(file = stderr(), "data type  -> SP", "\n")
    protein_out <- protein %>% dplyr::select(contains("ProteinAccessions"), contains("ProteinDescriptions"), 
                                             contains("ProteinNames"), contains("Genes"), contains("Quantity"))
    precursor_col <- protein %>% dplyr::select(contains("Precursors"))
    precursor_col$average <- round(rowMeans(precursor_col), 1)
    
    protein_out <- protein_out %>% add_column(precursor_col$average, .after = "PG.Genes")
    colnames(protein_out)[1:5] <- c("Accession", "Description", "ProteinName", "Gene", "PrecursorsAvg")

    #in case missing values reported as NaN
    protein_out[, 5:ncol(protein_out)] <- sapply(protein_out[, 5:ncol(protein_out)], as.numeric)
    protein_out[5:ncol(protein_out)][protein_out[5:ncol(protein_out)] == "NaN"] <- 0

  }else {
    cat(file = stderr(), "protein_to_protein data source not recognized", "\n")
  }
  Simple_Excel(protein_out, "Protein_Protein_Raw", str_c(dpmsr_set$file$extra_prefix, "_Protein_Protein_Raw", "_Protein_to_Protein_Raw.xlsx", collapse = " "))
return(protein_out)
}

#----------------------------------------------------------------------------------------
peptide_to_peptide <- function(){
  cat(file = stderr(), "peptide_to_peptide", "\n")
  peptide_groups <- dpmsr_set$data$data_raw_peptide
  
  if (dpmsr_set$x$data_source == "PD") {
    cat(file = stderr(), "peptide_to_peptide, PD data", "\n")
    peptide_out <- peptide_groups %>% dplyr::select(Confidence, Master.Protein.Accessions, Master.Protein.Descriptions, 
                                                  Sequence, Modifications, 
                                                  (starts_with("Positions.in.") & ends_with("Proteins")), 
                                                  (starts_with("Modifications.in.") & ends_with("Proteins")), 
                                                  contains('RT.in.min.by.Search.Engine.'), 
                                                  contains('Percolator.SVM'),  
                                                  contains("Percolator.q.Value"), contains("Abundance.F"))
    
    if (ncol(peptide_out) != (10 + dpmsr_set$y$sample_number))
    {
      shinyalert("Oops!", "Number of columns extracted is not as expected", type = "error")  
    }
    
    colnames(peptide_out)[1:10] <- c("Confidence", "Accession", "Description", "Sequence", "Modifications", "PositionMaster", "ModificationMaster",
                                    "Retention.Time", "SVM.Score", "q-Value")
    peptide_out <- subset(peptide_out, Confidence %in% ("High"))
  
  }else if (dpmsr_set$x$data_source == "SP") {
    cat(file = stderr(), "peptide_to_peptide, SP data", "\n")
    
    peptide_out <- peptide_groups %>% dplyr::select(contains('ProteinAccessions'), contains('ProteinDescriptions'), contains('Genes'), 
                                                    contains('ModifiedSequence'), contains('PeptidePosition'),
                                                    contains('ProteinPTMLocations'),
                                                    contains("TotalQuantity"))
    
    if (ncol(peptide_out) != (6 + dpmsr_set$y$sample_number))
    {
      shinyalert("Oops!", "Number of columns extracted is not as expected", type = "error")  
    }
    
    colnames(peptide_out)[1:6] <- c("Accession", "Description", "Genes", "Sequence", "PeptidePosition", "PTMLocations")  
  }
  
  # set "Filtered" in TotalQuantity to NA
  peptide_out[peptide_out ==  "Filtered"] <- NA
  peptide_out[8:ncol(peptide_out)] <- as.data.frame(lapply(peptide_out[8:ncol(peptide_out)], as.numeric))
  
  Simple_Excel(peptide_out, "Peptide_Peptide_Raw",  str_c(dpmsr_set$file$extra_prefix, "_Peptide_to_Peptide_Raw.xlsx", collapse = " "))
  cat(file = stderr(), "peptide_to_peptide complete", "\n")
  return(peptide_out)
}

#----------------------------------------------------------------------------------------
precursor_to_precursor <- function(){
  cat(file = stderr(), "precursor_to_precursor", "\n")
  precursor_groups <- dpmsr_set$data$data_raw_precursor
  
  precursor_colnames <- c("Accession", "Description", "Genes", "Organisms", "Sequence", "PrecursorId", "PeptidePosition")  
  n_col <- length(precursor_colnames)
  
  precursor_out <- precursor_groups %>% dplyr::select(contains('ProteinAccessions'), contains('ProteinDescriptions'), contains('Genes'), contains('Organisms'),
                                                    contains('ModifiedSequence'), contains('PrecursorId'), contains('PeptidePosition'),
                                                    contains("TotalQuantity"))
    
  if (ncol(precursor_out) != (n_col + dpmsr_set$y$sample_number))
  {
    shinyalert("Oops!", "Number of columns extracted is not as expected", type = "error")  
  }
    
  colnames(precursor_out)[1:n_col] <- precursor_colnames  
  
  # set "Filtered" in TotalQuantity to NA
  precursor_out[precursor_out ==  "Filtered"] <- NA
  precursor_out[(n_col + 1):ncol(precursor_out)] <- as.data.frame(lapply(precursor_out[(n_col + 1):ncol(precursor_out)], as.numeric))
  
  precursor_out$Description <- str_c(precursor_out$Description, ", org=", precursor_out$Organisms) 
  precursor_out$Organisms <- NULL
  
  Simple_Excel(precursor_out, "precursor_precursor_Raw",  str_c(dpmsr_set$file$extra_prefix, "_Precursor_to_Precursor_Raw.xlsx", collapse = " "))
  cat(file = stderr(), "precursor_to_precursor complete", "\n")
  return(precursor_out)
}

#----------------------------------------------------------------------------------------
precursor_PTM_to_precursor_PTM <- function(){
  cat(file = stderr(), "precursor_PTM_to_precursor_PTM", "\n")
  
  precursor_groups <- dpmsr_set$data$data_raw_precursor
  precursor_colnames <- c("Accession", "Description", "Genes", "Organisms", "Stripped_Seq", "Sequence", "PrecursorId", "PeptidePosition", "PTMLocations") 
  n_col <- length(precursor_colnames)

  precursor_out <- precursor_groups %>% dplyr::select(contains('ProteinAccessions'), contains('ProteinDescriptions'), contains('Genes'), contains('Organisms'),
                                                      contains('StrippedSequence'), contains('ModifiedSequence'), contains('PrecursorId'), contains('PeptidePosition'),
                                                      contains('ProteinPTMLocations'), contains("TotalQuantity"))
  
  if (ncol(precursor_out) != (n_col + dpmsr_set$y$sample_number))
  {
    shinyalert("Oops!", "Number of columns extracted is not as expected", type = "error")  
  }
  
  colnames(precursor_out)[1:n_col] <- precursor_colnames  
  
  # set "Filtered" in TotalQuantity to NA
  precursor_out[precursor_out ==  "Filtered"] <- NA
  precursor_out[(n_col + 1):ncol(precursor_out)] <- as.data.frame(lapply(precursor_out[(n_col + 1):ncol(precursor_out)], as.numeric))
  
  precursor_out$Description <- str_c(precursor_out$Description, ", org=", precursor_out$Organisms) 
  precursor_out$Organisms <- NULL
  
  Simple_Excel(precursor_out, "precursor_precursor_Raw",  str_c(dpmsr_set$file$extra_prefix, "_Precursor_to_Precursor_Raw.xlsx", collapse = " "))
  cat(file = stderr(), "precursor_to_precursor complete", "\n")
  return(precursor_out)
}


#Top.Apex.RT.in.min,
#----------------------------------------------------------------------------------------
isoform_to_isoform <- function(){
  cat(file = stderr(), "isoform_to_isoform", "\n")

  if (is.null(dpmsr_set$data$data_raw_isoform)) {
    cat(file = stderr(), "isoform text file NOT found", "\n")
    shinyalert("Oops!", "Isoform data not imported.  TMT datasets do not automatically export isoform data.", type = "error")
    }
    else {
      cat(file = stderr(), "isoform text file found", "\n")
      peptide_groups <- dpmsr_set$data$data_raw_isoform
      peptide_out <- try(peptide_groups %>% dplyr::select(contains("Confidence.by"), Master.Protein.Accessions, Master.Protein.Descriptions, 
                                                        Sequence, Modifications, 
                                                        (starts_with("Positions.in.") & ends_with("Proteins")), 
                                                        (starts_with("Modifications.in.") & ends_with("Proteins")), 
                                                        Top.Apex.RT.in.min, 
                                                        contains('Percolator.SVM'),  
                                                        contains("Percolator.q.Value"), contains("Abundance.F")))
      if (class(peptide_out) == 'try-error') {
        cat(file = stderr(), "column select error - retry", "\n")
        peptide_out <- peptide_groups %>% dplyr::select(contains("Confidence.by"), Master.Protein.Accessions, Master.Protein.Descriptions,
                                                        Sequence, Modifications, 
                                                        (starts_with("Positions.in.") & ends_with("Proteins")), 
                                                        (starts_with("Modifications.in.") & ends_with("Proteins")), 
                                                        contains("Positions."),
                                                        contains('RT.in.min.by.'), 
                                                        contains('Percolator.SVM'), 
                                                        contains("Percolator.q.Value"), contains("Abundance.F"))
      }
      
      
      cat(file = stderr(), str_c("There are ", ncol(peptide_out) - dpmsr_set$y$sample_number, "/10 info columns"), "\n")
      
      if ((ncol(peptide_out) - dpmsr_set$y$sample_number) < 10) {
        cat(file = stderr(), "If this is TMT phos you will need to manually export the isoform text file, load the correct layout file before export", "\n")
      }
      
      
      if (ncol(peptide_out) != (10 + dpmsr_set$y$sample_number))
      {
        shinyalert("Oops!", "Number of isoform columns extracted is not as expected", type = "error")  
      }
    
      colnames(peptide_out)[1:10] <- c("Confidence", "Accession", "Description", "Sequence", "Modifications", "PositionMaster", "ModificationMaster",
                                       "Retention.Time", "SVM.Score", "q-Value")
      
      peptide_out <- subset(peptide_out, Confidence %in% ("High"))
      Simple_Excel_bg(peptide_out, "Protein_Peptide_Raw", str_c(dpmsr_set$file$extra_prefix, "_Isoform_to_Isoform_Raw.xlsx", collapse = " "))
      cat(file = stderr(), "isoform_to_isoform complete", "\n")
      return(peptide_out)
    }
}






#----------------------------------------------------------------------------------------

psm_set_fdr <- function(){
  psmfdr_dir <- create_dir(str_c(data_dir,"//PSM_FDR"))
  psm_prefix <- str_c(psmfdr_dir, file_prefix)
  
  forward_psm <- dpmsr_set$data$data_raw_psm
  decoy_psm <- dpmsr_set$data$data_raw_decoypsm

  forward_psm$fdr <- rep("forward", nrow(forward_psm))
  decoy_psm$fdr <- rep("decoy", nrow(decoy_psm))
  
  isv1 <- forward_psm$Ions.Score
  isv2 <- decoy_psm$Ions.Score
  isv3 <- c(isv1, isv2)
  
  fdr1 <- forward_psm$fdr
  fdr2 <- decoy_psm$fdr
  fdr3 <- c(fdr1, fdr2)
  
  test <- data.frame(cbind(isv3, fdr3), stringsAsFactors = FALSE)
  colnames(test) <- c("Ions_Score", "FDR")
  test <- test[order(test$Ions_Score),]
  rcount <- nrow(combo_psm)
  
  for (i in 1:rcount) {
    testthis <- data.frame(table(test$FDR[i:rcount]))
    test_fdr <- testthis$Freq[1] / testthis$Freq[2] * 100
    if (test_fdr <= 1.0000) {
      break
    }
  }
  
  ion_score_cutoff <- min(test$Ions_Score[i:rcount])
  
  return(ion_score_cutoff)
}

#----------------------------------------------------------------------------------------
add_imputed_column <- function(df){
  cat(file = stderr(), "add_imputed_column_start...", "\n")
  #check to see if this was already completed, if so skip step
  if ("Detected_Imputed" %in% colnames(df)) {
    return(df)
  }else{
      #imputed column for protein output
      cat(file = stderr(), "add_imputed_column protein output...", "\n")
      if ((dpmsr_set$x$raw_data_input == "Protein_Peptide" || dpmsr_set$x$raw_data_input == "Peptide" || dpmsr_set$x$raw_data_input == "Precursor") 
          && dpmsr_set$x$final_data_output == "Protein") {
        peptide_data <- df
        peptide_annotate <- peptide_data[,1:(dpmsr_set$y$info_columns)]
        peptide_data <- peptide_data[,(dpmsr_set$y$info_columns + 1):ncol(peptide_data)]
        peptide_data[is.na(peptide_data)] <- 0
        peptide_data[peptide_data > 0] <- 1
        #force to numeric for SP data
        #peptide_data <- as.data.frame(lapply(peptide_data, as.numeric))
        peptide_annotate <- peptide_annotate[, c("Accession", "Description")]
        test1 <- cbind(peptide_annotate, peptide_data)
        test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(list(sum))
        test3 <- test2
        test3[3:ncol(test2)][test2[3:ncol(test2)] > 1] <- 1
        dpmsr_set$data$protein_imputed_df <<- test3
        test2 <- data.frame(ungroup(test2))
        test2 <- test2[3:ncol(test2)]
        test2[test2 == 0] <- "-"
        test2 <- test2 %>% mutate_all(as.character)
        while (ncol(test2) > 1) {
          test2[,1] <- str_c(test2[,1], ".", test2[,2])
          test2[,2] <- NULL
        }
        dpmsr_set$data$protein_missing <<- test2[,1]
      }
    
    #imputed column for peptide output when rolling up from precursor
    cat(file = stderr(), "add_imputed_column peptide output from precursor...", "\n")
    if (dpmsr_set$x$raw_data_input == "Precursor" || dpmsr_set$x$raw_data_input == "Precursor_PTM") {
      precursor_data <- df
      precursor_annotate <- precursor_data[,1:(dpmsr_set$y$info_columns)]
      precursor_data <- precursor_data[,(dpmsr_set$y$info_columns + 1):ncol(precursor_data)]
      precursor_data[is.na(precursor_data)] <- 0
      precursor_data[precursor_data > 0] <- 1
      #force to numeric for SP data
      #precursor_data <- as.data.frame(lapply(precursor_data, as.numeric))
      precursor_annotate <- precursor_annotate[, c("Accession", "Description", "Sequence")]
      test1 <- cbind(precursor_annotate, precursor_data)
      test2 <- test1 %>% group_by(Accession, Description, Sequence) %>% summarise_all(list(sum))
      test3 <- test2
      test3[4:ncol(test2)][test2[4:ncol(test2)] > 1] <- 1
      dpmsr_set$data$precursor_imputed_df <<- test3
      test2 <- data.frame(ungroup(test2))
      test2 <- test2[4:ncol(test2)]
      test2[test2 == 0] <- "-"
      test2 <- test2 %>% mutate_all(as.character)
      while (ncol(test2) > 1) {
        test2[,1] <- str_c(test2[,1], ".", test2[,2])
        test2[,2] <- NULL
      }
      dpmsr_set$data$peptide_missing <<- test2[,1]
    }    
    
      #imputed column for peptide output
      cat(file = stderr(), "add_imputed_column peptide...", "\n")
      df_annotation <- df[,1:dpmsr_set$y$info_columns]
      df <- df[,(dpmsr_set$y$info_columns + 1):ncol(df)]
      df_data <- df
      df[df >  0] <- "1"
      imputed_df <- df
      imputed_df[is.na(imputed_df)] <- 0
      dpmsr_set$data$peptide_imputed_df <<- cbind(df_annotation, imputed_df)
      df[is.na(df)] <- "-"
      df <- df %>% mutate_all(as.character)
      while (ncol(df) > 1 ) {
        df[,1] <- str_c(df[,1], ".", df[,2])
        df[,2] <- NULL
      }
      
      df_annotation$Detected_Imputed <- df[,1]
      df <- cbind(df_annotation,df_data)
      #add another column to info columns
      cat(file = stderr(), "add_imputed_column_end...", "\n")
      return(df)
  }
}

#----------------------------------------------------------------------------------------
add_imputed_column_protein <- function(df){
  cat(file = stderr(), "add_imputed_column_protein_start...", "\n")
  #check to see if this was already completed, if so skip step
  #df<-dpmsr_set$data$data_protein
  if ("Detected_Proteins" %in% colnames(df)) {
    return(df)
  }else{
    #imputed column for protein output
    cat(file = stderr(), "add_imputed_column_protein protein output...", "\n")
      protein_data <- df
      protein_annotate <- protein_data[,1:(dpmsr_set$y$info_columns)]
      protein_data <- protein_data[,(dpmsr_set$y$info_columns + 1):ncol(protein_data)]
      test1 <- protein_data
      test1[is.na(test1)] <- 0
      test1[test1 > 0] <- 1
      test1 <- cbind(protein_annotate, test1)
      dpmsr_set$data$protein_imputed_df <<- test1
      test1 <- data.frame(ungroup(test1))
      test1 <- test1[,(dpmsr_set$y$info_columns + 1):ncol(protein_data)]
      test1[test1 == 0] <- "-"
      test1 <- test1 %>% mutate_all(as.character)
      while (ncol(test1) > 1 ) {
        test1[,1] <- str_c(test1[,1], ".", test1[,2])
        test1[,2] <- NULL
      }
      dpmsr_set$data$protein_missing <<- test1[,1]

      df <- cbind(protein_annotate, dpmsr_set$data$protein_missing, protein_data)
      colnames(df)[ncol(protein_annotate) + 1] <- "Protein_Imputed"
      cat(file = stderr(), "add_imputed_column_protein_end...", "\n")
    return(df)
  }
}

#----------------------------------------------------------------------------------------
set_final_info_columns <- function() {
  cat(file = stderr(), "function set_final_info_columns...", "\n")
  df <- dpmsr_set$data$final$impute
  df <- df %>% dplyr::select((!ends_with("_CV")))
  dpmsr_set$y$info_columns_final <<- ncol(df) - dpmsr_set$y$sample_number
  
}

