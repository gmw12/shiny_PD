protein_to_peptide <- function(){
  cat(file=stderr(), "protein_to_peptide", "\n")
  protein <- dpmsr_set$data$data_raw_protein
  peptide_groups <- dpmsr_set$data$data_raw_peptide
  
  protein_master <- subset(protein, Master %in% ("IsMasterProtein"))
  protein_high_master <- subset(protein_master, Protein.FDR.Confidence.Combined %in% ("High"))
  protein_razor <- subset(protein, Number.of.Razor.Peptides>0)
  
  razor_accessions <- protein_razor$Accession
  master_accessions <- protein_high_master$Accession  
  
  peptide_shared <- subset(peptide_groups,  Quan.Info %in% ("NotUnique"))
  peptide_noquan <- subset(peptide_groups,  Quan.Info %in% ("NoQuanValues"))
  peptide_unique <- peptide_groups[peptide_groups$Quan.Info=="",]
  
  peptide_shared_expand  <- peptide_shared %>% 
    mutate(Master.Protein.Accessions = strsplit(as.character(Master.Protein.Accessions), "; ", fixed = TRUE)) %>% 
    unnest(Master.Protein.Accessions)
  
  if(dpmsr_set$x$peptides_to_use=="Razor"){
    peptide_shared_expand <- subset(peptide_shared_expand, Master.Protein.Accessions %in% razor_accessions )
    protein_razor_lookup <- protein_razor %>% dplyr::select(Accession, Description, Number.of.Peptides, 
                                                            Coverage.in.Percent, Number.of.Unique.Peptides, Number.of.Razor.Peptides)
    peptide_shared_expand <- merge(peptide_shared_expand, protein_razor_lookup, by.x="Master.Protein.Accessions", by.y="Accession")
    peptide_shared_expand$unique <- str_c(peptide_shared_expand$Annotated.Sequence, peptide_shared_expand$Modifications)
    peptide_shared_expand <- peptide_shared_expand[order(peptide_shared_expand$unique, -peptide_shared_expand$Number.of.Peptides, 
                                                         peptide_shared_expand$Coverage.in.Percent, 
                                                         -peptide_shared_expand$Number.of.Razor.Peptides),]
    peptide_final <- peptide_shared_expand[!duplicated(peptide_shared_expand$unique),]
    peptide_final$Master.Protein.Descriptions <- peptide_final$Description
    peptide_final <- peptide_final[1:(ncol(peptide_groups))]
    peptide_final <- rbind(peptide_unique, peptide_final)
  }else if (dpmsr_set$x$peptides_to_use=="Shared"){
    peptide_final <- rbind(peptide_unique, peptide_shared_expand)
  }else{
    peptide_final <- peptide_unique
  }
  
  peptide_final <- peptide_final[order(peptide_final$Master.Protein.Accessions, peptide_final$Sequence),]
  peptide_out <- peptide_final %>% dplyr::select(Confidence, Master.Protein.Accessions, Master.Protein.Descriptions, 
                                                 Sequence, Modifications,
                                                 contains('RT.in.min.by.Search.Engine.'), 
                                                 starts_with('mz.in.Da.by.Search.Engine.'), 
                                                 contains('Charge.by.Search.Engine.'), 
                                                 contains('Percolator.SVM'), 
                                                 contains("Percolator.q.Value"), contains("Abundance.F"))
  
  
  if(ncol(peptide_out) != (10 + dpmsr_set$y$sample_number))
  {
    shinyalert("Oops!", str_c("Number of columns extracted is not as expected ", ncol(peptide_out), "/", (10+dpmsr_set$y$sample_number)), type = "error")  
  }
  
  colnames(peptide_out)[1:10] <- c("Confidence", "Accession", "Description", "Sequence", "Modifications", "Retention.Time","Da","mz", "Ion.Score", "q-Value")
  peptide_out <- subset(peptide_out, Accession %in% master_accessions )
  Simple_Excel(peptide_out, str_c(dpmsr_set$file$extra_prefix,"_ProteinPeptide_to_Peptide_Raw.xlsx", collapse = " "))
  return(peptide_out)
}


#----------------------------------------------------------------------------------------
protein_to_protein <- function(){
  cat(file=stderr(), "protein_to_protein", "\n")
  protein <- dpmsr_set$data$data_raw_protein
  protein <- subset(protein, Master %in% ("IsMasterProtein"))
  protein <- subset(protein, Protein.FDR.Confidence.Combined %in% ("High"))
  
  protein_out <- protein %>% dplyr::select(Accession, Description, Number.of.Protein.Unique.Peptides, 
                                           contains("Abundance"), -contains("Abundance.Count"))
  colnames(protein_out)[1:3] <- c("Accession", "Description", "Unique.Peptides")
  Simple_Excel(protein_out, str_c(dpmsr_set$file$extra_prefix,"_Protein_to_Protein_Raw.xlsx", collapse = " "))
return(protein_out)
}

#----------------------------------------------------------------------------------------
peptide_to_peptide <- function(){
  cat(file=stderr(), "peptide_to_peptide", "\n")
  peptide_groups <- dpmsr_set$data$data_raw_peptide
  peptide_out <- peptide_groups %>% dplyr::select(Confidence, Master.Protein.Accessions, Master.Protein.Descriptions, 
                                                Sequence, Modifications, Positions.in.Master.Proteins, Modifications.in.Master.Proteins,
                                                contains('RT.in.min.by.Search.Engine.'), 
                                                contains('Percolator.SVM'),  
                                                contains("Percolator.q.Value"), contains("Abundance.F"))
  
  if(ncol(peptide_out) != (10 + dpmsr_set$y$sample_number))
  {
    shinyalert("Oops!", "Number of columns extracted is not as expected", type = "error")  
  }
  
  colnames(peptide_out)[1:10] <- c("Confidence", "Accession", "Description", "Sequence", "Modifications", "PositionMaster", "ModificationMaster",
                                  "Retention.Time", "SVM.Score", "q-Value")
  peptide_out <- subset(peptide_out, Confidence %in% ("High"))
  Simple_Excel(peptide_out, str_c(dpmsr_set$file$extra_prefix,"_Peptide_to_Peptide_Raw.xlsx", collapse = " "))
  cat(file=stderr(), "peptide_to_peptide complete", "\n")
  return(peptide_out)
}

#Top.Apex.RT.in.min,
#----------------------------------------------------------------------------------------
isoform_to_isoform <- function(){
  cat(file=stderr(), "isoform_to_isoform", "\n")
  peptide_groups <- dpmsr_set$data$data_raw_isoform
  peptide_out <- try(peptide_groups %>% dplyr::select(contains("Confidence.by"), Master.Protein.Accessions, Master.Protein.Descriptions, 
                                                    Sequence, Modifications, Positions.in.Master.Proteins, Modifications.in.Master.Proteins,
                                                    Top.Apex.RT.in.min, 
                                                    contains('Percolator.SVM'),  
                                                    contains("Percolator.q.Value"), contains("Abundance.F")))
  if (class(peptide_out) == 'try-error') {
    cat(file=stderr(), "column select error - retry", "\n")
    peptide_out <- peptide_groups %>% dplyr::select(contains("Confidence.by"), Master.Protein.Accessions, Master.Protein.Descriptions,
                                                    Sequence, Modifications, Positions.in.Master.Proteins, Modifications.in.Master.Proteins,
                                                    contains('Positions.'),
                                                    contains('RT.in.min.by.'), 
                                                    contains('Percolator.SVM'), 
                                                    contains("Percolator.q.Value"), contains("Abundance.F"))
  }
  
  if(ncol(peptide_out) != (10 + dpmsr_set$y$sample_number))
  {
    shinyalert("Oops!", "Number of columns extracted is not as expected", type = "error")  
  }

  colnames(peptide_out)[1:10] <- c("Confidence", "Accession", "Description", "Sequence", "Modifications", "PositionMaster", "ModificationMaster",
                                   "Retention.Time", "SVM.Score", "q-Value")
  
  peptide_out <- subset(peptide_out, Confidence %in% ("High"))
  Simple_Excel(peptide_out, str_c(dpmsr_set$file$extra_prefix,"_Isoform_to_Isoform_Raw.xlsx", collapse = " "))
  cat(file=stderr(), "isoform_to_isoform complete", "\n")
  return(peptide_out)
}


#peptide_data <- xdf
#--- collapse peptide to protein-------------------------------------------------------------
collapse_peptide <- function(peptide_data){
  info_columns <- ncol(peptide_data) - dpmsr_set$y$sample_number
  peptide_annotate <- peptide_data[1:(info_columns)]
  peptide_data <- peptide_data[(info_columns+1):ncol(peptide_data)]
  peptide_data[is.na(peptide_data)] <- 0
  peptide_annotate <- peptide_annotate[, c("Accession", "Description")]
  peptide_annotate$Peptides <- 1
  peptide_annotate$Peptides <- as.numeric(peptide_annotate$Peptides)
  test1 <- cbind(peptide_annotate, peptide_data)
  #test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(funs(sum))
  test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(list(sum))
  test2 <- data.frame(ungroup(test2))
  #add imputed column info
  if ((dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") 
      && dpmsr_set$x$final_data_output == "Protein" && !as.logical(dpmsr_set$x$tmt_spqc_norm)   ){
            test2 <- add_column(test2, dpmsr_set$data$protein_missing, .after = "Peptides")
            dpmsr_set$y$info_columns_final <<- ncol(test2)-dpmsr_set$y$sample_number
            names(test2)[4] <- "PD_Detected_Peptides"
  }
  return(test2)
} 

#--- collapse peptide to protein-------------------------------------------------------------
collapse_peptide_stats <- function(peptide_data, info_columns){
  peptide_annotate <- peptide_data[1:info_columns]
  peptide_data <- peptide_data[(info_columns+1):ncol(peptide_data)]
  peptide_data[is.na(peptide_data)] <- 0
  peptide_annotate <- peptide_annotate[, c("Accession", "Description")]
  peptide_annotate$Peptides <- 1
  peptide_annotate$Peptides <- as.numeric(peptide_annotate$Peptides)
  test1 <- cbind(peptide_annotate, peptide_data)
  #test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(funs(sum))
  test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(list(sum))
  test2 <- data.frame(ungroup(test2))
  return(test2)
} 

#----------------------------------------------------------------------------------------

psm_set_fdr <- function(){
  psmfdr_dir <- create_dir(str_c(data_dir,"//PSM_FDR"))
  psm_prefix <- str_c(psmfdr_dir, file_prefix)
  
  forward_psm <- dpmsr_set$data$data_raw_psm
  decoy_psm<- dpmsr_set$data$data_raw_decoypsm

  forward_psm$fdr <- rep("forward", nrow(forward_psm))
  decoy_psm$fdr <- rep("decoy", nrow(decoy_psm))
  
  isv1 <- forward_psm$Ions.Score
  isv2 <- decoy_psm$Ions.Score
  isv3 <- c(isv1, isv2)
  
  fdr1 <- forward_psm$fdr
  fdr2 <- decoy_psm$fdr
  fdr3 <- c(fdr1, fdr2)
  
  test <- data.frame(cbind(isv3, fdr3), stringsAsFactors = FALSE)
  colnames(test)<-c("Ions_Score", "FDR")
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
  #check to see if this was already completed, if so skip step
  if("PD_Detected_Peptides" %in% colnames(df)){
    return(df)
  }else{
      #imputed column for protein output
      if ((dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Peptide") 
          && dpmsr_set$x$final_data_output == "Protein"){
        peptide_data <- df
        peptide_annotate <- peptide_data[1:(dpmsr_set$y$info_columns)]
        peptide_data <- peptide_data[(dpmsr_set$y$info_columns+1):ncol(peptide_data)]
        peptide_data[is.na(peptide_data)] <- 0
        peptide_data[peptide_data>0] <- 1
        peptide_annotate <- peptide_annotate[, c("Accession", "Description")]
        test1 <- cbind(peptide_annotate, peptide_data)
        test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(list(sum))
        test3 <- test2
        test3[3:ncol(test2)][test2[3:ncol(test2)] > 1] <- 1
        dpmsr_set$data$Protein_imputed_df <<- test3
        test2 <- data.frame(ungroup(test2))
        test2 <- test2[3:ncol(test2)]
        test2[test2==0] <- "-"
        test2 <- test2 %>% mutate_all(as.character)
        while (ncol(test2)>1) {
          test2[,1] <- str_c(test2[,1], ".", test2[,2])
          test2[,2] <- NULL
        }
        dpmsr_set$data$protein_missing <<- test2[,1]
      }
      
      #imputed column for peptide output
      df_annotation <- df[1:dpmsr_set$y$info_columns]
      df <- df[(dpmsr_set$y$info_columns+1):ncol(df)]
      df_data <- df
      df[df>0] <- "1"
      imputed_df<- df
      imputed_df[is.na(imputed_df)] <- 0
      dpmsr_set$data$peptide_imputed_df <<- cbind(df_annotation, imputed_df)
      df[is.na(df)] <- "-"
      df <- df %>% mutate_all(as.character)
      while (ncol(df)>1) {
        df[,1] <- str_c(df[,1], ".", df[,2])
        df[,2] <- NULL
      }
      
      df_annotation$PD_Detected_Peptides <- df[,1]
      df<- cbind(df_annotation,df_data)
      #add another column to info columns
      return(df)
  }
}



# df<-dpmsr_set$data$data_to_norm
# df$PD_Detected_Peptides<-NULL
# dpmsr_set$y$info_columns<-ncol(df) - dpmsr_set$y$sample_number
