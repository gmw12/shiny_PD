  
run_motifx <- function(input, output, data_in){
  require(stringr)
  require(stringi)
  require(PTMphinder)
  require(rmotifx)

  plot_dir <- str_c(dpmsr_set$file$output_dir, input$select_final_data_motif, "//")
  
  fasta_txt_file <- input$motif_fasta
  fasta_txt_file <- unlist(fasta_txt_file$files[[as.character(0)]][2])
  #---------Create background data ----------------------------------------
  # Set parameters for extractBackground (s - sequence list; c - central character; w - width of motifs
  parsed_ref <- read.table(fasta_txt_file, header=FALSE, row.names=NULL, sep="\t")
  s <- unlist(parsed_ref[,2])
  w <- 15
  pval_motif <- input$pval_motif * 1e-5
  
  if(input$pval_filter != 0 & input$fc_filter != 0) {
    filter_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_motif),"_Pval")] <= as.numeric(input$pval_filter) &
                          (data_in[ ,str_c(as.character(input$select_data_comp_motif),"_FC")] >= as.numeric(input$fc_filter) ))
    FC <- "Up"
    ptm_data <- create_motifx_input(filter_df, parsed_ref)
    motifx_S <- motifx_calc(s, "S", w, "Up", ptm_data, parsed_ref, pval_motif, input$pval_filter, input$fc_filter)
    motifx_T <- motifx_calc(s, "T", w, "Up", ptm_data, parsed_ref, pval_motif, input$pval_filter, input$fc_filter)
    motifx_all <- rbind(motifx_S, motifx_T)
    
    filter_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_motif),"_Pval")] <= as.numeric(input$pval_filter) &
                          (data_in[ ,str_c(as.character(input$select_data_comp_motif),"_FC")] <= -as.numeric(input$fc_filter)))
    FC <- "Down"
    ptm_data <- create_motifx_input(filter_df, parsed_ref)
    motifx_S <- motifx_calc(s, "S", w, "Down", ptm_data, parsed_ref, pval_motif, input$pval_filter, input$fc_filter)
    motifx_T <- motifx_calc(s, "T", w, "Down", ptm_data, parsed_ref, pval_motif, input$pval_filter, input$fc_filter)
    motifx_all <- rbind(motifx_all, motifx_S, motifx_T)
    
    Simple_Excel(motifx_all, str_c(plot_dir, "MotifX_", input$select_data_comp_motif, ".xlsx"))
    }
  
  return(motifx_all)
}





#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

motifx_calc <- function(s, c, w, FC, ptm_data, parsed_ref, pval_motif, pval_filter, fc_filter){

  extractBack_Example <- extractBackground(s, c, w)
  
  # Locate PTMs within full-length proteins and extract neighboring motifs
  phindPTMs_Example <- phindPTMs(ptm_data, parsed_ref)
  
  # Reformat foreground sequences for motif-x
  foreground_Seqs <- unlist(strsplit(phindPTMs_Example[,"Flank_Seq"], split = "[.]"))
  foreground_Seqs_Filtered <- foreground_Seqs[which(lapply(foreground_Seqs, nchar)==15)]
  
  # Run motif-x using foreground and background sequences from PTMphinder functions above
  motifx_data <- motifx(foreground_Seqs_Filtered, extractBack_Example, central.res = c, min.seqs = 20, pval.cutoff = 1e-5)
  
  if (!is.null(motifx_data)){
    motifx_data <- add_column(motifx_data, FC, .before = 1)
    motifx_data <- add_column(motifx_data, pval_filter, .before = 1)
    motifx_data <- add_column(motifx_data, fc_filter, .before = 1)
    motifx_data <- add_column(motifx_data, c, .before = 1)
    motifx_data <- add_column(motifx_data, dpmsr_set$y$comp_groups$comp_name[1], .before = 1)
    names(motifx_data)[1:5] <- c("comparison", "central.res", "foldchange", "pval", "direction")
  }
  
  return(motifx_data)
}



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

create_motifx_input <- function(filter_df, parsed_ref){
  #---------Create input file for PTM data ----------------------------------------

  MotifPhos <- data.frame(cbind(filter_df$Accession, filter_df$Sequence, filter_df$Modifications))
  colnames(MotifPhos) <- c("Accession", "Sequence", "Modifications")
  MotifPhos$Accession <- as.character(MotifPhos$Accession)
  MotifPhos$Sequence <- as.character(MotifPhos$Sequence)
  MotifPhos$Modifications <- as.character(MotifPhos$Modifications)
  
  MotifPhos$Accession <- gsub(";.*", "", MotifPhos$Accession)
  MotifPhos <- subset(MotifPhos, Accession != 1000 )
  
  MotifPhos$PhosOnly <- str_extract(MotifPhos$Modifications, "\\dxP.+\\]")
  MotifPhos$Total_Sites <- substr(MotifPhos$PhosOnly,0,1)
  MotifPhos$PTM_Loc <- str_extract_all(MotifPhos$PhosOnly, "[STY]\\d+")
  MotifPhos$PTM_Loc <- gsub("[(]", "", MotifPhos$PTM_Loc)
  MotifPhos$PTM_Loc <- gsub("[)]", "", MotifPhos$PTM_Loc)
  MotifPhos$PTM_Loc <- gsub("[\"]", "", MotifPhos$PTM_Loc)
  MotifPhos$PTM_Loc <- gsub(",", ";", MotifPhos$PTM_Loc)
  MotifPhos <- subset(MotifPhos, PTM_Loc != 'character0')
  MotifPhos$PTM_Loc <- gsub("c", "", MotifPhos$PTM_Loc)
  MotifPhos$PTM_Loc <- gsub(" ", "", MotifPhos$PTM_Loc)
  
  MotifPhos$PTM_Score <- str_extract_all(MotifPhos$PhosOnly, "[(]\\d+\\.*\\d*")
  MotifPhos$PTM_Score <- gsub("[(]", "", MotifPhos$PTM_Score)
  MotifPhos$PTM_Score <- gsub("[)]", "", MotifPhos$PTM_Score)
  MotifPhos$PTM_Score <- gsub("[\"]", "", MotifPhos$PTM_Score)
  MotifPhos$PTM_Score <- gsub(",", ";", MotifPhos$PTM_Score)
  MotifPhos <- subset(MotifPhos, PTM_Score != 'character0')
  MotifPhos$PTM_Score <- gsub("c", "", MotifPhos$PTM_Score)
  MotifPhos$PTM_Score <- gsub(" ", "", MotifPhos$PTM_Score)
  
  MotifPhos$Identifier <- str_c("PhosPeptide", seq.int(nrow(MotifPhos)))
  
  df <- data.frame(cbind(MotifPhos$Identifier, MotifPhos$Accession, MotifPhos$Sequence, 
                         MotifPhos$Total_Sites, MotifPhos$PTM_Loc, MotifPhos$PTM_Score))
  colnames(df) <- c("Identifier", "Protein_ID", "Peptide_Seq", "Total_Sites", "PTM_Loc", "PTM_Score")
  
  #---------Create background data ----------------------------------------
  ptm_data <- df
  ptm_data$Protein_ID <- gsub( " .*$", "", ptm_data$Protein_ID)
  ptm_data <-subset(ptm_data, Protein_ID %in% parsed_ref$V1)

  return(ptm_data)
}



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

create_phos_database <- function(){
  raw_fasta <- read.csv2("4402_AFF293_051717.fasta", header=FALSE)
  colnames(raw_fasta) <- "fasta"
  raw_fasta$fasta <- as.character(raw_fasta$fasta)
  test_s <- ""
  test_a <- ""
  test_d <- ""
  find_split <- str_locate(raw_fasta[1,], "\\|")
  test_a <- c(test_a, substring(raw_fasta[1,], 2, find_split[1]-1))
  test_d <- c(test_d, substring(raw_fasta[1,], find_split[1]+1))
  temp_seq <-""
  for (i in 2:nrow(raw_fasta)) {
    testa <- substring(raw_fasta[i,],1,1)
    if (identical(testa, ">" )) {
      find_split <- str_locate(raw_fasta[i,], "\\|")
      test_a <- c(test_a, substring(raw_fasta[i,], 2, find_split[1]-1) )
      test_d <- c(test_d, substring(raw_fasta[i,], find_split[1]+1))
      test_s <- c(test_s, temp_seq)
      temp_seq <- ""
    }else{
      temp_seq <- str_c(temp_seq, raw_fasta[i,])
    }
  }
  test_s <- c(test_s, temp_seq)
  
  test_a <- as.data.frame(test_a)
  test_d <- as.data.frame(test_d)
  test_s <- as.data.frame(test_s)
  new_fasta<-(cbind(test_a, test_d, test_s))
  new_fasta <- as.data.frame(new_fasta[-1,])
  colnames(new_fasta) <- c("Accession", "Description", "Sequence")
  new_fasta$Accession <- as.character(new_fasta$Accession)
  new_fasta$Description <- as.character(new_fasta$Description)
  new_fasta$Sequence <- as.character(new_fasta$Sequence)
  Simple_Excel(new_fasta, "4402_AFF293_051717.xlsx")
  

}


