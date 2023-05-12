  
run_motifx <- function(input, output, data_in){
  cat(file=stderr(), "start motifx..." , "\n") 
  require(stringr)
  require(stringi)
  require(PTMphinder)
  require(rmotifx)

  #---------Create background data ----------------------------------------
  # Set parameters for extractBackground (s - sequence list; c - central character; w - width of motifs
  parsed_ref <- dpmsr_set$data$phos$background
  s <- parsed_ref$Sequence
  w <- 15
  pval_motif <- input$pval_motif * 1e-5
  min_seq <- as.integer(input$motif_min_seq)
  
   #parsed_refx <<- parsed_ref
   #sx <<- s
   #wx <<- w
   ##pval_motifx <<- pval_motif
   #data_inx <<- data_in
   pvalue_cutoff <<- input$pvalue_cutoff
   foldchange_cutoff <<- input$foldchange_cutoff
   select_data_comp_motif <<- input$select_data_comp_motif
  
  
  
  
  motifx_up <- future({
  cat(file=stderr(), "motifx up..." , "\n") 
    filter_df <- subset(data_in, data_in$Stats == "Up" ) 
      if (nrow(filter_df > 0)){
        FC <- "Up"
        ptm_data <- create_motifx_input(filter_df, parsed_ref, select_data_comp_motif, "Up")
        cat(file=stderr(), str_c("ptm_data has ", nrow(ptm_data), " entries") , "\n") 
        motifx_S <- motifx_calc(s, "S", w, "Up", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
        cat(file=stderr(), str_c("motifx_S has ", nrow(motifx_S), " entries") , "\n") 
        motifx_T <- motifx_calc(s, "T", w, "Up", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
        cat(file=stderr(), str_c("motifx_T has ", nrow(motifx_T), " entries") , "\n") 
        motifx_Y <- motifx_calc(s, "Y", w, "Up", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
        cat(file=stderr(), str_c("motifx_Y has ", nrow(motifx_Y), " entries") , "\n") 
        motifx_up <- rbind(motifx_S, motifx_T, motifx_Y)
        cat(file=stderr(), str_c("motifx_up has ", nrow(motifx_up), " entries") , "\n") 
      }else{
        motifx_up <- NULL
        shinyalert("Oops!", "No peptides passed stats for Up", type = "error")
      }
      motifx_up
  },seed=TRUE)
   
  motifx_down <- future({
  cat(file=stderr(), "motifx down..." , "\n")     
    filter_df <- subset(data_in, data_in$Stats == "Down" )
    if (nrow(filter_df > 0)){
      FC <<- "Down"
      ptm_data <- create_motifx_input(filter_df, parsed_ref, select_data_comp_motif, "Down")
      motifx_S <- motifx_calc(s, "S", w, "Down", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
      motifx_T <- motifx_calc(s, "T", w, "Down", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
      motifx_Y <- motifx_calc(s, "Y", w, "Down", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
      motifx_down <- rbind(motifx_S, motifx_T, motifx_Y)
    }else{
      motifx_down <- NULL
      shinyalert("Oops!", "No peptides passed stats for Down", type = "error")
    }
    motifx_down
  },seed=TRUE)
  
  motifx_updown <- future({  
  cat(file=stderr(), "motifx updown..." , "\n")  
    filter_df<- subset(data_in, data_in$Stats == "Up" | data_in$Stats == "Down") 
    if (nrow(filter_df > 0)){
      FC <- "UpDown"
      ptm_data <- create_motifx_input(filter_df, parsed_ref, select_data_comp_motif, "UpDown")
      motifx_S <- motifx_calc(s, "S", w, "UpDown", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
      motifx_T <- motifx_calc(s, "T", w, "UpDown", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
      motifx_Y <- motifx_calc(s, "Y", w, "UpDown", ptm_data, parsed_ref, pval_motif, min_seq, pvalue_cutoff, foldchange_cutoff, select_data_comp_motif)
      motifx_updown <- rbind(motifx_S, motifx_T, motifx_Y)
    }else{
      motifx_updown <- NULL
      shinyalert("Oops!", "No peptides passed stats for UpDown", type = "error")
    }
    motifx_updown
  },seed=TRUE)
  
  cat(file=stderr(), str_c("compiling motifs... ", nrow(value(motifx_up)), "... ", nrow(value(motifx_down)), "... ", nrow(value(motifx_updown)) ) , "\n")   

  
    motifx_all <- rbind(value(motifx_up), value(motifx_down), value(motifx_updown))
    
    if(is.null(motifx_all)){
      shinyalert("Oops!", "No motifx found", type = "error")
    }else {
      cat(file=stderr(), "write motifx excel..." , "\n")  
      # save filename to dpmsr_set so downloadhandler can grab it
      dpmsr_set$data$phos$filename <<- str_c("MotifX_", input$select_data_comp_motif, ".xlsx")
      if (site_user=="dpmsr"){
        Simple_Excel_bg(motifx_all, "MotifX", str_c(dpmsr_set$file$phos, "MotifX_", input$select_data_comp_motif, ".xlsx"))
      }
    }
    
    cat(file=stderr(), "motifx... complete" , "\n") 
    gc()
  return(motifx_all)
}





#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

motifx_calc <- function(s, c, w, FC, ptm_data, parsed_ref, pval_motif, min_seq, pval_filter, fc_filter, comparison){
  
  extractBack_Example <- extractBackground(s, c, w)
  
  # Locate PTMs within full-length proteins and extract neighboring motifs
  phindPTMs_Example <- phindPTMs(ptm_data, parsed_ref)

  
  # Reformat foreground sequences for motif-x
  foreground_Seqs <- unlist(strsplit(phindPTMs_Example[,"Flank_Seq"], split = "[.]"))
  
  foreground_Seqs_Filtered <- foreground_Seqs[which(lapply(foreground_Seqs, nchar)==15)]
  
  cat(file=stderr(), str_c("Number of foreground sequences for = ", length(foreground_Seqs_Filtered)), "\n")  
  
  # Run motif-x using foreground and background sequences from PTMphinder functions above
  motifx_data <- try(motifx(foreground_Seqs_Filtered, extractBack_Example, central.res = c, min.seqs = min_seq, pval.cutoff = pval_motif))
  
  if ((!is.null(motifx_data)) & (class(motifx_data) != "try-error")){
    motifx_data <- add_column(motifx_data, FC, .before = 1)
    motifx_data <- add_column(motifx_data, pval_filter, .before = 1)
    motifx_data <- add_column(motifx_data, fc_filter, .before = 1)
    motifx_data <- add_column(motifx_data, c, .before = 1)
    motifx_data <- add_column(motifx_data, comparison, .before = 1)
    names(motifx_data)[1:5] <- c("comparison", "central.res", "foldchange", "pval", "direction")
  } else {
    motifx_data <- NULL
  }
  
  cat(file=stderr(), str_c("Number of motifs for ", c," = ", nrow(motifx_data)), "\n")
  
  if (site_user=="dpmsr"){
    Simple_Excel_bg(phindPTMs_Example, "MotifX", str_c(dpmsr_set$file$phos, "MotifX_phindPTM_", comparison, "_", FC, ".xlsx"))
  }
  gc()
  return(motifx_data)
} 



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

create_motifx_input <- function(filter_df, parsed_ref, comparison, direction_filter){
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
  
  #--------- Subset to insure proteins are in database ----------------------------------------
  ptm_data <- df
  ptm_data$Protein_ID <- gsub( " .*$", "", ptm_data$Protein_ID)
  ptm_data <-subset(ptm_data, Protein_ID %in% parsed_ref$Accession)
  
  if (site_user=="dpmsr"){
    cat(file=stderr(), "write ptm_data to excel..." , "\n")
    Simple_Excel_bg(ptm_data, "MotifX", str_c(dpmsr_set$file$phos, "MotifX_ptmdata_", comparison, "_", direction_filter, ".xlsx"))
  }
  
  return(ptm_data)
}



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

create_phos_database <- function(session, input, output){
  cat(file=stderr(), "Create MotifX sequences start...", "\n")
  dpmsr_set$x$motif_set <<- 1
  
  require(stringr)

  fasta_path <- parseFilePaths(volumes, input$motif_fasta)
  cat(file=stderr(), str_c("fasta_path... ",fasta_path), "\n")
  
  #direct read for legacy dpmsr_set files
  #fasta_txt_file <- "/home/dpmsr/shared/mouse.fasta"
  fasta_txt_file <- fasta_path$datapath
  
  #raw_fasta <- read.csv2(fasta_txt_file, header=FALSE)
  raw_fasta <- data.frame(read_lines(fasta_txt_file)  )
  
  raw_fasta <- data.frame(lapply(raw_fasta, function(x) gsub(">sp\\|", ">", x)), stringsAsFactors=F)
             
  if(input$accession_split == "Bar"){
    split_char <- "\\|"
  }else  if(input$accession_split == "Space"){
    split_char <- " "
  }
    
  # \\| space
  
  colnames(raw_fasta) <- "fasta"
  raw_fasta$fasta <- as.character(raw_fasta$fasta)
  test_s <- ""
  test_a <- ""
  find_split <- str_locate(raw_fasta[1,], split_char)
  test_a <- c(test_a, substring(raw_fasta[1,], 2, find_split[1]-1))
  temp_seq <-""
  for (i in 2:nrow(raw_fasta)) {
    testa <- substring(raw_fasta[i,],1,1)
    if (identical(testa, ">" )) {
      find_split <- str_locate(raw_fasta[i,], split_char)
      test_a <- c(test_a, substring(raw_fasta[i,], 2, find_split[1]-1) )
      test_s <- c(test_s, temp_seq)
      temp_seq <- ""
    }else{
      temp_seq <- str_c(temp_seq, raw_fasta[i,])
    }
  }
  test_s <- c(test_s, temp_seq)
  
  test_a <- as.data.frame(test_a)
  test_s <- as.data.frame(test_s)
  new_fasta<-(cbind(test_a, test_s))
  new_fasta <- as.data.frame(new_fasta[-1,])
  colnames(new_fasta) <- c("Accession", "Sequence")
  new_fasta$Sequence <- as.character(new_fasta$Sequence)
  
  file_name <- fasta_path$name
  file_name <- str_replace(file_name, ".fasta", ".txt")
  
  dpmsr_set$data$phos$background <<- new_fasta
  
  if (site_user=="dpmsr"){
    write.csv2(new_fasta, str_c(dpmsr_set$file$phos, file_name), row.names = FALSE)
  }
  
  cat(file=stderr(), "Create MotifX sequences end...", "\n")
}

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

create_phos_database_custom <- function(){
  require(stringr)
  raw_fasta <- read.csv2("/Users/gregwaitt/Documents/Data/5560_crov2_Sep2020.fasta", header=FALSE)
  colnames(raw_fasta) <- "fasta"
  raw_fasta$fasta <- as.character(raw_fasta$fasta)
  test_s <- ""
  test_a <- ""
  test_d <- ""
  find_split <- str_locate(raw_fasta[1,], " ")
  test_a <- c(test_a, substring(raw_fasta[1,], 2, find_split[1]-1))
  test_d <- c(test_d, substring(raw_fasta[1,], find_split[1]+1))
  temp_seq <-""
  for (i in 2:nrow(raw_fasta)) {
    testa <- substring(raw_fasta[i,],1,1)
    if (identical(testa, ">" )) {
      find_split <- str_locate(raw_fasta[i,], " ")
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
  if (site_user=="dpmsr"){
    Simple_Excel_bg(new_fasta, "MotifX", "5560_crov2_Sep2020.xlsx")
  }
}

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
testonly <- function(){
  parsed_ref <- parsed_refx
  s <- sx
  w <- wx
  pval_motif <- pval_motifx
  data_in <- data_inx
  
  parsed_ref <- dpmsr_set$data$phos$background
  names(parsed_ref) <- c("V1", "V2")
  s <- list(parsed_ref$Sequence)
  s <- unlist(s)
  
  s <- unlist(parsed_ref[,2])
  extractBack_Example <- extractBackground(s, "S", 15)
  
  filter_df <- subset(data_in, data_in$Stats == "Up" ) 
  FC <- "Up"
  ptm_data <- create_motifx_input(filter_df, parsed_ref)
  motifx_S <- motifx_calc(s, "S", w, "Up", ptm_data, parsed_ref, pval_motif, input$pvalue_cutoff, input$foldchange_cutoff, input$select_data_comp_motif)
  
  test_parsed_ref <- read.table("Mmusculus_110419.txt", header=FALSE, row.names=NULL, sep="\t")                     
  test_s <- unlist(test_parsed_ref[,2])
  extractBack_Example <- extractBackground(test_s, "S", 15)
  
  parsed_ref <- dpmsr_set$data$phos$background
  names(parsed_ref) <- c("V1", "V2")
  parsed_ref <- rbind(c("Accession", "Sequence"), parsed_ref)
  s <- unlist(parsed_ref[,2])
  extractBack_Example <- extractBackground(s, "S", 15)
  
  s<-parsed_ref$Sequence
  
  
  
}
#--------------------------------------------------------------------------------------------

create_padded_seq <- function(){
  
  filter_df <- dpmsr_set$data$final$sltmm
  
  parsed_ref <- dpmsr_set$data$phos$background
  
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
  
  #--------- Subset to insure proteins are in database ----------------------------------------
  ptm_data <- df
  ptm_data$Protein_ID <- gsub( " .*$", "", ptm_data$Protein_ID)
  ptm_data <-subset(ptm_data, Protein_ID %in% parsed_ref$Accession)
  
  if (site_user=="dpmsr"){
    cat(file=stderr(), "write ptm_data to excel..." , "\n") 
    Simple_Excel_bg(ptm_data, "MotifX", str_c(dpmsr_set$file$phos, "Padded_Seq.xlsx"))
  }
  
  phindPTMs_Example <- phindPTMs(ptm_data, parsed_ref)
  df_PTMs <- data.frame(phindPTMs_Example)
  
  ptm_data$Prot_Loc <- df_PTMs$Prot_Loc
  ptm_data$Flank_Seq <- df_PTMs$Flank_Seq
  ptm_data$Ambiguity <- df_PTMs$Ambiguity
  ptm_data$Prot_Seq <- df_PTMs$Prot_Seq
  
  if (site_user=="dpmsr"){
    Simple_Excel_bg(ptm_data, "MotifX", str_c(dpmsr_set$file$phos, "Padded_Seq.xlsx"))
  }
  
  return(ptm_data)
}





