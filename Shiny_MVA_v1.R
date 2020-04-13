  #create data frame for comparisons
set_mva_groups <- function(session, input, output){
  mva_data <- dpmsr_set$data$impute[[input$select_final_data_mva]]
  mva_data <- mva_data[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
  if(input$select_final_data_mva == "impute"){
    colnames(mva_data) <- dpmsr_set$design$Header3
  }else{
    colnames(mva_data) <- dpmsr_set$design$Header2
  }
  comp_groups <- data.frame(seq(from=1, to=as.numeric(input$mva_comp)))
  colnames(comp_groups) <- "CompNumber"
  color_choices <- distinctColorPalette(as.numeric(input$mva_comp)*2)
  color_choices_N <- color_choices[c(TRUE, FALSE)]
  color_choices_D <- color_choices[c(FALSE, TRUE)]
  
  for (i in 1:input$mva_comp){
    mva_data_N <- mva_data
    mva_data_D <- mva_data
    
    for(mva_group in input[[str_c('var_',i,'N')]]){mva_data_N <- mva_data_N %>% dplyr::select(contains(mva_group))   }
    for(mva_group in input[[str_c('var_',i,'D')]]){mva_data_D <- mva_data_D %>% dplyr::select(contains(mva_group))   }
    
    comp_groups$comp_N[i] <- paste(input[[str_c('var_',i,'N')]], collapse = "_")
    comp_groups$comp_D[i] <- paste(input[[str_c('var_',i,'D')]], collapse = "_")
    comp_groups$comp_name[i] <- str_c(comp_groups$comp_N[i], "_v_", comp_groups$comp_D[i])
    comp_groups$com_N_headers[i] <- list(colnames(mva_data_N))
    comp_groups$com_D_headers[i] <- list(colnames(mva_data_D))
    comp_groups$fc[i] <- str_c(comp_groups$comp_name[i], "_FC")
    comp_groups$fc2[i] <- str_c(comp_groups$comp_name[i], "_FC2")
    comp_groups$pval[i] <- str_c(comp_groups$comp_name[i], "_pval")
    comp_groups$limma_pval[i] <- str_c(comp_groups$comp_name[i], "_limma_pval")
    comp_groups$exactTest[i] <- str_c(comp_groups$comp_name[i], "_exactTest")
    comp_groups$mf[i] <- str_c(comp_groups$comp_name[i], "_MF")
    comp_groups$N_count[i] <- ncol(mva_data_N)
    comp_groups$D_count[i] <- ncol(mva_data_D)
    comp_groups$N_color[1] <- color_choices_N[i]
    comp_groups$D_color[1] <- color_choices_D[i]
    updatePickerInput(session, inputId = str_c('var_',i,'N'), label=str_c("Comp",i,"_N  selected -> ", ncol(mva_data_N)) )
    updatePickerInput(session, inputId = str_c('var_',i,'D'), label=str_c("Comp",i,"_D  selected -> ", ncol(mva_data_D)) )
    dpmsr_set$y$mva[[str_c("comp", i, "_N")]] <<- input[[str_c('var_',i,'N')]]
    dpmsr_set$y$mva[[str_c("comp", i, "_D")]] <<- input[[str_c('var_',i,'D')]]
    comp_groups$sample_numbers_N[i] <- list(match(colnames(mva_data_N), names(mva_data)))
    comp_groups$sample_numbers_D[i] <- list(match(colnames(mva_data_D), names(mva_data)))
  }
  
  dpmsr_set$y$mva$groups <<- comp_groups
  dpmsr_set$y$mva$comp_number <<- input$mva_comp
  
}


#--------------------------------------------------------------------------------------------------------------------------------
#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
mva_stat_calc <- function(session, input, output) {
  data_in <- dpmsr_set$data$impute[[input$select_final_data_mva]]
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    data_in <- collapse_peptide(data_in)
  }
  if (as.logical(dpmsr_set$x$peptide_ptm_out)){
    data_in <- data_in [grep(dpmsr_set$x$peptide_report_grep, data_in$Modifications),]
  }
  
  annotate_in <- data_in[1:dpmsr_set$y$info_columns_final]
  data_in <- data_in[(dpmsr_set$y$info_columns_final+1):ncol(data_in)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  #generate %CV's for each group
  for(i in 1:nrow(dpmsr_set$y$mva$groups)) 
  {
    stat_df[ , str_c(dpmsr_set$y$mva$groups$comp_N[i], "_CV")] <- percentCV_gw(data_in %>% dplyr::select(contains(dpmsr_set$y$mva$groups$comp_N[i])))
    stat_df[ , str_c(dpmsr_set$y$mva$groups$comp_D[i], "_CV")] <- percentCV_gw(data_in %>% dplyr::select(contains(dpmsr_set$y$mva$groups$comp_D[i])))
  } 

  #generate pvalue and FC for each comparison
  for(i in 1:nrow(dpmsr_set$y$mva$groups))
  {
    comp_N_data <- data_in %>% dplyr::select(contains(dpmsr_set$y$mva$groups$comp_N[i]))
    comp_D_data <- data_in %>% dplyr::select(contains(dpmsr_set$y$mva$groups$comp_D[i]))
    stat_df[ ,dpmsr_set$y$mva$groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$mva$groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$mva$groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
    #stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
    #stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
  } 
  
  # Create tables for excel--------------------------------------------------
  if(input$select_final_data_mva == "impute"){
    colnames(data_in) <- dpmsr_set$design$Header3
  }else{
    colnames(data_in) <- dpmsr_set$design$Header2
  }
  data_table <- cbind(annotate_in, data_in, stat_df[2:ncol(stat_df)])
  dpmsr_set$data$mva$final <<- data_table
  dpmsr_set$y$mva$pvalue_cutoff <<- input$mva_pvalue_cutoff
  dpmsr_set$y$mva$foldchange_cutoff <<- input$mva_foldchange_cutoff
  return()
}

#----------------------------------------------------------------------------------------
# create final excel documents
MVA_Final_Excel <- function() {
  require(openxlsx)
    
    filename <- str_c(dpmsr_set$file$extra_prefix, "_MVA_final.xlsx")
    df <- dpmsr_set$data$mva$final
    #remove FC2 from df for excel
    df <- df[,-grep(pattern="_FC2", colnames(df))]
    
    if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
      df_raw <- dpmsr_set$data$finalraw
    }else{
      df_raw <- dpmsr_set$data$finalraw
      #if PTM out need to reduce raw data frame for excel
      if (as.logical(dpmsr_set$x$peptide_ptm_out)){
        df_raw <- df_raw[grep(dpmsr_set$x$peptide_report_grep, df_raw$Modifications),]
      }
    }
    
    # save the option to save normalized table with raw
    #df2 <- cbind(df_raw, df[(dpmsr_set$y$info_columns_final+1):ncol(df)])
    
    df2 <- df
    
    if (as.logical(dpmsr_set$x$accession_report_out)){
      df <-subset(df, Accession %in% dpmsr_set$x$accession_report_list )
      df2 <-subset(df2, Accession %in% dpmsr_set$x$accession_report_list )
    }
    
    nextsheet <- 1
    
    wb <- createWorkbook()
    
    if(dpmsr_set$x$raw_data_input=="Protein_Peptide" || dpmsr_set$x$raw_data_input=="Protein"){
      
      raw_protein <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Protein_to_Protein_Raw.xlsx"))
      raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_ProteinPeptide_to_Peptide_Raw.xlsx"))
      
      if(as.logical(dpmsr_set$x$peptide_isoform) && dpmsr_set$x$raw_data_input=="Peptide"){
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Isoform_to_Isoform_Raw.xlsx"))
      }
      
      if(!as.logical(dpmsr_set$x$peptide_isoform) && dpmsr_set$x$raw_data_input=="Peptide"){
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Peptide_to_Peptide_Raw.xlsx"))
      }
      
      addWorksheet(wb, "Raw Peptide Data")
      writeData(wb, sheet = nextsheet, raw_peptide)
      nextsheet <- nextsheet +1
      addWorksheet(wb, "Raw Protein Data")
      writeData(wb, sheet = nextsheet, raw_protein) 
      nextsheet <- nextsheet +1
    }else if (dpmsr_set$x$raw_data_input=="Peptide"){
      if(as.logical(dpmsr_set$x$peptide_isoform)){
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Isoform_to_Isoform_Raw.xlsx"))
      }else{
        raw_peptide <- read_excel(str_c(dpmsr_set$file$extra_prefix,"_Peptide_to_Peptide_Raw.xlsx"))
      }
      addWorksheet(wb, "Raw Peptide Data")
      writeData(wb, sheet = nextsheet, raw_peptide) 
      nextsheet <- nextsheet +1
    }
    
    #addWorksheet(wb, deparse(substitute(df2)))
    addWorksheet(wb, "Normalized Data")
    writeData(wb, sheet = nextsheet, df2)  
    nextsheet <- nextsheet +1
    z=1
    for(i in 1:as.numeric(dpmsr_set$y$mva$comp_number))  {
      comp_string <- dpmsr_set$y$mva$groups$comp_name[i]
      df_N <- df %>% dplyr::select(unlist(dpmsr_set$y$mva$groups$com_N_headers[i]))  
      df_D <- df %>% dplyr::select(unlist(dpmsr_set$y$mva$groups$com_D_headers[i]) ) 
      df_Comp <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$comp_name [i]) ) 
      df_final <- cbind(df[1:(dpmsr_set$y$info_columns_final)], df_N, df_D, df_Comp)
      
      filtered_df <- subset(df_final, df_final[ , dpmsr_set$y$mva$groups$pval[i]] <= as.numeric(dpmsr_set$y$mva$pvalue_cutoff) &  
                              (df_final[ , dpmsr_set$y$mva$groups$fc[i]] >= as.numeric(dpmsr_set$y$mva$foldchange_cutoff) | 
                                 df_final[ , dpmsr_set$y$mva$groups$fc[i]] <= -as.numeric(dpmsr_set$y$mva$foldchange_cutoff))) 
      addWorksheet(wb, comp_string)
      writeData(wb, sheet = nextsheet, filtered_df)
      nextsheet <- nextsheet +1
      z <- z+2
    }
    saveWorkbook(wb, filename, overwrite = TRUE)
    
}






