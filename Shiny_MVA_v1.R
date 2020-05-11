
#--------------------------------------------------------------------------------------
stat_prep <- function(){
  ncores <- detectCores()
  if (is.na(ncores)) {ncores <- 1}
  data_list <- dpmsr_set$data$impute
  dpmsr_set$data$final <<- mclapply(data_list, stat_prep_parallel, mc.cores = ncores)
  
  #check that the parallel processing went through, if not do it again one at a time
  check_stats_prep_parallel(data_list)
  
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    dpmsr_set$data$finalraw <<- collapse_peptide(dpmsr_set$data$normalized$impute)
    Simple_Excel(dpmsr_set$data$finalraw, str_c(dpmsr_set$file$extra_prefix, "_final_raw_protein.xlsx", collapse = " "))
  } else{
    dpmsr_set$data$finalraw <<- dpmsr_set$data$normalized$impute
    Simple_Excel(dpmsr_set$data$finalraw, str_c(dpmsr_set$file$extra_prefix, "_final_raw_peptide.xlsx", collapse = " "))
  }
}

#--------------------------------------------------------------------------------------
save_final_nostats <- function(){
  #save to excel final without stats
  for (name in names(dpmsr_set$data$final)){
    cat(file=stderr(), str_c("Saving final excel w/o stats....", name), "\n")
    filename <- str_c(dpmsr_set$file$output_dir, name, "//Final_", name, "_noStats.xlsx")
    Simple_Excel_name(dpmsr_set$data$final[[name]], filename, name)
  }
  
  for (name in names(dpmsr_set$data$impute)){
    cat(file=stderr(), str_c("Saving imputed peptide excel w/o stats....", name), "\n")
    filename <- str_c(dpmsr_set$file$output_dir, name, "//Imputed_Peptide_", name, ".xlsx")
    Simple_Excel_name(dpmsr_set$data$impute[[name]], filename, name)
  }
}

#--------------------------------------------------------------------------------------
stat_prep_parallel <- function(data_list){
  stat_prep_collapse(data_list)
}

#--------------------------------------------------------------------------------------
stat_prep_collapse <- function(data_in) {
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    data_in <- collapse_peptide(data_in)
  }
  
  if (as.logical(dpmsr_set$x$peptide_ptm_out)){
    data_in <- data_in [grep(dpmsr_set$x$peptide_report_grep, data_in$Modifications),]
  }
  
  info_columns_final <- ncol(data_in)-dpmsr_set$y$sample_number
  annotate_in <- data_in[1:info_columns_final]
  data_in <- data_in[(info_columns_final+1):ncol(data_in)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  #generate %CV's for each group
  for(i in 1:dpmsr_set$y$group_number) 
  {
    stat_df[ , str_c(dpmsr_set$y$sample_groups$Group[i], "_CV")] <- percentCV_gw(data_in[dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i]])
  } 
  
  stat_df <- stat_df[2:ncol(stat_df)]
  colnames(data_in) <- dpmsr_set$design$Header2
  data_table <- cbind(annotate_in, data_in, stat_df)
  
  
  return(data_table)
}

#--------------------------------------------------------------------------------------
check_stats_prep_parallel <- function(data_list){
  for(data_name in names(data_list)){
    if(is.null(dpmsr_set$data$final[[data_name]]    )){
      cat(file=stderr(), str_c("Stats Prep Parallel function...", data_name), "\n")
      dpmsr_set$data$final[[data_name]] <<- stat_prep_collapse(dpmsr_set$data$impute[[data_name]])
    }
  }
}

#--------------------------------------------------------------------------------------

 #create data frame for comparisons
set_stat_groups <- function(session, input, output){
  cat(file=stderr(), "set_stat_groups....1", "\n")
  stats_data <- dpmsr_set$data$final[[input$select_final_data_stats]]
  stats_data <- stats_data[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
  
  if(input$select_final_data_stats == "impute"){
    colnames(stats_data) <- dpmsr_set$design$Header3
  }else{
    colnames(stats_data) <- dpmsr_set$design$Header2
  }
  
  #save spqc information
  dpmsr_set$y$stats$comp_spqc <<- input$comp_spqc
  spqc_df <- stats_data %>% dplyr::select(contains(input$comp_spqc))
  dpmsr_set$y$stats$comp_spqc_sample_numbers <<- list(match(colnames(spqc_df), names(stats_data)))
  dpmsr_set$y$stats$comp_spqc_sample_headers <<- list(colnames(spqc_df))
  
  comp_groups <- data.frame(seq(from=1, to=as.numeric(input$comp_number)))
  colnames(comp_groups) <- "CompNumber"
  color_choices <- distinctColorPalette(as.numeric(input$comp_number)*2)
  color_choices_N <- color_choices[c(TRUE, FALSE)]
  color_choices_D <- color_choices[c(FALSE, TRUE)]
  cat(file=stderr(), "set_stats_groups....2", "\n")
  
  for (i in 1:input$comp_number){
    cat(file=stderr(), str_c("stats_comp ", i, " of ", input$comp_number), "\n")
    stats_data_N <- stats_data
    stats_data_D <- stats_data
    
    for(stats_group in input[[str_c('comp_',i,'N')]]){stats_data_N <- stats_data_N %>% dplyr::select(contains(stats_group))   }
    for(stats_group in input[[str_c('comp_',i,'D')]]){stats_data_D <- stats_data_D %>% dplyr::select(contains(stats_group))   }
    
    comp_groups$comp_N[i] <- paste(input[[str_c('comp_',i,'N')]], collapse = "_")
    comp_groups$comp_D[i] <- paste(input[[str_c('comp_',i,'D')]], collapse = "_")
    comp_groups$comp_name[i] <- str_c(comp_groups$comp_N[i], "_v_", comp_groups$comp_D[i])
    comp_groups$com_N_headers[i] <- list(colnames(stats_data_N))
    comp_groups$com_D_headers[i] <- list(colnames(stats_data_D))
    comp_groups$fc[i] <- str_c(comp_groups$comp_name[i], "_FC")
    comp_groups$fc2[i] <- str_c(comp_groups$comp_name[i], "_FC2")
    comp_groups$pval[i] <- str_c(comp_groups$comp_name[i], "_pval")
    comp_groups$limma_pval[i] <- str_c(comp_groups$comp_name[i], "_limma_pval")
    comp_groups$exactTest[i] <- str_c(comp_groups$comp_name[i], "_exactTest")
    comp_groups$mf[i] <- str_c(comp_groups$comp_name[i], "_MF")
    comp_groups$cv_N[i] <- str_c(comp_groups$comp_N[i], "_CV")  
    comp_groups$cv_D[i] <- str_c(comp_groups$comp_D[i], "_CV")  
    comp_groups$N_count[i] <- ncol(stats_data_N)
    comp_groups$D_count[i] <- ncol(stats_data_D)
    comp_groups$N_color[i] <- color_choices_N[i] #first i was 1
    comp_groups$D_color[i] <- color_choices_D[i] #first i was 1
    updatePickerInput(session, inputId = str_c('comp_',i,'N'), label=str_c("Comp",i,"_N  selected -> ", ncol(stats_data_N)) )
    updatePickerInput(session, inputId = str_c('comp_',i,'D'), label=str_c("Comp",i,"_D  selected -> ", ncol(stats_data_D)) )
    updateTextInput(session, inputId = str_c("volcano",i,"_stats_plot_title"),  value = str_c("Volcano ", comp_groups$comp_name[i]))
    dpmsr_set$y$stats[[str_c("comp", i, "_N")]] <<- input[[str_c('comp_',i,'N')]]
    dpmsr_set$y$stats[[str_c("comp", i, "_D")]] <<- input[[str_c('comp_',i,'D')]]
    comp_groups$sample_numbers_N[i] <- list(match(colnames(stats_data_N), names(stats_data)))
    comp_groups$sample_numbers_D[i] <- list(match(colnames(stats_data_D), names(stats_data)))
  }
  
  cat(file=stderr(), "set_stats_groups....4", "\n")
  dpmsr_set$y$stats$groups <<- comp_groups
  dpmsr_set$y$stats$comp_number <<- input$comp_number
  
  updatePickerInput(session, "stats_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name)
  updatePickerInput(session, "stats_oneprotein_plot_comp", choices = dpmsr_set$y$stats$groups$comp_name, selected = dpmsr_set$y$stats$groups$comp_name[1])
  updateTextInput(session, "stats_barplot_title", value = str_c("Barplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  updateTextInput(session, "stats_boxplot_title", value = str_c("Boxplot ", dpmsr_set$y$stats$groups$comp_name[input$mva_plot_comp]))
  updateSelectInput(session, "stats_select_data_comp", choices = dpmsr_set$y$stats$groups$comp_name, 
                    selected = dpmsr_set$y$stats$groups$comp_name[1])

}




#--------------------------------------------------------------------------------------------------------------------------------
#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_calc <- function(session, input, output) {
  data_in <- dpmsr_set$data$final[[input$select_final_data_stats]]

  annotate_in <- data_in[1:dpmsr_set$y$info_columns_final]
  data_in <- data_in[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  imputed_df <- dpmsr_set$data$Protein_imputed_df[3:ncol(dpmsr_set$data$Protein_imputed_df)]  
  
  for(i in 1:nrow(dpmsr_set$y$stats$groups)) 
  {
    stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_N[i], "_CV")] <- percentCV_gw(data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ])
    stat_df[ , str_c(dpmsr_set$y$stats$groups$comp_D[i], "_CV")] <- percentCV_gw(data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ])
  } 

  stat_df[ , str_c(dpmsr_set$y$stats$comp_spqc, "_CV")] <- percentCV_gw(data_in[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)])
  
  #generate pvalue and FC for each comparison
  for(i in 1:nrow(dpmsr_set$y$stats$groups))
  {
    comp_N_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_data <- data_in[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    stat_df[ ,dpmsr_set$y$stats$groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$stats$groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
    #stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
    #stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
    comp_N_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_N[i]) ]
    comp_D_imputed <- imputed_df[,unlist(dpmsr_set$y$stats$groups$sample_numbers_D[i]) ]
    stat_df[ ,dpmsr_set$y$stats$groups$mf[i]] <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
  } 
  
  # Create final tables--------------------------------------------------
  if(input$select_final_data_stats == "impute"){
    colnames(data_in) <- dpmsr_set$design$Header3
  }else{
    colnames(data_in) <- dpmsr_set$design$Header2
  }
  data_table <- cbind(annotate_in, data_in, stat_df[2:ncol(stat_df)])
  dpmsr_set$data$stats$final <<- data_table
  
  #single shot observers
  dpmsr_set$data$stats$final_comp <<- input$select_final_data_stats
  dpmsr_set$y$stats$pvalue_cutoff <<- input$stats_pvalue_cutoff
  dpmsr_set$y$stats$foldchange_cutoff <<- input$stats_foldchange_cutoff
  dpmsr_set$y$stats$stats_spqc_cv_filter <<- input$stats_spqc_cv_filter
  dpmsr_set$y$stats$stats_spqc_cv_filter_factor <<- input$stats_spqc_cv_filter_factor
  dpmsr_set$y$stats$stats_comp_cv_filter <<- input$stats_spqc_cv_filter
  dpmsr_set$y$stats$stats_comp_cv_filter_factor <<- input$stats_spqc_cv_filter_factor
  return()
}








#----------------------------------------------------------------------------------------
# create final excel documents
stats_Final_Excel <- function(session, input, output) {
  require(openxlsx)
    
    filename <- str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp, "//", input$final_stats_name)
    df <- dpmsr_set$data$stats$final
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
      
      addWorksheet(wb, "Sample Info")
      writeData(wb, sheet = nextsheet, dpmsr_set$protocol)
      nextsheet <- nextsheet +1
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
    for(i in 1:as.numeric(dpmsr_set$y$stats$comp_number))  {
      comp_string <- dpmsr_set$y$stats$groups$comp_name[i]
      
      filtered_df <- stats_filter(dpmsr_set$data$stats$final, i)  
  
      cat(file=stderr(), str_c("Saving excel...", comp_string), "\n")
      addWorksheet(wb, comp_string)
      writeData(wb, sheet = nextsheet, filtered_df)
      nextsheet <- nextsheet +1
      z <- z+2
    }
    cat(file=stderr(), "writting excel to disk...", "\n")
    saveWorkbook(wb, filename, overwrite = TRUE)
    
}

# data table filter ------------------------------------------------------

stats_data_table_filter <- function(session, input, output) {
  
  comp_string <- input$stats_select_data_comp
  comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
  
  if(input$stats_include_all) {
      df_filter <- dpmsr_set$data$stats$final
      sample_number <- dpmsr_set$y$sample_number
    }else{
      df <- dpmsr_set$data$stats$final
      df_N <- df %>% dplyr::select(unlist(dpmsr_set$y$stats$groups$com_N_headers[comp_number]))  
      df_D <- df %>% dplyr::select(unlist(dpmsr_set$y$stats$groups$com_D_headers[comp_number]) ) 
      df_cv_N <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$cv_N[comp_number]) ) 
      df_cv_D <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$cv_D[comp_number]) ) 
      df_spqc <- df %>% dplyr::select(unlist(dpmsr_set$y$stats$comp_spqc_sample_headers ))  
      df_spqc_cv <- df %>% dplyr::select(contains(str_c(dpmsr_set$y$stats$comp_spqc, "_CV") )) 
      df_Comp <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$comp_name [comp_number]) ) 
      df_filter <- cbind(df[1:(dpmsr_set$y$info_columns_final)], df_N, df_D, df_spqc, df_cv_N, df_cv_D, df_spqc_cv, df_Comp)
      sample_number <- ncol(df_N) + ncol(df_D) + ncol(df_spqc)
    }
  
  if(input$stats_add_filters){
    df_filter <- stats_filter(df_filter, comp_number)
  }
  
  
  if(input$stats_data_topn != 0 ){
    df_filter$sum <- rowSums(df_filter[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+sample_number)])
    df_filter <- df_filter[order(-df_filter$sum),]                      
    df_filter <- df_filter[1:input$stats_data_topn,]
    df_filter$sum <- NULL
  }
  
  if(input$stats_data_accession != "0" ){
    df_filter <-subset(df_filter, Accession %in% as.character(input$stats_data_accession)  )
  }
  
  if(input$stats_data_description != "0") {
    df_filter <-df_filter[grep(as.character(input$stats_data_description), df_filter$Description), ]
  }
  
  return(df_filter)
}



# peptide zscore ------------------------------------------------------

peptide_zscore <- function(df_peptide, info_columns) {
  info_df <- df_peptide[1:info_columns]
  df <- df_peptide[(info_columns+1):ncol(df_peptide)]
  df<- log(df, 2)
  z_mean <- rowMeans(df)
  z_stdev  <- apply(df[1:ncol(df)], 1, sd)
  
  df <- apply(df, MARGIN = 2, function(x) (x-z_mean)/z_stdev  )
  df <- data.frame(cbind(info_df, df), stringsAsFactors = FALSE)
  
  return(df)
  
}


# one protein data ------------------------------------------------------

oneprotein_data <- function(session, input, output) {

    df <- dpmsr_set$data$stats$final
    df <-subset(df, Accession %in% as.character(input$stats_oneprotein_accession)  )
    df <- df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
    
    comp_string <- input$stats_oneprotein_plot_comp
    cat(file=stderr(), "Create stats plots" , "\n")
    
    for (i in 1:length(comp_string)){
      comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string[i])
      if (i==1){
        comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
      }else{
        comp_rows <- c(comp_rows, dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
      }
    }
    #add spqc to plots
    if(input$stats_oneprotein_plot_spqc){
      comp_rows <- c(comp_rows, dpmsr_set$y$stats$comp_spqc_sample_numbers)
    }
    
    comp_rows <- sort(unique(unlist(comp_rows)), decreasing = FALSE)
    df <- df[,comp_rows]
    namex <- dpmsr_set$design$Label[comp_rows]
    color_list <- dpmsr_set$design$colorlist[comp_rows]
    groupx <- dpmsr_set$design$Group[comp_rows]
    
    df_peptide <- dpmsr_set$data$impute[[dpmsr_set$data$stats$final_comp]]
    df_peptide <- subset(df_peptide, Accession %in% as.character(input$stats_oneprotein_accession)  )
    peptide_info_columns <- ncol(df_peptide) - dpmsr_set$y$sample_number
    colnames(df_peptide)[(peptide_info_columns+1):ncol(df_peptide)] <- dpmsr_set$design$Header1
    testcn <<- colnames(df_peptide)
    #df_peptide <- df_peptide[(peptide_info_columns+1):ncol(df_peptide)]
    df_peptide <- df_peptide[,(c(1:peptide_info_columns,(comp_rows+peptide_info_columns))) ]
    #sort peptides by intensity, keep highest abundant peptide
    df_peptide$sum <- rowSums(df_peptide[(peptide_info_columns+1):ncol(df_peptide)])
    df_peptide <- df_peptide[order(df_peptide$sum), ]
    df_peptide$sum <- NULL
    df_peptide <- df_peptide[!duplicated(df_peptide$Sequence),]
    
    if(input$stats_use_zscore){df_peptide <- peptide_zscore(df_peptide, peptide_info_columns)}


    return(list("df"=df, "df_peptide"=df_peptide, "namex"=namex, "color_list"=color_list, "comp_string"=comp_string, "peptide_info_columns"=peptide_info_columns))
}
    



# stat filter ------------------------------------------------------  
    stats_filter <- function(df, comp_number) { 
      
      filtered_df <- filter(df, df[[dpmsr_set$y$stats$groups$pval[comp_number]]] <= as.numeric(dpmsr_set$x$pvalue_cutoff) &  
                              (df[[dpmsr_set$y$stats$groups$fc[comp_number]]] >= as.numeric(dpmsr_set$x$foldchange_cutoff) | 
                                 df[[dpmsr_set$y$stats$groups$fc[comp_number]]] <= -as.numeric(dpmsr_set$x$foldchange_cutoff))) 
      
      filtered_df <- filter(filtered_df, filtered_df[[dpmsr_set$y$stats$groups$mf[comp_number] ]] >= as.numeric(dpmsr_set$x$missing_factor))
      
      
      if(dpmsr_set$y$stats$stats_spqc_cv_filter){
        filtered_df <- filter(filtered_df, filtered_df[[str_c(dpmsr_set$y$stats$comp_spqc, "_CV")]] <= as.numeric(dpmsr_set$y$stats$stats_spqc_cv_filter_factor ) )
      }
      
      if(dpmsr_set$y$stats$stats_comp_cv_filter){
        filtered_df <- filter(filtered_df, filtered_df[[dpmsr_set$y$stats$groups$cv_N[comp_number]]] <= as.numeric(dpmsr_set$y$stats$stats_comp_cv_filter_factor) | 
                                filtered_df[[dpmsr_set$y$stats$groups$cv_D[comp_number]]] <= as.numeric(dpmsr_set$y$stats$stats_comp_cv_filter_factor) )
      }
     
      return(filtered_df)
       
    }      
    
    
  peptide_position_lookup <- function(session, input, output)  {
    # create peptide lookup table
    peptide_pos_lookup <- dpmsr_set$data$data_raw_peptide %>% dplyr::select(Sequence, Positions.in.Master.Proteins)
    colnames(peptide_pos_lookup) <- c("Sequence", "Position")
    peptide_pos_lookup$Position <- gsub("\\]; \\[", "xxx",  peptide_pos_lookup$Position)
    s <- strsplit(peptide_pos_lookup$Position, split = "; ")
    peptide_pos_lookup <- data.frame(Sequence = rep(peptide_pos_lookup$Sequence, sapply(s, length)), Position = unlist(s))
    peptide_pos_lookup$Count <- str_count(peptide_pos_lookup$Position, "xxx")
    peptide_pos_lookup <- peptide_pos_lookup %>% separate(Position, c("Accession", "Start_Stop"), sep = " ")
    peptide_pos_lookup$Start_Stop <- gsub("\\[", "",  peptide_pos_lookup$Start_Stop)
    peptide_pos_lookup$Start_Stop <- gsub("\\]", "",  peptide_pos_lookup$Start_Stop)
    peptide_pos_lookup$Start_Stop <- gsub("xxx", ", ",  peptide_pos_lookup$Start_Stop)
    peptide_pos_lookup$AS <- str_c(peptide_pos_lookup$Accession, " ", peptide_pos_lookup$Sequence)
    s <- strsplit(peptide_pos_lookup$Start_Stop, split = ", ")
    peptide_pos_lookup <- data.frame(AS = rep(peptide_pos_lookup$AS, sapply(s, length)), Position = unlist(s))
    peptide_pos_lookup <- peptide_pos_lookup %>% separate(AS, c("Accession", "Sequence"), sep = " ")
    peptide_pos_lookup <- peptide_pos_lookup %>% separate(Position, c("Start", "Stop"), sep = "-")
    peptide_pos_lookup$dup <- str_c(peptide_pos_lookup$Accession, peptide_pos_lookup$Sequence, peptide_pos_lookup$Start, peptide_pos_lookup$Stop)
    peptide_pos_lookup <- peptide_pos_lookup[!duplicated(peptide_pos_lookup$dup),]
    peptide_pos_lookup$dup <- NULL
    
    peptide_pos_lookup <- subset(peptide_pos_lookup, Accession %in% as.character(input$stats_oneprotein_accession)  )
    #peptide_pos_lookup <- subset(peptide_pos_lookup, Accession %in% "P60202"  )  
  
    return(peptide_pos_lookup)
    
    }