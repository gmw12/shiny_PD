

#----------------------------------------------------------------------------------------
# create final excel documents
Simple_Excel <- function(df, filename) {
  require(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, deparse(substitute(df)))
  writeData(wb, sheet =1, df)  
  saveWorkbook(wb, filename, overwrite = TRUE)
}
#----------------------------------------------------------------------------------------
# create final excel documents
Simple_Excel_name <- function(df, filename, sheet_name) {
  require(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet =1, df)  
  saveWorkbook(wb, filename, overwrite = TRUE)
}
#----------------------------------------------------------------------------------------
# create final excel documents
Final_Excel <- function() {
  require(openxlsx)
  
  #colnames(dpmsr_set$data$finalraw)[(dpmsr_set$y$info_columns_final+1):ncol(dpmsr_set$data$finalraw)] <- dpmsr_set$design$Header1
  colnames(dpmsr_set$data$finalraw[(dpmsr_set$y$info_columns_final+1):ncol(dpmsr_set$data$finalraw)]) <- dpmsr_set$design$Header1
  
  for(df_name in names(dpmsr_set$data$final)){
    filename <- str_c(dpmsr_set$file$output_dir, df_name, "//", dpmsr_set$x$file_prefix, "_", df_name, "_final.xlsx")
    df <- dpmsr_set$data$final[[df_name]]
    #remove FC2 from df for excel
    df <- df[,-grep(pattern="_FC2$", colnames(df))]
   
    if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
      df_raw <- dpmsr_set$data$finalraw
    }else{
        df_raw <- dpmsr_set$data$finalraw
         #if PTM out need to reduce raw data frame for excel
        if (as.logical(dpmsr_set$x$peptide_report_ptm)){
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
    for(i in 1:as.numeric(dpmsr_set$x$comp_number))  {
      comp_string <- dpmsr_set$y$comp_groups$comp_name[i]
      #assign(comp_string, subset(df, get(comp_pval_groups[i])<=pvalue_cutoff & (get(comp_fc_groups[i])>=fc_cutoff | get(comp_fc_groups[i])<= -fc_cutoff)) )
      filtered_df <- subset(df, df[ , dpmsr_set$y$comp_groups$pval[i]] <= as.numeric(dpmsr_set$x$pvalue_cutoff) &  
                              (df[ , dpmsr_set$y$comp_groups$fc[i]] >= as.numeric(dpmsr_set$x$foldchange_cutoff) | 
                                 df[ , dpmsr_set$y$comp_groups$fc[i]] <= -as.numeric(dpmsr_set$x$foldchange_cutoff))) 
      addWorksheet(wb, comp_string)
      writeData(wb, sheet = nextsheet, filtered_df)
      nextsheet <- nextsheet +1
      z <- z+2
    }
    saveWorkbook(wb, filename, overwrite = TRUE)
    
  }
}


#----------------------------------------------------------------------------------------
save_design <- function(session, input){
  design_data <- input$design_file
  design_name <- unlist(design_data$files[[as.character(0)]][2])
  file.copy(design_name, str_c(dpmsr_set$file$backup_dir, basename(design_name)))
}

#----------------------------------------------------------------------------------------
save_data <- function(data_file){
  file.copy(data_file, str_c(dpmsr_set$file$backup_dir, basename(data_file)))
}

#----------------------------------------------------------------------------------------
create_dir <- function(name){
  if(is_dir(name)) {
    #added file delete, dir delete not working on customer shiny server
    do.call(file.remove, list(list.files(name, full.names = TRUE)))
    dir_delete(name)
    dir_create(name)
    }else{
    dir_create(name)
    }
  name <- str_replace_all(name, "/", "//")
  name <- str_c(name, "//")
  return(name)
}

create_dir_old <- function(name){
  if(is_dir(name)) { 
    unlink(file.path(".", name), recursive = TRUE, force=TRUE)} 
    ifelse (!dir.exists(file.path(".", name)), dir_create(name))
  name <- str_replace_all(name, "/", "//")
  name <- str_c(".", name, "//")
  return(name)
}

# name<-dpmsr_set$file$data_dir
# dir.exists(file.path(name))
# dir.create(file.path(".", name))
# name<-substring(name,2)
# name <- str_replace_all(name, "//", "/")
# name <- dpmsr_set$x$file_prefix
# dir_create(test)
# is_dir(test)
# dir_create("testdir/testdir")
# dir_delete("testdir")
# 
# 
# np <- dpmsr_set$x$new_path
# np<-substring(np,2)
# np<-str_c(np,"/test")
# dir_create(dpmsr_set$file$data_dir)
#----------------------------------------------------------------------------------------
file_set <- function(){
  dpmsr_set$file$data_path <<- dpmsr_set$x$new_path
  dpmsr_set$file$data_dir <<- str_c(dpmsr_set$file$data_path, dpmsr_set$x$file_prefix)
  dpmsr_set$file$output_dir <<- create_dir(dpmsr_set$file$data_dir)
  dpmsr_set$file$backup_dir <<- create_dir(str_c(dpmsr_set$file$data_dir,"/Backup"))
  dpmsr_set$file$extra_dir <<- create_dir(str_c(dpmsr_set$file$data_dir,"/Extra"))
  dpmsr_set$file$qc_dir <<- create_dir(str_c(dpmsr_set$file$data_dir,"/QC"))
  dpmsr_set$file$extra_prefix <<- str_c(dpmsr_set$file$extra_dir,dpmsr_set$x$file_prefix) 
  dpmsr_set$file$string <<- create_dir(str_c(dpmsr_set$file$data_dir,"/String"))
  dpmsr_set$file$phos <<- create_dir(str_c(dpmsr_set$file$data_dir,"/Phos"))
  dpmsr_set$file$app_dir <<- create_dir(str_c(dpmsr_set$file$data_dir,"/Backup/App"))
  
  r_files <- list.files()
  for (i in 1:length(r_files)){
    if(grepl(".R$", r_files[i])){
      file.copy(r_files[i], str_c(dpmsr_set$file$app_dir, r_files[i]))
    }
  }
  

}







