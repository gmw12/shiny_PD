#----------------------------------------------------------------------------------------
# clear shiny files and memory
clear_memory <- function(input, output, session) {
  cat(file=stderr(), "clean old data", "\n")
  rm(list = ls(envir = .GlobalEnv), pos = .GlobalEnv, inherits = FALSE)
  cat(file=stderr(), "reload libraries and functions", "\n")
  rm(list = ls())
  gc()
  #source("Shiny_Functions_v1.R")
  #set_user(session, input, output)
  #app_startup(session, input, output)
}

#----------------------------------------------------------------------------------------
# create final excel documents

Simple_Excel <- function(df, sheetname, filename) {
  cat(file=stderr(), str_c("Simple_Excel -> ", filename), "\n")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetname)
  openxlsx::writeData(wb, sheet=1, df)  
  openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)
}

Simple_Excel_bg <- function(df, sheetname, filename) {
  cat(file=stderr(), str_c("Simple_Excel_bg -> ", filename), "\n")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetname)
  openxlsx::writeData(wb, sheet=1, df)  
  test <- callr::r_bg(openxlsx::saveWorkbook, args = list(wb, filename, overwrite = TRUE), supervise = TRUE)
}


#----------------------------------------------------------------------------------------
# create final excel documents
Simple_read_Excel <- function(path, sheet_name) {
  data <- readxl::read_excel(path, sheet_name)
  return(data)
}

#----------------------------------------------
# read excel in new backgroup process
Simple_read_Excel <- function(path, sheet_name, col_names, skip) {
  readxl::read_excel(path=path, sheet=sheet_name, col_names = col_names, skip=skip)
} 

bg_read_Excel <- function(path, sheet_name, col_names, skip) {
  output <- callr::r_bg(func=Simple_read_Excel, args=list(path=path, sheet=sheet_name, col_names = col_names, skip=skip), supervise = TRUE)
  output$wait()
  return(output$get_result())
}

bg_read_Excel2 <- function(count) {
  cat(file=stderr(), str_c("printing this"), "\n")
  print("printing this")
  temp_list <- file_args[[count]]
  output <- callr::r_bg(func=Simple_read_Excel, args=list(path=temp_list$path, sheet=temp_list$sheet, col_names = temp_list$col_names, skip=temp_list$skip), supervise = TRUE)
  output$wait()
  return(output$get_result())
}

parallel_read_Excel <- function(file_count){
  ncores <- (detectCores()/2)
  cat(file=stderr(), str_c("ncores = ", ncores), "\n")
  if (is.na(ncores) | ncores < 1) {ncores <- 1}
  erase2 <<- mclapply(file_count, bg_read_Excel2, mc.cores = ncores)
}

Simple_fread <- function(file) {
  data.table::fread(file=file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
} 

#----------------------------------------------------------------------------  
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
save_data_bg <- function(file1, dir1){
  file2 <- paste(dir1, basename(file1))
  file.copy(file1, file2)
}

#----------------------------------------------------------------------------------------
create_dir <- function(name){
  cat(file=stderr(), "create_dir...", "\n")
  if(is_dir(name)) {
    #added file delete, dir delete not working on customer shiny server
    cat(file=stderr(), "dir exists, deleting...", "\n")
    do.call(file.remove, list(list.files(name, full.names = TRUE)))
    dir_delete(name)
    dir_create(name)
    }else{
    dir_create(name)
    }
  name <- str_replace_all(name, "/", "//")
  name <- str_c(name, "//")
  cat(file=stderr(), str_c(name, " created...", "\n"))
  return(name)
}

#----------------------------------------------------------------------------------------
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


#----------------------------------------------------------------------------------------
# create final excel documents
save_object <- function(df, file_cust_name) {
  dpmsr_set <- df
  save(dpmsr_set, file=file_cust_name)
}


#----------------------------------------------------------------------------------------
# get file extention
getExtension <- function(file){ 
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
} 


#----------------------------------------------------------------------------------------
# log app usuage
app_log <- function() {
  cat(file=stderr(), "app_log...", "\n")
  log_file <- "/data/ShinyData/app_log/app_log.csv"
  
  df <- data.frame(Sys.time(), dpmsr_set$x$file_prefix, dpmsr_set$x$raw_data_input, dpmsr_set$x$final_data_output, dpmsr_set$x$peptide_report_ptm,  dpmsr_set$x$tmt_spqc_norm)
  colnames(df) <- c("Date", "File", "Input", "Output", "PTM", "TMT_norm")
  
  if(!is_dir("/data/ShinyData/app_log/")) {
    dir_create("/data/ShinyData/app_log/")
  }
  
  if(file.exists(log_file)){
    temp_df <- data.table::fread(file=log_file, header = TRUE, stringsAsFactors = FALSE, sep = ";")
    df <- rbind(temp_df, df)
  }

  write.csv2(df, log_file, row.names = FALSE)
  cat(file=stderr(), "app_log...finished", "\n")
}