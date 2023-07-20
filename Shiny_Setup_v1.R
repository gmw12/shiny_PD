set_user <- function(session, input, output){
  
  #set user to unkown to force app to find correct usr
  site_user <<- "unknown"
  volumes <<- "unknown"
  
  cat(file=stderr(), str_c("site_user default --> ", site_user), "\n")
  cat(file=stderr(), str_c("volumes --> ", volumes), "\n")
  
  cat(file=stderr(), "setting user...", "\n")
  while (site_user == "unknown"){
    if(Sys.info()["sysname"]=="Darwin" ){
      volumes <<- c(dd='/Users/gregwaitt/Documents/Data', wd='.', Home = fs::path_home(),  getVolumes()())
      #version determines website content
      site_user <<- "dpmsr"
      #volumes <<- c(wd='.', Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
    }else if (Sys.info()["nodename"] == "titanshinyu20"){
      #for titan_black VM
      volumes <<- c(dd='/home/dpmsr/shared/h_drive', dd2='/home/dpmsr/shared/other_black', RawData='/home/dpmsr/shared/RawData', wd='.', Home = fs::path_home(), getVolumes()())
      site_user <<- "dpmsr"
    }else if (Sys.info()["nodename"] == "greg-GS63VR-7RF"){
      #for greg linux laptop
      volumes <<- c(dd='/home/dpmsr/shared', wd='.', Home = fs::path_home(), getVolumes()())
      site_user <<- "dpmsr"
    }else if (Sys.info()["nodename"] == "greg-ThinkPad-W550s"){
      #for greg linux laptop
      volumes <<- c(dd='/home/dpmsr/shared', wd='.', Home = fs::path_home(), getVolumes()())
      site_user <<- "dpmsr"
    }else if (Sys.info()["nodename"] == "mascot"){
      #for mascot search pc
      volumes <<- c(h1='/mnt/h_black1', h2='/mnt/h_black2', dc='/mnt/DataCommons', wd='.', Home = fs::path_home(), getVolumes()())
      site_user <<- "dpmsr"
    }else{
      #for public website
      volumes <<- c(dd='/data', wd='.', Home = fs::path_home(), getVolumes()())
      site_user <<- "not_dpmsr"
    }
  }
  cat(file=stderr(), "setting user...end", "\n")
  
  cat(file=stderr(), str_c("site_user set to -->  ", site_user), "\n")
  cat(file=stderr(), str_c("volumes --> ", volumes), "\n")
  
  return(site_user)
}


#---------------------------------------------------------------------
load_dpmsr_set <- function(session, input, output){
  cat(file=stderr(), "creating dpmsr_set...", "\n")
  design_list <- c("General", "QC", "Fill_Norm", "Filters", "Protein", "TMT_PTM", "Report")
  design_data <- parseFilePaths(volumes, input$design_file)
  new_path <- str_extract(design_data$datapath, "^/.*/")
  cat(file=stderr(), str_c("loading design file from ", new_path), "\n")
  
  #design <- read_excel(design_data$datapath, sheet="SampleList")
  #design <- bg_read_Excel(path=design_data$datapath, sheet="SampleList", col_names = TRUE, skip=0)
  #protocol <- bg_read_Excel(path=design_data$datapath, sheet="Protocol", col_names = TRUE, skip = 0)
  
  #practicing using parallel processing
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  data <- foreach (sheet=c("SampleList", "Protocol")) %dopar% {
    bg_read_Excel(path=design_data$datapath, sheet=sheet, col_names = TRUE, skip=0)
  }
  design <- data[[1]]
  protocol <- data[[2]]
  
  dpmsr_set <<- list(design = design)
  dpmsr_set$protocol <<- protocol
  temp_list<-list()
  
  for(i in design_list){
    design_tab<-bg_read_Excel(design_data$datapath, sheet=i, col_names = FALSE, skip = 1)
    for(j in 1:nrow(design_tab)){
      lname <- design_tab[j,1]
      ldata <- design_tab[j,2]
      temp_list[[as.character(lname)]] <- ldata
    }
  }
  #set default values
  dpmsr_set$x <<- temp_list
  dpmsr_set$x$new_path <<- new_path
  dpmsr_set$x$volumes <<- volumes
  dpmsr_set$x$design_name <<- design_data$datapath
  dpmsr_set$x$default_root <<- 'wd'
  dpmsr_set$x$default_path <<- ''
  
  #set default volume
  if(grepl("h_drive", dpmsr_set$x$new_path)){
    dpmsr_set$x$default_root <<- 'dd'
    dpmsr_set$x$default_path <<- substring(dpmsr_set$x$new_path, nchar(volumes[[dpmsr_set$x$default_root]])+1)
  }else if (grepl("other_black", dpmsr_set$x$new_path)){
    dpmsr_set$x$default_root <<- 'dd2'
    dpmsr_set$x$default_path <<- substring(dpmsr_set$x$new_path, nchar(volumes[[dpmsr_set$x$default_root]])+1)
  }
  
  substring(dpmsr_set$x$new_path, nchar(volumes[[dpmsr_set$x$default_root]])+1)
  
  #set default values for parameters not in sample list file
  dpmsr_set$x$rollup_method <<- "Sum"
  dpmsr_set$y$rollup_topN_count <<- 3
  dpmsr_set$x$checkbox_norm_include <<- FALSE
  dpmsr_set$x$checkbox_norm_exclude <<- FALSE
  dpmsr_set$x$directlfq <<- FALSE
  dpmsr_set$x$data_source <<- 1
  dpmsr_set$x$pathway_set <<- 0 
  dpmsr_set$x$motif_set <<- 0
  dpmsr_set$app_version <<- app_version
  if (input$primary_group){
    dpmsr_set$x$primary_group <<- TRUE
  }else{
    dpmsr_set$x$primary_group <<- FALSE
  }
  
  cat(file=stderr(), "dpmsr_set created...", "\n")

  #check design file for errors
  if(any(duplicated(dpmsr_set$design$ID))) {
    shinyalert("Oops!", "Duplicated ID's in Sample List File", type = "error")
  }
  cat(file=stderr(), "dpmsr_set checked for duplicate ID's...", "\n")

  
  #check if comparisons too long
  if(check_comp_name_length() ) {
    shinyalert("Oops!", "Comparison names too long for excel sheet names", type = "error")
  } 
  cat(file=stderr(), "dpmsr_set checked for comparison names too long", "\n")
}


#------------------------

check_comp_name_length <- function(){
  cat(file=stderr(), "check_comp_name_length...", "\n")
  names <- unique(dpmsr_set$design$Group)
  names_len <- sort(nchar(names), decreasing = TRUE)
  if (length(names) > 1){
    longest_comp <- names_len[1] + names_len[2]
  }else{
    longest_comp <- names_len[1] + names_len[1]
  }
  if (longest_comp > 27) {
    return(TRUE)
  }else {
      return(FALSE)
    }
}


#------------------------
#create initial dpmsr_set S3 data class
dpmsr_set_create <- function(design){
  dpmsr_set <<- list(design=design)
  class(dpmsr_set) <<- "dpmsr_set"
}

#----------------------------------------------------------------------------------------
#load design elements into dpmsr_set
load_design <- function(session, input){
  cat(file=stderr(), "load_design...", "\n")
  design_list <- c("General", "QC", "Fill_Norm", "Filters", "Protein", "TMT_PTM", "Report")
  for(x in design_list){
    for(i in 1:nrow(dpmsr_set$design[[x]]) ){
      dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]] <<- dpmsr_set$design[[x]][i,2]
      if(dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]] == "TRUE" || dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]] == "FALSE" ){
        dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]]<<-as.logical(dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]])
      }else if(!grepl("\\D", dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]] )){
        dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]]<<-as.numeric(dpmsr_set$x[[dpmsr_set$design[[x]][i,1]]])
      }
    }
  }

  #dpmsr_set$x$test <<- "test"
  #dpmsr_set$x$peptide_norm_grep <<- str_replace_all(dpmsr_set$x$peptide_norm_grep, "/", "\\\\")
  #dpmsr_set$x$peptide_impute_grep <<- str_replace_all(dpmsr_set$x$peptide_impute_grep, "/", "\\\\")
  #dpmsr_set$x$peptide_report_grep <<- str_replace_all(dpmsr_set$x$peptide_report_grep, "/", "\\\\")
  #dpmsr_set$x$int_cutoff <- as.numeric(dpmsr_set$x$int_cutoff)
}

#----------------------------------------------------------------------------------------
load_data <- function(session, input, volumes){
  start_time <- Sys.time()
  cat(file=stderr(), "load_data (Proteome Discoverer)...", "\n")
  raw_data <- parseFilePaths(dpmsr_set$x$volumes, input$raw_files)
  
  #Loop through and load data in background
  for (i in 1:nrow(raw_data) ){
    raw_name <- raw_data$datapath[i]
    if (grepl("_PeptideGroups.txt", raw_name)){
      cat(file=stderr(), "loading raw peptide data...", "\n")
      temp_df1 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("_Proteins.txt", raw_name)){
      cat(file=stderr(), "loading raw protein data...", "\n")
      temp_df2 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("_PSMs.txt", raw_name)){
      cat(file=stderr(), "loading raw psm data...", "\n")
      temp_df3 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("MSMSSpectrumInfo.txt", raw_name)){
      cat(file=stderr(), "loading raw msms data...", "\n")
      temp_df4 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("InputFiles.txt", raw_name)){
      cat(file=stderr(), "loading raw inputfile data...", "\n")
      temp_df5 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("_DecoyPeptideGroups.txt", raw_name)){
      cat(file=stderr(), "loading raw decoy peptide data...", "\n")
      temp_df6 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("_DecoyProteins.txt", raw_name)){
      cat(file=stderr(), "loading raw decoy protein data...", "\n")
      temp_df7 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("_DecoyPSMs.txt", raw_name)){
      cat(file=stderr(), "loading raw decoy psm data...", "\n")
      temp_df8 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    } else if (grepl("_PeptideIsoforms.txt", raw_name)){
      cat(file=stderr(), "loading raw peptide isoform data...", "\n")
      temp_df9 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    }else if (grepl("_LCMSFeatures.txt", raw_name)){
      cat(file=stderr(), "loading feature data...", "\n")
      temp_df10 <- callr::r_bg(func=Simple_fread, args=list(file=raw_name), supervise = TRUE)
    }
  }  
  #loop back through and assign dataframe from backgroup process  
  for (i in 1:nrow(raw_data) ){
    raw_name <- raw_data$datapath[i]
    if (grepl("_PeptideGroups.txt", raw_name)){
      cat(file=stderr(), "assign raw peptide data...", "\n")
      temp_df1$wait()
      dpmsr_set$data$data_raw_peptide <<- temp_df1$get_result() %>% dplyr::select(contains(
        c("Confidence", "Accession", "Description", "Sequence", "Modifications", "Positions", "Abundance.F", "Retention.Time", "Ion.Score", "Percolator.SVM",
          "q.Value", "RT.in.min", "mz.in.Da.by.Search.Engine.", "Charge.by.Search.Engine.", "Quan.Info")
      ))
      dpmsr_set$data$data_raw_peptide <<- as.data.frame(dpmsr_set$data$data_raw_peptide)
    } else if (grepl("_Proteins.txt", raw_name)){
      cat(file=stderr(), "assign raw protein data...", "\n")
      temp_df2$wait()
      dpmsr_set$data$data_raw_protein <<- temp_df2$get_result() %>% dplyr::select(contains(
        c("Master", "Protein.FDR.Confidence.Combined", "Number.of.Razor.Peptides", "Accession", "Description",
          "Number.of.Peptides", "Coverage.in.Percent", "Number.of.Unique.Peptides", "Number.of.Razor.Peptides",
          "Number.of.Protein.Unique.Peptides", "Abundance.F")
      ))
      dpmsr_set$data$data_raw_protein <<- as.data.frame(dpmsr_set$data$data_raw_protein)
    } else if (grepl("_PSMs.txt", raw_name)){
      cat(file=stderr(), "assign raw psm data...", "\n")
      temp_df3$wait()
      dpmsr_set$data$data_raw_psm <<- as.data.frame(temp_df3$get_result())
    } else if (grepl("MSMSSpectrumInfo.txt", raw_name)){
      cat(file=stderr(), "assign raw msms data...", "\n")
      temp_df4$wait()
      dpmsr_set$data$data_raw_msms <<- as.data.frame(temp_df4$get_result())
    } else if (grepl("InputFiles.txt", raw_name)){
      cat(file=stderr(), "assign raw inputfile data...", "\n")
      temp_df5$wait()
      dpmsr_set$data$data_raw_inputfiles <<- as.data.frame(temp_df5$get_result())
    } else if (grepl("_DecoyPeptideGroups.txt", raw_name)){
      cat(file=stderr(), "assign raw decoy peptide data...", "\n")
      temp_df6$wait()
      dpmsr_set$data$data_raw_decoypeptide <<- as.data.frame(temp_df6$get_result())
    } else if (grepl("_DecoyProteins.txt", raw_name)){
      cat(file=stderr(), "assign raw decoy protein data...", "\n")
      temp_df7$wait()
      dpmsr_set$data$data_raw_decoyprotein <<- as.data.frame(temp_df7$get_result())
    } else if (grepl("_DecoyPSMs.txt", raw_name)){
      cat(file=stderr(), "assign raw decoy psm data...", "\n")
      temp_df8$wait()
      dpmsr_set$data$data_raw_decoypsm <<- as.data.frame(temp_df8$get_result())
    } else if (grepl("_PeptideIsoforms.txt", raw_name)){
      cat(file=stderr(), "assign raw peptide isoform data...", "\n")
      temp_df9$wait()
      dpmsr_set$data$data_raw_isoform <<- temp_df9$get_result() %>% dplyr::select(contains(
        c("Confidence", "Accession", "Description", "Sequence", "Modifications", "Positions", "Abundance.F", "Retention.Time", "Ion.Score", "Percolator.SVM",
          "q.Value", "RT.in.min", "mz.in.Da.by.Search.Engine.", "Charge.by.Search.Engine.", "Quan.Info")
      ))
      dpmsr_set$data$data_raw_isoform <<- as.data.frame(dpmsr_set$data$data_raw_isoform)
    }else if (grepl("_LCMSFeatures.txt", raw_name)){
      cat(file=stderr(), "assign feature data...", "\n")
      temp_df10$wait()
      dpmsr_set$data$data_features <<- as.data.frame(temp_df10$get_result())
    }
  }
  
  #backup text files
  cat(file=stderr(), "backing up text files...", "\n")
  for (i in 1:nrow(raw_data) ){
    data_file = raw_data$datapath[i]
    testme <- callr::r_bg(func=save_data_bg, args=list(file1=data_file, dir1=dpmsr_set$file$backup_dir), supervise=TRUE)
  }
  save_data(dpmsr_set$x$design_name)
  
  #reassign to data.frame as data.table has subtle differences
  dpmsr_set$data$data_raw_peptide <<- data.frame(dpmsr_set$data$data_raw_peptide)
  dpmsr_set$data$data_raw_protein <<- data.frame(dpmsr_set$data$data_raw_protein)
  
  gc()
  cat(file=stderr(), str_c("Data load time ---> ", Sys.time()-start_time), "\n")
}




#----------------------------------------------------------------------------------------
load_data_sp <- function(session, input, volumes){
  cat(file=stderr(), "load_data (Spectronaut)...", "\n")
  raw_data <- parseFilePaths(dpmsr_set$x$volumes, input$raw_files)
  cat(file=stderr(), str_c("file path = ", raw_data), "\n")
  
  #design_data <- parseFilePaths(volumes, input$design_file)
  new_path <- str_extract(raw_data$datapath, "^/.*/")
  volumes["wd"] <- new_path
  dpmsr_set$data$data_raw_protein <<- read_excel(raw_data$datapath)
  save_data(raw_data$datapath)
  save_data(dpmsr_set$x$design_name)
  gc()
  cat(file=stderr(), "Spectronaut data loaded...", "\n")
}

#----------------------------------------------------------------------------------------
prepare_data <- function(session, input) {  #function(data_type, data_file_path){
  cat(file=stderr(), "prepare_data...", "\n")
   data_type <- input$radio_input
   if(as.numeric(data_type) == 1){
    cat(file=stderr(), "prepare data_type 1", "\n")
    dpmsr_set$data$data_peptide_start <<- protein_to_peptide()
    dpmsr_set$data$data_protein_start <<- protein_to_protein()
    dpmsr_set$y$state <<- "Peptide"
  }else if(data_type ==2){
    cat(file=stderr(), "prepare data_type 2", "\n")
    dpmsr_set$data$data_protein_start <<- protein_to_protein()
    dpmsr_set$y$state <<- "Protein"
  }else if(data_type ==3){
    cat(file=stderr(), "prepare data_type 3", "\n")
    dpmsr_set$data$data_peptide_start <<- peptide_to_peptide()
    dpmsr_set$y$state <<- "Peptide"
  }else{msgBox <- tkmessageBox(title = "Whoops",
                           message = "Invalid input output in design file", icon = "info", type = "ok")}

   if(input$checkbox_isoform){
     dpmsr_set$data$data_peptide_isoform_start <<- isoform_to_isoform()
   }
   gc()
}


#----------------------------------------------------------------------------------------
set_sample_groups <- function(session, input, output){
  cat(file=stderr(), "set_sample_groups...", "\n")
  
  sample_info <- dpmsr_set$design
  #sample_info <- dpmsr_set$design[,1:7]
  
  #check if using primary group for filter and impute
  if (input$primary_group){
    group_type <- "PrimaryGroup"
    sample_info$PrimaryGroup <- sapply(sample_info$Group, function(string) str_split(string, "_")[[1]][1])
  }else{
    group_type <- "Group"
  }
  
  #check to see if designfile/samplelist is sorted, if not then sort it
  sample_info <- check_design_sort(group_type, sample_info)
  
  # number of groups, comparisons are not hardcoded, R will extract information from excel SampleData file
  sample_number <- nrow(sample_info)

  #count number of samples in each group - add column
  sample_info$Count <- sapply(sample_info[[group_type]], function(string) sum(string==sample_info[[group_type]]))
  
  # create unqiue group dataframe with sample number
  sample_groups <- sample_info[7:ncol(sample_info)]
  sample_groups <- sample_groups %>% distinct(sample_groups[[group_type]], .keep_all = TRUE)
  #replace Group names with itself or the PrimaryGroup Name
  sample_groups$Group <- sample_groups[[group_type]]

  #count number of groups
  group_number <- nrow(sample_groups)
  
  #assign start and end columns for each group in pd output
  sample_groups$start <- 1
  sample_groups$end <- sample_groups$Count[1]
  for(i in 2:(group_number)) {
    sample_groups$start[i] <- sample_groups$start[i-1] + sample_groups$Count[i-1]
    sample_groups$end[i] <- sample_groups$start[i] + sample_groups$Count[i] - 1
  }

  #assign colors to groups
  color_choices <-distinctColorPalette(group_number)
  group_color <- color_choices[1:group_number]
  sample_groups$colorlist <- color_choices[1:group_number]
  sample_info$colorlist <- with(sample_groups, colorlist[match(sample_info[[group_type]], sample_groups[[group_type]])])
  sample_groups$title <- sample_groups[[group_type]]
    
  #organize column headers for final output
  sample_info$Header1 <- str_c(sample_info$ID, " ", sample_info$Group)
  sample_info$Header2 <- str_c(sample_info$ID, " ", sample_info$Group, " Normalized")
  sample_info$Header3 <- str_c(sample_info$ID, " ", sample_info$Group, " Imputed")


  #erase if no crash
  #group_title <- paste(sample_groups$title, collapse=" ")
  #group_title <- gsub("\\)", "\\),",  group_title)
  

  #create factor vector
  group_factor <- rep(1, sample_groups$Count[1])
  for(i in 2:nrow(sample_groups)) {
    group_factor <- c(group_factor, rep(i, sample_groups$Count[i]))
  }
  group_factor <-factor(group_factor)
  
  #save to dpmsr_set
  dpmsr_set$y$sample_number <<- sample_number
  dpmsr_set$y$sample_groups <<- sample_groups
  dpmsr_set$y$group_number <<- group_number
  dpmsr_set$y$group_factor <<- group_factor
  dpmsr_set$design <<- sample_info
  
  #create dataframe to hold cv summaries for normalization strategies
  dpmsr_set$data$summary_cv <<- data.frame(sample_groups[[group_type]])
  
  #delete later 5/4/23
  #create list of all sample variables for Multivariate Analysis - in order from left to right
  #groups <- list()
  #for(i in 1:length(sample_info$Group)){
  #  groups <-  c(groups, unlist(str_split(sample_info$Group[i], "_")))
  #}

  #dpmsr_set$y$uniquegroups <- c(unique(groups), "Min.Samples", "All.Samples")
  #dpmsr_set$y$uniquegroups <<- dpmsr_set$y$uniquegroups[!is.na(dpmsr_set$y$uniquegroups)]

  #save unique group names delete this too
  dpmsr_set$y$uniquegroups <<- unique(sample_info$Group)
  
}



#----------------------------------------------------------------------------------------
isoform_set <- function(){
  cat(file=stderr(), "isoform_set...", "\n")
  norm_dir <- create_dir(str_c(dpmsr_set$file$data_dir,"//Norm"))
  norm_prefix <- str_c(norm_dir, file_prefix)
  if(ptm_peptide_only){norm_data <- special_norm(data_raw) }
  norm_data <- order_columns(norm_data)
  norm_data <- filter_data(norm_data)
  Simple_Excel(norm_data, "Norm_Data", str_c(norm_prefix, "_NormData.xlsx"))
  norm_data <- norm_data[(ncol(data_raw)-sample_number+1):ncol(data_raw)]
  norm_data[is.na(norm_data)] <- 0.0
  data_raw <- peptide_isoform
}
#----------------------------------------------------------------------------------------

set_input <- function(input_type){
  if(input_type ==1) {dpmsr_set$x$raw_data_input <<- "Protein_Peptide"}
    else if(input_type ==2) {dpmsr_set$x$raw_data_input <<- "Protein"}
    else if(input_type ==3 || input_type ==4) {dpmsr_set$x$raw_data_input <<- "Peptide"}
    else if(input_type ==5 || input_type ==6) {dpmsr_set$x$raw_data_input <<- "PSM"}
}
#----------------------------------------------------------------------------------------

set_output <- function(output_type){
  if(output_type ==1) {dpmsr_set$x$final_data_output <<- "Protein"}
    else if(output_type ==2) {dpmsr_set$x$final_data_output <<- "Peptide"} 
}

#--------------------------------------------------------------------------------------
set_tmt_sqpc_norm <- function(spqc_tmt){
  if(spqc_tmt ==0) {dpmsr_set$x$tmt_spqc_norm <<- FALSE}
    else if(spqc_tmt == 1) {dpmsr_set$x$tmt_spqc_norm <<- TRUE}  
}

#--------------------------------------------------------------------------------------
set_out_ptm <- function(ptm_out){
  if(ptm_out ==0) {dpmsr_set$x$peptide_ptm_out <<- FALSE}
    else if(ptm_out== 1) {dpmsr_set$x$peptide_ptm_out <<- TRUE}  
}


#----------------------------------------------------------------------------------------
#check if sample list is sorted (groups are together)

check_design_sort <- function(group_type, sample_info){
  sample_info$check <- 0
  
  if (sample_info[[group_type]][1] == sample_info[[group_type]][2]) {
    sample_info$check[1] <- 1
  }
  for (i in (2:nrow(sample_info))){
    if (sample_info[[group_type]][i] == sample_info[[group_type]][i-1] || sample_info[[group_type]][i] == sample_info[[group_type]][i+1]){
      sample_info$check[i] <- 1
    } 
  }
  
  if(sum(sample_info$check) >= nrow(sample_info)){
    cat(file=stderr(), "The sample design file is properly sorted", "\n")
  }else{
    cat(file=stderr(), "The sample design file is NOT properly sorted", "\n")
    # sort dataframe, if SPQC exists put it at bottom
    sample_info_spqc <- sample_info[sample_info[[group_type]]=='SPQC',]
    sample_info <- sample_info[sample_info[[group_type]]!='SPQC',]
    
    sample_info <- arrange(sample_info, sample_info[[group_type]], Replicate)
    sample_info <- rbind(sample_info, sample_info_spqc)
    cat(file=stderr(), "The sample design was sorted alphabetically with SPQC at bottom", "\n")
  }
  
  sample_info$check <- NULL
  return(sample_info)
}
