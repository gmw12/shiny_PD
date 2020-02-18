
load_dpmsr_set <- function(session, input){
  design_list <- c("General", "QC", "Fill_Norm", "Filters", "Protein", "TMT_PTM")
  design_data <- input$design_file
  design_name <- unlist(design_data$files[[as.character(0)]][2])
  design<-read_excel(design_name, sheet="SampleList")
  dpmsr_set <<- list(design = design)
  temp_list<-list()
  for(i in design_list){
    design_tab<-read_excel(design_name, sheet=i, col_names = FALSE, skip = 1)
    for(j in 1:nrow(design_tab)){
      lname <- design_tab[j,1]
      ldata <- design_tab[j,2]
      temp_list[[as.character(lname)]] <- ldata
    }
  }
  dpmsr_set$x <<- temp_list
}


#create initial dpmsr_set S3 data class
dpmsr_set_create <- function(design){
  dpmsr_set <<- list(design=design)
  class(dpmsr_set) <<- "dpmsr_set"
}

#----------------------------------------------------------------------------------------
#load design elements into dpmsr_set
load_design <- function(session, input){
  design_list <- c("General", "QC", "Fill_Norm", "Filters", "Protein", "TMT_PTM")
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
  dpmsr_set$x$test <<- "test"
  dpmsr_set$x$peptide_grep <<- str_replace_all(dpmsr_set$x$peptide_grep, "/", "\\\\")
  dpmsr_set$x$peptide_report_grep <<- str_replace_all(dpmsr_set$x$peptide_report_grep, "/", "\\\\")
  dpmsr_set$x$int_cutoff <- as.numeric(dpmsr_set$x$int_cutoff)
}


#----------------------------------------------------------------------------------------
load_data <- function(session, input){
  raw_data <- input$raw_files
  for (i in 0:(length(raw_data$files)-1)   ){
    raw_name <- unlist(raw_data$files[[as.character(i)]][2])
    if (grepl("PeptideGroups.txt", raw_name)){
        dpmsr_set$data$data_raw_peptide <<- read.delim(raw_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        save_data(raw_name)
      } else if (grepl("Proteins.txt", raw_name)){
        dpmsr_set$data$data_raw_protein <<- read.delim(raw_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        save_data(raw_name)
      } else if (grepl("PSMs.txt", raw_name)){
        dpmsr_set$data$data_raw_psm <<- read.delim(raw_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        save_data(raw_name)
      } else if (grepl("MSMSSpectrumInfo.txt", raw_name)){
        dpmsr_set$data$data_raw_msms <<- read.delim(raw_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        save_data(raw_name)
      } else if (grepl("InputFiles.txt", raw_name)){
        dpmsr_set$data$data_raw_inputfiles <<- read.delim(raw_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        save_data(raw_name)
      }
  }
}

#----------------------------------------------------------------------------------------
prepare_data <- function(session, input) {  #function(data_type, data_file_path){
   data_type <- input$radio_input
   if(as.numeric(data_type) == 1){
    save_data(input$file_raw_protein_peptide)
    raw_data <<- input$file_raw_protein_peptide
    dpmsr_set$data$data_protein_peptide <<- read_excel(raw_data$datapath, 1)
    dpmsr_set$data$data_peptide <<- protein_to_peptide(dpmsr_set$data$data_protein_peptide)
    dpmsr_set$y$state <<- "Peptide"
  }else if(data_type ==2){
    save_data(input$file_raw_protein)
    raw_data <- input$file_raw_protein
    dpmsr_set$data$data_protein <<- read_excel(raw_data$datapath, 1)
    dpmsr_set$y$state <<- "Protein"
  }else if(data_type ==3){
    save_data(input$file_raw_peptide)
    raw_data <- input$file_raw_peptide
    dpmsr_set$data$data_peptide <<- read_excel(raw_data$datapath, 1)
    dpmsr_set$y$state <<- "Peptide"
  }else if(data_type ==4){
    save_data(input$file_raw_peptide)
    save_data(input$file_raw_peptide_decoy)
    save_data(input$file_raw_psm)
    save_data(input$file_raw_psm_decoy)
    raw_data <- input$file_raw_peptide
    dpmsr_set$data$data_peptide <<- read_excel(raw_data$datapath, 1)
    raw_data <- input$file_raw_peptide_decoy
    dpmsr_set$data$data_peptide_decoy <<- read_excel(raw_data$datapath, 1)
    raw_data <- input$file_raw_psm
    dpmsr_set$data$data_psm <<- read_excel(raw_data$datapath, 1)
    raw_data <- input$file_raw_psm_decoy
    dpmsr_set$data$data_psm_decoy <<- read_excel(raw_data$datapath, 1)
    dpmsr_set$data$data_peptide <<- peptide_psm_set_fdr()
    dpmsr_set$y$state <<- "Peptide"
  }else{msgBox <- tkmessageBox(title = "Whoops",
                           message = "Invalid input output in design file", icon = "info", type = "ok")}

   if(input$checkbox_isoform){
     save_data(input$file_raw_peptide_isoform)
     raw_data <- input$file_raw_peptide_isoform
     dpmsr_set$data$data_peptide_isoform <<- read_excel(raw_data$datapath, 1)
   }
}


#----------------------------------------------------------------------------------------
set_sample_groups<- function(){
  sample_info <- dpmsr_set$design
  
  # number of groups, comparisons are not hardcoded, R will extract information from excel SampleData file
  sample_number <- nrow(sample_info)
  #comp_number <- length(grep(x = colnames(sample_info), pattern = "^Comp"))
  #count number of samples in each group - add column
  sample_info$Count <- sapply(sample_info$Group, function(string) sum(string==sample_info$Group))
  
  # create unqiue group dataframe with sample number
  sample_groups <- sample_info[7:ncol(sample_info)]
  sample_groups <- sample_groups %>% distinct(Group, .keep_all = TRUE)
  group_number <- nrow(sample_groups)
  
  #assign start and end columns for each group in pd output
  sample_groups$start <- 1
  sample_groups$end <- sample_groups$Count[1]
  for(i in 2:(group_number)) {
    sample_groups$start[i] <- sample_groups$start[i-1] + sample_groups$Count[i-1]
    sample_groups$end[i] <- sample_groups$start[i] + sample_groups$Count[i] - 1
  }
  
  # #create data frame for comparisons
  # comp_groups <- data.frame(seq(from=1, to=comp_number))
  # colnames(comp_groups) <- "CompNumber"
  # for(i in 1:comp_number){
  #   comp_N <- grep("N", sample_groups[[i+1]], ignore.case=TRUE)
  #   comp_D <- grep("D", sample_groups[[i+1]], ignore.case=TRUE)
  #   comp_groups$comp_N[i] <- comp_N 
  #   comp_groups$comp_D[i] <- comp_D
  #   comp_groups$N_start[i] <- sample_groups$start[comp_N]
  #   comp_groups$N_end[i] <- sample_groups$end[comp_N]
  #   comp_groups$D_start[i] <- sample_groups$start[comp_D]
  #   comp_groups$D_end[i] <- sample_groups$end[comp_D]
  #   comp_groups$comp_name[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D])
  #   comp_groups$fc[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_FC")
  #   comp_groups$fc2[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_FC2")
  #   comp_groups$pval[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_Pval")
  #   comp_groups$limma_pval[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_LimmaPval")
  #   comp_groups$exactTest[i] <- str_c(sample_groups$Group[comp_N], "_v_", sample_groups$Group[comp_D], "_ExactTest")  
  #   comp_groups$Ncount <- sample_groups$end[comp_N] - sample_groups$start[comp_N] +1
  #   comp_groups$Dcount <- sample_groups$end[comp_D] - sample_groups$start[comp_D] +1
  # }
  #color_choices <- rainbow(group_number, s = 1, v = 1, start = 0, end = max(1, group_number - 1)/group_number, alpha = 1)
  color_choices <-distinctColorPalette(group_number)
  
  group_color <- color_choices[1:group_number]
  sample_groups$colorlist <- color_choices[1:group_number]
  sample_info$colorlist <- with(sample_groups, colorlist[match(sample_info$Group, Group)])
  
  #sample_groups$title <- str_c(sample_groups$Group,"(",sample_groups$colorlist,")")
  sample_groups$title <- sample_groups$Group

  group_title <- paste(sample_groups$title, collapse=" ")
  group_title <- gsub("\\)", "\\),",  group_title)
  
  #organize column headers for final output
  sample_info$Header1 <- str_c(sample_info$ID, " ", sample_info$Group)
  sample_info$Header2 <- str_c(sample_info$ID, " ", sample_info$Group, " Normalized")
  sample_info$Header3 <- str_c(sample_info$ID, " ", sample_info$Group, " Imputed")
  
  #create factor vector
  group_factor <- rep(1, sample_groups$Count[1])
  for(i in 2:nrow(sample_groups)) {
    group_factor <- c(group_factor, rep(i, sample_groups$Count[i]))
  }
  group_factor <-factor(group_factor)
  
  #save to dpmsr_set
  dpmsr_set$y$sample_number <<- sample_number
  dpmsr_set$y$sample_groups <<- sample_groups
  # dpmsr_set$y$comp_groups <<- comp_groups
  # dpmsr_set$y$comp_number <<- comp_number
  dpmsr_set$y$group_number <<- group_number
  dpmsr_set$y$group_factor <<- group_factor
  dpmsr_set$design <<- sample_info
  
  #create dataframe to hold cv summaries for normalization strategies
  dpmsr_set$data$summary_cv <<- data.frame(sample_groups$Group)
}



#----------------------------------------------------------------------------------------
isoform_set <- function(){
  norm_dir <- create_dir(str_c(dpmsr_set$file$data_dir,"//Norm"))
  norm_prefix <- str_c(norm_dir, file_prefix)
  if(ptm_peptide_only){norm_data <- special_norm(data_raw) }
  norm_data <- order_columns(norm_data)
  norm_data <- filter_data(norm_data)
  Simple_Excel(norm_data, str_c(norm_prefix, "_NormData.xlsx"))
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
dupliatedstuff <- function(){
  excel_order <- sample_info$PD_Order
  group_list <- sample_groups$Group
  treatment_groups <- sample_info$Group
  color_list<- sample_info$colorlist
}


