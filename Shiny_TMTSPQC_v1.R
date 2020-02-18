
#-------------------------------------------------------------------------------------
# TMT SPQC - normalize TMT, SL, IRS
#-------------------------------------------------------------------------------------
tmt_spqc_normalize <- function(){
  
  sets_tmt <- as.numeric(design_info[36,2])
  tmt_kit <- as.numeric(design_info[37,2])
  filter_peptides <- as.logical(design_info[38,2])
  sd_factor <- as.numeric(design_info[39,2])
  
  all_data <- forward_data
  
  # create subdirectory to store csv and plots
  TMT_dir <- create_dir(str_c(data_dir,"//TMT_IRS"))
  TMT_prefix <- TMT_dir #str_c(TMT_dir, "_")
  testthis <- str_c(TMT_dir, basename(input_data))
  #file.copy(input_data, str_c(TMT_dir, basename(input_data)))
  #file.copy(study_design, str_c(TMT_dir, basename(study_design)))
  
  #function Bar plot-------------------------------------------------
  barplot_TMT <- function(x,y) {
    x[is.na(x)] <- 0.0
    x <- colSums(x)
    png(filename=str_c(TMT_dir, y, "_barplot.png"), width = 888, height = 571)  
    barplot(x, 
            col = color_list,
            names.arg = treatment_groups,
            cex.names = .8,
            las=2,
            main = y)
    dev.off()
  }
  
  # function create column to flag and delete rows with no data
  delete_empty_row <- function(data_in){
    data_in$na_count <- apply(data_in[(info_count+1):ncol(data_in)], 1, function(x) sum(is.na(x)))
    data_in <- subset(data_in, data_in$na_count != tmt_kit)
    data_in <- data_in[1:(info_count + tmt_kit)]
    return(data_in)
  }
  
  
  #--- collapse peptide to protein-------------------------------------------------------------
  collapse_peptide_tmt <- function(peptide_data){
    peptide_annotate <- peptide_data[, c("Accession", "Description")]
    peptide_data <- peptide_data[(info_count+1):ncol(peptide_data)]
    peptide_annotate$Peptides <- 1
    peptide_annotate$Peptides <- as.numeric(peptide_annotate$Peptides)
    test1 <- cbind(peptide_annotate, peptide_data)
    test2 <- test1 %>% group_by(Accession, Description) %>% summarise_all(list(sum))
    test2 <- ungroup(test2)
    test2 %>% mutate_if(is.numeric, ~round(., 1))
    return(test2)
  } 
  
  #break down PD output to peptide only with protein assignments
  if(protein_peptide_input)
  {
    tmt_data <- protein_to_peptide(all_data)
  }else{
    tmt_data <- all_data
  }
  
  excel_order <- sample_info$PD_Order
  info_count <- ncol(tmt_data) - sample_number
  info_columns <- tmt_data[1:info_count]
  original_info_headers <- names(info_columns)
  all_peptide <- tmt_data #temp
  #all_peptide <-  all_peptide[(info_count+1):ncol(all_peptide)][, (excel_order)]
  tmt_data <-  tmt_data[(info_count+1):ncol(tmt_data)]
  
  
  if(protein_peptide_input)
  {
    input_type <- "Peptide"
    output_type <- "Protein"
    raw_peptide <- cbind(info_columns, all_peptide)
    Simple_Excel(raw_peptide, str_c(TMT_prefix, "TMT_raw_", input_type, ".xlsx"))
  }else if (data_input_type =="TMT SPQC Phos"){
    input_type <- "Peptide"
    output_type <- "Peptide"
  }else{
    input_type <- "Protein"
    output_type <- "Protein"
  }
  
  barplot_TMT(tmt_data, str_c("Raw_", input_type))
  
  
  # split combined data from PD into individual sets
  start_col<-1
  for(i in 1:sets_tmt){
    irs_set <- subset(sample_info, sample_info$Set==i)
    irs_set$setorder <- seq(start_col, (nrow(irs_set)+start_col-1))
    irs_set$sequence <- seq(1, nrow(irs_set))
    temp_data <- cbind(info_columns, tmt_data[ ,irs_set$PD_Order])
    temp_data <- delete_empty_row(temp_data)
    temp_data[is.na(temp_data)] <- 0.0
    assign(str_c("Raw",i,"_input"), temp_data)
    assign(str_c("IRS",i,"_set"), irs_set)
    treatment_groups<-irs_set$Group
    color_list<-irs_set$colorlist
    barplot_TMT(temp_data[(info_count+1):ncol(temp_data)], str_c("Raw_", input_type, "_",i))
    start_col <- start_col + tmt_kit
  }
  
  
  # SL normalize sets, using target created from all TMT sets
  #sum each channel and then take mean of all channels
  temp_data <- get(str_c("Raw",1,"_input"))
  target <- colSums(temp_data[(info_count+1):ncol(temp_data)])
  for(i in 2:sets_tmt){
    temp_data <- get(str_c("Raw",i,"_input"))
    target <- c(target, colSums(temp_data[(info_count+1):ncol(temp_data)]))
  }
  target <- mean(target)
  
  for(i in 1:sets_tmt){
    temp_data <- get(str_c("Raw",i,"_input"))
    irs_info <- temp_data[1:info_count]
    temp_data <- temp_data[(info_count+1):ncol(temp_data)]
    norm_facs <- target / colSums(temp_data)
    temp_sl <- sweep(temp_data, 2, norm_facs, FUN = "*")
    color_list<-irs_set$colorlist
    barplot_TMT(temp_sl, str_c("SL_" , input_type, "_",i))
    temp_sl <- cbind(irs_info, temp_sl)
    assign(str_c("Norm",i,"_SL_input"), temp_sl)
    Simple_Excel(temp_sl, str_c(TMT_prefix, "Norm",i, "_SL_", input_type, ".xlsx"))
  }
  
  #Filter data based on CV of SPQC if requested, collapse peptide to protein
  for(i in 1:sets_tmt){
    temp_data <- get(str_c("Norm",i,"_SL_input"))
    if (filter_peptides){
      irs_set <- get(str_c("IRS",i,"_set"))
      spqc_set <- subset(irs_set, irs_set$Group == "SPQC")
      spqc_count <- nrow(spqc_set)
      #this assumes the spqc are orded in the front of the set
      if (spqc_count >1){
        temp_data$average <- apply(temp_data[(info_count+1):(info_count+spqc_count)], 1, FUN = mean)
        temp_data$stddev <- apply(temp_data[(info_count+1):(info_count+spqc_count)], 1, FUN = sd)
        temp_data$cv <- temp_data$stddev / temp_data$average * 100
        temp_data$bin <- ntile(temp_data$average, 5)  
        #histogram_tmt(temp_data$cv,"Set_CV","Log2 CV Distribution") 
        sd_info <- subset(temp_data) %>% group_by(bin) %>% 
          summarize(min = min(average), max = max(average), min_cv= min(cv), max_cv = max(cv), av_cv = mean(cv), sd_cv = sd(cv))
        temp_data$maxcv <- NA  
        for (j in 1:nrow(temp_data)){
          temp_data$maxcv[j] <- sd_info$av_cv[temp_data$bin[j]] + (sd_factor * sd_info$sd_cv[temp_data$bin[j]])
        }
        temp_data <- subset(temp_data, temp_data$cv < temp_data$maxcv)
      }else{
        temp_data <- subset(temp_data, temp_data[spqc_set$sequence]>0)
      }
      
      temp_data <- temp_data[1:(info_count+tmt_kit)]
      color_list<-irs_set$colorlist
      barplot_TMT(temp_data[(info_count+1):ncol(temp_data)], str_c("Norm_Peptide_Filtered_",i))
      assign(str_c("Norm",i,"_SL_filtered_peptide"), temp_data)
      Simple_Excel(temp_data, str_c(TMT_prefix, 'Norm_filtered',i,'_peptide.xlsx'))
    }
    
    #collapse to protein level
    if (protein_peptide_input)  {
      temp_protein <- collapse_peptide_tmt(temp_data)
    }else{
      temp_protein <- temp_data
    }
    #barplot_TMT(temp_protein[], str_c("SL_Protein_Filtered_",i))
    assign(str_c("IRS",i,"_output"), temp_protein)
    if (filter_peptides){
      barplot_TMT(temp_protein[(ncol(temp_protein)-tmt_kit+1):ncol(temp_protein)], str_c("SL_", output_type, "_Filtered_",i))
      Simple_Excel(temp_protein, str_c(TMT_prefix, "SL_Filtered",i, "_", output_type, ".xlsx"))
    }else {
      barplot_TMT(temp_protein[(ncol(temp_protein)-tmt_kit+1):ncol(temp_protein)], str_c("SL_", output_type, "_",i))
      Simple_Excel(temp_protein, str_c(TMT_prefix, "SL_",i, "_", output_type, ".xlsx"))
    }
  }
  
  
  #identify common proteins or peptides
  if (output_type == "Protein")
  {
    common_data <- unlist(data.frame(get("IRS1_output")[ , "Accession"]))
    for(i in 2:sets_tmt){
      common_data2 <- unlist(data.frame(get(str_c("IRS",i,"_output"))[ , "Accession"]))
      common_data <- intersect(common_data, common_data2)
    }
    for(i in 1:sets_tmt){
      temp_data <- get(str_c("IRS",i,"_output"))
      temp_data <- subset(temp_data, temp_data$Accession %in% common_data)
      temp_data <- temp_data[order(temp_data$Accession),]
      barplot_TMT(temp_data[(ncol(temp_data)-tmt_kit+1):ncol(temp_data)], str_c("SL_Common", output_type, "_Filtered_",i))
      assign(str_c("IRS",i,"_common"), subset(temp_data, temp_data$Accession %in% common_data))
    }
  }else{
    IRS1_output$uniqueID <- paste(IRS1_output$`Master Protein Accessions`, IRS1_output$`Annotated Sequence`, IRS1_output$Modifications)
    common_data <- unlist(data.frame(get("IRS1_output")[ , "uniqueID"]))
    for(i in 2:sets_tmt){
      temp_IRS <- get(str_c("IRS",i,"_output"))
      temp_IRS$uniqueID <- paste(temp_IRS$`Master Protein Accessions`, temp_IRS$`Annotated Sequence`, temp_IRS$Modifications)
      assign(str_c("IRS",i,"_output"), temp_IRS)
      common_data2 <- unlist(data.frame(get(str_c("IRS",i,"_output"))[ , "uniqueID"]))
      common_data <- intersect(common_data, common_data2)
    }
    for(i in 1:sets_tmt){
      temp_data <- get(str_c("IRS",i,"_output"))
      temp_data <- subset(temp_data, temp_data$uniqueID %in% common_data)
      temp_data <- temp_data[order(temp_data$uniqueID),]
      temp_data$uniqueID <- NULL
      barplot_TMT(temp_data[(ncol(temp_data)-tmt_kit+1):ncol(temp_data)], str_c("SL_Common", output_type, "_Filtered_",i))
      assign(str_c("IRS",i,"_common"), temp_data)
      #assign(str_c("IRS",i,"_common"), subset(temp_data, temp_data$Accession %in% common_data))
    }
  } 
  
  
  
  
  #reset info columns for protein level data
  irs_info <- data.frame(IRS1_common[1:(ncol(IRS1_common)-tmt_kit)])
  info_count <- ncol(irs_info)
  
  #create SPQC data frames for each set, with sums
  for(i in 1:sets_tmt){
    spqc_data <- get(str_c("IRS",i,"_common"))
    irs_set <- get(str_c("IRS",i,"_set"))
    spqc_set <- subset(irs_set, irs_set$Group == "SPQC")
    spqc_set$sequence <- spqc_set$sequence + info_count
    spqc_data <- spqc_data[, (spqc_set$sequence)]
    spqc_data$sums <- rowSums(spqc_data)
    assign(str_c("SPQC",i,"_data"), spqc_data)
  }
  
  #create dataframe with SPQC sums, add column with geometric mean (of sums of each protein per set)
  spqc_all <- data.frame(SPQC1_data$sums)
  for(i in 2:sets_tmt){
    spqc_data <- get(str_c("SPQC",i,"_data"))
    spqc_all <- cbind(spqc_all, spqc_data$sums)
  }
  spqc_all$geoavg <- apply(spqc_all, 1, function(x) exp(mean(log(x))))
  
  
  for(i in 1:sets_tmt){
    spqc_data <- get(str_c("SPQC",i,"_data"))
    spqc_data$fact <- spqc_data$sums/spqc_all$geoavg
    assign(str_c("SPQC",i,"_data"), spqc_data)
  }
  
  #apply IRS norm to each set
  for(i in 1:sets_tmt){
    temp_data <- get(str_c('IRS',i,'_common'))[(info_count+1):(info_count+tmt_kit)]
    spqc_data <- get(str_c('SPQC',i,'_data'))
    temp_data <- temp_data/spqc_data$fact
    barplot_TMT(temp_data[(ncol(temp_data)-tmt_kit+1):ncol(temp_data)], str_c("SL_IRS_Common", output_type, "_Filtered_",i))
    assign(str_c("IRS",i,"_final"), temp_data)
    Simple_Excel(temp_data, str_c(TMT_prefix, "IRS_SL_Filtered",i, "_", output_type, ".xlsx"))
  }
  
  #recombine sets for output
  final_SL_IRS <- IRS1_final
  for(i in 2:sets_tmt){
    final_SL_IRS <- cbind(final_SL_IRS, get(str_c("IRS",i,"_final")))
  }
  #reorder combined set in original order
  set_original <- IRS1_set
  for(i in 2:sets_tmt){
    irs_set <- get(str_c('IRS',i,'_set'))
    set_original <- rbind(set_original, irs_set) 
  }
  set_original$reorder <- seq(1,(tmt_kit*sets_tmt))
  set_original <- set_original[order(set_original$PD_Order),]
  final_SL_IRS <- final_SL_IRS[ ,set_original$reorder]
  #colnames(irs_info)<-original_info_headers
  final_SL_IRS <- cbind(irs_info, final_SL_IRS)
  
  Simple_Excel(final_SL_IRS, str_c(TMT_prefix, "TMT_SL_IRS_Norm.xlsx"))
  
  treatment_groups<-sample_info$Group
  color_list<-sample_info$colorlist
  barplot_TMT(final_SL_IRS[(info_count+1):ncol(final_SL_IRS)], str_c("Final_SL_IRS_Common", output_type, "_Filtered_"))
  
  return(final_SL_IRS)
}

