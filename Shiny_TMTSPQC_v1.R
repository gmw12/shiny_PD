tmt_spqc_prep <- function(){
  # create subdirectory to store csv and plots
  dpmsr_set$file$TMT_dir <<- create_dir(str_c(dpmsr_set$file$data_dir,"/TMT_IRS"))
  bar_plot_TMT(dpmsr_set$data$normalized$sl, str_c("TMT_IRS_Step0"), dpmsr_set$file$TMT_dir, dpmsr_set$design)
  
}
#-------------------------------------------------------------------------------------
# TMT SPQC - normalize TMT, SL, IRS
#-------------------------------------------------------------------------------------
tmt_spqc_normalize <- function(data_in){
  #Reset incase reanlayzing in same session
  if(dpmsr_set$x$raw_data_input == "Protein_Peptide" || dpmsr_set$x$raw_data_input == "Peptide")
  {
    dpmsr_set$y$state <<- "Peptide" 
  }
  sets_tmt <- as.numeric(dpmsr_set$x$tmt_spqc_sets)
  tmt_kit <- as.numeric(dpmsr_set$x$tmt_spqc_channels)
  sd_factor <- as.numeric(dpmsr_set$x$tmt_spqc_sd_factor_filter)

  dpmsr_set$y$info_columns <<- ncol(data_in) - dpmsr_set$y$sample_number
#--------------------------------------------------------------------------------------------------------------  
  
#---Step 1-2, split combined data into indiv sets, remove empty rows------------------------------------------------------------
  cat(file=stderr(), "TMT IRS step1,2...", "\n")
  for(i in 1:sets_tmt){
    #create design for each TMT set
    irs_design <- dpmsr_set$design[dpmsr_set$design$Set==i,]
    assign(str_c("TMT",i,"_design"), irs_design)
    #assign(str_c("Test",i,"_design"), irs_design, envir = .GlobalEnv)
    #create df for each TMT set
    irs_set <- data_in[,(dpmsr_set$y$info_columns + dpmsr_set$design[dpmsr_set$design$Set==i,]$PD_Order)]
    irs_set <- cbind(data_in[1:dpmsr_set$y$info_columns], irs_set)
    assign(str_c("Test_",i,"_Step1_data"), irs_set, envir = .GlobalEnv)
    assign(str_c("IRS_",i,"_Step1_data"), irs_set)
    Simple_Excel(irs_set, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step1.xlsx", collapse = " "))
    bar_plot_TMT(irs_set, str_c("TMT", i, "_IRS_Step1"),dpmsr_set$file$TMT_dir,irs_design)
    #remove empty rows 
    irs_set$na_count <- apply(irs_set[(dpmsr_set$y$info_columns+1):ncol(irs_set)], 1, function(x) sum(!is.na(x)))
    irs_set <- subset(irs_set, irs_set$na_count > 1)
    irs_set$na_count <- NULL
    #irs_set[is.na(irs_set)] <- 0.0
    #assign(str_c("Test_",i,"_Step2_data"), irs_set, envir = .GlobalEnv)
    assign(str_c("IRS_",i,"_Step2_data"), irs_set)
    Simple_Excel(irs_set, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step2.xlsx", collapse = " "))
    #bar_plot_TMT(irs_set, str_c("TMT", i, "_IRS_Step2"), dpmsr_set$file$TMT_dir,irs_design)
  } 
  
#----------------------------------------------------------------------------------------------------    
  
#---Step3, Impute missing values with random value from bottom 5% of measured values
  cat(file=stderr(), "TMT IRS step3...", "\n")
  data_dist <- as.vector(t(data_in[(dpmsr_set$y$info_columns+1):ncol(data_in)]))
  data_dist <- data_dist[!is.na(data_dist)]
  data_dist <- data_dist[data_dist>0]
  data_dist <- data.frame(data_dist)
  data_dist <- log2(data_dist)
  data_dist$bin <- ntile(data_dist, 20)  
  data_dist <- data_dist[data_dist$bin==1,]
  bottom5_max <- max(data_dist$data_dist)
  bottom5_min <- min(data_dist$data_dist)
  
  for(i in 1:sets_tmt){
    irs_df_temp <- get(str_c("IRS_",i,"_Step2_data"))
    irs_info <- irs_df_temp[1:dpmsr_set$y$info_columns]
    irs_df_temp <- irs_df_temp[(dpmsr_set$y$info_columns+1):ncol(irs_df_temp)]
    irs_df_temp <- log2(irs_df_temp)
  
    for (j in 1:nrow(irs_df_temp)){
      for (k in 1:ncol(irs_df_temp)){
        if ( is.na(irs_df_temp[j,k])) {irs_df_temp[j,k] <- runif(1, bottom5_min, bottom5_max) } 
      }
    }
    irs_df_temp <- 2^irs_df_temp
    irs_df_temp <- cbind(irs_info, irs_df_temp)
    assign(str_c("IRS_",i,"_Step3_data"), irs_df_temp)
    Simple_Excel(irs_set, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step3.xlsx", collapse = " "))
    #bar_plot_TMT(irs_set, str_c("TMT", i, "_IRS_Step3"), dpmsr_set$file$TMT_dir,irs_design)
}


#----------------------------------------------------------------------------------------------------  

#---Step4------Filter data based on CV of SPQC if requested, collapse peptide to protein
  cat(file=stderr(), "TMT IRS step4...", "\n")
  for(i in 1:sets_tmt){
    if (as.logical(dpmsr_set$x$tmt_filter_peptide)){
      filter_data <- get(str_c("IRS_",i,"_Step3_data"))
      irs_design <- get(str_c("TMT",i,"_design"))
      spqc_set <- subset(irs_set, irs_set$Group == "SPQC")
      spqc_temp <- filter_data %>% dplyr::select(contains("SPQC"))
      spqc_count <- ncol(spqc_temp)
      
      if (spqc_count>1){
        spqc_temp$average <- apply(spqc_temp[1:spqc_count], 1, FUN = mean)
        spqc_temp$stddev <- apply(spqc_temp[1:spqc_count], 1, FUN = sd)
        spqc_temp$cv <- spqc_temp$stddev / spqc_temp$average * 100
        spqc_temp$bin <- ntile(spqc_temp$average, 5)  
        sd_info <- subset(spqc_temp) %>% group_by(bin) %>% 
          summarize(min = min(average), max = max(average), min_cv= min(cv), max_cv = max(cv), av_cv = mean(cv), sd_cv = sd(cv))
        spqc_temp$maxcv <- lapply(spqc_temp$bin, function(x) (sd_info$av_cv[x] + (as.numeric(dpmsr_set$x$tmt_spqc_sd_factor_filter) * sd_info$sd_cv[x] )))
        spqc_temp$maxcv <- as.numeric(spqc_temp$maxcv)
        spqc_temp$pass <- spqc_temp$maxcv - spqc_temp$cv
        filter_data <- cbind(filter_data, spqc_temp$pass)
        filter_data <- subset(filter_data, filter_data$`spqc_temp$pass` >0)
        filter_data$`spqc_temp$pass` <- NULL
      }
      
      assign(str_c("IRS_",i,"_Step4_data"), filter_data)
      Simple_Excel(filter_data, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step4.xlsx", collapse = " "))
      bar_plot_TMT(filter_data, str_c("TMT", i, "_IRS_Step4"), dpmsr_set$file$TMT_dir, get(str_c("TMT",i,"_design")))
    }else{
      assign(str_c("Test_",i,"_Step4_data"), get(str_c("IRS_",i,"_Step3_data")), envir = .GlobalEnv) 
      assign(str_c("IRS_",i,"_Step4_data"), get(str_c("IRS_",i,"_Step3_data")) )
    }
  }


#----------------------------------------------------------------------------------------------------  

#---Step5----collapse peptide to protein  
  cat(file=stderr(), "TMT IRS step5...", "\n")
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    for(i in 1:sets_tmt){
        if (exists(str_c("IRS_",i,"_Step4_data"))){
          peptide_data <- get(str_c("IRS_",i,"_Step4_data"))
        }else{
          peptide_data <- get(str_c("IRS_",i,"_Step3_data"))
        }
        peptide_annotate <- peptide_data[, c("Accession", "Description")]
        peptide_annotate$Peptides <- 1
        peptide_annotate$Peptides <- as.numeric(peptide_annotate$Peptides)
        peptide_data <- peptide_data[(dpmsr_set$y$info_columns+1):ncol(peptide_data)]
        peptide_data[is.na(peptide_data)] <- 0
        peptide_data <- cbind(peptide_annotate, peptide_data)
        protein_data <- peptide_data %>% group_by(Accession, Description) %>% summarise_all(list(sum))
        protein_data <- data.frame(ungroup(protein_data))
        assign(str_c("IRS_",i,"_Step5_data"), protein_data)
        assign(str_c("Test_",i,"_Step5_data"), protein_data, envir = .GlobalEnv)
        Simple_Excel(protein_data, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step5.xlsx", collapse = " "))
        #bar_plot_TMT(protein_data, str_c("TMT", i, "_IRS_Step5"), dpmsr_set$file$TMT_dir, get(str_c("TMT",i,"_design")))
    }
    #set info column for protein data
    dpmsr_set$y$info_columns<<-3
  }else{
    assign(str_c("IRS_",i,"_Step5_data"), get(str_c("IRS_",i,"_Step4_data")))
  }
 
#----------------------------------------------------------------------------------------------------    

#---Step6--identify common proteins or peptides
  cat(file=stderr(), "TMT IRS step6...", "\n")
  if (dpmsr_set$x$final_data_output == "Protein")
  {
    # create vector of common accession for protein dataframes
    common_data <- unlist(data.frame(get("IRS_1_Step5_data"))[ , "Accession"])
    for(i in 2:sets_tmt){
      common_data2 <- unlist(data.frame(get(str_c("IRS_",i,"_Step5_data"))[ , "Accession"]))
      common_data <- intersect(common_data, common_data2)
    }
    for(i in 1:sets_tmt){
      temp_data <- get(str_c("IRS_",i,"_Step5_data"))
      temp_data <- subset(temp_data, temp_data$Accession %in% common_data)
      temp_data <- temp_data[order(temp_data$Accession),]
      assign(str_c("IRS_",i,"_Step6_data"), temp_data)
      Simple_Excel(temp_data, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step6.xlsx", collapse = " "))
      #bar_plot_TMT(temp_data, str_c("TMT", i, "_IRS_Step6"), dpmsr_set$file$TMT_dir, get(str_c("TMT",i,"_design")))
      #combine peptide count from multiple sets into one column for output later
      if(i==1){
        peptide_count <- temp_data$Peptides
      }else{
        peptide_count <- str_c(peptide_count, ", ", temp_data$Peptides)
      }
    }
  }else{
    # for peptide - NOT TESTED
    IRS_1_Step5_data$uniqueID <- paste(IRS_1_Step5_data$Accession, IRS_1_Step5_data$Sequence, IRS_1_Step5_data$Modifications)
    common_data <- unlist(IRS_1_Step5_data$uniqueID)
    for(i in 2:sets_tmt){
      temp_IRS <- get(str_c("IRS_",i,"_Step5_data"))
      temp_IRS$uniqueID <- paste(temp_IRS$Accession, temp_IRS$Sequence, temp_IRS$Modifications)
      assign(str_c("IRS_",i,"_Step5_data"), temp_IRS)
      common_data2 <- unlist(temp_IRS$uniqueID)
      common_data <- intersect(common_data, common_data2)
    }
    for(i in 1:sets_tmt){
      temp_data <- get(str_c("IRS_",i,"_Step5_data"))
      temp_data <- subset(temp_data, temp_data$uniqueID %in% common_data)
      temp_data <- temp_data[order(temp_data$uniqueID),]
      temp_data$uniqueID <- NULL
      assign(str_c("IRS_",i,"_Step6_data"), temp_data)
      Simple_Excel(temp_data, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step6.xlsx", collapse = " "))
      #bar_plot_TMT(temp_data, str_c("TMT", i, "_IRS_Step6"), dpmsr_set$file$TMT_dir, get(str_c("TMT",i,"_design")))
    }
  } 
  

#----------------------------------------------------------------------------------------------------    
  
#---Step7--IRS normalization-------------------------------------------------------------------
  cat(file=stderr(), "TMT IRS step7...", "\n")
  #create SPQC data frames for each set, with sums
  for(i in 1:sets_tmt){
    irs_data <- get(str_c("IRS_",i,"_Step6_data"))
    spqc_temp <- irs_data %>% dplyr::select(contains("SPQC"))
    spqc_temp$sums <- rowSums(spqc_temp)
    assign(str_c("SPQC",i,"_data"), spqc_temp)
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
    temp_data <- get(str_c("IRS_",i,"_Step6_data"))
    #temp_data <- get(str_c("IRS_",i,"_Step6_data"))[(dpmsr_set$y$info_columns+1):(dpmsr_set$y$info_columns+tmt_kit)]
    spqc_data <- get(str_c('SPQC',i,'_data'))
    temp_data[(dpmsr_set$y$info_columns+1):(dpmsr_set$y$info_columns+tmt_kit)] <- temp_data[(dpmsr_set$y$info_columns+1):(dpmsr_set$y$info_columns+tmt_kit)]/spqc_data$fact
    assign(str_c("IRS_",i,"_Step7_data"), temp_data)
    Simple_Excel(temp_data, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Step7.xlsx", collapse = " "))
    #bar_plot_TMT(temp_data, str_c("TMT", i, "_IRS_Step7"), dpmsr_set$file$TMT_dir, get(str_c("TMT",i,"_design")))
  }
  
  
  #----------------------------------------------------------------------------------------------------    
  
  #---Step8--Recombine Data------------------------------------------------------------------- 
  
  #recombine sets for output
  final_SL_IRS <- IRS_1_Step7_data
  for(i in 2:sets_tmt){
    final_SL_IRS <- cbind(final_SL_IRS, get(str_c("IRS_",i,"_Step7_data"))[(dpmsr_set$y$info_columns+1):(dpmsr_set$y$info_columns+tmt_kit)] )
  }
  final_SL_IRS$Peptides <- peptide_count
  
  #reorder combined set in original order
  set_original <- TMT1_design
  for(i in 2:sets_tmt){
    irs_set <- get(str_c('TMT',i,'_design'))
    set_original <- rbind(set_original, irs_set) 
  }
  set_original$reorder <- seq(1,(tmt_kit*sets_tmt))
  set_original <- set_original[order(set_original$PD_Order),]
  final_info_columns <- final_SL_IRS[1:dpmsr_set$y$info_columns] 
  final_SL_IRS<-final_SL_IRS[(dpmsr_set$y$info_columns+1):ncol(final_SL_IRS)] 
  final_SL_IRS <- final_SL_IRS[ ,set_original$reorder]
  final_SL_IRS <- cbind(final_info_columns, final_SL_IRS)
  
  Simple_Excel(final_SL_IRS, str_c(dpmsr_set$file$TMT_dir, "TMT", i , "_Final.xlsx", collapse = " "))
  bar_plot_TMT(final_SL_IRS, str_c("TMT_IRS_Final"), dpmsr_set$file$TMT_dir, dpmsr_set$design)
  


#--------------------------------------------------------------------
#Combined Bar plot-------------------------------------------------
for(j in 1:7){
    step <- str_c("Step", j)
    for(i in 1:sets_tmt){
     plot_data <- get(str_c("IRS_",i,"_",step,"_data"))
     plot_design <- get(str_c("TMT",i,"_design"))
     info_columns <- ncol(plot_data) - nrow(plot_design)
     plot_data<- plot_data[(info_columns+1):ncol(plot_data)]
     namex <- plot_design$Label
     datay <- colSums(plot_data, na.rm = TRUE)
     if(i==1){
        df2 <- data.frame(namex)
        df2$Total_Intensity <- datay
        plot_design_all <- plot_design
     }else{
        df3 <- data.frame(namex)
        df3$Total_Intensity <- datay
        df2 <- rbind(df2,df3)
        plot_design_all <- rbind(plot_design_all, plot_design)
      }
    }
    colnames(df2) <- c("Sample", "Total_Intensity")
    df2$Sample <- factor(df2$Sample, levels = df2$Sample)
    file_name <- str_c(dpmsr_set$file$TMT_dir, "TMT_IRS_", step, "_barplot.png")
    ymax <- max(datay)
    ggplot(data=df2, aes(x=Sample, y=Total_Intensity)) +
      geom_bar(stat="identity", fill=plot_design_all$colorlist)+ theme_classic() + 
      ggtitle(str_c("TMT_IRS_", step)   ) + 
      xlab(NULL)+
      ylab(NULL)+
      #scale_y_discrete(labels = NULL) +
      coord_cartesian(ylim=NULL, expand = TRUE) +
      theme(plot.title = element_text(hjust = 0.5), 
            axis.text.x = element_text(size=5, angle = 90,  color="black"),
            axis.text.y = element_text(size=5,  color="black"),
      ) 
    ggsave(file_name, width=5, height=4)
}
  #save files
  dpmsr_set$data$impute$impute <<- final_SL_IRS 
  dpmsr_set$data$impute$sl <<- final_SL_IRS 
  dpmsr_set$y$state <<- "Protein"
  dpmsr_set$y$info_columns_final <<- 3
  return()
}    
    
    

#-----------------------------
# TNT specific functions
#-----------------------------


#Bar plot-------------------------------------------------
bar_plot_TMT <- function(df,plot_title,plot_dir,irs_design) {
  info_columns <- ncol(df) - nrow(irs_design)
  df <- df[(info_columns+1):ncol(df)]
  namex <- irs_design$Label
  datay <- colSums(df, na.rm = TRUE)
  df2 <- data.frame(namex)
  df2$Total_Intensity <- datay
  colnames(df2) <- c("Sample", "Total_Intensity")
  df2$Sample <- factor(df2$Sample, levels = df2$Sample)
  file_name <- str_c(plot_dir, plot_title, "_barplot.png")
  ymax <- max(datay)
  ggplot(data=df2, aes(x=Sample, y=Total_Intensity)) +
    geom_bar(stat="identity", fill=irs_design$colorlist)+ theme_classic() + 
    ggtitle(plot_title) + 
    xlab(NULL)+
    ylab(NULL)+
    #scale_y_discrete(labels = NULL) +
    coord_cartesian(ylim=NULL, expand = TRUE) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(size=5, angle = 90,  color="black"),
          axis.text.y = element_text(size=5,  color="black"),
    ) 
  ggsave(file_name, width=5, height=4)
  return("done")
}