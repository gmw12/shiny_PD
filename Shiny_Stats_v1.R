  #create data frame for comparisons
set_comp_groups <- function(){
  comp_groups <- data.frame(seq(from=1, to=as.numeric(dpmsr_set$x$comp_number)))
  colnames(comp_groups) <- "CompNumber"
  sample_groups<-data.frame(dpmsr_set$y$sample_groups)
  rownames(sample_groups)<-sample_groups[,1]
  for(i in 1:as.numeric(dpmsr_set$x$comp_number)){
    comp_N <- as.character(dpmsr_set$x[[str_c("comp",i,"N")]])
    comp_D <- as.character(dpmsr_set$x[[str_c("comp",i,"D")]])
    comp_groups$comp_N[i] <- comp_N 
    comp_groups$comp_D[i] <- comp_D
    comp_groups$N_start[i] <- sample_groups$start[which(sample_groups$Group == comp_N)]
    comp_groups$N_end[i] <- sample_groups$end[which(sample_groups$Group == comp_N)]
    comp_groups$D_start[i] <- sample_groups$start[which(sample_groups$Group == comp_D)]
    comp_groups$D_end[i] <- sample_groups$end[which(sample_groups$Group == comp_D)]
    comp_groups$comp_name[i] <- str_c(sample_groups$Group[which(sample_groups$Group == comp_N)], "_v_", sample_groups$Group[which(sample_groups$Group == comp_D)])
    comp_groups$fc[i] <- str_c(sample_groups$Group[which(sample_groups$Group == comp_N)], "_v_", sample_groups$Group[which(sample_groups$Group == comp_D)], "_FC")
    comp_groups$fc2[i] <- str_c(sample_groups$Group[which(sample_groups$Group == comp_N)], "_v_", sample_groups$Group[which(sample_groups$Group == comp_D)], "_FC2")
    comp_groups$pval[i] <- str_c(sample_groups$Group[which(sample_groups$Group == comp_N)], "_v_", sample_groups$Group[which(sample_groups$Group == comp_D)], "_Pval")
    comp_groups$limma_pval[i] <- str_c(sample_groups$Group[which(sample_groups$Group == comp_N)], "_v_", sample_groups$Group[which(sample_groups$Group == comp_D)], "_LimmaPval")
    comp_groups$exactTest[i] <- str_c(sample_groups$Group[which(sample_groups$Group == comp_N)], "_v_", sample_groups$Group[which(sample_groups$Group == comp_D)], "_ExactTest")  
    comp_groups$mf[i] <- str_c(sample_groups$Group[which(sample_groups$Group == comp_N)], "_v_", sample_groups$Group[which(sample_groups$Group == comp_D)], "_MF") 
    comp_groups$Ncount <- sample_groups$end[which(sample_groups$Group == comp_N)] - sample_groups$start[which(sample_groups$Group == comp_N)] +1
    comp_groups$Dcount <- sample_groups$end[which(sample_groups$Group == comp_D)] - sample_groups$start[which(sample_groups$Group == comp_D)] +1
  }
  
  dpmsr_set$y$comp_groups <<- comp_groups
}

#--------------------------------------------------------------------------------------
apply_stats <- function(){
  #parallel stat calc
  ncores <- detectCores()
  data_list <- dpmsr_set$data$impute
  dpmsr_set$data$final <<- mclapply(data_list, stats_parallel, mc.cores = ncores)
  
  #check that the parallel processing went through, if not do it again one at a time
  check_stats_parallel(data_list)
  
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    dpmsr_set$data$finalraw <<- collapse_peptide(dpmsr_set$data$normalized$impute)
    Simple_Excel(dpmsr_set$data$finalraw, str_c(dpmsr_set$file$extra_prefix, "_final_raw_protein.xlsx", collapse = " "))
  } else{
    dpmsr_set$data$finalraw <<- dpmsr_set$data$normalized$impute
    Simple_Excel(dpmsr_set$data$finalraw, str_c(dpmsr_set$file$extra_prefix, "_final_raw_peptide.xlsx", collapse = " "))
  }
}

#--------------------------------------------------------------------------------------
stats_parallel <- function(data_list){
  stat_calc(data_list)
}

#--------------------------------------------------------------------------------------
check_stats_parallel <- function(data_list){
  for(data_name in names(data_list)){
    if(is.null(dpmsr_set$data$final[[data_name]]    )){
      cat(file=stderr(), str_c("stats Parallel function...", data_name), "\n")
      dpmsr_set$data$final[[data_name]] <<- stat_calc(dpmsr_set$data$impute[[data_name]])
    }
  }
}

#--------------------------------------------------------------------------------------------------------------------------------
#Fold change, pvalue, export volcano, return organized table for output-------------------------------------------------
stat_calc <- function(data_in) {
  if (dpmsr_set$y$state=="Peptide" && dpmsr_set$x$final_data_output == "Protein"){
    data_in <- collapse_peptide(data_in)
  }
  
  if (as.logical(dpmsr_set$x$peptide_ptm_out)){
    data_in <- data_in [grep(dpmsr_set$x$peptide_report_grep, data_in$Modifications),]
  }
  
  info_columns_final <- ncol(data_in)-dpmsr_set$y$sample_number
  annotate_in <- data_in[1:info_columns_final]
  data_in <- data_in[(info_columns_final+1):ncol(data_in)]
  #add imputed data frame for stats, data is 1 or 0
  imputed_df <- dpmsr_set$data$Protein_imputed_df[3:ncol(dpmsr_set$data$Protein_imputed_df)]
  #start df for stats
  stat_df <- annotate_in[1:1]
  #generate %CV's for each group
  for(i in 1:dpmsr_set$y$group_number) 
  {
    stat_df[ , str_c(dpmsr_set$y$sample_groups$Group[i], "_CV")] <- percentCV_gw(data_in[dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i]])
  } 

  #generate pvalue and FC for each comparison
  for(i in 1:as.numeric(dpmsr_set$x$comp_number)){
    comp_N_data <- data_in[dpmsr_set$y$comp_groups$N_start[i]:dpmsr_set$y$comp_groups$N_end[i]]
    comp_D_data <- data_in[dpmsr_set$y$comp_groups$D_start[i]:dpmsr_set$y$comp_groups$D_end[i]]
    stat_df[ ,dpmsr_set$y$comp_groups$fc[i]] <- foldchange_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$comp_groups$fc2[i]] <- foldchange_decimal_gw(comp_N_data, comp_D_data)
    stat_df[ ,dpmsr_set$y$comp_groups$pval[i]] <- pvalue_gw(comp_N_data, comp_D_data)
    #stat_df[ , comp_groups$limma_pval[i]] <- limma_gw(comp_N_data, comp_D_data, comp_groups$comp_name[i], plot_dir)
    #stat_df[ , comp_groups$exactTest[i]] <- exactTest_gw(comp_N_data, comp_D_data)
    
    #add column for imputation factor
    comp_N_imputed <- imputed_df[dpmsr_set$y$comp_groups$N_start[i]:dpmsr_set$y$comp_groups$N_end[i]]
    comp_D_imputed <- imputed_df[dpmsr_set$y$comp_groups$D_start[i]:dpmsr_set$y$comp_groups$D_end[i]]
    stat_df[ ,dpmsr_set$y$comp_groups$mf[i]] <- missing_factor_gw(comp_N_imputed, comp_D_imputed)
  } 
  
  # Create tables for excel--------------------------------------------------
  stat_df <- stat_df[2:ncol(stat_df)]
  colnames(data_in) <- dpmsr_set$design$Header2
  data_table <- cbind(annotate_in, data_in, stat_df)
  return(data_table)
}



#--------------------------------------------------------------------------------------------------------------------------------
# dostats <- function(data_in, data_title) {
#   # reduce to phos only
#   if (ptm_peptide_only){
#     data_in<- data_in[grep(peptide_report_grep, data_in$Modifications),]  
#     Plot_All_gw(data_in[(info_columns+1):ncol(data_in)], str_c(data_title, "_PTM"))
#   }
#   #generate plots for proteins of interest
#   plot_dir <- str_c(output_dir, data_title, "//")
#   BioID_normalize_gw(data_in, data_title, plot_dir)
# 
#   data_out <- stat_test_gw(data_in, data_title, plot_dir)
#   cv_stats_gw(data_out, data_title)
#   if (qc_spike){qc_spike_gw(data_in, data_title)}
#   #final excel output
#   Final_Excel_gw(data_out, str_c(file_prefix1, "_", data_title, "_final.xlsx"))
#   return(data_out)
# }


#qc spike metrics ---------------------------------
# qc_spike_gw <- function(data_in, data_title) {
#   spike_protein <-subset(data_in, Accession %in% qc_spike_id)  
#   spike_protein <- spike_protein[(info_columns_final+1):ncol(spike_protein)]
#   spike_total <- colSums(spike_protein)
#   
#   total_qc_spike <-get("total_qc_spike",.GlobalEnv)
#   total_qc_spike <- cbind(total_qc_spike, spike_total)
#   colnames(total_qc_spike)[colnames(total_qc_spike) == 'spike_total'] <- data_title
#   assign("total_qc_spike",total_qc_spike,.GlobalEnv)
# }


#qc spike metrics finalize---------------------------------
# qc_spike_final_gw <- function(data_in) {
#   data_in <- aggregate(data_in, by=list(Category=data_in$SpikeLevel), FUN=mean)
#   data_in$Category <- NULL
#   row_count <- nrow(data_in)
#   for(i in 2:row_count) {
#     data_in <- rbind(data_in, (data_in[i,]/data_in[1,]) )
#   }
#   Simple_Excel(data_in, str_c(output_dir, "QC_Spike.xlsx"))
#   return(data_in)
# }



#cohensD ---------------------------------
cohend_gw <- function(x,y) {
  cohend_est = try(cohen.d(as.numeric(x),as.numeric(y), 
                           na.rm = TRUE, pooled=FALSE, paired=FALSE))
  if (is(cohend_est, "try-error")) return(NA) else return(signif((cohend_est$estimate), digits = 3))
}

#missing factor ---------------------------------
missing_factor_gw <- function(x,y) {
  mf_x <- rowSums(x)/ncol(x)
  mf_y <- rowSums(y)/ncol(y)
  df_mf <- data.frame(cbind(mf_x, mf_y), stringsAsFactors = FALSE)
  df_mf$max <- apply(df_mf, 1, FUN = function(x) {max(x, na.rm=TRUE)})
  return(signif(df_mf$max, digits = 3))
}


#fold change ---------------------------------
foldchange_gw <- function(x,y) {
  if(!as.logical(dpmsr_set$x$pair_comp)){
    ave_x = rowMeans(x)
    ave_y = rowMeans(y)
    test = ave_x/ave_y
  }else{
    sn <- ncol(x)
    indiv_fc <- x
    for (i in 1:sn){
      indiv_fc[i] <- (x[i]/y[i])
    }
    test <- rowMeans(indiv_fc)
  }
  fc <- ifelse ((test >= 1), test, -1/test)
  return(signif(fc, digits = 3))
}

#fold change pair---------------------------------
foldchange_pair_gw <- function(x,y) {
  sn <- ncol(x)
  indiv_fc <- x
  for (i in 1:sn){
    indiv_fc[i] <- (x[i]/y[i])
  }
  test <- rowMeans(indiv_fc)
  fc <- ifelse ((test >= 1), test, -1/test)
  return(signif(fc, digits = 3))
}



#fold change decimal ---------------------------------
foldchange_decimal_gw <- function(x,y) {
  if(!as.logical(dpmsr_set$x$pair_comp)){
    ave_x = rowMeans(x)
    ave_y = rowMeans(y)
    test = ave_x/ave_y
  }else{
    sn <- ncol(x)
    indiv_fc <- x
    for (i in 1:sn){
      indiv_fc[i] <- (x[i]/y[i])
    }
    test <- rowMeans(indiv_fc) 
  }
  fc <- test
  return(signif(fc, digits = 3))
}

#fold change pair decimal---------------------------------
foldchange_pair_decimal_gw <- function(x,y) {
  sn <- ncol(x)
  indiv_fc <- x
  for (i in 1:sn){
    indiv_fc[i] <- (x[i]/y[i])
  }
  test <- rowMeans(indiv_fc)
  fc <- test
  return(signif(fc, digits = 3))
}


#Percent CV ---------------------------------
percentCV_gw <- function(x) {
  ave <- rowMeans(x)
  n <- ncol(x)
  sd <- apply(x[1:n], 1, sd)
  cv <- (100 * sd / ave)
  return(signif(cv, digits = 3))
}


#t.test ---------------------------------
ttest_gw <- function(x,y) {
  if(!as.logical(dpmsr_set$x$pair_comp)){
    ttest_pvalue = try(t.test(x,y, 
                              alternative="two.sided",
                              var.equal = FALSE,
                              paired = FALSE), silent=TRUE)
  }else{
    ttest_pvalue = t.test(x,y, paired = TRUE)
  }
  if (is(ttest_pvalue, "try-error")) return(NA) else return(signif((ttest_pvalue$p.value), digits = 3))
}


pvalue_gw <- function(x,y){
  x <- log2(x)
  y <- log2(y)
  temp_pval <- rep(NA, nrow(x))
  for(i in 1:nrow(x)) 
  {
    temp_pval[i] <- ttest_gw(as.numeric(x[i,]), as.numeric(y[i,]) )   
  }
  return(temp_pval)
}

exactTest_gw <- function(x,y){
  #x <- log2(x)
  #y <- log2(y)
  et <- exactTestDoubleTail(x,y) 
  return(et)
} 
  


#x<-comp_N_data
#y<-comp_D_data
#comp_name <- comp_groups$comp_name[1]
limma_gw <- function(x,y, comp_name, plot_dir){
  xy <- cbind(x,y)
  xy <- log2(xy)
  n <- ncol(x)
  d <- ncol(y)
  design <- model.matrix(~ 0+factor(c(rep(1,n), rep(0,d))))
  colnames(design) <- c("group1", "group2")
  contrast.matrix <- makeContrasts(group2-group1, levels=design)
  fit <- lmFit(xy,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  topfit <- topTable(fit2,coef=1, sort="none", number=Inf)
  data_out <-topfit$P.Value
  try(limma_qq(fit2$t, comp_name, plot_dir), silent = TRUE)
  try(limma_ma(topfit, comp_name, plot_dir), silent = TRUE)
  try(limma_volcano(fit2, comp_name, plot_dir), silent = TRUE)
  return(data_out)
}


old_limma_gw <- function(x,y){
  xy <- cbind(x,y)
  xy <- log2(xy)
  n <- ncol(x)
  d <- ncol(y)
  design <- cbind(Grp1=1,Grp2vs1=c(rep(1,n), rep(0,d)))
  #design <- c(0,0,0,1,1,1)
  # Ordinary fit
  fit <- lmFit(xy,design)
  fit <- eBayes(fit)
  topfit <- topTable(fit,coef=2, adjust.method="BH", sort="none", number=Inf)
  data_out <-topfit$P.Value
  return(data_out)
}






