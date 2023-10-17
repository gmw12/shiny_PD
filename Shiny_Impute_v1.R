
apply_impute <- function(session, input, output){
  cat(file = stderr(), "apply_impute function...1", "\n")

  # save set of random numbers to be used for impute, if reimpute numbers the same
  set.seed(123)
  #get number of missing values
  total_missing <- table(is.na(dpmsr_set$data$normalized$impute))
  cat(file = stderr(), str_c(total_missing[2], " missing values"), "\n")
  dpmsr_set$y$rand_impute <<- runif(total_missing[2]*2, min = -1, max = 1)
  cat(file = stderr(), "created random number dataframe", "\n")
  dpmsr_set$y$rand_count <<- 1
  
  cat(file = stderr(), "apply_impute function...2", "\n")
  ncores <- detectCores()
  cat(file = stderr(), str_c("ncores = ", ncores), "\n")
  if (is.na(ncores) | ncores < 1) {ncores <- 1}
  norm_list <- dpmsr_set$y$norm_list
  dpmsr_set$data$impute <<- mclapply(norm_list, impute_parallel, mc.cores = ncores)
  #need to complete TMM norms after imputation ()
  
  cat(file = stderr(), "apply_impute function...3", "\n")
  #check that the parallel processing went through, if not do it again one at a time
  check_impute_parallel(norm_list)
  
  #check if directlfq norm was run and impute the protein level data
  check_impute_protein_directlfq()
  
  norm_list2 <- dpmsr_set$y$norm_list2
  
  info_columns <- ncol(dpmsr_set$data$impute$impute) - dpmsr_set$y$sample_number
  
  cat(file = stderr(), "apply_impute function...4", "\n")
  if (2 %in% dpmsr_set$y$norm_list2)
    {dpmsr_set$data$impute$tmm <<- tmm_normalize(dpmsr_set$data$impute$impute, "TMM_Norm", info_columns)
  }
  cat(file = stderr(), "apply_impute function...5", "\n")
  if (3 %in% dpmsr_set$y$norm_list2)
    {dpmsr_set$data$impute$sltmm <<- tmm_normalize(dpmsr_set$data$impute$sl, "SLTMM_Norm", info_columns)
  }
  cat(file = stderr(), "apply_impute function...6", "\n")
  if (11 %in% dpmsr_set$y$norm_list2)
    {dpmsr_set$data$impute$protein <<- protein_normalize(dpmsr_set$data$impute$impute, "Protein_Norm", info_columns)
  }
  cat(file = stderr(), "apply_impute function...end", "\n")
}


#--------------------------------------------------------------------------------
impute_parallel <- function(norm_type) {
  if (norm_type == 99) {data_raw_impute <- impute_only(dpmsr_set$data$normalized$impute, "impute")}
    else if (norm_type == 98) {data_raw_impute_bottom <- impute_bottom(dpmsr_set$data$normalized$impute, "impute_test")}
    else if (norm_type == 1) {data_sl_impute <- impute_only(dpmsr_set$data$normalized$sl, "sl")}
    else if (norm_type == 4) {data_quantile_impute <- impute_only(dpmsr_set$data$normalized$quantile, "quantile")}
    else if (norm_type == 5) {data_lr_impute <- impute_only(dpmsr_set$data$normalized$lr, "lr")}
    else if (norm_type == 6) {data_loess_impute <- impute_only(dpmsr_set$data$normalized$loess, "loess")}
    else if (norm_type == 7) {data_vsn_impute <- impute_only(dpmsr_set$data$normalized$vsn, "vsn")}
    else if (norm_type == 8) {data_ti_impute <- impute_only(dpmsr_set$data$normalized$ti, "ti")}
    else if (norm_type == 9) {data_mi_impute <- impute_only(dpmsr_set$data$normalized$mi, "mi")}
    else if (norm_type == 10) {data_ai_impute <- impute_only(dpmsr_set$data$normalized$ai, "ai")}
    else if (norm_type == 13) {data_directlfq_impute <- impute_only(dpmsr_set$data$normalized$directlfq, "directlfq")}
}

#--------------------------------------------------------------------------------------
check_impute_parallel <- function(norm_list){
  for (norm_name in names(norm_list)) {
    if (is.null(dpmsr_set$data$impute[[norm_name]]    )) {
      cat(file = stderr(), str_c("check parallel, apply_impute function...",norm_name), "\n")
      dpmsr_set$data$impute[[norm_name]] <<- impute_only(dpmsr_set$data$normalized[[norm_name]], norm_name  )
    }else{
      cat(file = stderr(), str_c("check parallel, found --> ",norm_name), "\n")
    }
  }
}

#--------------------------------------------------------------------------------------
check_impute_protein_directlfq <- function(){
  
  if (is.null(dpmsr_set$data$directlfq$directlfq_protein)) {
    cat(file = stderr(), str_c("DirectLFQ not run, no protein level impute necesary"), "\n")
  }else{
    df <- impute_only(dpmsr_set$data$directlfq$directlfq_protein, "directlfq_protein")
    dpmsr_set$data$directlfq$directlfq_protein_impute <<- df
  }
}

#--------------------------------------------------------------------------------------
# data_raw_impute <- impute_only(dpmsr_set$data$normalized$impute, "impute")
# norm_name  <- "impute"
# data_out <- dpmsr_set$data$normalized$impute
#data_in <- dpmsr_set$data$normalized$sl
#norm_name <-

#--------------------------------------------------------------------------------
impute_only <-  function(data_in, norm_name){
  cat(file = stderr(), str_c("impute_only... ", norm_name), "\n")
  cat(file = stderr(), str_c("impute method... ", dpmsr_set$x$impute_method), "\n")
  
  #data_in <- dpmsr_set$data$normalized$sl
  
  info_columns <- ncol(data_in) - dpmsr_set$y$sample_number
  distribution_data <- data_in
  annotation_data <- data_in[1:info_columns]
  data_out <- data_in[(info_columns + 1):ncol(data_in)]
  
  if (dpmsr_set$x$impute_method == "Duke" ||
      dpmsr_set$x$impute_method == "KNN" || dpmsr_set$x$impute_method == "LocalLeastSquares") {
    data_out <- impute_multi(data_out, distribution_data, info_columns)
  } else if (dpmsr_set$x$impute_method == "Floor") {
    data_out[is.na(data_out)] <- dpmsr_set$x$area_floor
  } else if (dpmsr_set$x$impute_method == "Average/Group") {
    data_out <- impute_average_group(data_out)
  } else if (dpmsr_set$x$impute_method == "Minimium") {
    data_out <- impute_minimum(data_out)
  } else if (dpmsr_set$x$impute_method == "MLE") {
    data_out <- impute_mle(data_out)  
  } else if (dpmsr_set$x$impute_method == "BottomX") {
    data_out <- impute_bottomx(data_out, distribution_data, info_columns)  
  } else if (dpmsr_set$x$impute_method == "Average/Global") {
    data_out <- impute_average_global(data_out)  
  } else {
    data_out[is.na(data_out)] <- 0.0}
  data_out <- data.frame(lapply(data_out, as.numeric))
  data_out <- cbind(annotation_data, data_out)
  return(data_out)
}

#--------------------------------------------------------------------------------
# forcing option of bottom x without misalignment filter, allows check for signfig because of misaligned filter
#--------------------------------------------------------------------------------
impute_bottom <-  function(data_in, norm_name){
  cat(file = stderr(), str_c("impute_bottom... ", norm_name), "\n")
  
  info_columns <- ncol(data_in) - dpmsr_set$y$sample_number
  distribution_data <- data_in
  annotation_data <- data_in[1:info_columns]
  data_out <- data_in[(info_columns + 1):ncol(data_in)]
  data_out <- impute_bottomx(data_out, distribution_data, info_columns) 
  data_out <- data.frame(lapply(data_out, as.numeric))
  data_out <- cbind(annotation_data, data_out)
  return(data_out)
}


#--------------------------------------------------------------------------------
# imputation of missing data
impute_multi <- function(data_in, distribution_in, info_columns){
  #distribution_in <- distribution_data
  #data_in <- data_out
  
  #Use all data for distribution or only ptm
  if (as.logical(dpmsr_set$x$peptide_ptm_impute)) {
    if ("Modifications" %in% colnames(distribution_in)) {
      distribution_data <- distribution_in[grep(dpmsr_set$x$peptide_impute_grep, distribution_in$Modifications),]
    } else {
      distribution_data <- distribution_in[grep(dpmsr_set$x$peptide_impute_grep, distribution_in$Sequence),]
    } 
  }else{
    distribution_data <- distribution_in
  }
  
  distribution_data <- distribution_data[(info_columns + 1):ncol(distribution_data)] 
  distribution_data <- log(distribution_data,2)
  #distribution_data[is.na(distribution_data)] <- 0.0
  
  # reset rand_count every time it starts to impute a new normalization
  rand_count <- 1
  
  data_in <- log(data_in,2)
  #data_in[is.na(data_in)] <- 0.0

  for (i in 1:dpmsr_set$y$group_number) {
    # calculate stats for each sample group
    #assign(dpmsr_set$y$sample_groups$Group[i], data.frame(data_in[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])]))
    #df <- get(dpmsr_set$y$sample_groups$Group[i])
    df <- data.frame(data_in[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])])
    
    # adding if statement incase there are non quant groups (1 sample)
    if (ncol(df) > 1) {
      
      df$sum <- rowSums(df, na.rm = TRUE)
      df$rep <- dpmsr_set$y$sample_groups$Count[i]
      #df$min <- dpmsr_set$y$sample_groups$Count[i]/2
      df$max_missing <- dpmsr_set$y$sample_groups$Count[i]*((100 - as.numeric(dpmsr_set$x$missing_cutoff))/100)
      df$max_misaligned <- dpmsr_set$y$sample_groups$Count[i]*(as.numeric(dpmsr_set$x$misaligned_cutoff)/100)
      df$missings <- rowSums(is.na(df[1:dpmsr_set$y$sample_groups$Count[i]]))
      df$average <- apply(df[1:dpmsr_set$y$sample_groups$Count[i]], 1, FUN = function(x) {mean(x, na.rm = TRUE)})
      
      #separate calc for distribution data - for low ptm sets where impute cold be skewed if keep the high abundance data
      df2 <- data.frame(distribution_data[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])])
      #df2$sum <- rowSums(df2, na.rm = TRUE)
      #df2$rep <- dpmsr_set$y$sample_groups$Count[i]
      #df2$min <- dpmsr_set$y$sample_groups$Count[i]/2
      #df2$missings <- rowSums(is.na(df2[1:dpmsr_set$y$sample_groups$Count[i]]))
      df2$average <- apply(df2[1:dpmsr_set$y$sample_groups$Count[i]], 1, FUN = function(x) {mean(x, na.rm = TRUE)})
      df2$sd <- apply(df2[1:dpmsr_set$y$sample_groups$Count[i]], 1, FUN = function(x) {sd(x, na.rm = TRUE)})
      df2$bin <- ntile(df2$average, 20)  
      #sd_info <- subset(df2, missings ==0) %>% group_by(bin) %>% summarize(min = min(average), max = max(average), sd = mean(sd))
      sd_info <- subset(df2, !is.na(sd)) %>% group_by(bin) %>% summarize(min = min(average), max = max(average), sd = mean(sd))
      for (x in 1:(nrow(sd_info) - 1)) {sd_info$max[x] <- sd_info$min[x + 1]}
      sd_info$max[nrow(sd_info)] <- 100
      sd_info$min2 <- sd_info$min
      sd_info$min2[1] <- 0
      sd_info <- sd_info[-21,]
      
      
      # if the number of missing values <= minimum then will impute based on normal dist of measured values
      if (dpmsr_set$x$impute_method == "Duke") {  
        find_rows <- which(df$missings > 0 & df$missings <= df$max_missing)
          for (j in find_rows) {
            findsd <- sd_info %>% filter(df$average[j] >= min2, df$average[j] <= max)
            for (k in 1:dpmsr_set$y$sample_groups$Count[i]) {
              if (is.na(df[j,k])) {
                #nf <-  rnorm(1, 0, 1)
                #nf <- mean(runif(4, min=-1, max=1))
                df[j,k] = df$average[j] + (dpmsr_set$y$rand_impute[rand_count] * findsd$sd[1])
                rand_count <- rand_count + 1
              }
            }
          }
        }
      
      
      # if number of missing greater than minimum and measured value is above intensity cuttoff then remove measured value
      df$missings <- rowSums(is.na(df[1:dpmsr_set$y$sample_groups$Count[i]]))  #recalc if Duke filled, filters may overlap
      if (as.logical(dpmsr_set$x$duke_misaligned)) {
        cat(file = stderr(), str_c("finding and removing misaligned values"), "\n")
        find_rows <- which(df$missings > df$max_misaligned  & df$average >= log(dpmsr_set$x$int_cutoff,2) )
        for (j in find_rows) {
          for (k in 1:dpmsr_set$y$sample_groups$Count[i]) {
                df[j,k] <- NA
           }
          }
      }
    } # end of if ncol(df)>1
    
    #save group dataframe
    assign(dpmsr_set$y$sample_groups$Group[i], df[1:dpmsr_set$y$sample_groups$Count[i]])
    gc()
  }
  

  #get first group
  df3 <- get(dpmsr_set$y$sample_groups$Group[1])
  #add remaining groups
  for (i in 2:dpmsr_set$y$group_number)  {
    df3 <- cbind(df3, get(dpmsr_set$y$sample_groups$Group[i]))
    }
  
  if (dpmsr_set$x$impute_method == "LocalLeastSquares") {df3 <- impute_lls(df3)}
  if (dpmsr_set$x$impute_method == "KNN") {df3 <- impute_knn(df3)}  
  
  
  df3 <- data.frame(2^df3)

  if (dpmsr_set$x$impute_method == "Duke") {
    df3 <- impute_bottomx(df3, distribution_in, info_columns)
    }
  
  return(df3)
}


#--------------------------------------------------------------------------------
# imputation of missing data
impute_bottomx <- function(data_in, distribution_data, info_columns){
  cat(file = stderr(), "apply_bottomx function...", "\n")
  
  # reset rand_count every time it starts to impute a new normalization
  rand_count <- 1
  
  #Use all data for distribution or only ptm
  if (as.logical(dpmsr_set$x$peptide_ptm_impute)) {
    if ("Modifications" %in% colnames(distribution_data)) {
      distribution_data <- distribution_data[grep(dpmsr_set$x$peptide_impute_grep, distribution_data$Modifications),]
    } else {
      distribution_data <- distribution_data[grep(dpmsr_set$x$peptide_impute_grep, distribution_data$Sequence),]
    } 
  }
  
  distribution_data <- distribution_data[(info_columns + 1):ncol(distribution_data)] 
  distribution_data <- log(distribution_data,2)
  
  
  #calc 100 bins for Bottom X%
  data_dist <- as.vector(t(distribution_data))
  data_dist <- data_dist[!is.na(data_dist)]
  data_dist <- data_dist[data_dist > 0]
  data_dist <- data.frame(data_dist)
  data_dist$bin <- ntile(data_dist, 100)  
  bottomx_min <- min(data_dist[data_dist$bin == 1,]$data_dist)
  bottomx_max <- max(data_dist[data_dist$bin == as.numeric(dpmsr_set$x$bottom_x),]$data_dist)
  
  data_in <- log(data_in,2)
  test <- apply(is.na(data_in), 2, which)
  test <- test[lapply(test,length) > 0]

  for (n in names(test)) {
    for (l in (test[[n]])) {
      #data_in[[n]][l] <- mean(runif(4, bottomx_min, bottomx_max))
      # uses stored random numbers from -1 to 1
      data_in[[n]][l] <- bottomx_min + (abs(as.numeric(dpmsr_set$y$rand_impute[rand_count])) * (bottomx_max - bottomx_min))
      rand_count <- rand_count + 1
  }
  }
  
  data_in <- data.frame(2^data_in)
  return(data_in)
}

#--------------------------------------------------------------------------------
# Local least squares imputation (lls)
impute_lls <- function(data_in){
  column_names <- colnames(data_in)
  data_in[data_in==0] <- NA
  tdata_in <-as.data.frame(t(data_in))
  tdata_in_impute <- llsImpute(tdata_in, k=150, allVariables = TRUE)
  data_out <- tdata_in_impute@completeObs
  data_out <- as.data.frame(t(data_out))
  colnames(data_out) <- column_names
  return(data_out)
}

#--------------------------------------------------------------------------------
# Local least squares imputation (lls)
impute_knn <- function(data_in){
  require(impute)
  knndata <- log2(data_in)
  knndata[knndata==-Inf] = NA
  knndata[knndata==0] <- NA
  knndata <- data.matrix(knndata)
  knnimpute <- impute.knn(knndata, k=10)
  data_out <- knnimpute$data
  return(data_out)
}


#--------------------------------------------------------------------------------
impute_average_group <- function(data_in){
  data_in[data_in==0] <- NA
  for(i in 1:dpmsr_set$y$group_number) 
  {
    assign(dpmsr_set$y$sample_groups$Group[i], data.frame(data_in[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])]))
    df <- get(dpmsr_set$y$sample_groups$Group[i])
    df$sum <- rowSums(df, na.rm=TRUE)
    df$rep <- dpmsr_set$y$sample_groups$Count[i]  
    df$min <- dpmsr_set$y$sample_groups$Count[i]/2
    df$missings <- rowSums(is.na (df[1:dpmsr_set$y$sample_groups$Count[i]]))
    df$average <- df$sum / (df$rep - df$missings)
    
    for (j in 1:nrow(df)){
      for (k in 1:dpmsr_set$y$sample_groups$Count[i]){
        if (is.na (df[j,k]) ) { df[j,k] <- df$average[j]}
      }
    }
    assign(dpmsr_set$y$sample_groups$Group[i], df[1:dpmsr_set$y$sample_groups$Count[i]])
  }
  data_out <- get(dpmsr_set$y$sample_groups$Group[1])
  for(i in 2:dpmsr_set$y$group_number)  {data_out <- cbind(data_out, get(dpmsr_set$y$sample_groups$Group[i]))}
  return(data_out)
}

#--------------------------------------------------------------------------------
impute_average_global <- function(data_in){
  data_in[data_in==0] <- NA
  df <- data_in
  df$average <- apply(df, 1, FUN = function(x) {mean(x[x > 0], na.rm=TRUE)})
    for (j in 1:nrow(df)){
      for (k in 1:dpmsr_set$y$sample_number){
        if (is.na (df[j,k]) ) { df[j,k] <- df$average[j]}
      }
    }

  data_out <- df[1:dpmsr_set$y$sample_number]
  return(data_out)
}



#--------------------------------------------------------------------------------

impute_minimum <- function(df){
  df[df==0] <- NA
  df$minimum <- apply(df, 1, FUN = function(x) {min(x[x > 0], na.rm=TRUE)})
  for (j in 1:nrow(df)){
    for (k in 1:dpmsr_set$y$sample_number){
      if (is.na(df[j,k])) {df[j,k] = df$minimum[j]}
    }
  }
  return(df[1:dpmsr_set$y$sample_number])
}

#--------------------------------------------------------------------------------
# Imputing missing values using the EM algorithm proposed in section 5.4.1 of Schafer (1997).
impute_mle <- function(df){
  require(imp4p)
  df[df==0] <- NA
  df <- log2(df)
  df_mle <- impute.mle(df, dpmsr_set$y$group_factor) 
  df_mle <- data.frame(df_mle)
  return(df_mle)
}



#--------------------------------------------------------------------------------
adjust_intensity_cutoff <- function(session, input, output){
  cat(file=stderr(), "Calclulate new intensity cutoff", "\n")
  if(as.numeric(dpmsr_set$x$int_cutoff_sd) < 0) {
    dpmsr_set$x$int_cutoff <<- dpmsr_set$y$intensity_mean + (as.numeric(dpmsr_set$x$int_cutoff_sd) * dpmsr_set$y$intensity_sd )
    cat(file=stderr(), str_c("std < 0    ",  dpmsr_set$x$int_cutoff ), "\n")
  }else{
    dpmsr_set$x$int_cutoff <<- dpmsr_set$y$intensity_mean + (dpmsr_set$x$int_cutoff_sd * dpmsr_set$y$intensity_sd )
    cat(file=stderr(), str_c("std > 0    ",  dpmsr_set$x$int_cutoff ), "\n")
  }
  dpmsr_set$x$int_cutoff <<- trunc(2^dpmsr_set$x$int_cutoff,0)
  output$text_i2 <- renderText(str_c("Intensity cutoff:  ", as.character(dpmsr_set$x$int_cutoff)))
  
}



