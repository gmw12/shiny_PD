
apply_impute <- function(){
  cat(file=stderr(), "apply_impute function...", "\n")
  ncores <- detectCores()
  if (is.na(ncores)) {ncores <- 1}
  norm_list <- dpmsr_set$y$norm_list
  dpmsr_set$data$impute <<- mclapply(norm_list, impute_parallel, mc.cores = ncores)
  #need to complete TMM norms after imputation ()
  
  #check that the parallel processing went through, if not do it again one at a time
  check_impute_parallel(norm_list)
  
  norm_list2 <- dpmsr_set$y$norm_list2
  if (2 %in% dpmsr_set$y$norm_list2)
  {dpmsr_set$data$impute$tmm <<- tmm_normalize(dpmsr_set$data$impute$impute, dpmsr_set$data$impute$impute, "TMM_Norm")
    Simple_Excel(dpmsr_set$data$impute$tmm, str_c(dpmsr_set$file$extra_prefix, "_TMM_", dpmsr_set$x$impute_method, "_impute.xlsx", collapse = " "))
    }
  if (3 %in% dpmsr_set$y$norm_list2)
  {dpmsr_set$data$impute$sltmm <<- tmm_normalize(dpmsr_set$data$impute$sl,dpmsr_set$data$impute$sl, "SLTMM_Norm")
    Simple_Excel(dpmsr_set$data$impute$sltmm, str_c(dpmsr_set$file$extra_prefix, "_SLTMM_", dpmsr_set$x$impute_method, "_impute.xlsx", collapse = " "))
  }
  if (11 %in% dpmsr_set$y$norm_list2)
  {dpmsr_set$data$impute$protein <<- protein_normalize(dpmsr_set$data$impute$impute, "Protein_Norm")
    Simple_Excel(dpmsr_set$data$impute$protein, str_c(dpmsr_set$file$extra_prefix, "_Protein_", dpmsr_set$x$impute_method, "_impute.xlsx", collapse = " "))
    }
}


#--------------------------------------------------------------------------------
impute_parallel <- function(norm_type){
  if(norm_type==99){data_raw_impute <- impute_only(dpmsr_set$data$normalized$impute, "impute")}
    else if(norm_type==1){data_sl_impute <- impute_only(dpmsr_set$data$normalized$sl, "sl")}
    else if(norm_type==4){data_quantile_impute <- impute_only(dpmsr_set$data$normalized$quantile, "quantile")}
    else if(norm_type==5){data_lr_impute <- impute_only(dpmsr_set$data$normalized$lr, "lr")}
    else if(norm_type==6){data_loess_impute <- impute_only(dpmsr_set$data$normalized$loess, "loess")}
    else if(norm_type==7){data_vsn_impute <- impute_only(dpmsr_set$data$normalized$vsn, "vsn")}
    else if(norm_type==8){data_ti_impute <- impute_only(dpmsr_set$data$normalized$ti, "ti")}
    else if(norm_type==9){data_mi_impute <- impute_only(dpmsr_set$data$normalized$mi, "mi")}
    else if(norm_type==10){data_ai_impute <- impute_only(dpmsr_set$data$normalized$ai, "ai")}
}

#--------------------------------------------------------------------------------------
check_impute_parallel <- function(norm_list){
  for(norm_name in names(norm_list)){
    if(is.null(dpmsr_set$data$impute[[norm_name]]    )){
      cat(file=stderr(), str_c("apply_impute function...",norm_name), "\n")
      dpmsr_set$data$impute[[norm_name]]<<-impute_only(dpmsr_set$data$normalized[[norm_name]], norm_name  )
    }
  }
}

#--------------------------------------------------------------------------------
impute_only <-  function(data_out, norm_name){
  distribution_data <- data_out
  annotation_data <- data_out[1:dpmsr_set$y$info_columns]
  data_out <- data_out[(dpmsr_set$y$info_columns+1):ncol(data_out)]
  if (dpmsr_set$x$impute_method == "Duke" ||
      dpmsr_set$x$impute_method == "KNN" || dpmsr_set$x$impute_method =="LocalLeastSquares"){
    data_out <- impute_multi(data_out, distribution_data)
  } else if (dpmsr_set$x$impute_method== "Floor") {
    data_out[is.na(data_out)] <- dpmsr_set$x$area_floor
  } else if (dpmsr_set$x$impute_method == "Average") {
    data_out <- impute_average(data_out)
  } else if (dpmsr_set$x$impute_method == "Minimium") {
    data_out <- impute_minimum(data_out)
  } else if (dpmsr_set$x$impute_method == "MLE") {
    data_out <- impute_mle(data_out)  
  } else if (dpmsr_set$x$impute_method == "BottomX") {
    data_out <- impute_bottomx(data_out, distribution_data)  
  } else {
    data_out[is.na(data_out)] <- 0.0}
  data_out<-data.frame(lapply(data_out, as.numeric))
  data_out<-cbind(annotation_data, data_out)
  Simple_Excel(data_out, str_c(dpmsr_set$file$extra_prefix, "_", norm_name, "_", dpmsr_set$x$impute_method, "_impute.xlsx", collapse = " "))
  return(data_out)
}

#--------------------------------------------------------------------------------
# imputation of missing data
impute_multi <- function(data_in, distribution_in){
  #Use all data for distribution or only ptm
  if (as.logical(dpmsr_set$x$peptide_ptm_impute)){
    distribution_data <- distribution_in[grep(dpmsr_set$x$peptide_grep, distribution_data$Modifications),]
  }else{
    distribution_data <- distribution_in
  }
  distribution_data <- distribution_data[(dpmsr_set$y$info_columns+1):ncol(distribution_data)] 
  distribution_data <- log(distribution_data,2)
  #distribution_data[is.na(distribution_data)] <- 0.0
  
  data_in <- log(data_in,2)
  #data_in[is.na(data_in)] <- 0.0

  for(i in 1:dpmsr_set$y$group_number){
    # calculate stats for each sample group
    #assign(dpmsr_set$y$sample_groups$Group[i], data.frame(data_in[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])]))
    #df <- get(dpmsr_set$y$sample_groups$Group[i])
    df <- data.frame(data_in[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])])
    df$sum <- rowSums(df, na.rm = TRUE)
    df$rep <- dpmsr_set$y$sample_groups$Count[i]
    df$min <- dpmsr_set$y$sample_groups$Count[i]/2
    df$missings <- rowSums(is.na(df[1:dpmsr_set$y$sample_groups$Count[i]]))
    df$average <- apply(df[1:dpmsr_set$y$sample_groups$Count[i]], 1, FUN = function(x) {mean(x, na.rm=TRUE)})
    
    #separate calc for distribution data - for low ptm sets where impute cold be skewed if keep the high abundance data
    df2 <- data.frame(distribution_data[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])])
    #df2$sum <- rowSums(df2, na.rm = TRUE)
    #df2$rep <- dpmsr_set$y$sample_groups$Count[i]
    #df2$min <- dpmsr_set$y$sample_groups$Count[i]/2
    #df2$missings <- rowSums(is.na(df2[1:dpmsr_set$y$sample_groups$Count[i]]))
    df2$average <- apply(df2[1:dpmsr_set$y$sample_groups$Count[i]], 1, FUN = function(x) {mean(x, na.rm=TRUE)})
    df2$sd <- apply(df2[1:dpmsr_set$y$sample_groups$Count[i]], 1, FUN = function(x) {sd(x, na.rm=TRUE)})
    df2$bin <- ntile(df2$average, 20)  
    #sd_info <- subset(df2, missings ==0) %>% group_by(bin) %>% summarize(min = min(average), max = max(average), sd = mean(sd))
    sd_info <- subset(df2, !is.na(sd)) %>% group_by(bin) %>% summarize(min = min(average), max = max(average), sd = mean(sd))
    for (x in 1:19){sd_info$max[x] <- sd_info$min[x+1]}
    sd_info$max[nrow(sd_info)] <- 100
    sd_info$min2 <- sd_info$min
    sd_info$min2[1] <- 0
    sd_info <- sd_info[-21,]
    
    
    # if the number of missing values <= minimum then will impute based on normal dist of measured values
    if (dpmsr_set$x$impute_method == "Duke"){  
      find_rows <- which(df$missings >0 & df$missings<=df$min)
        for (j in find_rows){
          findsd <- sd_info %>% filter(df$average[j] >= min2, df$average[j]<= max)
          for (k in 1:dpmsr_set$y$sample_groups$Count[i]){
            if (is.na(df[j,k])) {
              nf <-  rnorm(1, 0, 1)
              df[j,k] = df$average[j] + (nf*findsd$sd[1])
            }
          }
        }
      }
    

    # if number of missing greater than minimum and measured value is above intensity cuttoff then remove measured value
    if (dpmsr_set$x$missing_50){
      find_rows <- which(df$missings>df$min)
      for (j in find_rows){
        for (k in 1:dpmsr_set$y$sample_groups$Count[i]){
            if (!is.na(df[j,k]) && df[j,k] > log(dpmsr_set$x$int_cutoff,2)) {
              df[j,k] <- NA
            }
         }
        }
      }
   

    #save group dataframe
    assign(dpmsr_set$y$sample_groups$Group[i], df[1:dpmsr_set$y$sample_groups$Count[i]])
  }
  
  #saving original data with imputed data frame for return
  df3 <- get(dpmsr_set$y$sample_groups$Group[1])
  for(i in 2:dpmsr_set$y$group_number)  {df3 <- cbind(df3, get(dpmsr_set$y$sample_groups$Group[i]))}
  
  if (dpmsr_set$x$impute_method == "LocalLeastSquares"){df3 <- impute_lls(df3)}
  if (dpmsr_set$x$impute_method == "KNN"){df3 <- impute_knn(df3)}  
  
  
  df3 <- data.frame(2^df3)

  if (dpmsr_set$x$impute_method == "Duke"){df3 <- impute_bottomx(df3, distribution_in)}
  
  return(df3)
}


#--------------------------------------------------------------------------------
# imputation of missing data
impute_bottomx <- function(data_in, distribution_data){
  cat(file=stderr(), "apply_bottomx function...", "\n")
  #Use all data for distribution or only ptm
  if (as.logical(dpmsr_set$x$peptide_ptm_impute)){
    distribution_data <- distribution_data[grep(dpmsr_set$x$peptide_grep, distribution_data$Modifications),]
  }  
  distribution_data <- distribution_data[(dpmsr_set$y$info_columns+1):ncol(distribution_data)] 
  distribution_data <- log(distribution_data,2)
  
  
  #calc 100 bins for Bottom X%
  data_dist <- as.vector(t(distribution_data))
  data_dist <- data_dist[!is.na(data_dist)]
  data_dist <- data_dist[data_dist>0]
  data_dist <- data.frame(data_dist)
  data_dist$bin <- ntile(data_dist, 100)  
  bottomx_min <- min(data_dist[data_dist$bin==1,]$data_dist)
  bottomx_max <- max(data_dist[data_dist$bin==as.numeric(dpmsr_set$x$bottom_x),]$data_dist)
  
  data_in <- log(data_in,2)
  test <- apply(is.na(data_in), 2, which)
  
  for(n in names(test)){
    for(l in (test[[n]]))
      data_in[[n]][l] <- runif(1, bottomx_min, bottomx_max)
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
#if less than half values exist all values go to area_floor
impute_average <- function(data_in){
  data_in[data_in==0] <- NA
  for(i in 1:dpmsr_set$y$group_number) 
  {
    assign(dpmsr_set$y$sample_groups$Group[i], data.frame(data_in[c(dpmsr_set$y$sample_groups$start[i]:dpmsr_set$y$sample_groups$end[i])]))
    df <- get(dpmsr_set$y$sample_groups$Group[i])
    df$sum <- rowSums(df, na.rm=TRUE)
    df$rep <- dpmsr_set$y$sample_groups$Count[i]  
    df$min <- dpmsr_set$y$sample_groups$Count[i]/2
    df$missings <- rowSums(df[1:dpmsr_set$y$sample_groups$Count[i]] == "0")
    df$missings <- rowSums(is.na (df[1:dpmsr_set$y$sample_groups$Count[i]]))
    df$average <- df$sum / (df$rep - df$missings)
    
    for (j in 1:nrow(data_in)){
      for (k in 1:dpmsr_set$y$sample_groups$Count[i]){
        if (is.na (df[j,k]) && df$missings[j] <= df$min[j]) { df[j,k] <- df$average[j]}
        else if (df$missings[j] > df$min[j]) { df[j,k] = area_floor}
      }
    }
    assign(dpmsr_set$y$sample_groups$Group[i], df[1:dpmsr_set$y$sample_groups$Count[i]])
  }
  data_out <- get(dpmsr_set$y$sample_groups$Group[i])
  for(i in 2:dpmsr_set$y$group_number)  {data_out <- cbind(data_out, get(dpmsr_set$y$sample_groups$Group[i]))}
  return(data_out)
}

#--------------------------------------------------------------------------------
# fix missings - by replicate group if more than half of the values exist, replace missings with average, if less than half values exist all values go to area_floor
impute_minimum <- function(df){
  df[df==0] <- NA
  df$minimum <- apply(df, 1, FUN = function(x) {min(x[x > 0])})
  for (j in 1:nrow(df)){
    for (k in 1:dpmsr_set$y$sample_number){
      if (df[j,k] == "0") {df[j,k] = df$minimum[j]}
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
