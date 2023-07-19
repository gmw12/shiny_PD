#----------------------------------------------------------------------------------------
# create column to flag and delete rows with no data or only one data value
#require row to contain at least one group that has at least 2 values
filter_data <- function(session, input, output, df){
  #test_df <<- df
  
  cat(file=stderr(), "start filter_data...", "\n")
  start <- Sys.time()
  excel_list <- list()
  
  # Step 1 remove peptides/proteins below minumum count requirement overall 
  cat(file=stderr(), "step 1, remove below minimum...", "\n")
  info_columns <- ncol(df) - dpmsr_set$y$sample_number
  df[df==0] <- NA
  df$measured_count <- apply(df[, (info_columns+1):ncol(df)], 1, function(x) sum(!is.na(x)))
  df <- subset(df, df$measured_count >= input$minimum_measured_all)
  df <- df[,1:dpmsr_set$y$total_columns]
  excel_list[["Raw_Filter"]] <- df

  # Step 2 remove peptides/proteins below minumum count requirement in groups 
  cat(file=stderr(), "step 2, remove below minimum group requirement...", "\n")
  if(as.logical(dpmsr_set$x$require_x)){
    for(i in 1:dpmsr_set$y$group_number) {
      df$na_count <- apply(df[, (dpmsr_set$y$sample_groups$start[i]+dpmsr_set$y$info_columns):
                                  (dpmsr_set$y$sample_groups$end[i]+dpmsr_set$y$info_columns)], 1, function(x) sum(!is.na(x)))
      df$na_count <- df$na_count / dpmsr_set$y$sample_groups$Count[i]
      colnames(df)[colnames(df)=="na_count"] <- dpmsr_set$y$sample_groups$Group[i]
    }
    df$max <- apply(df[, (dpmsr_set$y$total_columns+1):(ncol(df))], 1, function(x) max(x))
    df <- subset(df, df$max >= dpmsr_set$x$require_x_cutoff )
    df <- df[1:dpmsr_set$y$total_columns]
    
    excel_list[["Raw_Filter_Req2"]] <- df
  }
  
  #Step 3 optional filter by cv of specific sample group
  cat(file=stderr(), "step 3, cv minimum...", "\n")
  if(as.logical(dpmsr_set$x$filter_cv)){
    df<-filter_by_cv(df, info_columns)
    excel_list[["Raw_Filter_CV"]] <- df
  }
  
  cat(file=stderr(), str_c("filter_data completed in ", Sys.time()-start), "\n")
  
  # write excel in parallel in background
  x <- foreach(i=1:length(excel_list)) %dopar% {
    excel_name <- names(excel_list[i])
    Simple_Excel_bg(excel_list[[excel_name]], excel_name, str_c(dpmsr_set$file$extra_prefix, "_", excel_name, ".xlsx", collapse = " "))
  }
  
  gc()
  return(df)
}

#----------------------------------------------------------------------------------------

#filter by CV of a specific group (QCPool)
filter_by_cv <- function(df, info_columns){
  cat(file=stderr(), str_c("filter by cv"), "\n")
  start <- info_columns + dpmsr_set$y$sample_groups$start[dpmsr_set$y$sample_groups$Group == dpmsr_set$x$filter_cv_group]
  end <- info_columns + dpmsr_set$y$sample_groups$end[dpmsr_set$y$sample_groups$Group == dpmsr_set$x$filter_cv_group]
  df$filterCV <- percentCV_gw(df[start:end])
  df <- subset(df, df$filterCV < dpmsr_set$x$filter_cv_cutoff)
  df <- df[-ncol(df)]
  return(df)
}