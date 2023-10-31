run_go_analysis <- function(session, input, output){
  require(topGO)
  
  comp_string <- input$select_data_comp_go
  comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
  
  cat(file = stderr(), str_c("Run go analysis...1" ), "\n")

  data_in <- dpmsr_set$data$stats[[comp_string]]
  
  if (input$go_direction == 'Up') {
    go_df <- subset(data_in, data_in$Stats == "Up" ) 
  }else if (input$go_direction == 'Down') {
    go_df <- subset(data_in, data_in$Stats == "Down" ) 
  }else{
    go_df <- subset(data_in, data_in$Stats == "Up" | data_in$Stats == "Down") 
  }
  

  atest <- go_df$Accession
  selection <- atest
  background <- data_in$Accession
  ont_go <- input$select_ont_go
  
    #selection <<- atest; background <<- data_in$Accession; ont_go <<- input$select_ont_go
    
    
    #-----------------------------

    cat(file=stderr(), str_c("Run go analysis...2" ), "\n")
    
    ont_go <- input$select_ont_go
    
    fun_topgodata <- function(selection, background, myGENE2GO, ont_go) {
      topgodata <- ViSEAGO::create_topGOdata(
        geneSel = selection,
        allGenes = background,
        gene2GO = myGENE2GO, 
        ont = ont_go,
        nodeSize = 5
      )
      return(topgodata)
    }
    
    topgodata <- callr::r_bg(fun_topgodata, args = list(selection, background, myGENE2GO, ont_go), supervise = TRUE)
    topgodata$wait()
    topgodata <- topgodata$get_result()
    
    
    cat(file = stderr(), str_c("Run go analysis...3" ), "\n")  
    fun_resultFisher <- function(topgodata) {
      resultFisher <- topGO::runTest(
        topgodata,
        algorithm = "classic",
        statistic = "fisher"
      )
      return(resultFisher)
    }

    
    cat(file = stderr(), str_c("Run go analysis...4" ), "\n")
    fun_resultKS <- function(topgodata){
      resultKS <- topGO::runTest(
        topgodata,
        algorithm = "classic",
        statistic = "ks"
      )
      return(resultKS)
    }
    
    cat(file = stderr(), str_c("Run go analysis...5" ), "\n")
    fun_resultKS.elim <- function(topgodata){
      resultKS.elim <- topGO::runTest(
        topgodata,
        algorithm = "elim",
        statistic = "ks"
      )
      return(resultKS.elim)
    }
     
    resultFisher <- callr::r_bg(fun_resultFisher, args = list(topgodata), supervise = TRUE)
    resultKS <- callr::r_bg(fun_resultKS, args = list(topgodata), supervise = TRUE)
    resultKS.elim <- callr::r_bg(fun_resultKS.elim, args = list(topgodata), supervise = TRUE)
    
    resultFisher$wait()
    resultKS$wait()
    resultKS.elim$wait()
    
    resultFisher <- resultFisher$get_result() 
    resultKS <- resultKS$get_result()   
    resultKS.elim <- resultKS.elim$get_result()  
    
    
    cat(file = stderr(), str_c("Run go analysis...6" ), "\n")
    
    # merge results from topGO
    
    allRes <- GenTable(topgodata, classicFisher = resultFisher,
                       classicKS = resultKS, elimKS = resultKS.elim,
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10) 
    
  Simple_Excel_bg(allRes, "Go", str_c(dpmsr_set$file$string, "GO_", input$select_data_comp_go, "_", 
                          input$select_ont_go, ".xlsx", collapse = " "))
  gc()
  return(allRes)
}


  #----create go lookup data for volcano plots ---------------------
setup_go_volcano <- function(session, input, output){  
  comp_string <- input$select_data_comp_go
  comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
  data_in <- dpmsr_set$data$stats[[comp_string]]
  
  cat(file = stderr(), str_c("Setup go volcano...1" ), "\n")
  volcano_df <- cbind(data_in$Accession,
                      data_in[ ,dpmsr_set$y$stats$groups$pval[comp_number]],
                      data_in[ ,dpmsr_set$y$stats$groups$fc2[comp_number]]) 
  volcano_df <- data.frame(volcano_df, stringsAsFactors = FALSE)
  
  names(volcano_df) <- c("Accession", "pvalue", "foldchange")
  selection <- as.character(volcano_df$Accession)
  
  cat(file = stderr(), str_c("Setup go volcano...2" ), "\n")
  
  myGENE2GO_df <- myGENE2GO@BP
  
  myGENE2GO_lookup <- data.frame(names(myGENE2GO_df),stringsAsFactors = FALSE)
  myGENE2GO_lookup$Go <- 1
  names(myGENE2GO_lookup)[1] <- "Accession"
  
  for (i in 1:nrow(myGENE2GO_lookup)) {
    test <- myGENE2GO_df[[myGENE2GO_lookup$Accession[i]]]
    test_str <- paste(test, collapse = ", ")
    myGENE2GO_lookup$Go[i] <- test_str
  }
  
  #Simple_Excel(myGENE2GO_lookup, "myGENE2GO_lookup.xlsx")
  mergedf <- merge(x = volcano_df, y = myGENE2GO_lookup, by.x = "Accession", by.y = "Accession")
  #Simple_Excel(dpmsr_set$data$pathway$mergedf, "merged_df.xlsx")
  cat(file = stderr(), str_c("Setup go volcano...end" ), "\n")
  mergedf
}



#-------------------------------------------------------------------------------------------
create_go_volcano <- function(session, input, output){
  cat(file=stderr(), str_c("create_go_volcano...1" ), "\n")
  comp_string <- input$select_data_comp_go
  comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
  data_in <- dpmsr_set$data$stats[[comp_string]]
  
  sub_df <- dpmsr_set$data$pathway$mergedf[grep(as.character(input$go_volcano_id), dpmsr_set$data$pathway$mergedf$Go),]
  
  cat(file=stderr(), str_c("create_go_volcano...2" ), "\n")
  sub_df$log_pvalue <- -log(as.numeric(sub_df$pvalue), 10)
  sub_df$log_fc <- log(as.numeric(sub_df$foldchange), 2)

  cat(file=stderr(), str_c("create_go_volcano...3" ), "\n")
  testdf <- data.frame(cbind(data_in$Accession, data_in$Description))
  colnames(testdf) <- c("Accession", "Description")
  volcano_data <- merge(x=sub_df, y=testdf, by.x="Accession", by.y="Accession")
  cat(file=stderr(), str_c("create_go_volcano...end" ), "\n")
  gc()
  return(volcano_data)
}





#-------------------------------------------------------------------------------------------


test_e <- function()
{
  
  envir = environment()
  print(envir)
  
}