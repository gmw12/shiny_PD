run_go_analysis <- function(input, output, data_in){

    if(input$pval_filter_go != 0 & input$fc_filter_go > 0) {
      go_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_go),"_Pval")] <= as.numeric(input$pval_filter) &
                        (data_in[ ,str_c(as.character(input$select_data_comp_go),"_FC")] >= as.numeric(input$fc_filter) ))
      atest <- go_df$Accession
    }
    
    if(input$pval_filter_go != 0 & input$fc_filter_go < 0) {
      go_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_go),"_Pval")] <= as.numeric(input$pval_filter) &
                        (data_in[ ,str_c(as.character(input$select_data_comp_go),"_FC")] <= as.numeric(input$fc_filter) ))
      atest <- go_df$Accession    
    } 
    
    selection <- atest
    background <- data_in$Accession
    #-----------------------------


    topgodata <- ViSEAGO::create_topGOdata(
      geneSel=selection,
      allGenes=background,
      gene2GO=dpmsr_set$pathway$myGENE2GO, 
      ont=input$select_ont_go,
      nodeSize=5
    )

    resultFisher<-topGO::runTest(
      topgodata,
      algorithm ="classic",
      statistic = "fisher"
    )
    
    resultKS<-topGO::runTest(
      topgodata,
      algorithm ="classic",
      statistic = "ks"
    )
    
    resultKS.elim<-topGO::runTest(
      topgodata,
      algorithm ="elim",
      statistic = "ks"
    )
    
    allRes <- GenTable(topgodata, classicFisher = resultFisher,
                       classicKS = resultKS, elimKS = resultKS.elim,
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10) 
    
  Simple_Excel(allRes, str_c(dpmsr_set$file$string, "GO_", input$select_data_comp_go, "_", 
                          input$select_ont_go, ".xlsx", collapse = " "))
  
  return(allRes)
}


  #----create go lookup data for volcano plots ---------------------
setup_go_volcano <- function(input, output, data_in){  
  
  volcano_df <- cbind(data_in$Accession,  
                 data_in[ ,str_c(as.character(input$select_data_comp_go),"_Pval")], 
                 data_in[ ,str_c(as.character(input$select_data_comp_go),"_FC2")]) 
  volcano_df <- data.frame(volcano_df, stringsAsFactors = FALSE)
  names(volcano_df) <- c("Accession", "pvalue", "foldchange")
  selection <- as.character(volcano_df$Accession)
  
  
  myGENE2GO_df <- dpmsr_set$pathway$myGENE2GO@BP
  
  myGENE2GO_lookup <- data.frame(names(myGENE2GO_df),stringsAsFactors = FALSE)
  myGENE2GO_lookup$Go <- 1
  names(myGENE2GO_lookup)[1] <- "Accession"
  
  for(i in 1:nrow(myGENE2GO_lookup)){
    test <- myGENE2GO_df[[myGENE2GO_lookup$Accession[i]]]
    test_str <- paste(test, collapse = ", ")
    myGENE2GO_lookup$Go[i] <- test_str
  }
  
  #Simple_Excel(myGENE2GO_lookup, "myGENE2GO_lookup.xlsx")
  dpmsr_set$data$pathway$mergedf <<- merge(x=volcano_df, y=myGENE2GO_lookup, by.x="Accession", by.y="Accession")
  #Simple_Excel(dpmsr_set$data$pathway$mergedf, "merged_df.xlsx")
  
}


#-------------------------------------------------------------------------------------------
create_go_volcano <- function(input, output, session){
  
  sub_df <- dpmsr_set$data$pathway$mergedf[grep(as.character(input$go_volcano_id), dpmsr_set$data$pathway$mergedf$Go),]

  sub_df$log_pvalue <- -log(as.numeric(sub_df$pvalue), 10)
  sub_df$log_fc <- log(as.numeric(sub_df$foldchange), 2)
  
  #Simple_Excel(mergedf, "mergedf.xlsx")
  #Simple_Excel(sub_df, "sub_df.xlsx")
  
  return(sub_df)
}


test_e <- function()
{
  
  envir = environment()
  print(envir)
  
}