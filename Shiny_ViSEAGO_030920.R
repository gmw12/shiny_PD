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

    selection_gene <- bitr(selection, fromType = "UNIPROT",
                     toType = c("ENTREZID", "ENSEMBL", "SYMBOL","GENENAME","PATH", "ONTOLOGY"),
                     OrgDb = dpmsr_set$pathway$tax_db)
    background_gene <- bitr(background, fromType = "UNIPROT",
                            toType = c("ENTREZID", "ENSEMBL", "SYMBOL","GENENAME","PATH", "ONTOLOGY"),
                            OrgDb = dpmsr_set$pathway$tax_db)
    
    # require(ensembldb)
    # require(EnsDb.Hsapiens.v86)
    # edb <- EnsDb.Hsapiens.v86
    # selection_gene <<- ensembldb::select(edb, keys = selection, keytype = "UNIPROTID",
    #                          columns = "ENTREZID")
    # background_gene <<- ensembldb::select(edb, keys = background, keytype = "UNIPROTID",
    #                           columns = "ENTREZID")
    
    
    ###################
    # create topGOdata

    topgodata <- ViSEAGO::create_topGOdata(
      geneSel=selection_gene$ENTREZID,
      allGenes=background_gene$ENTREZID,
      gene2GO=dpmsr_set$pathway$myGENE2GO, 
      ont="BP", #input$select_ont_go,
      nodeSize=5
    )

    test_topgo <- topGO::runTest(
      topgodata,
      algorithm = "classic",  # input$select_algorithm,
      statistic = "fisher"
    )
    
    
    ###################
    # merge results
prayer <- function(x,y) { 
    go_sResults <- ViSEAGO::merge_enrich_terms(
      Input=list(
        comp=c("x","y")
      ), envir = environment()
    )
    return(go_sResults)
}
    
    go_sResults <- prayer(topgodata, test_topgo)  
    
  goresult <- go_sResults@data
  colnames(goresult) <- c("GO.ID", "Term", "Definition", "Gene_Freq", "pvalue", "log10_pvalue", "Signif_Gene", "Signif_Gene_Symbol")  
  
  Simple_Excel(goresult, str_c(dpmsr_set$file$string, "GO_", input$select_data_comp_go, "_", 
                          input$select_ont_go,"_",input$select_algorithm, ".xlsx", collapse = " "))
  
  return(goresult)
}


  #----create go lookup data for volcano plots ---------------------
setup_go_volcano <- function(input, output, data_in){  
  
  volcano_df <- cbind(data_in$Accession,  
                 data_in[ ,str_c(as.character(input$select_data_comp_go),"_Pval")], 
                 data_in[ ,str_c(as.character(input$select_data_comp_go),"_FC2")]) 
  volcano_df <- data.frame(volcano_df, stringsAsFactors = FALSE)
  names(volcano_df) <- c("Accession", "pvalue", "foldchange")
  selection <- as.character(volcano_df$Accession)
  
  testallgene <- AnnotationDbi::select(dpmsr_set$pathway$tax_db, keys=selection, 
                                       columns = c("ENTREZID","GO", "GENENAME"), keytype = "UNIPROT")
  
  
  testallgene <- testallgene[testallgene$ONTOLOGY==input$select_ont_go,]
  golookup <- data.frame(cbind(testallgene$UNIPROT, testallgene$GO))
  goolookup <- golookup %>% distinct()
  goolookup <- na.omit(goolookup)
  
  # v_results_df <- v_results_df[order(v_results_df$condition.pvalue  )]
  # firstgo <- v_results_df$GO.ID[1]
  # goproteins <- goolookup[goolookup$X2==firstgo,]
  
  test1 <- goolookup %>% group_by(X1) %>% mutate(Y1 = paste0(X2, collapse = ",")) 
  test1$X2 <- NULL
  testtme <- unique(test1)
  testtme$X1 <-as.character(testtme$X1)
  
  dpmsr_set$data$pathway$mergedf <<- merge(x=volcano_df, y=testtme, by.x="Accession", by.y="X1")
  Simple_Excel(dpmsr_set$data$pathway$mergedf, "merged_df.xlsx")
  
}


#-------------------------------------------------------------------------------------------
create_go_volcano <- function(input, output, session){
  
  sub_df <- dpmsr_set$data$pathway$mergedf[grep(as.character(input$go_volcano_id), dpmsr_set$data$pathway$mergedf$Y1),]

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