run_wiki <- function(input, output, data_in){
  
  if(input$pval_filter_wiki != 0 & input$fc_filter_wiki > 0) {
    go_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_wiki),"_Pval")] <= as.numeric(input$pval_filter) &
                          (data_in[ ,str_c(as.character(input$select_data_comp_wiki),"_FC")] >= as.numeric(input$fc_filter) ))
    atest <<- go_df$Accession
  }
   
  if(input$pval_filter_wiki != 0 & input$fc_filter_wiki < 0) {
    go_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_wiki),"_Pval")] <= as.numeric(input$pval_filter) &
                     (data_in[ ,str_c(as.character(input$select_data_comp_wiki),"_FC")] <= as.numeric(input$fc_filter) ))
    atest <- go_df$Accession    
  } 

  
  test.df <- bitr(atest, fromType = "UNIPROT",
                  toType = c("ENTREZID", "ENSEMBL", "SYMBOL","GENENAME","PATH", "ONTOLOGY"),
                  OrgDb = dpmsr_set$pathway$tax_db)
  
  universe.df <- bitr(data_in$Accession, fromType = "UNIPROT",
                   toType = c("ENTREZID", "ENSEMBL", "SYMBOL","GENENAME","PATH", "ONTOLOGY"),
                   OrgDb = dpmsr_set$pathway$tax_db)
  
 # wp2gene <<- gmt_get(tax_choice)
  
  #wp2gene <- read.gmt(wpgmtfile)
  #dpmsr_set$pathway$wp2gene <<-   dpmsr_set$pathway$wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  wpid2gene <-   dpmsr_set$pathway$wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name <-   dpmsr_set$pathway$wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  
  ewp <- enricher(gene=test.df$ENTREZID, universe=universe.df$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  ewp <- setReadable(ewp, dpmsr_set$pathway$tax_db, keyType = "ENTREZID")
  
  Simple_Excel(ewp@result, str_c(dpmsr_set$file$string, "Wiki_", input$select_data_comp_wiki, "_", 
                          input$select_ont_wiki,input$select_level_wiki, ".xlsx", collapse = " "))
  #ggo_result <<- ggo@result[1:5]
  #ggo_result <- ggo_result[order(-ggo_result$Count),]
  return(ewp@result)
}



run_profile <- function(input, output, data_in){
  
  if(input$pval_filter_profile != 0 & input$fc_filter_profile > 0) {
    go_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_profile),"_Pval")] <= as.numeric(input$pval_filter) &
                      (data_in[ ,str_c(as.character(input$select_data_comp_profile),"_FC")] >= as.numeric(input$fc_filter) ))
    atest <<- go_df$Accession
  }
  
  if(input$pval_filter_profile != 0 & input$fc_filter_profile < 0) {
    go_df <- subset(data_in, data_in[ ,str_c(as.character(input$select_data_comp_profile),"_Pval")] <= as.numeric(input$pval_filter) &
                      (data_in[ ,str_c(as.character(input$select_data_comp_profile),"_FC")] <= as.numeric(input$fc_filter) ))
    atest <- go_df$Accession    
  } 
  
  
  test.df <- bitr(atest, fromType = "UNIPROT",
                  toType = c("ENTREZID", "ENSEMBL", "SYMBOL","GENENAME","PATH", "ONTOLOGY"),
                  OrgDb = dpmsr_set$pathway$tax_db)
  
  ggo <- groupGO(gene     = test.df$ENTREZID,
                 OrgDb    = dpmsr_set$pathway$tax_db,
                 ont      = input$select_ont_profile,
                 level    = as.numeric(input$select_level_profile),
                 readable = TRUE) 
  
  Simple_Excel(ggo, str_c(dpmsr_set$file$string, "GO_Profile_", input$select_data_comp_profile, "_", 
                          input$select_ont_profile,input$select_level_profile, ".xlsx", collapse = " "))
  #ggo_result <<- ggo@result[1:5]
  #ggo_result <- ggo_result[order(-ggo_result$Count),]
  return(ggo)
}