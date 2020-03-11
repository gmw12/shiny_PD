setup_string <- function(input, output, data_in){
  
  dpmsr_set$string$string_db <<- NULL
  
  dpmsr_set$string$string_db <<- STRINGdb$new( version="10", species=9606,
                                               score_threshold=0, input_directory="")
  
  for (i in 1:dpmsr_set$x$comp_number){
    comp_name <- dpmsr_set$y$comp_groups$comp_name[i]
    pval_col <- str_c(comp_name, "_Pval")
    fc2_col <- str_c(comp_name, "_FC2")
    df <- data.frame(cbind(data_in[pval_col], data_in[fc2_col], data_in$Accession), stringsAsFactors = FALSE)
    names(df) <- c("pvalue", "logFC", "Uniprot")
    df$logFC <- as.numeric(df$logFC)
    df$pvalue <- as.numeric(df$pvalue)
    df$logFC <- log(df$logFC, 2)
    assign(str_c("dpmsr_set$string$comp", i), df, envir=.GlobalEnv)
    dpmsr_set$string[[comp_name]]  <<- dpmsr_set$string$string_db$map(df, "Uniprot", removeUnmappedRows = TRUE )
  }
  
  backgroundV <- dpmsr_set$string[[dpmsr_set$y$comp_groups$comp_name[1]]]$STRING_id 
  dpmsr_set$string$string_db$set_background(backgroundV)
}  
  
#-------------------------------------------------------------------
  
run_string <- function(input, output){ 
  input_fc <- log(input$fc_filter_string, 2)
  input_pval <- input$pval_filter_string
  input_comp <- input$select_data_comp_string
  
  df <- dpmsr_set$string[[input_comp]]
  df <- subset(df, pvalue < input_pval)
  df <- subset(df, logFC > input_fc)
  df <- df[order(-df$logFC),]
  
  if (nrow(df) > input$protein_number){
    hits <- df$STRING_id[1:input$protein_number]
  }else{
    hits <- df$STRING_id
  }
  
  dpmsr_set$string$string_db$plot_network(hits)
  
  string_file_name <- str_c(dpmsr_set$file$string, input_comp, ".png")
  
  dpmsr_set$string$string_db$get_png(hits, required_score=NULL, network_flavor="evidence", file=string_file_name, payload_id=NULL)
  
  linkthis <- dpmsr_set$string$string_db$get_link(hits, required_score=NULL, network_flavor="evidence", payload_id=NULL)
  
  return(list("string_file_name" = string_file_name, "linkthis"=linkthis))
}


#--------------------------------------------------------------------

run_string_enrich <- function(input, output, data_in){

  input_fc <- log(input$fc_filter_string_enrich, 2)
  input_pval <- input$pval_filter_string_enrich
  input_comp <- input$select_data_comp_string_enrich
  
  df <- dpmsr_set$string[[input_comp]]
  checkdf<<- df
  df <- subset(df, pvalue < input_pval)
  df <- subset(df, logFC > input_fc)
  df <- df[order(-df$logFC),]
  
  hits <- df$STRING_id
  
  enrichment <- dpmsr_set$string$string_db$get_enrichment(hits, category = input$select_string_enrich, methodMT = input$select_methodMT, iea = TRUE )
  
  Simple_Excel(enrichment, str_c(dpmsr_set$file$string, input_comp, "_", input$select_string_enrich, ".xlsx", collapse = " "))
  
return(enrichment)


}


#enrichment <- dpmsr_set$string$string_db$get_enrichment(checkdf, category = "Process", methodMT = "fdr", iea = TRUE )
