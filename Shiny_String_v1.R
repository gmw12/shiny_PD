setup_string <- function(session, input, output){
  
  dpmsr_set$string$string_db <<- NULL
  tax_choice <- input$select_organism
  
  #string_species <- get_STRING_species(version=10)
  
  if (version$major < 4){
    string_version <- "10"
  }else{
      string_version <- "11.0"
    }
  
  
  if (site_user == "dpmsr"){
    string_input_dir <- ""
  }else{
    string_input_dir <- "/data/ShinyData"
  }
  
  
  
  if(input$select_organism=="Human"){
          dpmsr_set$string$string_db <<- STRINGdb$new(version=string_version, species=9606,
                                               score_threshold=0, input_directory=string_input_dir)
  }
  
  if(input$select_organism=="Mouse"){
    dpmsr_set$string$string_db <<- STRINGdb$new( version=string_version, species=10090,
                                                 score_threshold=0, input_directory=string_input_dir)
  } 
  
  for (i in 1:dpmsr_set$x$comp_number){
    comp_name <- dpmsr_set$y$stats$groups$comp_name[i]
    data_in <- dpmsr_set$data$stats[[comp_name]]
    pval_col <- dpmsr_set$y$stats$groups$pval[i]
    fc2_col <- dpmsr_set$y$stats$groups$fc2[i]
    df <- data.frame(cbind(data_in[pval_col], data_in[fc2_col], data_in$Accession), stringsAsFactors = FALSE)
    names(df) <- c("pvalue", "logFC", "Uniprot")
    df$logFC <- as.numeric(df$logFC)
    df$pvalue <- as.numeric(df$pvalue)
    df$logFC <- log(df$logFC, 2)
    assign(str_c("dpmsr_set$string$comp", i), df, envir=.GlobalEnv)
    dpmsr_set$string[[comp_name]]  <<- dpmsr_set$string$string_db$map(df, "Uniprot", removeUnmappedRows = TRUE )
  }
  
  backgroundV <- dpmsr_set$string[[dpmsr_set$y$stats$groups$comp_name[1] ]]$STRING_id 
  dpmsr_set$string$string_db$set_background(backgroundV)
  cat(file=stderr(), "function setup_string complete...", "\n")
}  
  
#-------------------------------------------------------------------
  
run_string <- function(session, input, output){
  cat(file=stderr(), "run string step 1", "\n")
  input_fc_up <- log(input$foldchange_cutoff, 2)
  input_fc_down <- log(1/input$foldchange_cutoff, 2)
  
  input_pval <- input$pvalue_cutoff
  input_comp <- input$select_data_comp_string
  
  cat(file=stderr(), "run string step 2", "\n")
  
  df <- dpmsr_set$string[[input_comp]]
  df <- subset(df, pvalue < input_pval)
  
  cat(file=stderr(), "run string step 3", "\n")
  
  if (input$string_direction == "Up"){
    df <- subset(df, logFC >= input_fc_up)
  }else if (input$string_direction == "Down"){
    df <- subset(df, logFC <= input_fc_down)
  }else {
    df <- subset(df, logFC >= input_fc_up | logFC <= input_fc_down )
  }
  
  df <- df[order(-df$logFC),]

  cat(file=stderr(), "run string step 4", "\n")
  
  if (nrow(df) > input$protein_number){
    hits <- df$STRING_id[1:input$protein_number]
  }else{
    hits <- df$STRING_id
  }
  
  
  dpmsr_set$string$string_db$plot_network(hits)
  
  cat(file=stderr(), "run string step 5", "\n")
  string_file_name <- str_c(dpmsr_set$file$string, input_comp, ".png")
  
  cat(file=stderr(), "run string step 6", "\n")
  dpmsr_set$string$string_db$get_png(hits, required_score=NULL, network_flavor="evidence", file=string_file_name, payload_id=NULL)
  
  cat(file=stderr(), "run string step 7", "\n")
  linkthis <- dpmsr_set$string$string_db$get_link(hits, required_score=NULL, network_flavor="evidence", payload_id=NULL)
  
  return(list("string_file_name" = string_file_name, "linkthis"=linkthis))
}


#--------------------------------------------------------------------

run_string_enrich <- function(input, output, data_in){
  
  input_fc_up <- log(input$foldchange_cutoff, 2)
  input_fc_down <- log(1/input$foldchange_cutoff, 2)
  
  input_pval <- input$pvalue_cutoff
  
  input_comp <- input$select_data_comp_string_enrich
  
  df <- dpmsr_set$string[[input_comp]]
  df <- subset(df, pvalue < input_pval)
  
  if (input$string_enrich_direction == "Up"){
    df <- subset(df, logFC >= input_fc_up)
  }else if (input$string_enrich_direction == "Down"){
    df <- subset(df, logFC <= input_fc_down)
  }else {
    df <- subset(df, logFC >= input_fc_up | logFC <= input_fc_down )
  }

  df <- df[order(-df$logFC),]
  
  hits <- df$STRING_id
  
  enrichment <- dpmsr_set$string$string_db$get_enrichment(hits, category = input$select_string_enrich, methodMT = input$select_methodMT, iea = TRUE )
  
  Simple_Excel(enrichment, str_c(dpmsr_set$file$string, input_comp, "_", input$select_string_enrich, ".xlsx", collapse = " "))
  
return(enrichment)


}


#enrichment <- dpmsr_set$string$string_db$get_enrichment(checkdf, category = "Process", methodMT = "fdr", iea = TRUE )
