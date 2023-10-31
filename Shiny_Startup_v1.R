
app_startup <- function(session, input, output){
  
  # check if dpmsr_set file is loaded
  dpmsr_present <<- reactiveValues()
  dpmsr_present$test <<- exists("dpmsr_set")
  
  cat(file = stderr(), "Shiny Server started ...2", "\n")
  
  shinyjs::disable("action_load_design")
  shinyjs::disable("action_load_dpmsr_set")
  shinyjs::disable("load_data")
  
  cat(file = stderr(), "Shiny Server started ...3", "\n")
  
  output$app_version_text <- renderText({str_c('app version = ', app_version)})
  output$app_version_text2 <- renderText({str_c('app version = ', app_version)})
  
  output$text_n1 <- renderText("Check all Normalization Strategies")
  output$text_f1 <- renderText("Normalization Filters")
  output$text_i1 <- renderText("Select Imputation Method")
  output$text_r1 <- renderText("Select Protein Rollup Method")
  
  output$text_impute_ptm <- renderText("Impute Distributions using Impute PTM grep (Load Data)")
  
  cat(file = stderr(), "Shiny Server started ...4", "\n")
  
  shinyFileChoose(input, 'design_file', session = session, roots = volumes, filetypes = c('', 'xlsx'))

  shinyFileChoose(input,'dpmsr_set_file', roots = volumes, session = session, 
                  filetypes = c('','dpmsr_set'), defaultPath = '', defaultRoot = 'wd')
  
  shinyFileChoose(input,'motif_fasta', roots = volumes, session = session, 
                  filetypes = c('','fasta'), defaultPath = '', defaultRoot = 'wd')
  
  shinyFileChoose(input,'test', roots = volumes, session = session, 
                  filetypes = c('' , 'xlsx', 'jpg', 'png'), defaultPath = '', defaultRoot = 'wd')
  
  # shinyDirChoose(input,'new_save_directory', roots=volumes, session=session, 
  #                 defaultPath='', defaultRoot='wd')
  
  cat(file = stderr(), "Shiny Server started ...6", "\n")
  
  # test to see if dpmsr_set exisits - if so will update widgets with defaults
  observe({
    if (dpmsr_present$test) {
      cat(file = stderr(), "dpmsr_set found...", "\n")
      update_widget_all(session, input, output)
      inputloaddata_render(session, input, output)
      inputfilterapply_render(session, input, output)
      inputnorm_render(session, input, output)
      qc_render(session, input, output)
      inputproteinselect_render(session, input, output)
      update_dpmsr_set_from_widgets(session, input, output)
      cat(file = stderr(), "dpmsr_set loaded...", "\n")
    }else{
      cat(file = stderr(), "dpmsr_set doest not exist...", "\n")
    }
    
  })
  
  
  # function home to shinyjs hide/enable observers
  hide_enable(session, input, output)
  
  cat(file = stderr(), "Shiny Server started ...7", "\n")
  
  hideTab(inputId = "nlp1", target = "load")
  
  if (site_user == "dpmsr") {
    updateTabsetPanel(session, "nlp1", selected = "tp_load_design")
  }else{
    updateTabsetPanel(session, "nlp1", selected = "tp_customer_load")
  }
}


#----------------------------------------------------------------------------------------
# update dpmsr_set file (older versions missing variables)

update_version <- function(session, input, output){
  
  if (is.null(dpmsr_set$x$primary_group)) {
    cat(file = stderr(), ("Older dpmsr_set file, adding primary group field"), "\n")
    dpmsr_set$x$primary_group <<- FALSE
  }
  #does not have this field
  if (is.null(dpmsr_set$x$rollup_method)) {
    cat(file = stderr(), ("Older dpmsr_set file, adding rollup method"), "\n")
    dpmsr_set$x$rollup_method <<- "Sum"
  }
  
  if (is.null(dpmsr_set$x$directlfq)) {dpmsr_set$x$directlfq <<- 0}
  if (is.null(dpmsr_set$x$pathway_set)) {dpmsr_set$x$pathway_set <<- 0}
  
  #update sample groups
  if (ncol(dpmsr_set$y$sample_groups) < 7) {
    cat(file = stderr(), ("Older dpmsr_set file, rerunning function set_sample_groups"), "\n")
    set_sample_groups(session, input, output)
  }
 
  #change imputed column name
  for (name in names(dpmsr_set$data$impute)) {
    df <- dpmsr_set$data$impute[[name]]
    names(df)[names(df) == 'PD_Detected_Peptides'] <- 'Detected_Imputed'
    dpmsr_set$data$impute[[name]] <<- df
  }
  
  
  
  if ("Genes" %in% colnames(dpmsr_set$data$impute$impute)) {
    cat(file = stderr(), ("Gene column found in impute"), "\n")
  }else{
    cat(file = stderr(), ("Gene column NOT found in impute"), "\n")
    #loop through data and add Gene column
    for (name in names(dpmsr_set$data$impute)) {
      df <- dpmsr_set$data$impute[[name]]
      add_gene <- str_extract(df$Description, "GN=\\w*")
      add_gene <- gsub("GN=", "", add_gene)
      dpmsr_set$data$impute[[name]] <<- df %>% add_column(Genes = add_gene, .before = 3) 
    }
    dpmsr_set$y$info_columns <<- dpmsr_set$y$info_columns + 1
    dpmsr_set$y$info_columns_final <<- dpmsr_set$y$info_columns_final + 1
  }
  
  if ("Genes" %in% colnames(dpmsr_set$data$final$impute)) {
    cat(file = stderr(), ("Gene column found in final"), "\n")
  }else{
    cat(file = stderr(), ("Gene column NOT found in final"), "\n")
    #loop through data and add Gene column
    for (name in names(dpmsr_set$data$final)) {
      df <- dpmsr_set$data$final[[name]]
      add_gene <- str_extract(df$Description, "GN=\\w*")
      add_gene <- gsub("GN=", "", add_gene)
      dpmsr_set$data$final[[name]] <<- df %>% add_column(Genes = add_gene, .before = 3) 
    }
  }  
  
  
  

}