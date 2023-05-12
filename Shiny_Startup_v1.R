
app_startup <- function(session, input, output){
  
  # check if dpmsr_set file is loaded
  dpmsr_present <<- reactiveValues()
  dpmsr_present$test <<- exists("dpmsr_set")
  
  cat(file=stderr(), "Shiny Server started ...2", "\n")
  
  shinyjs::disable("action_load_design")
  shinyjs::disable("action_load_dpmsr_set")
  shinyjs::disable("load_data")
  
  cat(file=stderr(), "Shiny Server started ...3", "\n")
  
  output$app_version_text <- renderText({str_c('app version = ', app_version)})
  output$app_version_text2 <- renderText({str_c('app version = ', app_version)})
  
  output$text_n1 <- renderText("Check all Normalization Strategies")
  output$text_f1 <- renderText("Normalization Filters")
  output$text_i1 <- renderText("Select Imputation Method")
  output$text_impute_ptm <- renderText("Impute Distributions using Impute PTM grep (Load Data)")
  
  cat(file=stderr(), "Shiny Server started ...4", "\n")
  
  shinyFileChoose(input, 'design_file', session=session, roots=volumes, filetypes=c('', 'xlsx'))

  shinyFileChoose(input,'dpmsr_set_file', roots=volumes, session=session, 
                  filetypes=c('','dpmsr_set'), defaultPath='', defaultRoot='wd')
  
  shinyFileChoose(input,'motif_fasta', roots=volumes, session=session, 
                  filetypes=c('','fasta'), defaultPath='', defaultRoot='wd')
  
  shinyFileChoose(input,'test', roots=volumes, session=session, 
                  filetypes=c('' , 'xlsx', 'jpg', 'png'), defaultPath='', defaultRoot='wd')
  
  # shinyDirChoose(input,'new_save_directory', roots=volumes, session=session, 
  #                 defaultPath='', defaultRoot='wd')
  
  cat(file=stderr(), "Shiny Server started ...6", "\n")
  
  # test to see if dpmsr_set exisits - if so will update widgets with defaults
  observe({
    if (dpmsr_present$test){
      cat(file=stderr(), "dpmsr_set found...", "\n")
      update_widget_all(session, input, output)
      inputloaddata_render(session, input, output)
      inputfilterapply_render(session, input, output)
      inputnorm_render(session, input, output)
      qc_render(session, input, output)
      inputproteinselect_render(session, input, output)
      update_dpmsr_set_from_widgets(session, input, output)
      cat(file=stderr(), "dpmsr_set loaded...", "\n")
    }else{
      cat(file=stderr(), "dpmsr_set doest not exist...", "\n")
    }
    
  })
  
  
  # function home to shinyjs hide/enable observers
  hide_enable(session, input, output)
  
  cat(file=stderr(), "Shiny Server started ...7", "\n")
  
  hideTab(inputId = "nlp1", target = "load")
  
  if (site_user == "dpmsr"){
    updateTabsetPanel(session, "nlp1", selected = "tp_load_design")
  }else{
    updateTabsetPanel(session, "nlp1", selected = "tp_customer_load")
  }
}