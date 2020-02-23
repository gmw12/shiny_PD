rm(list = ls())

options(shiny.maxRequestSize=100*1024^2)

source("Shiny_Startup_v1.R")

function(input, output, session) {
    useShinyjs()
    dpmsr_present <- reactiveValues()
    dpmsr_present$test <- exists("dpmsr_set")
    
    shinyjs::disable("action_load_design")
    shinyjs::disable("load_data")
    
    output$text_n1 <- renderText("Check all Normalization Strategies")
    output$text_i1 <- renderText("Select Imputation Method")
    
    volumes <- c(wd='.', Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
    #uncomment for titan_black VM, and comment above line
    #volumes <- c(dd='/home/dpmsr/shared', wd='.', Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
    shinyFileChoose(input, 'design_file', session=session, roots=volumes, filetypes=c('', 'xlsx'))
    
    shinyFileChoose(input,'raw_files', roots=volumes, session=session, 
                    filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    
    shinyFileChoose(input,'dpmsr_set_file', roots=volumes, session=session, 
                    filetypes=c('','dpmsr_set'), defaultPath='', defaultRoot='wd')
    
    shinyFileChoose(input,'motif_fasta', roots=volumes, session=session, 
                     filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    
    # test to see if dpmsr_set exisits - if so will update widgets with defaults
    observe({
      if (dpmsr_present$test){
        update_shiny_defaults(session)
        update_dpmsr_set_from_widgets(session, input)
        inputloaddata_render(session, input, output)
        inputfilterapply_render(session, input, output)
        inputnorm_render(session, input, output)
        inputstartstats_render(session, input, output)
        inputproteinselect_render(session, input, output)
        selection_updates_new_start(session, input, output)
      }
    })
    
    # function home to shinyjs hide/enable observers
    hide_enable(session, input)
  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$action_load_design, {
    load_dpmsr_set(session, input, volumes)
    #reassign these to update volumes with path from design_file
    shinyFileChoose(input,'raw_files', roots=dpmsr_set$x$volumes, session=session,
                    filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    shinyFileChoose(input,'motif_fasta', roots=dpmsr_set$x$volumes, session=session, 
                    filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    set_sample_groups()
    file_set()
    save_design(session, input)
    updateTabsetPanel(session, "nlp1", selected = "tp_load_data")
    dpmsr_present$test <- exists("dpmsr_set")
    #update_shiny_defaults(session)
  })


#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------    
  observeEvent(input$load_data, {
    showModal(modalDialog("Working...", footer = NULL))
    #inFile2 <- input$file_raw_data
    #load_data(input$radio_input, inFile2$datapath)
    load_data(session, input, volumes)
    prepare_data(session, input)
    project_overview()
    inputloaddata_render(session, input, output)
    updateTabsetPanel(session, "nlp1", selected = "tp_overview")
    removeModal()
  })

#------------------------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$filter_apply, {
    showModal(modalDialog("Working...", footer = NULL))
    preprocess_data()

    removeModal()
    updateTabsetPanel(session, "nlp1", selected = "tp5")
    
    bar_plot(dpmsr_set$data$data_peptide[(dpmsr_set$y$info_columns+1):ncol(dpmsr_set$data$data_peptide)],"Raw", dpmsr_set$file$qc_dir)
    box_plot(dpmsr_set$data$data_peptide[(dpmsr_set$y$info_columns+1):ncol(dpmsr_set$data$data_peptide)],"Raw", dpmsr_set$file$qc_dir)
    
    inputfilterapply_render(session, input, output)
  })
  
#------------------------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$norm1, {
    showModal(modalDialog("Working...", footer = NULL))
    norm_prep()
    histogram_plot()
    
    inputnorm_render(session, input, output)
    
    if(input$checkbox_n1 == 1) {n1<-1}else{n1<-NULL}
    if(input$checkbox_n2 == 1) {n2<-2}else{n2<-NULL}
    if(input$checkbox_n3 == 1) {n3<-3}else{n3<-NULL}
    if(input$checkbox_n4 == 1) {n4<-4}else{n4<-NULL}
    if(input$checkbox_n5 == 1) {n5<-5}else{n5<-NULL}
    if(input$checkbox_n6 == 1) {n6<-6}else{n6<-NULL}
    if(input$checkbox_n7 == 1) {n7<-7}else{n7<-NULL}
    if(input$checkbox_n8 == 1) {n8<-8}else{n8<-NULL}
    if(input$checkbox_n9 == 1) {n9<-9}else{n9<-NULL}
    if(input$checkbox_n10 == 1) {n10<-10}else{n10<-NULL}
    if(input$checkbox_n11 == 1) {n11<-11}else{n11<-NULL}
    norm_list <-  list("impute"=99,"sl"=n1, "quantile"=n4,"lr"=n5,"loess"=n6,"vsn"=n7,"ti"=n8,"mi"=n9,"ai"=n10)
    dpmsr_set$y$norm_list <<- norm_list[lapply(norm_list, length) > 0]
    norm_list2 <-  list("tmm"=n2,"sltmm"=n3,"protein"=n11)
    dpmsr_set$y$norm_list2 <<- norm_list2[lapply(norm_list2, length) > 0]
    apply_norm()
    removeModal()
    updateTabsetPanel(session, "nlp1", selected = "tp6")
  })

#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$impute1, {
    showModal(modalDialog("Working...", footer = NULL))  
    apply_impute()
    removeModal()
    updateTabsetPanel(session, "nlp1", selected = "tp7")
  })
  
    
#------------------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$start_stats, {
    showModal(modalDialog("Working...", footer = NULL))  
    set_comp_groups()
    apply_stats()
    qc_apply()
    create_plots()
    
    inputstartstats_render(session, input, output)
    
    #update pulldown for proteinQC plots
    updateSelectInput(session, "norm_type", choices = names(dpmsr_set$data$final) , selected= "impute")
    updateSelectInput(session, "volcano_select", choices = names(dpmsr_set$data$final), selected = "impute")
    updateSelectInput(session, "select_data_comp", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
    updateSelectInput(session, "select_data_comp_motif", choices = dpmsr_set$y$comp_groups$comp_name, selected= dpmsr_set$y$comp_groups$comp_name[1])
    updateSelectInput(session, "select_final_data", choices = names(dpmsr_set$data$final), selected= "impute")
    updateSelectInput(session, "select_final_data_motif", choices = names(dpmsr_set$data$final), selected= "impute")
    updateSelectInput(session, "select_oneprotein_norm", choices = names(dpmsr_set$data$final), selected= "impute")
    updateSelectInput(session, "select_onepeptide_norm", choices = names(dpmsr_set$data$final), selected= "impute")
    
    updateTabsetPanel(session, "nlp1", selected = "tp8") 
    Final_Excel()
    removeModal()
    }) # end of observe stat button

#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------  

    
    observeEvent(input$protein_select_plots, {     
      showModal(modalDialog("Working...", footer = NULL))  
      create_select_plots()
    
      updateSelectInput(session, "protein_select", choices = dpmsr_set$y$protein_list , selected= "ADH")
      
      inputproteinselect_render(session, input, output)
      
      removeModal()
    })
 
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    
    observeEvent(input$oneprotein_show, { 
      filter_df <- dpmsr_set$data$final[[input$select_oneprotein_norm]]
      
      if(input$oneprotein_accession != "0" ){
        filter_df <-subset(filter_df, Accession %in% as.character(input$oneprotein_accession)  )
      }
      stats_df <- filter_df[(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number+1):ncol(filter_df)]
      stats_df <- stats_df[,-grep(pattern="_FC2$", colnames(stats_df))]
      stats_df <- t(stats_df)
      filter_df <- filter_df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
      plot_dir <- str_c(dpmsr_set$file$output_dir, input$select_oneprotein_norm, "//")
      bar_plot(filter_df, str_c(input$oneprotein_accession), plot_dir)
      
      output$oneprotein_plot_select <- renderImage({
        list(src=str_c(plot_dir, input$oneprotein_accession, "_barplot.png"),  
             contentType = 'image/png', width=500, height=400, alt="this is alt text")
      }, deleteFile = FALSE)
      
      output$oneprotein_stats<- renderRHandsontable({
        rhandsontable(stats_df, colHeaders = NULL, rowHeaderWidth = 300) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) })
      
    })
    
#-------------------------------------------------------------------------------------------------------------      
#-------------------------------------------------------------------------------------------------------------  
    
    
    observeEvent(input$onepeptide_show, { 
      filter_df <- dpmsr_set$data$final[[input$select_onepeptide_norm]]
      
      if(input$onepeptide_sequence != "0" ){
        filter_df <-subset(filter_df, Sequence %in% as.character(input$onepeptide_sequence)  )
      }
      stats_df <- filter_df[(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number+1):ncol(filter_df)]
      stats_df <- stats_df[,-grep(pattern="_FC2$", colnames(stats_df))]
      stats_df <- t(stats_df)
      filter_df <- filter_df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
      plot_dir <- str_c(dpmsr_set$file$output_dir, input$select_onepeptide_norm, "//")
      bar_plot(filter_df, str_c(input$onepeptide_sequence), plot_dir)
      
      output$onepeptide_plot_select <- renderImage({
        list(src=str_c(plot_dir, input$onepeptide_sequence, "_barplot.png"),  
             contentType = 'image/png', width=500, height=400, alt="this is alt text")
      }, deleteFile = FALSE)
      
      output$onepeptide_stats<- renderRHandsontable({
        rhandsontable(stats_df, colHeaders = NULL, rowHeaderWidth = 300) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) })
      
    })
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    
    
observeEvent(input$data_show, { 
    filter_df <- dpmsr_set$data$final[[input$select_final_data]]

    if(input$data_topn != 0 ){
      filter_df$sum <- rowSums(filter_df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)])
      filter_df <- filter_df[order(-filter_df$sum),]                      
      filter_df <- filter_df[1:input$data_topn,]
      filter_df$sum <- NULL
    }
    
    if(input$data_accession != "0" ){
      filter_df <-subset(filter_df, Accession %in% as.character(input$data_accession)  )
    }
      
    if(input$data_description != "0") {
      filter_df <-filter_df[grep(as.character(input$data_description), filter_df$Description), ]
    }

    if(input$data_pvalue != 0 & input$data_foldchange1 != 0 & input$data_foldchange2 != 0) {
      filter_df <- subset(filter_df, filter_df[ ,str_c(as.character(input$select_data_comp),"_Pval")] <= as.numeric(input$data_pvalue) &  
                              (filter_df[ ,str_c(as.character(input$select_data_comp),"_FC")] >= as.numeric(input$data_foldchange1) | 
                                 filter_df[ ,str_c(as.character(input$select_data_comp),"_FC")] <= -as.numeric(input$data_foldchange2))) 
    }
    
    output$data_final<- renderRHandsontable({
      rhandsontable(filter_df, rowHeaders = NULL, columnHeaderwidth = 200, width = 1200, height = 500) %>%
        hot_cols(fixedColumnsLeft = 1, colWidths = 80, halign = "htCenter" ) %>%
        #hot_col(col = "Peptides", format="0") %>%
        hot_col(col = dpmsr_set$y$info_columns_final:(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number), format="0", digits =0) %>%
        hot_col(col = (dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number+1):(ncol(dpmsr_set$data$final$impute)), format="0.00") %>%
        hot_col(col = "Description", halign = "htLeft", colWidths = 350) %>%
        hot_rows(rowHeights = 10)
    })
  })  
    
  
#-------------------------------------------------------------------------------------------------------------      
#-------------------------------------------------------------------------------------------------------------   
observeEvent(input$save_dpmsr_set, { 
  save(dpmsr_set, file=str_c(dpmsr_set$file$output_dir, dpmsr_set$x$file_prefix, ".dpmsr_set"))
})  

    
#-------------------------------------------------------------------------------------------------------------      
#------------------------------------------------------------------------------------------------------------- 
observeEvent(input$load_dpmsr_set, { 
  showModal(modalDialog("Working...", footer = NULL))  
  dpmsr_file <- parseFilePaths(volumes, input$dpmsr_set_file)
  load(file = dpmsr_file$datapath, envir = .GlobalEnv)
  
  update_shiny_defaults(session)
  update_dpmsr_set_from_widgets(session, input)
  inputloaddata_render(session, input, output)
  inputfilterapply_render(session, input, output)
  inputnorm_render(session, input, output)
  inputstartstats_render(session, input, output)
  inputproteinselect_render(session, input, output)
  selection_updates_new_start(session, input, output)
  removeModal()
})  

    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$motif_show, {
      showModal(modalDialog("Working...", footer = NULL))  

      filter_df <- dpmsr_set$data$final[[input$select_final_data_motif]]
      motif_data <- run_motifx(input, output, filter_df)
      
      output$motif_table<- renderRHandsontable({
        rhandsontable(motif_data, rowHeaders = NULL) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) %>%
          hot_col(col = "comparison", halign = "htCenter", colWidths = 150) %>%
          hot_col(col = "motif", halign = "htCenter", colWidths = 100) %>%
          hot_col(col = "fold.increase", halign = "htCenter", colWidths = 100)
      })
      removeModal()
    })
    
    
# 
#      observe({
#        if (!is.null(input$motif_fasta)){
#          fasta_txt_file <- input$motif_fasta
#          fasta_txt_file <- unlist(fasta_txt_file$files[[as.character(0)]][2])
#          output$fasta <- renderText(fasta_txt_file)
#        }
#      })
#      
    output$fasta <- renderText({
      if (req(typeof(input$motif_fasta)=="list")) {
        motif_path <- parseFilePaths(dpmsr_set$x$volumes, input$motif_fasta)
        basename(motif_path$datapath)
        #unlist(input$motif_fasta$files[[as.character(0)]][2])
      }
    })
  
}



