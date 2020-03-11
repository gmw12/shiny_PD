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
    selection_update_after_stats(session, input, output)
 
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
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$set_pathway, {
      showModal(modalDialog("Working...", footer = NULL))  
      set_pathway(input, output, session)
      removeModal()
      updateTabsetPanel(session, "nlp1", selected = "wiki")
    }
    )
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$wiki_show, {
      showModal(modalDialog("Working...", footer = NULL))  
      
      wiki_df <- dpmsr_set$data$final[[input$select_final_data_wiki]]
      wiki_data <<- run_wiki(input, output, wiki_df)
      wiki_data <<- wiki_list_data[[1]]
      
      output$wiki_table<- renderRHandsontable({
        rhandsontable(wiki_data, rowHeaders = NULL, readOnly = TRUE) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) %>%
          hot_col(col = "ID", halign = "htCenter", colWidths = 120) %>%
          hot_col(col = "Description", halign = "htCenter", colWidths = 250) %>%
          hot_col(col = "pvalue", halign = "htCenter", colWidths = 150, format = '0.000000') %>% 
          hot_col(col = "p.adjust", halign = "htCenter", colWidths = 150, format = '0.000000') %>% 
          hot_col(col = "qvalue", halign = "htCenter", colWidths = 100, format = '0.00000') %>% 
          hot_col(col = "geneID", halign = "htCenter", colWidths = 400)
      })
      removeModal()
      }
    )
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$profile_show, {
      showModal(modalDialog("Working...", footer = NULL))  
      
      profile_df <- dpmsr_set$data$final[[input$select_final_data_profile]]
      profile_data <- run_profile(input, output, profile_df)
      
      go_profile_result <<- profile_data@result[1:5]
      go_profile_result <<- go_profile_result[order(-go_profile_result$Count),]
      
      output$go_profile_table <- renderRHandsontable({
        rhandsontable(go_profile_result, rowHeaders = NULL, readOnly = TRUE) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) %>%
          hot_col(col = "ID", halign = "htCenter", colWidths = 120) %>%
          hot_col(col = "Description", halign = "htCenter", colWidths = 250) %>%
          hot_col(col = "Count", halign = "htCenter", colWidths = 50) %>%
          hot_col(col = "GeneRatio", halign = "htCenter", colWidths = 100) %>%
          hot_col(col = "geneID", halign = "htCenter", colWidths = 150)
      })
      
      output$go_profile_plot <- renderPlot({
        barplot(go_profile_data, title = str_c("Go Profile ", input$select_ont_profile, " level=", input$select_level_profile), 
                drop=TRUE, showCategory=12, order=TRUE)
      })
      
      removeModal()
    }
    )
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_show, {
      showModal(modalDialog("Working...", footer = NULL))  
      
      go_df <- dpmsr_set$data$final[[input$select_final_data_go]]
      #go_data <<- try(run_go_analysis(input, output, go_df), silent = FALSE)
      go_data <- run_go_analysis(input, output, go_df)
      #go_data <- goresult
      setup_go_volcano(input, output, go_df)
      #go_data <- try(go_data[order(go_data$pvalue),], silent = TRUE)
      
      output$go_table <- renderRHandsontable({
        rhandsontable(go_data, rowHeaders = NULL, readOnly = FALSE) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) %>%
          hot_col(col = "GO.ID", halign = "htCenter", colWidths = 100) %>%
          hot_col(col = "Term", halign = "htCenter", colWidths = 200) %>%
          hot_col(col = "Rank in classicFisher", halign = "htCenter", colWidths = 200)   %>%
          hot_col(col = "classicFisher", halign = "htCenter", colWidths = 100)  
      })
      
      # output$wiki_plot <- renderPlot({
      #   barplot(wiki_data, title = str_c("WikiPathways ", input$select_ont_wiki, " level=", input$select_level_wiki), 
      #           drop=TRUE, showCategory=12, order=TRUE)
      # })
      
      removeModal()
    }
    )   
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_volcano, {
      showModal(modalDialog("Working...", footer = NULL))  

      volcano_data <- create_go_volcano(input, output, session)
      removeModal()
      
      #output$scatterplot <- renderPlot({
      volcano_go_plot <- reactive({
        ggplot(volcano_data, aes(x = log_fc, y = log_pvalue)) +
          theme_minimal() +
          geom_point(alpha=0.4, size=input$plot_dot_size, color = input$volcano_dot_color) +
          xlab(input$plot_x_axis_label) + ylab(input$plot_y_axis_label) +
          scale_colour_gradient(low = input$volcano_dot_color, high = input$volcano_dot_color) +
          ggtitle(input$plot_title)+    
          xlim(-max(volcano_data$log_fc), max(volcano_data$log_fc)) +
          theme(plot.title = element_text(size=input$plot_title_size, hjust = 0.5),
                axis.title = element_text(size=input$plot_label_size, color="black"),
                axis.text.x = element_text(size=10, color="black"),
                axis.text.y = element_text(size=10,  color="black"),
                legend.position = "none")
      })
      
      output$volcano_go_plot <- renderPlot({
        req(volcano_go_plot())
        volcano_go_plot()
      })
      
      output$Download <- downloadHandler(
        filename = function(){
          str_c(dpmsr_set$file$string, "GoVolcano_", input$select_data_comp_go, "_", input$go_volcano_id, "_", 
                input$select_ont_go, ".png", collapse = " ")
        },
        content = function(file){
          req(volcano_go_plot())
          ggsave(file, plot = volcano_go_plot(), device = 'png')
        }
      )
      
      
      
      
      
      
      output$hover_info <- renderUI({
        hover <- input$plot_hover
        point <- nearPoints(volcano_data, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
        if (nrow(point) == 0) return(NULL)
        
        # calculate point position INSIDE the image as percent of total dimensions
        # from left (horizontal) and from top (vertical)
        left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
        top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
        
        # calculate distance from left and bottom side of the picture in pixels
        left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
        top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
        
        # create style property fot tooltip
        # background color is set so tooltip is a bit transparent
        # z-index is set so we are sure are tooltip will be on top
        style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                        "left:", left_px + 2, "px; top:", top_px + 2, "px;")
        
        # actual tooltip created as wellPanel
        wellPanel(
          style = style,
          p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                        "<b> FC: </b>", point$foldchange, "<br/>",
                        "<b> pvalue: </b>", point$pvalue, "<br/>")))
        )
      })
      
      
    }
    )   
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$setup_string, {
      showModal(modalDialog("Working...", footer = NULL))  
      setup_string(input, output, dpmsr_set$data$final[[input$select_final_data_string]])
      removeModal()
    }
    ) 
    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_string, {
      showModal(modalDialog("Working...", footer = NULL))  
      
      string_list <- run_string(input, output)
      
      output$string_plot <- renderImage({
        list(src=string_list$string_file_name,  
             contentType = 'image/png', width=2000, height=2000, alt="this is alt text")
      }, deleteFile = FALSE)
      
      output$string_link <- renderText({string_list$linkthis})
      
      removeModal()
    }
    ) 
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_string_enrich, {
      showModal(modalDialog("Working...", footer = NULL))  
      
      string_result <- run_string_enrich(input, output, dpmsr_set$data$final[[input$select_final_data_string]])
      
      output$string_table<- renderRHandsontable({
        rhandsontable(string_result, rowHeaders = NULL, readOnly = TRUE) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) %>%
          hot_col(col = "term_id", halign = "htCenter", colWidths = 100) %>%
          hot_col(col = "proteins", halign = "htCenter", colWidths = 100) %>%  
          hot_col(col = "hits", halign = "htCenter", colWidths = 100)  %>% 
          hot_col(col = "pvalue", halign = "htCenter", colWidths = 180, format = '0.0000000000000') %>% 
          hot_col(col = "pvalue_fdr", halign = "htCenter", colWidths = 180, format = '0.0000000000000') %>% 
          hot_col(col = "term_description", halign = "htCenter", colWidths = 600)
      })
      
      
      removeModal()
    }
    )    
  
}



