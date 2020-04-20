
options(shiny.maxRequestSize=100*1024^2)

source("Shiny_Startup_v1.R")

shinyServer(function(input, output, session) {
    cat(file=stderr(), "Shiny Server started ...1", "\n")
    useShinyjs()
    dpmsr_present <- reactiveValues()
    dpmsr_present$test <- exists("dpmsr_set")
    
    cat(file=stderr(), "Shiny Server started ...2", "\n")
    
    shinyjs::disable("action_load_design")
    shinyjs::disable("action_load_dpmsr_set")
    shinyjs::disable("load_data")
    
    cat(file=stderr(), "Shiny Server started ...3", "\n")
    
    output$text_n1 <- renderText("Check all Normalization Strategies")
    output$text_i1 <- renderText("Select Imputation Method")
    
    cat(file=stderr(), "Shiny Server started ...4", "\n")
    
    if(Sys.info()["sysname"]=="Darwin" ){
      volumes <- c(dd='/Users/gregwaitt/Documents/Data', wd='.', Home = fs::path_home(),  getVolumes()())
      #volumes <- c(wd='.', Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
    }else{
      #for titan_black VM
      volumes <- c(dd='/home/dpmsr/shared', wd='.', Home = fs::path_home(), getVolumes()())
    }
    
    cat(file=stderr(), "Shiny Server started ...5", "\n")
    
    shinyFileChoose(input, 'design_file', session=session, roots=volumes, filetypes=c('', 'xlsx'))
    
    shinyFileChoose(input,'raw_files', roots=volumes, session=session, 
                    filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    
    shinyFileChoose(input,'dpmsr_set_file', roots=volumes, session=session, 
                    filetypes=c('','dpmsr_set'), defaultPath='', defaultRoot='wd')
    
    shinyFileChoose(input,'motif_fasta', roots=volumes, session=session, 
                     filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    
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
        update_dpmsr_set_from_widgets(session, input)
        cat(file=stderr(), "dpmsr_set loaded...", "\n")
      }else{
        cat(file=stderr(), "dpmsr_set doest not exist...", "\n")
      }
      
    })
    
    # function home to shinyjs hide/enable observers
    hide_enable(session, input)
 
    cat(file=stderr(), "Shiny Server started ...7", "\n")
    
    hideTab(inputId = "nlp1", target = "load")
    updateTabsetPanel(session, "nlp1", selected = "tp1")
    
  #------------------------------------------------------------------------------------------------------  
    observeEvent(input$action_clear, {
      cat(file=stderr(), "clean old data", "\n")
      rm(list = ls(envir = .GlobalEnv), pos = .GlobalEnv, inherits = FALSE)
      cat(file=stderr(), "reload libraries and functions", "\n")
      source("Shiny_Startup_v1.R")
    })    
     
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$action_load_design, {
    cat(file=stderr(), "load design triggered", "\n")
    load_dpmsr_set(session, input, volumes)
    #reassign these to update volumes with path from design_file
    shinyFileChoose(input,'raw_files', roots=dpmsr_set$x$volumes, session=session,
                    filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    shinyFileChoose(input,'motif_fasta', roots=dpmsr_set$x$volumes, session=session, 
                    filetypes=c('','txt'), defaultPath='', defaultRoot='wd')
    set_sample_groups()
    #file_set()
    #save_design(session, input)
    updateTabsetPanel(session, "nlp1", selected = "tp_load_data")
    dpmsr_present$test <- exists("dpmsr_set")
    update_widget_startup(session, input, output)
    update_widget_filter(session, input, output)
    update_widget_norm(session, input, output)
    update_widget_impute(session, input, output)
    update_widget_stats(session, input, output)
  })

    cat(file=stderr(), "Shiny Server started ...8", "\n")
#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------    
  observeEvent(input$load_data, {
    cat(file=stderr(), "Load and save design", "\n")
    showModal(modalDialog("Loading Data...", footer = NULL))
    cat(file=stderr(), "file_set", "\n")
    file_set()
    cat(file=stderr(), "save design", "\n")
    save_design(session, input)
    cat(file=stderr(), "load data", "\n")
    load_data(session, input, volumes)
    cat(file=stderr(), "prepare data", "\n")
    prepare_data(session, input)
    removeModal()
    showModal(modalDialog("Creating Overview...", footer = NULL))
    cat(file=stderr(), "project overview", "\n")
    project_overview()
    inputloaddata_render(session, input, output)
    updateTabsetPanel(session, "nlp1", selected = "tp_overview")
    cat(file=stderr(), "order data", "\n")
    preprocess_order()
    removeModal()
  })

#------------------------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$filter_apply, {
    cat(file=stderr(), "filter apply triggered", "\n")
    showModal(modalDialog("Applying Filters...", footer = NULL))
    cat(file=stderr(), "preprocess filter", "\n")
    preprocess_filter()
    removeModal()
    updateTabsetPanel(session, "nlp1", selected = "tp5")
    
    bar_plot(dpmsr_set$data$data_peptide[(dpmsr_set$y$info_columns+1):ncol(dpmsr_set$data$data_peptide)],"Raw", dpmsr_set$file$qc_dir)
    box_plot(dpmsr_set$data$data_peptide[(dpmsr_set$y$info_columns+1):ncol(dpmsr_set$data$data_peptide)],"Raw", dpmsr_set$file$qc_dir)

    showModal(modalDialog("Preparing Data for Norm...", footer = NULL))
    cat(file=stderr(), "prepare data for normalization", "\n")
    norm_prep()
    removeModal()
    showModal(modalDialog("Preparing Data for Histogram Plot...", footer = NULL))
    histogram_plot()
    removeModal()
    
    inputfilterapply_render(session, input, output)
  })
  
#------------------------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$norm1, {
    cat(file=stderr(), "normalization triggered", "\n")
    showModal(modalDialog("Normalizing Data...", footer = NULL))
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
    if (input$checkbox_tmt){
      tmt_spqc_prep()
      updateTabsetPanel(session, "nlp1", selected = "tp_tmt")
    }else{
    updateTabsetPanel(session, "nlp1", selected = "tp6")
    }
  })

#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$impute1, {
    cat(file=stderr(), "impute triggered", "\n")
    showModal(modalDialog("Imputing Missing Data...", footer = NULL))  
    apply_impute()
    removeModal()
    # cat(file=stderr(), "Set comp groups...", "\n")
    # showModal(modalDialog("Setting groups...", footer = NULL))  
    # set_comp_groups()
    # removeModal()
    cat(file=stderr(), "stat prep...", "\n")
    showModal(modalDialog("Stat prep...", footer = NULL)) 
    stat_prep()
    removeModal()
    cat(file=stderr(), "qc_apply", "\n")
    showModal(modalDialog("QC stats...", footer = NULL)) 
    qc_apply()
    removeModal()
    cat(file=stderr(), "qc - create_plots", "\n")
    showModal(modalDialog("QC plots...", footer = NULL)) 
    create_qc_plots()
    qc_render(session, input, output)
    update_widget_post_processing(session, input, output)
    removeModal()
    updateTabsetPanel(session, "nlp1", selected = "tp_qc")
  })
  
#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------  
    observeEvent(input$tmt_irs_go, {
      cat(file=stderr(), "TMT IRS started...", "\n")
      showModal(modalDialog("IRS Normalization...", footer = NULL))  
      tmt_spqc_normalize(dpmsr_set$data$normalized$sl)
      removeModal()
    })
    

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
      stats_df <- t(stats_df)
      stats_df <- data.frame(stats_df)
      colnames(stats_df)<- "CV"
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
      stats_df <- t(stats_df)
      stats_df <- data.frame(stats_df)
      colnames(stats_df)<- "CV"
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
    observeEvent(input$check_stats, {
      showModal(modalDialog("Setting Stat groups...", footer = NULL))  
      cat(file=stderr(), "Setting Stat groups...", "\n")
      set_stat_groups(session, input, output)
      removeModal()
    }
    )  
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$start_stats, {
      showModal(modalDialog("Calculating stats...", footer = NULL))  
      cat(file=stderr(), "Calculating stats...", "\n")
      stat_calc(session, input, output)
      removeModal()
      showModal(modalDialog("Saving Final Excel...", footer = NULL))      
      cat(file=stderr(), "Saving Final Excel...", "\n")
      stats_Final_Excel()
      removeModal()
    }
    )      
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    
    
    observeEvent(input$create_stats_plots, {
      showModal(modalDialog("Working...", footer = NULL))  
      
      df <- dpmsr_set$data$stats$final[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)]
      comp_string <- input$stats_plot_comp
      cat(file=stderr(), "Create stats plots" , "\n")
      for (i in 1:length(comp_string)){
        comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string[i])
        if (i==1){
          comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
        }else{
          comp_rows <- c(comp_rows, dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
        }
      }
      #add spqc to plots
      if(input$stats_plot_spqc){
          comp_rows <- c(comp_rows, dpmsr_set$y$stats$comp_spqc_sample_numbers)
      }
      
      comp_rows <- sort(unique(unlist(comp_rows)), decreasing = FALSE)
      df <- df[,comp_rows]
      namex <- dpmsr_set$design$Label[comp_rows]
      color_list <- dpmsr_set$design$colorlist[comp_rows]
      groupx <- dpmsr_set$design$Group[comp_rows]
      
      interactive_barplot(session, input, output, df, namex, color_list)
      interactive_boxplot(session, input, output, df, namex, color_list)
      interactive_pca2d(session, input, output, df, namex, color_list, groupx)
      interactive_pca3d(session, input, output, df, namex, color_list, groupx)
      interactive_cluster(session, input, output, df, namex)
      interactive_heatmap(session, input, output, df, namex, groupx)
      removeModal()
    }
    ) 
    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$create_stats_volcano, {
      showModal(modalDialog("Working...", footer = NULL))  
      #interactive_stats_volcano1(session, input, output)
      #interactive_stats_volcano2(session, input, output)
      # for(i in 1:input$comp_number){
      #   do.call(str_c("interactive_stats_volcano",i), list(session, input, output) )
      # }
      for(i in 1:input$comp_number){
        do.call("interactive_stats_volcano", list(session, input, output, i) )
      }
      
      removeModal()
    }
    )  
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    
    
    observeEvent(input$stats_data_show, { 
      showModal(modalDialog("Working...", footer = NULL))  
      cat(file=stderr(), "stats data show triggered..." , "\n")
      filter_df <- dpmsr_set$data$stats$final
      
      if(input$stats_data_topn != 0 ){
        filter_df$sum <- rowSums(filter_df[(dpmsr_set$y$info_columns_final+1):(dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number)])
        filter_df <- filter_df[order(-filter_df$sum),]                      
        filter_df <- filter_df[1:input$stats_data_topn,]
        filter_df$sum <- NULL
      }
      
      if(input$stats_data_accession != "0" ){
        filter_df <-subset(filter_df, Accession %in% as.character(input$stats_data_accession)  )
      }
      
      if(input$stats_data_description != "0") {
        filter_df <-filter_df[grep(as.character(input$stats_data_description), filter_df$Description), ]
      }
      
      #if(input$stats_data_pvalue != 0 & input$stats_data_foldchange1 != 0 & input$stats_data_foldchange2 != 0) {
      if(input$stats_data_pvalue != 0) {
        df <- filter_df
        comp_string <- input$stats_select_data_comp
        comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
        df_N <- df %>% dplyr::select(unlist(dpmsr_set$y$stats$groups$com_N_headers[comp_number]))  
        df_D <- df %>% dplyr::select(unlist(dpmsr_set$y$stats$groups$com_D_headers[comp_number]) ) 
        df_spqc <- df[unlist(dpmsr_set$y$stats$comp_spqc_sample_numbers)]
        df_spqc_cv <- df %>% dplyr::select(contains(str_c(dpmsr_set$y$stats$comp_spqc, "_CV") )) 
        df_Comp <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$comp_name [comp_number]) ) 
        df_final <- cbind(df[1:(dpmsr_set$y$info_columns_final)], df_N, df_D, df_spqc, df_spqc_cv, df_Comp)
        
        filter_df <- subset(df_final, df_final[ , dpmsr_set$y$stats$groups$pval[comp_number]] <= input$stats_data_pvalue &  
                              (df_final[ , dpmsr_set$y$stats$groups$fc[comp_number]] >= input$stats_data_foldchange1 | 
                                 df_final[ , dpmsr_set$y$stats$groups$fc[comp_number]] <= -input$stats_data_foldchange2)) 
        filter_df <- subset(filter_df, filter_df[ , dpmsr_set$y$stats$groups$mf[comp_number]] >= input$stats_missing_factor )
        if(input$stats2_spqc_cv_filter_factor >0){
          filter_df <- subset(filter_df, filter_df[ , str_c(dpmsr_set$y$stats$comp_spqc, "_CV")] <= input$stats2_spqc_cv_filter_factor )
        }
        
        
        filter_df_colnames <- colnames(filter_df)
        filter_df_colnames <- gsub("_v_", " v ", filter_df_colnames)
        filter_df_colnames <- gsub("_FC", " FC", filter_df_colnames)
        filter_df_colnames <- gsub("_MF", " MF", filter_df_colnames)
        filter_df_colnames <- gsub("_pval", " pval", filter_df_colnames)
        filter_df_colnames <- gsub("_", "", filter_df_colnames)
        colnames(filter_df) <-  filter_df_colnames
        
      }
      
      output$stats_data_final<- renderRHandsontable({
        rhandsontable(filter_df, rowHeaders = NULL, columnHeaderwidth = 200, width = 1500, height = 1500) %>%
          hot_cols(fixedColumnsLeft = 1, colWidths = 80, halign = "htCenter" ) %>%
          hot_col(col = "Peptides", format="0")  %>%
          #hot_col(col = colnames(df_N), format="0", digits =0) %>%
          #hot_col(col = colnames(df_D), format="0", digits =0) %>%
          #hot_col(col = (dpmsr_set$y$info_columns_final+dpmsr_set$y$sample_number+1):(ncol(dpmsr_set$data$final$impute)), format="0.00") %>%
          hot_col(col = "Description", halign = "htLeft", colWidths = 350) %>%
          hot_rows(rowHeights = 10)
      })
      
      
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
      
      wiki_df <- dpmsr_set$data$stats$final
      wiki_data <- run_wiki(input, output, wiki_df)
      #wiki_data <- wiki_data[[1]]
      
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
      
      profile_df <- dpmsr_set$data$stats$final
      profile_data <- run_profile(input, output, profile_df)
      
      go_profile_result <- profile_data@result[1:5]
      go_profile_result <- go_profile_result[order(-go_profile_result$Count),]
      
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
        barplot(profile_data, title = str_c("Go Profile ", input$select_ont_profile, " level=", input$select_level_profile), 
                drop=TRUE, showCategory=12, order=TRUE)
      })
      
      removeModal()
    }
    )
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_show, {
      showModal(modalDialog("Working...", footer = NULL))  
      
      go_df <- dpmsr_set$data$stats$final
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
      interactive_go_volcano(session, input, output)
      removeModal()
    }
    )  
    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$setup_string, {
      cat(file=stderr(), "string_setup triggered", "\n")
      showModal(modalDialog("String setup...", footer = NULL))  
      setup_string(input, output, dpmsr_set$data$stats$final)
      removeModal()
    }
    ) 
    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_string, {
      cat(file=stderr(), "go string triggered", "\n")
      showModal(modalDialog("String Analysis...", footer = NULL))  
      
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
    

    
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------   
    observeEvent(input$save_dpmsr_set, { 
      save(dpmsr_set, file=str_c(dpmsr_set$file$output_dir, dpmsr_set$x$file_prefix, ".dpmsr_set"))
    })  
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #------------------------------------------------------------------------------------------------------------- 
    observeEvent(input$action_load_dpmsr_set, { 
      cat(file=stderr(), "load dpmsr_set triggered", "\n")
      showModal(modalDialog("Working...", footer = NULL))  
      dpmsr_file <<- parseFilePaths(volumes, input$dpmsr_set_file)
      load(file = dpmsr_file$datapath, envir = .GlobalEnv)
      #reset file locations
      dpmsr_set$file$data_dir <<- dirname(dpmsr_file$datapath)
      dpmsr_set$file$data_path <<- gsub("(.*)/.*","\\1",dpmsr_set$file$data_dir)
      dpmsr_set$file$output_dir <<- str_replace_all(dpmsr_set$file$data_dir, "/", "//")
      dpmsr_set$file$output_dir <<- str_c(dpmsr_set$file$output_dir, "//")
      dpmsr_set$file$backup_dir <<- str_c(dpmsr_set$file$output_dir, "Backup//")
      dpmsr_set$file$extra_dir <<- str_c(dpmsr_set$file$output_dir, "Extra//")
      dpmsr_set$file$qc_dir <<- str_c(dpmsr_set$file$output_dir, "QC//")
      dpmsr_set$file$string <<- str_c(dpmsr_set$file$output_dir, "String//")
      dpmsr_set$file$extra_prefix2 <<- str_c(dpmsr_set$file$extra_dir, dpmsr_set$x$file_prefix)
      #reload shiny 
      update_widget_all(session, input, output)
      update_dpmsr_set_from_widgets(session, input)
      inputloaddata_render(session, input, output)
      inputfilterapply_render(session, input, output)
      inputnorm_render(session, input, output)
      qc_render(session, input, output)
      inputproteinselect_render(session, input, output)
      updateTabsetPanel(session, "nlp1", selected = "tp8") 
      removeModal()
    })  
    
    
       
}
)

