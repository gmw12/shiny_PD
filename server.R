options(shiny.maxRequestSize = 4000*1024^2)

cat(file = stderr(), "Server.R started", "\n")

#set user
set_user(session, input, output)

#force setting when testing customer version
#volumes <<- c(dd='/data', wd='.', Home = fs::path_home(), getVolumes()())
#site_user <<- "not_dpmsr"


shinyServer(function(input, output, session) {
  app_version <<- '2023.10.31'
  cat(file = stderr(), "Shiny Server started ...1", "\n")
  useShinyjs()
  
  if (exists("site_user")) {
    cat(file = stderr(), str_c("Site user =  ", site_user), "\n")
    }else{
      cat(file = stderr(), "Site user does not exist...", "\n")
      source("Shiny_Libraries_v1.R")
      source("Shiny_Functions_v1.R")
      set_user(session, input, output)
    }

    if (site_user == "not_dpmsr") {
      shinyalert("Welcome!", "Thanks for using the DPMSR Visualization Tool.  Please contact greg.waitt@duke.edu with any questions or concerns.  We are currently updating the app and trying to resolve disconnect issues.  Please rerun 'Set Comparisons' and 'Start Analysis' on the Stats Tab", type = "info")
      #shinyalert("DPMSR Visualization Tool is currently off line.  We are upgrading the tool to work with our new Thermo Astral which primarily uses DIA workflows.  Sorry for the inconvience. Please contact greg.waitt@duke.edu with any questions or concerns. ", type = "info")
    }
    
    #find app current state and set new if clear
    app_startup(session, input, output)

#------------------------------------------------------------------------------------------------------
# Clear memory - cannot load dpmsr_set if one already loaded, button clears memory to start fresh  
  observeEvent(input$action_clear, {
    cat(file = stderr(), "Clear Memory", "\n")
    cat(file = stderr(), str_c("DPMSR exists? ", exists('dpmsr_set')), "\n")
    try(rm(list = ls(envir = .GlobalEnv), pos = .GlobalEnv, inherits = FALSE))
    #try(rm("dpmsr_set", pos = .GlobalEnv, inherits = FALSE))
    cat(file = stderr(), str_c("Deleting DPMSR, still exists? ", exists('dpmsr_set')), "\n")
    gc()
    session$reload()  
    })    
     
#------------------------------------------------------------------------------------------------------  
#Load design file
  observeEvent(input$action_load_design, {
    cat(file = stderr(), "Load design triggered", "\n")
    dpmsr_present$test <- exists("dpmsr_set")
    if (dpmsr_present$test) {
      shinyalert("Oops!", "design info already loaded - Please clear memory", type = "error")
    }else{
      showModal(modalDialog("Loading Sample Information...", footer = NULL))
      load_dpmsr_set(session, input, output)
      set_sample_groups(session, input, output)
      updateTabsetPanel(session, "nlp1", selected = "tp_load_data")
      dpmsr_present$test <- exists("dpmsr_set")
      update_widget_startup(session, input, output)
      update_widget_filter(session, input, output)
      update_widget_norm(session, input, output)
      update_widget_impute(session, input, output)
      update_widget_rollup(session, input, output)
      update_widget_stats(session, input, output)
      #load stat design table
      create_design_table(session, input, output)
      removeModal()
    }
    cat(file = stderr(), "Load design file complete", "\n")
  })

#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------    
  observeEvent(input$load_data, {
    if ((dpmsr_set$x$peptide_isoform == TRUE) & (dpmsr_set$x$raw_data_input != "Peptide")) {
      shinyalert("Oops!", "Peptide Isoform is selected but input is NOT Peptide?", type = "error")
    }else{
      cat(file = stderr(), "Load and save design", "\n")
      showModal(modalDialog("Loading Data...", footer = NULL))
      cat(file = stderr(), "file_set", "\n")
      file_set()
      cat(file = stderr(), "save design", "\n")
      save_design(session, input)
      
      cat(file = stderr(), "load data", "\n")
      if (dpmsr_set$x$data_source == "SP") {
        load_data_sp(session, input, volumes)
      }else{
        load_data(session, input, volumes)
      }
      cat(file = stderr(), "prepare data", "\n")
      prepare_data(session, input)
      removeModal()

      if (dpmsr_set$x$raw_data_input != "Protein" && dpmsr_set$x$data_source == "PD") { 
        showModal(modalDialog("Creating Overview...", footer = NULL))
        
        cat(file = stderr(), "project overview", "\n")
        project_overview(session, input, output)
        removeModal()
        showModal(modalDialog("Loading Plots...", footer = NULL))
        inputloaddata_render(session, input, output)
        updateTabsetPanel(session, "nlp1", selected = "tp_overview")
        removeModal()
      } else {
          inputloaddata_render(session, input, output)
          updateTabsetPanel(session, "nlp1", selected = "tp_filters")
      }
      
      
      cat(file = stderr(), "preprocess and order data", "\n")
      preprocess_order()
      cat(file = stderr(), "check_sample_id()", "\n")
      check_sample_id()
      cat(file = stderr(), "load data complete", "\n")
      #removeModal()
    }
  })

#------------------------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$filter_apply, {
    cat(file = stderr(), "filter apply triggered", "\n")

    if (input$checkbox_require_x & (input$require_x_cutoff < 0 | input$require_x_cutoff >= 1))
    {
      shinyalert("Oops!", "Check Require X filters", type = "error")  
      }else{
          showModal(modalDialog("Applying Filters...", footer = NULL))
          cat(file = stderr(), "preprocess filter", "\n")
          preprocess_filter(session, input, output)
          removeModal()
          updateTabsetPanel(session, "nlp1", selected = "tp_normalize")

          showModal(modalDialog("Preparing Data for Norm...", footer = NULL))
          cat(file = stderr(), "prepare data for normalization", "\n")
          norm_prep()

          if (dpmsr_set$x$raw_data_input == "Protein") {
            cat(file = stderr(), "raw protein plots...", "\n")
            #check info columns for rerun of filter (impute column could be added
            dpmsr_set$y$info_columns <<- ncol(dpmsr_set$data$data_protein) - dpmsr_set$y$sample_number

            #display raw protein data
            bar_plot(dpmsr_set$data$data_protein[(dpmsr_set$y$info_columns + 1):ncol(dpmsr_set$data$data_protein)],"Raw", dpmsr_set$file$qc_dir)
            box_plot(dpmsr_set$data$data_protein[(dpmsr_set$y$info_columns + 1):ncol(dpmsr_set$data$data_protein)],"Raw", dpmsr_set$file$qc_dir)
          }else if (dpmsr_set$x$raw_data_input == "Peptide") {
            cat(file = stderr(), "raw peptide plots...", "\n")
            #check info columns for rerun of filter (impute column could be added
            dpmsr_set$y$info_columns <<- ncol(dpmsr_set$data$data_peptide) - dpmsr_set$y$sample_number

            #display raw peptide data
            bar_plot(dpmsr_set$data$data_peptide[(dpmsr_set$y$info_columns + 1):ncol(dpmsr_set$data$data_peptide)],"Raw", dpmsr_set$file$qc_dir)
            box_plot(dpmsr_set$data$data_peptide[(dpmsr_set$y$info_columns + 1):ncol(dpmsr_set$data$data_peptide)],"Raw", dpmsr_set$file$qc_dir)
          }else {
            cat(file = stderr(), "raw precursor plots...", "\n")
            #check info columns for rerun of filter (impute column could be added
            dpmsr_set$y$info_columns <<- ncol(dpmsr_set$data$data_precursor) - dpmsr_set$y$sample_number
            
            #display raw precursor data
            bar_plot(dpmsr_set$data$data_precursor[(dpmsr_set$y$info_columns + 1):ncol(dpmsr_set$data$data_precursor)],"Raw", dpmsr_set$file$qc_dir)
            box_plot(dpmsr_set$data$data_precursor[(dpmsr_set$y$info_columns + 1 ):ncol(dpmsr_set$data$data_precursor)],"Raw", dpmsr_set$file$qc_dir)
          }

          #change to always display norm data
          cat(file = stderr(), "normalized plots...", "\n")
          info_columns <- ncol(dpmsr_set$data$norm_data) - dpmsr_set$y$sample_number
          bar_plot(dpmsr_set$data$norm_data[(info_columns + 1):ncol(dpmsr_set$data$norm_data)],"Normalization_Data", dpmsr_set$file$qc_dir)
          box_plot(dpmsr_set$data$norm_data[(info_columns + 1):ncol(dpmsr_set$data$norm_data)],"Normalization_Data", dpmsr_set$file$qc_dir)


          removeModal()

          showModal(modalDialog("Preparing Data for Histogram Plot...", footer = NULL))

          if (dpmsr_set$x$raw_data_input == "Protein") {
            histogram_plot(dpmsr_set$data$data_protein, "Intensity_Histogram")
          }else if (dpmsr_set$x$raw_data_input == "Peptide") {
            histogram_plot(dpmsr_set$data$data_peptide, "Intensity_Histogram")
          }else {
            histogram_plot(dpmsr_set$data$data_precursor, "Intensity_Histogram")
          }

          if (as.logical(dpmsr_set$x$peptide_ptm_norm) ) {
             histogram_plot(dpmsr_set$data$norm_data[, -which(names(dpmsr_set$data$norm_data) == 'Detected_Imputed')], "PTM_Only_Intensity_Histogram")
          }

          adjust_intensity_cutoff(session, input, output)

          removeModal()
          inputfilterapply_render(session, input, output)

      }
    
  })

    
#------------------------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$include_exclude, {
    cat(file = stderr(), "norm include/exclude filter triggered", "\n")    
    
    filter_norm()
    
    #recalculate these from above: change to always display norm data
    info_columns <- ncol(dpmsr_set$data$norm_data) - dpmsr_set$y$sample_number
    bar_plot(dpmsr_set$data$norm_data[(info_columns + 1):ncol(dpmsr_set$data$norm_data)],"Normalization_Data", dpmsr_set$file$qc_dir)
    box_plot(dpmsr_set$data$norm_data[(info_columns + 1):ncol(dpmsr_set$data$norm_data)],"Normalization_Data", dpmsr_set$file$qc_dir)         
    
    inputfilterapply_render(session, input, output)
    
  })
  
  
  
  
#------------------------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$norm1, {
    cat(file = stderr(), "normalization triggered", "\n")
    showModal(modalDialog("Normalizing Data...", footer = NULL))
    inputnorm_render(session, input, output)
    
    if (input$checkbox_n1 == 1) {n1 <- 1}else{n1 <- NULL}
    if (input$checkbox_n2 == 1) {n2 <- 2}else{n2 <- NULL}
    if (input$checkbox_n3 == 1) {n3 <- 3}else{n3 <- NULL}
    if (input$checkbox_n4 == 1) {n4 <- 4}else{n4 <- NULL}
    if (input$checkbox_n5 == 1) {n5 <- 5}else{n5 <- NULL}
    if (input$checkbox_n6 == 1) {n6 <- 6}else{n6 <- NULL}
    if (input$checkbox_n7 == 1) {n7 <- 7}else{n7 <- NULL}
    if (input$checkbox_n8 == 1) {n8 <- 8}else{n8 <- NULL}
    if (input$checkbox_n9 == 1) {n9 <- 9}else{n9 <- NULL}
    if (input$checkbox_n10 == 1) {n10 <- 10}else{n10 <- NULL}
    if (input$checkbox_n11 == 1) {n11 <- 11}else{n11 <- NULL}
    if (input$checkbox_n13 == 1) {n13 <- 13}else{n13 <- NULL}
    
    norm_list <-  list("impute" = 99, "impute_b2" = 98, "sl" = n1, "directlfq" = n13, "quantile" = n4,"lr" = n5,"loess" = n6,"vsn" = n7,"ti" = n8,"mi" = n9,"ai" = n10)
    dpmsr_set$y$norm_list <<- norm_list[lapply(norm_list, length) > 0]
    norm_list2 <-  list("tmm" = n2,"sltmm" = n3,"protein" = n11)
    dpmsr_set$y$norm_list2 <<- norm_list2[lapply(norm_list2, length) > 0]
    apply_norm()
    removeModal()
    if (input$checkbox_tmt) {
      tmt_spqc_prep()
      updateTabsetPanel(session, "nlp1", selected = "tp_tmt")
    }else{
    updateTabsetPanel(session, "nlp1", selected = "tp_impute")
    }
  })

#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------  
  observeEvent(input$impute1, {
    cat(file = stderr(), "impute triggered", "\n")
    showModal(modalDialog("Imputing Missing Data...", footer = NULL))  
    apply_impute(session, input, output)
    removeModal()
    # cat(file = stderr(), "Set comp groups...", "\n")
    # showModal(modalDialog("Setting groups...", footer = NULL))  
    # set_comp_groups()
    # removeModal()
    if (dpmsr_set$x$raw_data_input == "Protein" || dpmsr_set$x$final_data_output == "Peptide") {
      cat(file = stderr(), "impute 1.1,  Protein input or Peptide output, skipping rollup, stat prep...", "\n")
      stat_prep_rollup(session, input, output)
      updateTabsetPanel(session, "nlp1", selected = "tp_qc")
    }else{
      cat(file = stderr(), "impute 1.2 ...", "\n")
      updateTabsetPanel(session, "nlp1", selected = "tp_rollup") 
    }
    
  })
 
    #------------------------------------------------------------------------------------------------------   
    #------------------------------------------------------------------------------------------------------  
    observeEvent(input$rollup, {
      cat(file = stderr(), "action rollup triggered", "\n")
      stat_prep_rollup(session, input, output)
      updateTabsetPanel(session, "nlp1", selected = "tp_qc")
      cat(file = stderr(), "action rollup complete", "\n")
    })   
    
    #------------------------------------------------------------------------------------------------------   
    #------------------------------------------------------------------------------------------------------  
    observeEvent(input$adjust_intensity_cutoff, {

      showModal(modalDialog("Calclulate new intensity cutoff...", footer = NULL))  
      adjust_intensity_cutoff(session, input, output) 
      removeModal()     
    
    })
    
    #------------------------------------------------------------------------------------------------------   
    #------------------------------------------------------------------------------------------------------   
    observeEvent(input$calc_protein_spike, {
      cat(file = stderr(), "Calculating protein spike...", "\n")
      showModal(modalDialog("Calculating protein spike...", footer = NULL)) 
      qc_spike_start()
      qc_spike_render(session, input, output)
      removeModal()
    })
  
    
    #------------------------------------------------------------------------------------------------------   
    #------------------------------------------------------------------------------------------------------   
    observeEvent(input$update_qc_spike_levels, {
      dpmsr_set$design <- hot_to_r( input$protein_qc_spike_levels )
    })
    
#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------  
    observeEvent(input$tmt_irs_go, {
      cat(file = stderr(), "TMT IRS started...", "\n")
      showModal(modalDialog("IRS Normalization...", footer = NULL))  
      tmt_spqc_normalize(dpmsr_set$data$normalized$sl)
      removeModal()
    })
    
    
    #------------------------------------------------------------------------------------------------------   
    #------------------------------------------------------------------------------------------------------  
    observeEvent(input$tmt_irs_qc, { 
      cat(file = stderr(), "stat prep...", "\n")
      showModal(modalDialog("Stat prep...", footer = NULL)) 
      stat_prep()
      removeModal()
      cat(file = stderr(), "qc_apply", "\n")
      showModal(modalDialog("QC stats...", footer = NULL)) 
      qc_apply()
      removeModal()
      cat(file = stderr(), "qc - create_plots", "\n")
      showModal(modalDialog("QC plots...", footer = NULL)) 
      create_qc_plots()
      qc_render(session, input, output)
      update_widget_post_processing(session, input, output)
      save_final_nostats()
      removeModal()
      updateTabsetPanel(session, "nlp1", selected = "tp_qc")
    })   

#------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------  

    
    observeEvent(input$protein_select_plots, {     
      showModal(modalDialog("Working...", footer = NULL))  
      create_select_plots()
    
      updateSelectInput(session, "protein_select", choices = dpmsr_set$y$protein_list , selected = "ADH")
      
      inputproteinselect_render(session, input, output)
      
      removeModal()
    })
 
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    
    observeEvent(input$oneprotein_show, { 
      cat(file = stderr(), "observeEvent input$oneprotein_show...", "\n")
      filter_df <- dpmsr_set$data$final[[input$select_oneprotein_norm]]
      
      if (input$oneprotein_accession != "0" ) {
        cat(file = stderr(), str_c("filtering for accession... ", input$oneprotein_accession), "\n")
        filter_df <- subset(filter_df, Accession %in% as.character(input$oneprotein_accession)  )
      }else{
        cat(file = stderr(), str_c("filtering for accession... ", input$oneprotein_accession), "\n")
      }
      
      stats_df <- filter_df[(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number + 1):ncol(filter_df)]
      stats_df <- t(stats_df)
      stats_df <- data.frame(stats_df)
      colnames(stats_df) <- "CV"
      filter_df <- filter_df[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number)]
      plot_dir <- str_c(dpmsr_set$file$output_dir, input$select_oneprotein_norm, "//")
      bar_plot(filter_df, str_c(input$oneprotein_accession), plot_dir)
      
      output$oneprotein_plot_select <- renderImage({
        list(src = str_c(plot_dir, input$oneprotein_accession, "_barplot.png"),  
             contentType = 'image/png', width = 500, height = 400, alt = "this is alt text")
      }, deleteFile = FALSE)
      
      output$oneprotein_stats <- renderRHandsontable({
        rhandsontable(stats_df, colHeaders = NULL, rowHeaderWidth = 300) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) })
      
    })
    
#-------------------------------------------------------------------------------------------------------------      
#-------------------------------------------------------------------------------------------------------------  
    
    
    observeEvent(input$onepeptide_show, { 
      cat(file = stderr(), "observeEvent input$onepeptide_show...", "\n")
      filter_df <- dpmsr_set$data$final[[input$select_onepeptide_norm]]
      
      if (input$onepeptide_sequence != "0" ) {
        filter_df <- subset(filter_df, Sequence %in% as.character(input$onepeptide_sequence)  )
      }
      cat(file = stderr(), "observeEvent input$onepeptide_show...1", "\n")
      stats_df <- filter_df[(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number + 1):ncol(filter_df)]
      stats_df <- t(stats_df)
      stats_df <- data.frame(stats_df)
      colnames(stats_df) <- "CV"
      filter_df <- filter_df[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number)]
      plot_dir <- str_c(dpmsr_set$file$output_dir, input$select_onepeptide_norm, "//")
      cat(file = stderr(), "observeEvent input$onepeptide_show...2", "\n")
      bar_plot(filter_df, str_c(input$onepeptide_sequence), plot_dir)
      
      output$onepeptide_plot_select <- renderImage({
        list(src = str_c(plot_dir, input$onepeptide_sequence, "_barplot.png"),  
             contentType = 'image/png', width = 500, height = 400, alt = "this is alt text")
      }, deleteFile = FALSE)
      
      cat(file = stderr(), "observeEvent input$onepeptide_show...3", "\n")
      output$onepeptide_stats <- renderRHandsontable({
        rhandsontable(stats_df, colHeaders = NULL, rowHeaderWidth = 300) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) })
      
      cat(file = stderr(), "observeEvent input$onepeptide_show...end", "\n")
    })
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    
    
observeEvent(input$data_show, { 
    filter_df <- dpmsr_set$data$final[[input$select_final_data]]

    if (input$data_topn != 0 ) {
      filter_df$sum <- rowSums(filter_df[(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number)])
      filter_df <- filter_df[order(-filter_df$sum),]                      
      filter_df <- filter_df[1:input$data_topn,]
      filter_df$sum <- NULL
    }
    
    if (input$data_accession != "0" ) {
      filter_df <- subset(filter_df, Accession %in% as.character(input$data_accession)  )
    }
      
    if (input$data_description != "0") {
      filter_df <- filter_df[grep(as.character(input$data_description), filter_df$Description), ]
    }
    
    output$data_final <- renderRHandsontable({
      rhandsontable(filter_df, rowHeaders = NULL, columnHeaderwidth = 200, width = 1200, height = 500) %>%
        hot_cols(fixedColumnsLeft = 1, colWidths = 80, halign = "htCenter" ) %>%
        #hot_col(col = "Peptides", format="0") %>%
        hot_col(col = dpmsr_set$y$info_columns_final:(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number), format = "0", digits = 0) %>%
        hot_col(col = (dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number + 1):(ncol(dpmsr_set$data$final$impute)), format = "0.00") %>%
        hot_col(col = "Description", halign = "htLeft", colWidths = 350) %>%
        hot_rows(rowHeights = 10)
    })
  })  
    

    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$check_stats, {
      showModal(modalDialog("Setting Stat groups...", footer = NULL))  
      cat(file = stderr(), "Setting Stat groups...", "\n")
      try(set_stat_groups(session, input, output), silent = TRUE)
      #try(callr::r_bg(set_stat_groups, args = list(session, input, output), supervise = TRUE))
      if (is.null(input$comp_spqc)) {
        shinyalert("Oops!", "Please choose and SPQC group!", type = "error")
      }
      #update all of the dropdown stat group choices that are downstream from this point
      updateTextInput(session, "final_stats_name", value = str_c("Final_", input$select_final_data_stats,  "_stats.xlsx"))
      
      removeModal()
    })  
    
    observe({
      if (input$radio_output == 1) {
        output$stats_gear1 <- renderText({"Peptide Filters"})
        output$stats_gear2 <- renderText({"Protein Filters"})
      }
      else if (input$radio_output == 2) {
        output$stats_gear1 <- renderText({"Peptide Filters"})
        output$stats_gear2 <- renderText({"Protein Filters"})
      }
    
      
    }) 
    

    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$start_stats, {
      showModal(modalDialog("Calculating stats...", footer = NULL))  
      cat(file = stderr(), "Calculating stats...", "\n")
      #test_stat <- callr::r_bg(stat_calc2, args = list(session, input, output), supervise = TRUE)
      stat_calc2(session, input, output)
      removeModal()
    }
    )    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$save_stats, {
      showModal(modalDialog("Saving Final Excel...", footer = NULL))      
      cat(file = stderr(), "Saving Final Excel...", "\n")
      stats_Final_Excel(session, input, output)
      removeModal()
    }
    )   
    
    #-------------------------------------------------------------------------------------------------------------    
    output$download_stats_excel <- downloadHandler(
      file = function(){
        input$final_stats_name
      },
      content = function(file){
        fullname <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats, "//", input$final_stats_name)
        cat(file = stderr(), str_c("download_stats_excel fullname = ", fullname), "\n")
        file.copy(fullname, file)
      }
    )
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    
    
    observeEvent(input$create_stats_plots, {
      showModal(modalDialog("Creating Plots...", footer = NULL))  
      
      if (!is.null(dpmsr_set$data$final[[input$select_final_data_stats ]] )) {
      
        if (!is.null(input$stats_plot_comp)) {
   
            df <- dpmsr_set$data$final[[input$select_final_data_stats ]][(dpmsr_set$y$info_columns_final + 1):(dpmsr_set$y$info_columns_final + dpmsr_set$y$sample_number)]
            
            comp_string <- input$stats_plot_comp
            cat(file = stderr(), "Create stats plots" , "\n")
            for (i in 1:length(comp_string)) {
              comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string[i])
              if (i == 1) {
                comp_rows <- c(dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
              }else{
                comp_rows <- c(comp_rows, dpmsr_set$y$stats$groups$sample_numbers_N[comp_number],dpmsr_set$y$stats$groups$sample_numbers_D[comp_number] )
              }
            }
            
            #add spqc to plots
            cat(file = stderr(), "Stats Plots...1" , "\n")
            if (input$stats_plot_spqc) {
                comp_rows <- c(comp_rows, dpmsr_set$y$stats$comp_spqc_sample_numbers)
            }
            
            cat(file = stderr(), "Stats Plots...2" , "\n")
            comp_rows <- sort(unique(unlist(comp_rows)), decreasing = FALSE)
            df <- df[,comp_rows]
            
            #from protein normalization there will be standard deviation of 0 for the normalized protein - this will crash some of the pca/cluster/heatmap calcs
            df <- df[apply(df, 1, var) != 0, ]
            
            namex <- dpmsr_set$design$Label[comp_rows]
            color_list <- dpmsr_set$design$colorlist[comp_rows]
            if (dpmsr_set$x$primary_group) {
              groupx <- dpmsr_set$design$PrimaryGroup[comp_rows]
            }else{
              groupx <- dpmsr_set$design$Group[comp_rows]
            }
            
            cat(file = stderr(), "Stats Plots...3" , "\n")
            interactive_barplot(session, input, output, df, namex, color_list, "stats_barplot", input$stats_plot_comp)
            
            cat(file = stderr(), "Stats Plots...4" , "\n")
            interactive_boxplot(session, input, output, df, namex, color_list, input$stats_plot_comp)
            
            cat(file = stderr(), "Stats Plots...5" , "\n")
            interactive_pca2d(session, input, output, df, namex, color_list, groupx, input$stats_plot_comp)
            
            cat(file = stderr(), "Stats Plots...6" , "\n")
            interactive_pca3d(session, input, output, df, namex, color_list, groupx, input$stats_plot_comp)
            
            cat(file = stderr(), "Stats Plots...7" , "\n")
            interactive_cluster(session, input, output, df, namex, input$stats_plot_comp)
            
            cat(file = stderr(), "Stats Plots...8" , "\n")
            interactive_heatmap(session, input, output, df, namex, groupx, input$stats_plot_comp)
        
            
        }else{
          shinyalert("Oops!", "Please select comparison", type = "error")
        }
        
        
      }else{
        shinyalert("Oops!", "Data does not exist yet.  Did you impute?", type = "error")
      }
      
      removeModal()
      cat(file = stderr(), "Stats Plots...end" , "\n")
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
      
      if (!is.null(dpmsr_set$data$stats[[dpmsr_set$y$stats$groups$comp_name[1]]] )) {
      
        for (i in 1:input$comp_number) {
          do.call("interactive_stats_volcano", list(session, input, output, i) )
        }
      }else{
        shinyalert("Oops!", "Data does not exist yet.  Did you run stats?", type = "error")
      }
      removeModal()
    }
    )  
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  

    observeEvent(input$stats_data_show, { 
      showModal(modalDialog("Getting data...", footer = NULL))  
      cat(file = stderr(), "stats data show triggered..." , "\n")
      
      if (!is.null(dpmsr_set$data$stats[[dpmsr_set$y$stats$groups$comp_name[1]]] )  &  input$stats_select_data_comp != "Choice 1") {
        
        filter_df <- stats_data_table_filter(session, input, output)
        #test_df <<- filter_df
        
        filter_df_colnames <- colnames(filter_df)
        filter_df_colnames <- gsub("_v_", " v ", filter_df_colnames)
        filter_df_colnames <- gsub("_FC", " FC", filter_df_colnames)
        filter_df_colnames <- gsub("_CV", " CV", filter_df_colnames)
        filter_df_colnames <- gsub("_MF", " MF", filter_df_colnames)
        filter_df_colnames <- gsub("_pval", " pval", filter_df_colnames)
        filter_df_colnames <- gsub("_limmapval", " Limma pval", filter_df_colnames)
        filter_df_colnames <- gsub("_cohensd", " CohensD", filter_df_colnames)
        filter_df_colnames <- gsub("_adjpval", " adjpval", filter_df_colnames)
        filter_df_colnames <- gsub("_", ".", filter_df_colnames)
        colnames(filter_df) <-  filter_df_colnames

        if (dpmsr_set$x$final_data_output == "Protein") {
          stats_DT <- protein_table(session, input, output, filter_df)
        }else{
          stats_DT <- peptide_table(session, input, output, filter_df)
        }
  
        output$stats_data_final <-  DT::renderDataTable(stats_DT, selection = 'single' )
        #save data
        dpmsr_set$data$stats_DT <<- filter_df
        
        
        # get selections from data table for protein or peptide formats
        if (dpmsr_set$x$final_data_output == "Protein") {
            output$stats_data_final_protein <- renderPrint(stats_DT$x$data$Accession[as.numeric(unlist(input$stats_data_final_rows_selected)[1])] )
            observe({
              updateTextInput(session, "stats_oneprotein_accession", 
                              value = stats_DT$x$data$Accession[as.numeric(unlist(input$stats_data_final_rows_selected)[1])]  )
              updateSelectInput(session, "stats_oneprotein_plot_comp", selected = input$stats_select_data_comp)
              dpmsr_set$y$stats$accession_stat <<- stats_DT$x$data$Accession[as.numeric(unlist(input$stats_data_final_rows_selected)[1])] 
            })
        }else{
          output$stats_data_final_protein <- renderPrint(str_c(
            stats_DT$x$data$Accession[as.numeric(unlist(input$stats_data_final_rows_selected)[1])], "  ", 
            stats_DT$x$data$Sequence[as.numeric(unlist(input$stats_data_final_rows_selected)[1])], "  ", 
            stats_DT$x$data$Modification[as.numeric(unlist(input$stats_data_final_rows_selected)[1])] 
            ))
          observe({
            updateTextInput(session, "stats_onepeptide_accession", 
                            value = stats_DT$x$data$Accession[as.numeric(unlist(input$stats_data_final_rows_selected)[1])]  )
            updateTextInput(session, "stats_onepeptide_sequence", 
                            value = stats_DT$x$data$Sequence[as.numeric(unlist(input$stats_data_final_rows_selected)[1])] ) 
            updateTextInput(session, "stats_onepeptide_modification", 
                            value = stats_DT$x$data$Modification[as.numeric(unlist(input$stats_data_final_rows_selected)[1])] ) 
            updateSelectInput(session, "stats_onepeptide_plot_comp", selected = input$stats_select_data_comp)
            dpmsr_set$y$stats$accession_stat <<- stats_DT$x$data$Accession[as.numeric(unlist(input$stats_data_final_rows_selected)[1])] 
            dpmsr_set$y$stats$sequence_stat <<- stats_DT$x$data$Sequence[as.numeric(unlist(input$stats_data_final_rows_selected)[1])] 
            dpmsr_set$y$stats$modification_stat <<- stats_DT$x$data$Modifications[as.numeric(unlist(input$stats_data_final_rows_selected)[1])]
            #protein plot is still used for peptide output
            updateTextInput(session, "stats_oneprotein_accession", 
                            value = stats_DT$x$data$Accession[as.numeric(unlist(input$stats_data_final_rows_selected)[1])]  )
            updateSelectInput(session, "stats_oneprotein_plot_comp", selected = input$stats_select_data_comp)
          })
        }  
      }else{
        shinyalert("Oops!", "Data does not exist yet.  Did you run stats?", type = "error")
      }  
        
    removeModal()
    })  
    

    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    observeEvent(input$stats_data_update, {
      
      if (length(dpmsr_set$y$stats$accession_stat) != 0) {  
      
        if(dpmsr_set$x$final_data_output == "Protein"){  
          cat(file = stderr(), str_c("accesion to mark as not stat signifigant...   ", dpmsr_set$y$stats$accession_stat) , "\n") 
          row_number <- which(dpmsr_set$data$stats[[input$stats_select_data_comp]]$Accession == dpmsr_set$y$stats$accession_stat)
          cat(file = stderr(), str_c("row to mark as not stat signifigant...   ", row_number) , "\n") 
          dpmsr_set$data$stats[[input$stats_select_data_comp]]$Stats[row_number] <<- ""
        }else{
          cat(file = stderr(), str_c("peptide to mark as not stat signifigant...   ", dpmsr_set$y$stats$sequence_stat) , "\n") 
          row_number <- which(dpmsr_set$data$stats[[input$stats_select_data_comp]]$Sequence   == dpmsr_set$y$stats$sequence_stat &
                                dpmsr_set$data$stats[[input$stats_select_data_comp]]$Modifications   == dpmsr_set$y$stats$modification_stat)
          cat(file = stderr(), str_c("row to mark as not stat signifigant...   ", row_number) , "\n") 
          dpmsr_set$data$stats[[input$stats_select_data_comp]]$Stats[row_number] <<- ""
        }
      }else{
        
        shinyalert("Oops!", "Cannot update stat", type = "error")
        
      }
      
    })
    


    
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    observeEvent(input$stats_data_save, { 

        showModal(modalDialog("Saving Data...", footer = NULL))  
        cat(file = stderr(), "stats saving datatable to excel..." , "\n") 
        #Convert to R object
        #x <- hot_to_r(isolate(input$stats_data_final))
        #x <- stats_data_table_filter(session, input, output)
        filename <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats, "//", input$stats_data_filename)
        file_dir <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats) 
        
        if(!is_dir(file_dir)) {
          cat(file = stderr(), str_c("create_dir...", file_dir), "\n")
          dir_create(file_dir)
        }
        
        cat(file = stderr(), str_c("filename = ", filename) , "\n")
        Simple_Excel(dpmsr_set$data$stats_DT, "data", filename)
        removeModal()
    
    })

    
    #-------------------------------------------------------------------------------------------------------------    
    output$download_stats_data_save <- downloadHandler(
      file = function(){
        input$stats_data_filename
      },
      content = function(file){
        fullname <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats, "//", input$stats_data_filename)
        cat(file = stderr(), str_c("download_stats_data fullname = ", fullname), "\n")
        file.copy(fullname, file)
      }
    )
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    
    
    
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    observeEvent(input$create_stats_oneprotein_plots, { 
      cat(file = stderr(), str_c("Create stats oneprotein plots"), "\n")
      
      comp_number <- try(which(dpmsr_set$data$stats[[input$stats_oneprotein_plot_comp]] == input$stats_oneprotein_accession), silent = TRUE)
        
      if (length(comp_number) != 0) {  

        # do not run peptide section if TMT SPQC Norm
        if (!dpmsr_set$x$tmt_spqc_norm) {
          
          df_list <- oneprotein_data(session, input, output)
          #
          #test_df_list <<- df_list
          
          for (j in names(df_list)) {assign(j, df_list[[j]]) }
          
          interactive_barplot(session, input, output, df, namex, color_list, "stats_oneprotein_barplot", input$stats_oneprotein_plot_comp)
          
          peptide_pos_lookup <-  peptide_position_lookup(session, input, output, as.character(input$stats_oneprotein_accession))
          
          grouped_color <- unique(color_list)
          interactive_grouped_barplot(session, input, output, comp_string, df_peptide, peptide_info_columns, 
                                      input$stats_oneprotein_plot_comp, peptide_pos_lookup, grouped_color)
          
          #test_df_peptide <<- df_peptide
          #test_peptide_pos_lookup <<- peptide_pos_lookup
          
          df_peptide <- merge(df_peptide, peptide_pos_lookup, by = (c("Accession", "Sequence"))    )
          df_peptide$Start <- as.numeric(df_peptide$Start)
          df_peptide$Stop <- as.numeric(df_peptide$Stop)
          df_peptide <- df_peptide %>% dplyr::select(Stop, everything())
          df_peptide <- df_peptide %>% dplyr::select(Start, everything())
          df_peptide <- df_peptide[order(df_peptide$Start, df_peptide$Stop), ]
          
          sample_col_numbers <- seq(from = 14, to = ncol(df_peptide) )
          
          dpmsr_set$data$oneprotein_peptide_DT <<- df_peptide
          
          oneprotein_peptide_DT <-  DT::datatable(df_peptide,
                                     rownames = FALSE,
                                     extensions = c("FixedColumns"), #, "Buttons"),
                                     options = list(
                                       #dom = 'Bfrtipl',
                                       autoWidth = TRUE,
                                       scrollX = TRUE,
                                       scrollY = 500,
                                       scrollCollapse = TRUE,
                                       columnDefs = list(list(targets = c(0,1), visibile = TRUE, "width" = '30', className = 'dt-center'),
                                                         list(targets = c(2), visible = TRUE, "width" = '20', className = 'dt-center'),
                                                         list(
                                                           targets = c(5),
                                                           width = '250',
                                                           render = JS(
                                                             "function(data, type, row, meta) {",
                                                             "return type === 'display' && data.length > 35 ?",
                                                             "'<span title=\"' + data + '\">' + data.substr(0, 35) + '...</span>' : data;",
                                                             "}")
                                                         ),
                                                         list(
                                                           targets = c(6),
                                                           width = '150',
                                                           render = JS(
                                                             "function(data, type, row, meta) {",
                                                             "return type === 'display' && data.length > 35 ?",
                                                             "'<span title=\"' + data + '\">' + data.substr(0, 35) + '...</span>' : data;",
                                                             "}")
                                                         ),
                                                         list(
                                                           targets = c(12),
                                                           width = '100',
                                                           render = JS(
                                                             "function(data, type, row, meta) {",
                                                             "return type === 'display' && data.length > 20 ?",
                                                             "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                             "}")
                                                         )
                                       ),
                                       ordering = TRUE,
                                       orderClasses = TRUE,
                                       fixedColumns = list(leftColumns = 2),
                                       pageLength = 100, lengthMenu = c(10,50,100,200)),
                                     #buttons=c('copy', 'csv', 'excelHtml5', 'pdf')),
                                     callback = JS('table.page(3).draw(false);'
                                     ))
          
          oneprotein_peptide_DT <- oneprotein_peptide_DT %>%  formatRound(columns = c(sample_col_numbers), digits = 2)
          
          output$oneprotein_peptide_table <-  DT::renderDataTable({oneprotein_peptide_DT })
        }else {
          
          df_list <- TMT_IRS_protein_data(session, input, output)
          for (j in names(df_list)) {assign(j, df_list[[j]]) }

          interactive_barplot(session, input, output, df, namex, color_list, "stats_oneprotein_barplot", input$stats_oneprotein_plot_comp)

        } #end of TMT SPQC if 
        
      }else{
        shinyalert("Oops!", "No Accession...", type = "error")
      }
      removeModal()

    })
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    
    observeEvent(input$create_stats_onepeptide_plots, { 
      cat(file = stderr(), "create_stats_onepeptide_plots...", "\n")
      
      comp_test <- try(which(dpmsr_set$data$stats[[input$stats_onepeptide_plot_comp]] == input$stats_onepeptide_accession), silent =TRUE)
      comp_string <- input$stats_onepeptide_plot_comp
      comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string)
      
      cat(file = stderr(), "create_stats_onepeptide_plots...1", "\n")
      cat(file = stderr(), str_c("comp_number = ", comp_number, " -v- ",  length(dpmsr_set$y$stats$groups$comp_name)  ) , "\n")
      check_qc <- TRUE
      
      #test for all.samples v qc and adding qc
      # if(input$stats_onepeptide_plot_spqc & comp_number==length(dpmsr_set$y$stats$groups$comp_name)) {
      #   cat(file = stderr(), "Reset SPQC", "\n")
      #   updateCheckboxInput(session, "stats_onepeptide_plot_spqc",  value = FALSE)
      #   check_qc <- FALSE
      # }
      
      #test for SPQC group already in comparison, if so turn off and warn
      cat(file = stderr(), "create_stats_onepeptide_plots...2", "\n")
      if(input$stats_onepeptide_plot_spqc & grepl(dpmsr_set$y$stats$comp_spqc, dpmsr_set$y$stats$groups$comp_name[comp_number]) ) {
        cat(file = stderr(), "Reset SPQC", "\n")
        updateCheckboxInput(session, "stats_onepeptide_plot_spqc",  value = FALSE)
        check_qc <- FALSE
      }
      
      cat(file = stderr(), "create_stats_onepeptide_plots...3", "\n")
      if (length(comp_test)!=0 & check_qc){  
        
        df_list <- onepeptide_data(session, input, output)
        #df_peptide in df_list
        for(j in names(df_list)){assign(j, df_list[[j]]) }
        
        cat(file = stderr(), "create_stats_onepeptide_plots...4", "\n")
        interactive_barplot(session, input, output, df, namex, color_list, "stats_onepeptide_barplot", comp_string)
        
        peptide_pos_lookup <-  peptide_position_lookup(session, input, output, as.character(input$stats_onepeptide_accession))
        grouped_color <- unique(color_list)
        #interactive_grouped_peptide_barplot(session, input, output, comp_string, df_peptide, info_columns, comp_name, peptide_pos_lookup, grouped_color)
        interactive_grouped_peptide_barplot(session, input, output, comp_string, df_peptide, info_columns, input$stats_onepeptide_plot_comp, peptide_pos_lookup, grouped_color)
        
        cat(file = stderr(), "create_stats_onepeptide_plots...5", "\n")
        sample_col_numbers <- seq(from=12, to = ncol(df_peptide) )
        df_peptide <- cbind(df_peptide, df_peptide_stats)
        
        df_peptide <- merge(df_peptide, peptide_pos_lookup, by=(c("Accession", "Sequence"))    )
        df_peptide$Start <- as.numeric(df_peptide$Start)
        df_peptide$Stop <- as.numeric(df_peptide$Stop)
        df_peptide<- df_peptide %>% dplyr::select(Stop, everything())
        df_peptide <- df_peptide %>% dplyr::select(Start, everything())
        df_peptide <- df_peptide[order(df_peptide$Start, df_peptide$Stop), ]
        
        cat(file = stderr(), "create_stats_onepeptide_plots...6", "\n")
        onepeptide_peptide_DT <-  DT::datatable(df_peptide,
                                                rownames = FALSE,
                                                extensions = c("FixedColumns"), #, "Buttons"),
                                                options=list(
                                                  #dom = 'Bfrtipl',
                                                  autoWidth = TRUE,
                                                  scrollX = TRUE,
                                                  scrollY=500,
                                                  scrollCollapse=TRUE,
                                                  columnDefs = list(list(targets = c(0,1), visibile = TRUE, "width"='30', className = 'dt-center'),
                                                                    list(targets = c(2), visible = TRUE, "width"='20', className = 'dt-center'),
                                                                    list(
                                                                      targets = c(5),
                                                                      width = '250',
                                                                      render = JS(
                                                                        "function(data, type, row, meta) {",
                                                                        "return type === 'display' && data.length > 35 ?",
                                                                        "'<span title=\"' + data + '\">' + data.substr(0, 35) + '...</span>' : data;",
                                                                        "}")
                                                                    ),
                                                                    list(
                                                                      targets = c(6),
                                                                      width = '150',
                                                                      render = JS(
                                                                        "function(data, type, row, meta) {",
                                                                        "return type === 'display' && data.length > 35 ?",
                                                                        "'<span title=\"' + data + '\">' + data.substr(0, 35) + '...</span>' : data;",
                                                                        "}")
                                                                    ),
                                                                    list(
                                                                      targets = c(10),
                                                                      width = '100',
                                                                      render = JS(
                                                                        "function(data, type, row, meta) {",
                                                                        "return type === 'display' && data.length > 20 ?",
                                                                        "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                                        "}")
                                                                    )
                                                  ),
                                                  ordering = TRUE,
                                                  orderClasses= TRUE,
                                                  fixedColumns = list(leftColumns = 2),
                                                  pageLength = 100, lengthMenu = c(10,50,100,200)),
                                                #buttons=c('copy', 'csv', 'excelHtml5', 'pdf')),
                                                callback = JS('table.page(3).draw(false);'
                                                ))
        
        onepeptide_peptide_DT <- onepeptide_peptide_DT %>%  formatRound(columns=c(sample_col_numbers), digits=2)
        
        output$onepeptide_peptide_table<-  DT::renderDataTable({onepeptide_peptide_DT })
        
        
      }else{
        shinyalert("Oops!", "No Accession or SPQC add error", type = "error")
      }
      removeModal()
      cat(file = stderr(), "create_stats_onepeptide_plots...end", "\n")
      
    })
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  

    observeEvent(input$XXXXXX, { 
      
      showModal(modalDialog("Saving OneProtein Data...", footer = NULL))  
      cat(file = stderr(), "stats saving OneProtein datatable to excel..." , "\n") 
      
      df_list <- oneprotein_data(session, input, output)
      for(j in names(df_list)){
        assign(j, df_list[[j]]) 
        }
      
      filename <- str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp, "//", input$stats_oneprotein_data_filename)
      Simple_Excel(df_peptide, "data", filename)
      removeModal()
      
    })
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------  
    observeEvent(input$stats_oneprotein_data_save, { 
      
      showModal(modalDialog("Saving Data OneProtein...", footer = NULL))  
      cat(file = stderr(), "stats saving OneProtein datatable to excel..." , "\n") 
      
      filename <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats, "//", input$stats_oneprotein_data_filename)
      file_dir <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats) 

      if(!is_dir(file_dir)) {
        cat(file = stderr(), str_c("create_dir...", file_dir), "\n")
        dir_create(file_dir)
      }
      
      
      cat(file = stderr(), str_c("filename = ", filename) , "\n")
      df <- dpmsr_set$data$oneprotein_peptide_DT
      x <- unlist(df$Detected_Imputed)
      df$Detected_Imputed <- x
      Simple_Excel(df, "data", filename)
      removeModal()
      
    })


    #-------------------------------------------------------------------------------------------------------------    
    output$download_stats_oneprotein_data_save <- downloadHandler(
      file = function(){
        input$stats_oneprotein_data_filename
      },
      content = function(file){
        fullname <- str_c(dpmsr_set$file$output_dir, input$select_final_data_stats, "//", input$stats_oneprotein_data_filename)
        cat(file = stderr(), str_c("download_stats_data_oneprotein fullname = ", fullname), "\n")
        file.copy(fullname, file)
      }
    )
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$motif_show, {
      
      shinyalert("Hi", "Motif-X table will appear when analysis is complete.", type = "info")
      
      comp_string <- input$select_data_comp_motif
      comp_number <- which(dpmsr_set$y$stats$groups$comp_name == comp_string)
      
      filter_df <- dpmsr_set$data$stats[[comp_string]]
      motif_data <- run_motifx(input, output, filter_df)
      
      if (!is.null(motif_data)){
        output$motif_table<- renderRHandsontable({
          rhandsontable(motif_data, rowHeaders = NULL) %>%
            hot_cols(colWidths = 80, halign = "htCenter" ) %>%
            hot_col(col = "comparison", halign = "htCenter", colWidths = 150) %>%
            hot_col(col = "motif", halign = "htCenter", colWidths = 100) %>%
            hot_col(col = "fold.increase", halign = "htCenter", colWidths = 100)
        })
      }
      
    })

#      observe({
#        if (!is.null(input$motif_fasta)){
#          fasta_txt_file <- input$motif_fasta
#          fasta_txt_file <- unlist(fasta_txt_file$files[[as.character(0)]][2])
#          output$fasta <- renderText(fasta_txt_file)
#        }
#      })
      
 
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$parse_fasta, {
      showModal(modalDialog("Create sequence file for MotifX...", footer = NULL))  
      create_phos_database(session, input, output)
      
      output$fasta_example<- renderRHandsontable({
        rhandsontable(dpmsr_set$data$phos$background[1:10,], rowHeaders = NULL) %>%
          hot_cols(colWidths = 80, halign = "htCenter" ) %>%
          hot_col(col = "Accession", halign = "htCenter", colWidths = 200) %>%
          hot_col(col = "Sequence", halign = "htLeft", colWidths = 800)
      })
      
      removeModal()
    }
    )
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    
    output$fasta <- renderText({
      if (req(typeof(input$motif_fasta)=="list")) {
        motif_path <- parseFilePaths(volumes, input$motif_fasta)
        basename(motif_path$datapath)
        #unlist(input$motif_fasta$files[[as.character(0)]][2])
      }
    })
    
    #-------------------------------------------------------------------------------------------------------------    
    output$download_motifx_excel <- downloadHandler(
      file = function(){
        dpmsr_set$data$phos$filename
        #str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp, "//", input$final_stats_name)
      },
      content = function(file){
        fullName <- str_c(dpmsr_set$file$phos, dpmsr_set$data$phos$filename)
        file.copy(fullName, file)
      }
    )
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$set_pathway, {
      showModal(modalDialog("Downloading and Setting up databases...", footer = NULL))  
      set_pathway(input, output, session)
      removeModal()
      updateTabsetPanel(session, "nlp1", selected = "wiki")
      updateNavbarPage(session, "path", selected = "wiki")
    }
    )
    
    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$wiki_show, {
      showModal(modalDialog("Wiki Pathway Enrichment...", footer = NULL))  
      
      wiki_data <<- try(run_wiki(session, input, output), silent = TRUE)

      if ( class(wiki_data) != "try-error"){
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
      }else{
        shinyalert("Oops!", "Wiki Pathway enrichment failed due to insufficient gene IDs mapping to pathways", type = "error")
      }
      
      removeModal()
      
      
      fullName <- str_c(dpmsr_set$file$string, "Wiki_", input$select_data_comp_wiki, "_", 
                        input$select_ont_wiki,input$select_level_wiki, ".xlsx", collapse = " ")
      output$download_wiki_table <- downloadHandler(
        filename = function(){
          str_c("Wiki_", input$select_data_comp_wiki, "_", 
                input$select_ont_wiki,input$select_level_wiki, ".xlsx", collapse = " ")
        },
        content = function(file){
          file.copy(fullName, file)
        }
      )
      
      
      
    }
    )
    

    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$profile_show, {
      showModal(modalDialog("Creating Go Profile...", footer = NULL))  
      
      profile_data <- try(run_profile(session, input, output), silent=TRUE)
      
      if ( class(profile_data) != "try-error"){
      
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
        
        #diret plot of profile
         output$go_profile_plot <- renderPlot({
           barplot(profile_data, title = str_c("Go Profile ", input$select_ont_profile, " level=", input$select_level_profile), 
                   drop=TRUE, showCategory=12, order=TRUE)
         })
        
      }else{
        shinyalert("Oops!", "Go Profile enrichment failed...", type = "error")
      }
      
      removeModal()
      
      
      table_fullName <- str_c(dpmsr_set$file$string, "GO_Profile_", input$select_data_comp_profile, "_", 
                              input$select_ont_profile,input$select_level_profile, ".xlsx", collapse = " ")
      output$download_go_profile_table <- downloadHandler(
        filename = function(){
          str_c("GO_Profile_", input$select_data_comp_profile, "_", 
                input$select_ont_profile,input$select_level_profile, ".xlsx", collapse = " ")
        },
        content = function(file){
          file.copy(table_fullName, file)
        }
      )
      
      plot_fullName <- str_c(dpmsr_set$file$string, "GO_Profile_", input$select_data_comp_profile, "_", 
                             input$select_ont_profile,input$select_level_profile, ".png", collapse = " ")
      output$download_go_profile_plot <- downloadHandler(
        filename = function(){
          str_c("GO_Profile_", input$select_data_comp_profile, "_", 
                input$select_ont_profile,input$select_level_profile, ".png", collapse = " ")
        },
        content = function(file){
          file.copy(plot_fullName, file)
        }
      )
      
      
      
      
    }
    )
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_show, {
      showModal(modalDialog("Go Enrichment Analysis...", footer = NULL))  
      
      go_data <- try(run_go_analysis(session, input, output), silent = TRUE)
      #go_data <- goresult
      
      if ( class(go_data) != "try-error") {
        dpmsr_set$data$pathway$mergedf <<- setup_go_volcano(session, input, output)
        #go_data <- try(go_data[order(go_data$pvalue),], silent = TRUE)
        
        output$go_table <- renderRHandsontable({
          rhandsontable(go_data, rowHeaders = NULL, readOnly = FALSE) %>%
            hot_cols(colWidths = 80, halign = "htCenter" ) %>%
            hot_col(col = "GO.ID", halign = "htCenter", colWidths = 100) %>%
            hot_col(col = "Term", halign = "htCenter", colWidths = 200) %>%
            hot_col(col = "Rank in classicFisher", halign = "htCenter", colWidths = 200)   %>%
            hot_col(col = "classicFisher", halign = "htCenter", colWidths = 100)  
        })
      }else{
        shinyalert("Oops!", "Go Analysis failed...", type = "error")
      } 
        
      fullName <- str_c(dpmsr_set$file$string, "GO_", input$select_data_comp_go, "_", 
                        input$select_ont_go, ".xlsx", collapse = " ")
      output$download_go_table <- downloadHandler(
        filename = function(){
          str_c("GO_", input$select_data_comp_go, "_", 
                input$select_ont_go, ".xlsx", collapse = " ")
        },
        content = function(file){
          file.copy(fullName, file)
        }
      )
      
      removeModal()
      
    }
    )   
    
    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_volcano, {

      if (input$go_volcano_id != ""){
      
          showModal(modalDialog("Preparing Go Volcano...", footer = NULL))  
          
          comp_string <- input$select_data_comp_go
          comp_number <- which(grepl(comp_string, dpmsr_set$y$stats$groups$comp_name))
          data_in <- dpmsr_set$data$stats[[comp_string]]
          
          volcano_data <- interactive_go_volcano(session, input, output)
          volcano_go_data <- subset(data_in, Accession %in% volcano_data$Accession  )
          
          volcano_go_data_colnames <- colnames(volcano_go_data )
          volcano_go_data_colnames <- gsub("_v_", " v ", volcano_go_data_colnames)
          volcano_go_data_colnames <- gsub("_FC", " FC", volcano_go_data_colnames)
          volcano_go_data_colnames <- gsub("_CV", " CV", volcano_go_data_colnames)
          volcano_go_data_colnames <- gsub("_MF", " MF", volcano_go_data_colnames)
          volcano_go_data_colnames <- gsub("_pval", " pval", volcano_go_data_colnames)
          volcano_go_data_colnames <- gsub("_", ".", volcano_go_data_colnames)
          colnames(volcano_go_data) <-  volcano_go_data_colnames
          
          
          pval_cols <- colnames(volcano_go_data %>% dplyr::select(contains("pval") ) )
          sample_cols <- c(colnames(volcano_go_data %>% dplyr::select(contains("Normalized"))),
                           colnames(volcano_go_data %>% dplyr::select(contains("Imputed"))) )
          sample_col_numbers <- list(match(sample_cols, names(volcano_go_data)))
          sample_col_numbers <- unlist(sample_col_numbers)
          cv_cols <- colnames(volcano_go_data %>% dplyr::select(contains("CV") ) )
          mf_cols <- colnames(volcano_go_data %>% dplyr::select(contains("MF") ) )
          
          volcano_DT <-  DT::datatable(volcano_go_data,
                                     rownames = FALSE,
                                     extensions = c("FixedColumns"), #, "Buttons"),
                                     options=list(
                                       #dom = 'Bfrtipl',
                                       autoWidth = TRUE,
                                       scrollX = TRUE,
                                       scrollY=500,
                                       scrollCollapse=TRUE,
                                       columnDefs = list(list(targets = c(0), visibile = TRUE, "width"='30', className = 'dt-center'),
                                                         list(targets = c(2), visible = TRUE, "width"='20', className = 'dt-center'),
                                                         list(
                                                           targets = c(1),
                                                           width = '250',
                                                           render = JS(
                                                             "function(data, type, row, meta) {",
                                                             "return type === 'display' && data.length > 35 ?",
                                                             "'<span title=\"' + data + '\">' + data.substr(0, 35) + '...</span>' : data;",
                                                             "}")
                                                         ),
                                                         list(
                                                           targets = c(3),
                                                           width = '100',
                                                           render = JS(
                                                             "function(data, type, row, meta) {",
                                                             "return type === 'display' && data.length > 20 ?",
                                                             "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                             "}")
                                                         )
                                       ),
                                       ordering = TRUE,
                                       orderClasses= TRUE,
                                       fixedColumns = list(leftColumns = 1),
                                       pageLength = 100, lengthMenu = c(10,50,100,200)),
                                     #buttons=c('copy', 'csv', 'excelHtml5', 'pdf')),
                                     callback = JS('table.page(3).draw(false);'
                                     ))
          
          volcano_DT <- volcano_DT %>%  formatRound(columns=c(sample_col_numbers), digits=0)
          
          output$volcano_data_final <-  DT::renderDataTable({volcano_DT })
          
          Simple_Excel(volcano_go_data, "GoVolcano", str_c(dpmsr_set$file$string, "GoVolcano_", input$select_data_comp_go, "_", input$go_volcano_id, "_", 
                                     input$select_ont_go, ".xlsx", collapse = " "))
          
          
          fullName <- str_c(dpmsr_set$file$string, "GoVolcano_", input$select_data_comp_go, "_", input$go_volcano_id, "_", 
                            input$select_ont_go, ".xlsx", collapse = " ")
          output$download_go_volcano_table <- downloadHandler(
            filename = function(){
              str_c("GoVolcano_", input$select_data_comp_go, "_", input$go_volcano_id, "_", 
                    input$select_ont_go, ".xlsx", collapse = " ")
            },
            content = function(file){
              file.copy(fullName, file)
            }
          )
          
          
          removeModal()
          
      }else{
        shinyalert("Oops!", "GO ID is missing!", type = "error")
      }
    }
    )  

    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_string, {
      
      if(length(dpmsr_set$string) > dpmsr_set$x$comp_number) {
          cat(file = stderr(), "go string triggered", "\n")
          showModal(modalDialog("String Analysis...", footer = NULL))  
          
          string_list <- run_string(session, input, output)
          
          cat(file = stderr(), "string list complete, create plot", "\n")
          
          output$string_plot <- renderImage({
            list(src=string_list$string_file_name,  
                 contentType = 'image/png', width=1200, height=1200, alt="this is alt text")
          }, deleteFile = FALSE)
          
    
          fullName <- string_list$string_file_name
          output$download_string_plot <- downloadHandler(
            filename = function(){
              str_c(input$select_data_comp_string, ".png")
            },
            content = function(file){
              file.copy(fullName, file)
            }
          )
          
          # depreciated
           cat(file = stderr(), "create string link", "\n")
           #output$string_link <- renderText({dpmsr_set$string$link_network})
           #output$string_link <- renderText({string_list$linkthis})
           url <- a(dpmsr_set$string$link_network, href= dpmsr_set$string$link_network, target="_blank")
           output$string_link <- renderUI({
             tagList("URL link: ", url)
             })
           
          removeModal()
          
      }else{
        shinyalert("Oops!", "Need to run String Setup first...", type = "error")
      }    
    }
    ) 
    
    
    

    
    #-------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    observeEvent(input$go_string_enrich, {
      
      if(!is.null("dpmsr_set$string")) {
          showModal(modalDialog("Working...", footer = NULL))  
          string_result <- run_string_enrich(input, output)
          
          string_result <- dplyr::rename(string_result, genes = number_of_genes)
          string_result <- dplyr::rename(string_result, genes_background = number_of_genes_in_background)
          
          output$string_table<- renderRHandsontable({
            rhandsontable(string_result, rowHeaders = NULL, readOnly = TRUE) %>%
              hot_cols(colWidths = 80, halign = "htCenter" ) %>%
              hot_col(col = "term", halign = "htCenter", colWidths = 100) %>%
              hot_col(col = "genes_background", halign = "htCenter", colWidths = 100) %>%
              hot_col(col = "preferredNames",  halign = "htLeft", colWidths = 200) %>%
              hot_col(col = "inputGenes", halign = "htCenter", colWidths = 200) %>%
              hot_col(col = "p_value", halign = "htCenter", colWidths = 180, format = '0.0000000000000') %>% 
              hot_col(col = "fdr", halign = "htCenter", colWidths = 180, format = '0.0000000000000') %>% 
              hot_col(col = "description",  halign = "htLeft", colWidths = 400)
          })
          removeModal()
          
    
          
          fullName <- str_c(dpmsr_set$file$string, input$select_data_comp_string_enrich, "_", input$select_string_enrich, ".xlsx", collapse = " ")
          Simple_Excel_bg(string_result, "String", fullName)
          
          output$download_string_enrich_table <- downloadHandler(
            filename = function(){
              str_c(input$select_data_comp_string_enrich, "_", input$select_string_enrich, ".xlsx", collapse = " ")
            },
            content = function(file){
              file.copy(fullName, file)
            }
          )
      
      
      }else{
        
        shinyalert("Oops!", "Need to run String Setup first...", type = "error")
      }
      
    }
    )
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------   
    observeEvent(input$save_dpmsr_set, { 
      cat(file = stderr(), "save dpmsr_set triggered", "\n")
      
      if(!is.null(dpmsr_set$x$pathway_set) || !is.null(dpmsr_set$x$motif_set) ){
        if (dpmsr_set$x$pathway_set == 1 || dpmsr_set$x$motif_set == 1){
          shinyalert("FYI!", "Pathway or MotifX setup Complete!", type = "info")
          }else{
          shinyalert("FYI!", "Pathway or MotifX setup has not been run yet...", type = "info")
          }
        }else{
          shinyalert("FYI!", "Pathway or MotifX setup has not been run yet...", type = "info")
        }
      
      filename <- input$dpmsr_set_name
      future({ save(dpmsr_set, file=str_c(dpmsr_set$file$output_dir, filename, ".dpmsr_set")) })
   })  
    

    output$download_new_dpmsr_set <- downloadHandler(
      filename = function(){
        str_c(input$dpmsr_set_name, ".dpmsr_set")
      },
      content = function(file){
        fullName <- str_c(dpmsr_set$file$output_dir, input$dpmsr_set_name, ".dpmsr_set")
        file.copy(fullName, file)
      }
    )
    
    #-------------------------------------------------------------------------------------------------------------      
    #-------------------------------------------------------------------------------------------------------------   
    observeEvent(input$save_customer_dpmsr_set, { 
      cat(file = stderr(), "save customer dpmsr_set triggered", "\n")
      dpmsr_set_temp <- dpmsr_set
      
      dpmsr_set$data$data_raw_decoypsm <- NULL
      dpmsr_set$data$data_raw_decoypeptide <- NULL
      dpmsr_set$data$data_raw_decoyprotein <- NULL
      dpmsr_set$data$data_raw_inputfiles <- NULL
      dpmsr_set$data$data_raw_psm <- NULL
      #dpmsr_set$data$data_to_norm <- NULL
      #dpmsr_set$data$protein_imputed_df <- NULL
      #dpmsr_set$data$protein_missing <- NULL
      #dpmsr_set$data$normalized <- NULL
      #dpmsr_set$data$impute <- NULL
      dpmsr_set$overview <- NULL
      #dpmsr_set$protocol <- NULL
      #dpmsr_set$data$norm_data <- NULL
      
      save(dpmsr_set, file=str_c(dpmsr_set$file$output_dir, input$dpmsr_set_name_customer, ".dpmsr_set"))

      dpmsr_set <- dpmsr_set_temp
      dpmsr_set_temp <- NULL
    })     
    
    
    

    
    
    #-------------------------------------------------------------------------------------------------------------      
    #------------------------------------------------------------------------------------------------------------- 
    observeEvent(input$action_load_dpmsr_set, { 
      cat(file = stderr(), "load dpmsr_set triggered", "\n")
      showModal(modalDialog("Working...", footer = NULL))
      dpmsr_file <- parseFilePaths(volumes, input$dpmsr_set_file)

      load(file = dpmsr_file$datapath, envir = .GlobalEnv)  
      
      original_data_dir <- dpmsr_set$file$data_dir
      cat(file = stderr(), str_c("Original data directory --> ", original_data_dir), "\n")
      
      #reset file locations
      showModal(modalDialog("Working...", footer = NULL))
      dpmsr_set$file$data_dir <<- dirname(dpmsr_file$datapath)
      dpmsr_set$file$data_path <<- gsub("(.*)/.*","\\1",dpmsr_set$file$data_dir)
      dpmsr_set$file$output_dir <<- str_replace_all(dpmsr_set$file$data_dir, "/", "//")
      dpmsr_set$file$output_dir <<- str_c(dpmsr_set$file$output_dir, "//")        
      dpmsr_set$file$backup_dir <<- str_c(dpmsr_set$file$output_dir, "Backup//")
      dpmsr_set$file$extra_dir <<- str_c(dpmsr_set$file$output_dir, "Extra//")        
      dpmsr_set$file$qc_dir <<- str_c(dpmsr_set$file$output_dir, "QC//")
      dpmsr_set$file$string <<- str_c(dpmsr_set$file$output_dir, "String//")
      dpmsr_set$file$extra_prefix2 <<- str_c(dpmsr_set$file$extra_dir, dpmsr_set$x$file_prefix)
      dpmsr_set$file$phos <<- str_c(dpmsr_set$file$output_dir, "Phos//")
      
      #check dpmsr_set file, update if out of date
      update_version(session, input, output)

      #create dir for excel reports        
      if (!is_dir(str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp))) {
          create_dir(str_c(dpmsr_set$file$output_dir, dpmsr_set$data$stats$final_comp))
      }
      #reload shiny
      update_widget_all(session, input, output)
      update_dpmsr_set_from_widgets(session, input)
      inputloaddata_render(session, input, output)
      inputfilterapply_render(session, input, output)
      inputnorm_render(session, input, output)
      qc_render(session, input, output)
      inputproteinselect_render(session, input, output)
      updateTabsetPanel(session, "nlp1", selected = "tp8")
      hide_enable_after_dpmsr(session, input, output)
      removeModal()

      cat(file = stderr(), str_c("New data directory --> ", dpmsr_set$file$data_dir), "\n")
    })  
    
    

    # design_list <- c("General", "QC", "Fill_Norm", "Filters", "Protein", "TMT_PTM", "Report")
    # design_data <- parseFilePaths(volumes, input$design_file)
    # new_path <- str_extract(design_data$datapath, "^/.*/")
    # #new_path <- substr(new_path, 1, nchar(new_path)-1)
    # #new_path <- str_c(".", new_path)
    # volumes["wd"] <- new_path
    # design<-read_excel(design_data$datapath, sheet="SampleList")
    # 
    # 
    # output$downloadData <- downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #     copy_data <- parseFilePaths(volumes, input$test)
    #   },
    #   
    #   content <- function(file) {
    #     file.copy("out.zip", file)
    #   },
    #   contentType = "application/zip"
    # )   
    # 
    
    

    #-------------------------------------------------------------------------------------------------------------      
    #------------------------------------------------------------------------------------------------------------- 
    observeEvent(input$load_customer_dpmsr_set, { 
      cat(file = stderr(), "load customer dpmsr_set triggered", "\n")
      
      req(input$customer_dpmsr_set)
      
      if (!is.null(input$customer_dpmsr_set)) {fileUploaded <- TRUE  }  
      
      if (fileUploaded) {
        showModal(modalDialog("Working...", footer = NULL))  
        cat(file = stderr(), str_c("load file location:  ", input$customer_dpmsr_set$datapath), "\n")
        
        load(file = input$customer_dpmsr_set$datapath, envir = .GlobalEnv)
        
        original_data_dir <- dpmsr_set$file$data_dir
        cat(file = stderr(), str_c("Original data directory --> ", original_data_dir), "\n")
        
        update_try <- 1
        while (original_data_dir == dpmsr_set$file$data_dir) {
          cat(file = stderr(), str_c("dpmsr_set file update attempt --> ", update_try), "\n")
          if (update_try > 1) {
            load(file = input$customer_dpmsr_set$datapath, envir = .GlobalEnv)
          }
          update_try <- update_try + 1
          
          #check dpmsr_set file, update if out of date
          update_version(session, input, output)

          
          #reload shiny
          cat(file = stderr(), "update widgets", "\n")
          update_widget_all(session, input, output)
          
          cat(file = stderr(), "update dpmsr_set from widgets", "\n")
          update_dpmsr_set_from_widgets(session, input)
          
          
          #reset file locations
          tmp_dir <- format(Sys.time(), "%Y%m%d%H%M%S")
          cat(file = stderr(), "update file locations 1", "\n")
          #dpmsr_set$file$data_dir <<- dirname(dpmsr_file$datapath)
          
          dpmsr_set$file$tmp_dir <<- str_c("/data/ShinyData/", tmp_dir)
          dpmsr_set$file$data_dir <<- str_c(dpmsr_set$file$tmp_dir, "/", dpmsr_set$x$file_prefix)
          cat(file = stderr(), str_c("dpmsr_set$file_data_dir set to --->", dpmsr_set$file$data_dir), "\n")
          create_dir(dpmsr_set$file$data_dir)
          
          dpmsr_set$file$output_dir <<- str_replace_all(dpmsr_set$file$data_dir, "/", "//")
          dpmsr_set$file$output_dir <<- str_c(dpmsr_set$file$output_dir, "//")
          create_dir(dpmsr_set$file$output_dir)
          cat(file = stderr(), str_c("output_dir =", dpmsr_set$file$output_dir), "\n")
  
          dpmsr_set$file$string <<- str_c(dpmsr_set$file$output_dir, "String//")
          create_dir(dpmsr_set$file$string)
          cat(file = stderr(), str_c("string =", dpmsr_set$file$string), "\n")
          
          create_dir(str_c(dpmsr_set$file$output_dir,dpmsr_set$data$stats$final_comp))
          cat(file = stderr(), str_c("output_dir/final =", str_c(dpmsr_set$file$output_dir,dpmsr_set$data$stats$final_comp)), "\n")
          
          cat(file = stderr(), "update file locations... end", "\n")
        } 
        
        app_log()
        removeModal()
        updateTabsetPanel(session, "nlp1", selected = "tp_stats")
        
       } else{
          shinyalert("Oops!", "Please wait for file to upload...", type = "error")
        }
        
      
      cat(file = stderr(), str_c("Updated data directory --> ", dpmsr_set$file$data_dir), "\n")
    })  
    
    
    # This code will be run after the client has disconnected
    session$onSessionEnded(function() {
      if (site_user != "dpmsr") {
        cat(file = stderr(), "Running session end...", "\n")
        name <- dpmsr_set$file$tmp_dir
        cat(file = stderr(), str_c("Temp dir ", name,  " exists? ", dir.exists(name)), "\n")
        #do.call(file.remove, list(list.files(name, full.names = TRUE)))
        #dir_delete(name)
        system(str_c("rm -R ", name))
        cat(file = stderr(), str_c("Deleting ", name,  " still exists? ", dir.exists(name)), "\n")
        garbage_cleanup()
        cat(file = stderr(), str_c("DPMSR exists? ", exists('dpmsr_set')), "\n")
        try(rm(list = ls(envir = .GlobalEnv), pos = .GlobalEnv, inherits = FALSE))
        #try(rm("dpmsr_set", pos = .GlobalEnv, inherits = FALSE))
        cat(file = stderr(), str_c("Deleting DPMSR, still exists? ", exists('dpmsr_set')), "\n")
        gc()
        stopApp()
        #file_touch("restart.txt", access_time = Sys.time(), modification_time = Sys.time())
      }
    })
    
    
    
    
    #-------------------------------------------------------------------------------------------------------------      
    #------------------------------------------------------------------------------------------------------------- 
    
    
    cat(file = stderr(), "Shiny Server started ...end", "\n")  
}
)

