#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

inputloaddata_render <- function(session, input, output){
  
  output$text1 <- renderText({str_c("Original Data Format:  ",  as.character(dpmsr_set$x$raw_data_input)   )  })
  output$text2 <- renderText({str_c("Current Data Format:  ",  as.character(dpmsr_set$y$state)   )  })
  output$text3 <- renderText({str_c(as.character(dpmsr_set$y$state), " Count:  ",  nrow(dpmsr_set$data$data_peptide) )   })
  
  output$mass_accuracy<- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Mass_Accuracy.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$inj_summary <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"PSM_inj_summary.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  overview_df <- data.frame(names(dpmsr_set$overview), unlist(dpmsr_set$overview))
  output$project_overview <- renderRHandsontable({
    rhandsontable(overview_df, rowHeaders = NULL,  colHeaders = NULL, width = 800, height = 800) 
      })
  
  output$peptide_RT <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_RT.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$peptide_MZ <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_MZ.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$peptide_Charge <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_Charge.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$faims_psm_cv <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"FAIMS_PSM_CV.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE) 
  
  output$peptide_PI <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_PI.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$peptide_cruc <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_Cruc.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$peptide_aainfo <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_AAInfo.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$peptide_ai <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_AI.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$feature_width <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Feature_Peak_Width.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
}

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

inputfilterapply_render <- function(session, input, output){
  
  output$text4 <- renderText({str_c("Filtered ", as.character(dpmsr_set$y$state), " Count:  ",  nrow(dpmsr_set$data$data_peptide) )   })
  output$raw_bar <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Raw_barplot.png"), contentType = 'image/png', width=600, height=500, alt="this is alt text")
  }, deleteFile = FALSE)
  
}

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

inputnorm_render <- function(session, input, output){
  
  output$histogram <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Intensity_Histogram.png"), contentType = 'image/png', width=600, height=500, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$text_i2 <- renderText(str_c("Intensity cutoff:  ", as.character(dpmsr_set$x$int_cutoff)))
  
  output$tmt_sl <- renderImage({
    list(src=str_c(dpmsr_set$file$TMT_dir,"TMT_IRS_", input$tmt_step, "_barplot.png"), contentType = 'image/png', width=600, height=500, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$histogram_tmt <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Intensity_Histogram.png"), contentType = 'image/png', width=600, height=500, alt="this is alt text")
  }, deleteFile = FALSE)
  
}

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

qc_spike_render <- function(session, input, output){
  
  output$data_CV <- renderRHandsontable({
    cv_table <- dpmsr_set$data$summary_cv
    cv_table [,-1] <-trunc(round(cv_table [,-1],0))
    cv_table <- cv_table %>% mutate_all(as.character)
    rhandsontable(cv_table, readOnly = TRUE, rowHeaders = NULL, digits = 0)
  })
  
  output$data_proteinQC <- renderRHandsontable({
    rhandsontable(dpmsr_set$data$qc_spike_final, readOnly = TRUE, rowHeaders = NULL, digits = 0)
  })
  
  output$protein_qc_spike_levels <- renderRHandsontable({
    rhandsontable(dpmsr_set$design,  rowHeaders = NULL, digits = 0) %>%
      hot_col(col = "Set", halign = "htCenter", format="0", digits =0, readOnly = TRUE) %>%
      hot_col(col = "PD_Order", halign = "htCenter", format="0", digits =0, readOnly = TRUE) %>%
      hot_col(col = "Replicate", halign = "htCenter", format="0", digits =0, readOnly = TRUE) %>%
      hot_col(col = "QC Spike Level", halign = "htCenter", format="0", digits =0, readOnly = FALSE) %>%
      hot_col(col = "Header1", halign = "htLeft", width=200, readOnly = TRUE) %>%
      hot_col(col = "Header2", halign = "htLeft", width=250, readOnly = TRUE) %>%
      hot_col(col = "Header3", halign = "htLeft", width=250, readOnly = TRUE) 
  })
  
  
}

  #----------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------
  
  qc_render <- function(session, input, output){
    
  qc_spike_render(session, input, output)
    
  output$cv_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir, "CV_barplot.png"),  
         contentType = 'image/png', width=800, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$qc_spike_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir, input$norm_type, "_QC_Spike_barplot.png"),  
         contentType = 'image/png', width=800, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$adh_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir, "ADH_barplot.png"),  
         contentType = 'image/png', width=800, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  
  output$adh_boxplot <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir, "ADH_boxplot.png"),  
         contentType = 'image/png', width=800, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$impute_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "impute//impute_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$sl_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "sl//sl_", input$plot_select,".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$sltmm_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "sltmm//sltmm_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$tmm_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "tmm//tmm_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$quantile_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "quantile//quantile_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$loess_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "loess//loess_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$lr_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "lr//lr_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$ai_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "ai//ai_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$mi_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "mi//mi_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$ti_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "ti//ti_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$vsn_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "vsn//vsn_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$protein_plot <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "protein//protein_", input$plot_select, ".png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
}


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

inputproteinselect_render <- function(session, input, output){
  
  
  output$impute_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "impute//", input$protein_select, "_impute_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$sl_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "sl//", input$protein_select,"_sl_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$sltmm_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "sltmm//", input$protein_select, "_sltmm_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$tmm_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "tmm//", input$protein_select, "_tmm_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$quantile_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "quantile//", input$protein_select, "_quantile_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$loess_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "loess//", input$protein_select, "_loess_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$lr_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "lr//", input$protein_select, "_lr_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$ai_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "ai//", input$protein_select, "_ai_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$mi_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "mi//", input$protein_select, "_mi_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$ti_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "ti//", input$protein_select, "_ti_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$vsn_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "vsn//", input$protein_select, "_vsn_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$protein_plot_select <- renderImage({
    list(src=str_c(dpmsr_set$file$output_dir, "protein//", input$protein_select, "_protein_barplot.png"),  
         contentType = 'image/png', width=500, height=400, alt="this is alt text")
  }, deleteFile = FALSE)
  
  
}