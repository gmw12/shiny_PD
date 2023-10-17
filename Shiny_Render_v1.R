#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

inputloaddata_render <- function(session, input, output){
  cat(file = stderr(), "inputloaddata_render...", "\n")
  output$text1 <- renderText({str_c("Original Data Format:  ",  as.character(dpmsr_set$x$raw_data_input)   )  })
  output$text2 <- renderText({str_c("Current Data Format:  ",  as.character(dpmsr_set$y$state)   )  })
  
  if (dpmsr_set$x$raw_data_input == "Protein") {
    output$text3 <- renderText({str_c(as.character(dpmsr_set$y$state), " Count:  ",  nrow(dpmsr_set$data$data_protein_start) )   })  
  } else if (dpmsr_set$x$raw_data_input == "Protein_Peptide" ||  dpmsr_set$x$raw_data_input == "Peptide")  {
    output$text3 <- renderText({str_c(as.character(dpmsr_set$y$state), " Count:  ",  nrow(dpmsr_set$data$data_peptide_start) )   })
  } else if (dpmsr_set$x$raw_data_input == "Precursor" ||  dpmsr_set$x$raw_data_input == "Precursor_PTM")  {
    output$text3 <- renderText({str_c(as.character(dpmsr_set$y$state), " Count:  ",  nrow(dpmsr_set$data$data_precursor_start) )   })
  }
  
  output$mass_accuracy <- renderImage({
    list(src = str_c(dpmsr_set$file$qc_dir,"Mass_Accuracy.png"), contentType = 'image/png', width = 500, height = 300, alt = "this is alt text")
  }, deleteFile = FALSE)
  
  output$inj_summary <- renderImage({
    list(src = str_c(dpmsr_set$file$qc_dir,"PSM_inj_summary.png"), contentType = 'image/png', width = 500, height = 300, alt = "this is alt text")
  }, deleteFile = FALSE)
  
  overview_df <- data.frame(names(dpmsr_set$overview), unlist(dpmsr_set$overview))
  output$project_overview <- renderRHandsontable({
    rhandsontable(overview_df, rowHeaders = NULL,  colHeaders = NULL, width = 800, height = 800) 
      })
  
  observe({
    if (file.exists(str_c(dpmsr_set$file$qc_dir,"Peptide_RT.png"))) {
      output$peptide_RT <- renderImage({
        list(src = str_c(dpmsr_set$file$qc_dir,"Peptide_RT.png"), contentType = 'image/png', width = 500, height = 300, alt = "this is alt text")
      }, deleteFile = FALSE)
    }else{
      invalidateLater(50000, session)
    }
  }) 
  
  output$peptide_MZ <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_MZ.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$peptide_Charge <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_Charge.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$faims_psm_cv <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"FAIMS_PSM_CV.png"), contentType = 'image/png', width=500, height=300, alt="No FAIMS plot available")
  }, deleteFile = FALSE) 
  

    observe({
    if (file.exists(str_c(dpmsr_set$file$qc_dir,"Peptide_PI.png"))){
      output$peptide_PI <- renderImage({
        list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_PI.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
      }, deleteFile = FALSE)
    }else{
      invalidateLater(50000, session)
    }
  }) 
  
  observe({
    if (file.exists(str_c(dpmsr_set$file$qc_dir,"Peptide_Cruc.png"))){
      output$peptide_cruc <- renderImage({
        list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_Cruc.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
      }, deleteFile = FALSE)
    }else{
      invalidateLater(50000, session)
    }
  }) 
  
  observe({
    if (file.exists(str_c(dpmsr_set$file$qc_dir,"Peptide_AAInfo.png"))){
      output$peptide_aainfo <- renderImage({
        list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_AAInfo.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
      }, deleteFile = FALSE)
    }else{
      invalidateLater(1000, session)
    }
  }) 
  
  observe({
    if (file.exists(str_c(dpmsr_set$file$qc_dir,"Peptide_AI.png"))){
      output$peptide_ai <- renderImage({
        list(src=str_c(dpmsr_set$file$qc_dir,"Peptide_AI.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
      }, deleteFile = FALSE)
    }else{
      invalidateLater(1000, session)
    }
  }) 
  
  observe({
    if (file.exists(str_c(dpmsr_set$file$qc_dir,"Feature_Peak_Width.png"))){
      output$feature_width <- renderImage({
        list(src=str_c(dpmsr_set$file$qc_dir,"Feature_Peak_Width.png"), contentType = 'image/png', width=500, height=300, alt="this is alt text")
      }, deleteFile = FALSE)
    }else{
      invalidateLater(1000, session)
    }
  }) 
  
  output$comp1_text <- renderText({ "Comparison 1" })
  output$comp2_text <- renderText({ "Comparison 2" })
  output$comp3_text <- renderText({ "Comparison 3" })
  output$comp4_text <- renderText({ "Comparison 4" })
  output$comp5_text <- renderText({ "Comparison 5" })
  output$comp6_text <- renderText({ "Comparison 6" })
  output$comp7_text <- renderText({ "Comparison 7" })
  output$comp8_text <- renderText({ "Comparison 8" })
  output$comp9_text <- renderText({ "Comparison 9" })
  output$comp10_text <- renderText({ "Comparison 10" })
  output$comp11_text <- renderText({ "Comparison 11" })  
  output$comp12_text <- renderText({ "Comparison 12" })
  
}

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

inputfilterapply_render <- function(session, input, output){
  cat(file = stderr(), "inputfilterapply_render...", "\n")
  
  if (dpmsr_set$x$raw_data_input == "Protein_Peptide" || dpmsr_set$x$raw_data_input == "Peptide") {
      output$text4 <- renderText({str_c("Filtered ", as.character(dpmsr_set$y$state), " Count:  ",  nrow(dpmsr_set$data$data_peptide) ) })
  }else if (dpmsr_set$x$raw_data_input == "Precursor_PTM" || dpmsr_set$x$raw_data_input == "Precursor") {
      output$text4 <- renderText({str_c("Filtered ", as.character(dpmsr_set$y$state), " Count:  ",  nrow(dpmsr_set$data$data_precursor) )   })
  }
    
  output$raw_bar <- renderImage({
    list(src = str_c(dpmsr_set$file$qc_dir,"Raw_barplot.png"), contentType = 'image/png', width = 600, height = 500, alt = "this is alt text")
  }, deleteFile = FALSE)
  
  output$norm_data_bar <- renderImage({
    list(src = str_c(dpmsr_set$file$qc_dir,"Normalization_Data_barplot.png"), contentType = 'image/png', width = 600, height = 500, alt = "No PTM Only barplot!")
  }, deleteFile = FALSE)
  
  
  
}

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

inputnorm_render <- function(session, input, output){
  cat(file=stderr(), "inputnorm_render...", "\n")
  output$histogram <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"Intensity_Histogram.png"), contentType = 'image/png', width=600, height=500, alt="this is alt text")
  }, deleteFile = FALSE)
  
  output$ptm_histogram <- renderImage({
    list(src=str_c(dpmsr_set$file$qc_dir,"PTM_Only_Intensity_Histogram.png"), contentType = 'image/png', width=600, height=500, alt="No PTM Histogram!")
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
  cat(file=stderr(), "qc_spike_render...", "\n")
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
      hot_col(col = "QC Spike Level", halign = "htCenter", format="0.0", digits =0, readOnly = FALSE) %>%
      hot_col(col = "Header1", halign = "htLeft", width=200, readOnly = TRUE) %>%
      hot_col(col = "Header2", halign = "htLeft", width=250, readOnly = TRUE) %>%
      hot_col(col = "Header3", halign = "htLeft", width=250, readOnly = TRUE) 
  })
  
  
}

  #----------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------
  
  qc_render <- function(session, input, output){
    cat(file=stderr(), "qc_render...", "\n")
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
    
    output$directlfq_plot <- renderImage({
      list(src=str_c(dpmsr_set$file$output_dir, "directlfq//directlfq_", input$plot_select, ".png"),  
           contentType = 'image/png', width=500, height=400, alt="this is alt text")
    }, deleteFile = FALSE)
    
    output$directlfq_full_plot <- renderImage({
      list(src=str_c(dpmsr_set$file$output_dir, "directlfq_full//directlfq_full_", input$plot_select, ".png"),  
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
  
  cat(file=stderr(), "inputproteinselect_render...", "\n")
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
