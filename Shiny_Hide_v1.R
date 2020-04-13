
#----------------------------------------------------------------------
#----------------------------------------------------------------------
 hide_enable <- function(session, input){
   
   observeEvent(input$design_file, {
      shinyjs::enable("action_load_design")
  })
  
   observeEvent(input$dpmsr_set_file, {
     shinyjs::enable("action_load_dpmsr_set")
   })
   
   observeEvent(input$raw_files, {
     shinyjs::enable("load_data")
   })
   
  # observe({
  #   if (is.null(input$file_raw_data) || input$file_raw_data == "") {
  #     shinyjs::disable("load_data")
  #   } else {
  #     shinyjs::enable("load_data")
  #   }
  # })   
  
  observe({
    if ((input$radio_input==1 || input$radio_input==3) && input$radio_output==1 ){
      shinyjs::show("razor")
    } else {
      shinyjs::hide("razor")
    }
  })
  
  observe({
    if (input$radio_impute==2) {
      shinyjs::show("impute_floor")
    } else {
      shinyjs::hide("impute_floor")
    }
  })
  
  
  observe({
    if (is.null(input$dpmsr_set_file)) {               #input$impute1==0
      shinyjs::enable("start_stats")
    } else {
      shinyjs::enable("start_stats")
    }
  })
  
  observe({
    if (input$radio_input==1 || input$radio_input==2 || input$radio_output==1) {
      shinyjs::hide("checkbox_isoform")
    } else {
      shinyjs::show("checkbox_isoform")
    }
  })   
  
  
  
  observe({
    if (input$radio_output==1) {
      shinyjs::hide("checkbox_norm_ptm")
      shinyjs::hide("checkbox_norm_impute")
      shinyjs::hide("filter_grep")
      shinyjs::hide("checkbox_out_ptm")
      shinyjs::hide("report_grep")
      updateCheckboxInput(session, "checkbox_isoform", value=FALSE)
    } else {
      shinyjs::show("checkbox_norm_ptm")
      shinyjs::show("checkbox_norm_impute")
      shinyjs::show("filter_grep")
      shinyjs::show("checkbox_out_ptm")
      shinyjs::show("report_grep")
    }
  })   
  
  
  
  observe({
    if (!input$checkbox_filtercv) {
      shinyjs::hide("text_filtercvgroup")
      shinyjs::hide("number_filtercvcutoff")
    } else {
      shinyjs::show("text_filtercvgroup")
      shinyjs::show("number_filtercvcutoff")
    }
  })  
  
  
  observe({
    if (input$checkbox_tmt) {
      updateCheckboxInput(session, "checkbox_n2", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n3", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n4", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n5", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n6", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n7", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n8", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n9", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n10", value = FALSE) 
      updateCheckboxInput(session, "checkbox_n11", value = FALSE) 
      shinyjs::hide("checkbox_n2")
      shinyjs::hide("checkbox_n3")
      shinyjs::hide("checkbox_n4")
      shinyjs::hide("checkbox_n5")
      shinyjs::hide("checkbox_n6")
      shinyjs::hide("checkbox_n7")
      shinyjs::hide("checkbox_n8")
      shinyjs::hide("checkbox_n9")
      shinyjs::hide("checkbox_n10")
      shinyjs::hide("checkbox_n11")
    } else {
      shinyjs::show("checkbox_n2")
      shinyjs::show("checkbox_n3")
      shinyjs::show("checkbox_n4")
      shinyjs::show("checkbox_n5")
      shinyjs::show("checkbox_n6")
      shinyjs::show("checkbox_n7")
      shinyjs::show("checkbox_n8")
      shinyjs::show("checkbox_n9")
      shinyjs::show("checkbox_n10")
      shinyjs::show("checkbox_n11")
    }
  })  
  
  
  observe({
    if (input$checkbox_tmt){
    hideTab(inputId = "nlp1", target = "tp6")
    showTab(inputId = "nlp1", target = "tp_tmt")
  }else{
    showTab(inputId = "nlp1", target = "tp6")
    hideTab(inputId = "nlp1", target = "tp_tmt")
       }
  })
  
  
  observe({
    if (!input$checkbox_out_ptm){
      hideTab(inputId = "nlp1", target = "tp_phos")
    }else{
      showTab(inputId = "nlp1", target = "tp_phos")
    }
  })
  
  
  
  
  
  observe({
    if (input$checkbox_n3){
      updateCheckboxInput(session, "checkbox_n1", value = TRUE) 
      shinyjs::disable("checkbox_n1")
    }else{
      shinyjs::enable("checkbox_n1")
    }
  })
  
  
  observe({
    if (!input$checkbox_n11){
      shinyjs::hide("protein_norm_list")
    }else {
      shinyjs::show("protein_norm_list")
    }
  })
  
  observe({
    if (!input$checkbox_tmt_filter){
      shinyjs::hide("tmt_filter_sd")
    }else {
      shinyjs::show("tmt_filter_sd")
    }
  })
  
  
  
  
  observe({
    test_comp <- input$comp_number
    
    if (test_comp<6) {
      shinyjs::hide("comp_6N")
      shinyjs::hide("comp_6D")
    } else {
      shinyjs::show("comp_6N")
      shinyjs::show("comp_6D")
    }
    if (input$comp_number<5) {
      shinyjs::hide("comp_5N")
      shinyjs::hide("comp_5D")
    } else {
      shinyjs::show("comp_5N")
      shinyjs::show("comp_5D")
    }
    if (input$comp_number<4) {
      shinyjs::hide("comp_4N")
      shinyjs::hide("comp_4D")
    } else {
      shinyjs::show("comp_4N")
      shinyjs::show("comp_4D")
    }
    if (input$comp_number<3) {
      shinyjs::hide("comp_3N")
      shinyjs::hide("comp_3D")
    } else {
      shinyjs::show("comp_3N")
      shinyjs::show("comp_3D")
    }
    if (input$comp_number<2) {
      shinyjs::hide("comp_2N")
      shinyjs::hide("comp_2D")
    } else {
      shinyjs::show("comp_2N")
      shinyjs::show("comp_2D")
    }
    
    if (!input$checkbox_nc1){
      shinyjs::hide("sl_plot")
      shinyjs::hide("sl_plot_select")
    }else {
      shinyjs::show("sl_plot")
      shinyjs::show("sl_plot_select")
    }
    
    if (!input$checkbox_nc2){
      shinyjs::hide("tmm_plot")
      shinyjs::hide("tmm_plot_select")
    }else {
      shinyjs::show("tmm_plot")
      shinyjs::show("tmm_plot_select")
    }
    
    if (!input$checkbox_nc3){
      shinyjs::hide("sltmm_plot")
      shinyjs::hide("sltmm_plot_select")
    }else {
      shinyjs::show("sltmm_plot")
      shinyjs::show("sltmm_plot_select")
    }
    
    if (!input$checkbox_nc4){
      shinyjs::hide("quantile_plot")
      shinyjs::hide("quantile_plot_select")
    }else {
      shinyjs::show("quantile_plot")
      shinyjs::show("quantile_plot_select")
    }
    
    if (!input$checkbox_nc5){
      shinyjs::hide("lr_plot")
      shinyjs::hide("lr_plot_select")
    }else {
      shinyjs::show("lr_plot")
      shinyjs::show("lr_plot_select")
    }
    
    if (!input$checkbox_nc6){
      shinyjs::hide("loess_plot")
      shinyjs::hide("loess_plot_select")
    }else {
      shinyjs::show("loess_plot")
      shinyjs::show("loess_plot_select")
    }
    
    if (!input$checkbox_nc7){
      shinyjs::hide("vsn_plot")
      shinyjs::hide("vsn_plot_select")
    }else {
      shinyjs::show("vsn_plot")
      shinyjs::show("vsn_plot_select")
    }
    
    if (!input$checkbox_nc8){
      shinyjs::hide("ti_plot")
      shinyjs::hide("ti_plot_select")
    }else {
      shinyjs::show("ti_plot")
      shinyjs::show("ti_plot_select")
    }
    
    if (!input$checkbox_nc9){
      shinyjs::hide("mi_plot")
      shinyjs::hide("mi_plot_select")
    }else {
      shinyjs::show("mi_plot")
      shinyjs::show("mi_plot_select")
    }
    
    if (!input$checkbox_nc10){
      shinyjs::hide("ai_plot")
      shinyjs::hide("ai_plot_select")
    }else {
      shinyjs::show("ai_plot")
      shinyjs::show("ai_plot_select")
    }
    
    if (!input$checkbox_nc11){
      shinyjs::hide("protein_plot")
      shinyjs::hide("protein_plot_select")
    }else {
      shinyjs::show("protein_plot")
      shinyjs::show("protein_plot_select")
    }
    
    
    observe({
      
      #test_comp <- dpmsr_set$x$comp_number
      testcomp <- input$comp_number
      
      if (test_comp<6) {
        shinyjs::hide("volcano_plot6")
      } else {
        shinyjs::show("volcano_plot6")
      }
      
      if (test_comp<5) {
        shinyjs::hide("volcano_plot5")
      } else {
        shinyjs::show("volcano_plot5")
      }
      
      if (test_comp<4) {
        shinyjs::hide("volcano_plot4")
      } else {
        shinyjs::show("volcano_plot4")
      }
      
      if (test_comp<3) {
        shinyjs::hide("volcano_plot3")
      } else {
        shinyjs::show("volcano_plot3")
      }
      
      if (test_comp<2) {
        shinyjs::hide("volcano_plot2")
      } else {
        shinyjs::show("volcano_plot2")
      }
      
    })   
    
    
    
  })    
  
  
  
  observe({
    test_mva_comp <- as.numeric(input$mva_comp)
    if (test_mva_comp<12) {
      shinyjs::hide("var_12N")
      shinyjs::hide("var_12D")
    } else {
      shinyjs::show("var_12N")
      shinyjs::show("var_12D")
    }
    if (test_mva_comp<11) {
      shinyjs::hide("var_11N")
      shinyjs::hide("var_11D")
    } else {
      shinyjs::show("var_11N")
      shinyjs::show("var_11D")
    }   
    if (test_mva_comp<10) {
      shinyjs::hide("var_10N")
      shinyjs::hide("var_10D")
    } else {
      shinyjs::show("var_10N")
      shinyjs::show("var_10D")
    }     
    if (test_mva_comp<9) {
      shinyjs::hide("var_9N")
      shinyjs::hide("var_9D")
    } else {
      shinyjs::show("var_9N")
      shinyjs::show("var_9D")
    }   
    if (test_mva_comp<8) {
      shinyjs::hide("var_8N")
      shinyjs::hide("var_8D")
    } else {
      shinyjs::show("var_8N")
      shinyjs::show("var_8D")
    }  
    if (test_mva_comp<7) {
      shinyjs::hide("var_7N")
      shinyjs::hide("var_7D")
    } else {
      shinyjs::show("var_7N")
      shinyjs::show("var_7D")
    }  
    if (test_mva_comp<6) {
      shinyjs::hide("var_6N")
      shinyjs::hide("var_6D")
    } else {
      shinyjs::show("var_6N")
      shinyjs::show("var_6D")
    }  
    if (test_mva_comp<5) {
      shinyjs::hide("var_5N")
      shinyjs::hide("var_5D")
    } else {
      shinyjs::show("var_5N")
      shinyjs::show("var_5D")
    }  
    if (test_mva_comp<4) {
      shinyjs::hide("var_4N")
      shinyjs::hide("var_4D")
    } else {
      shinyjs::show("var_4N")
      shinyjs::show("var_4D")
    }  
    if (test_mva_comp<3) {
      shinyjs::hide("var_3N")
      shinyjs::hide("var_3D")
    } else {
      shinyjs::show("var_3N")
      shinyjs::show("var_3D")
    }  
    if (test_mva_comp<2) {
      shinyjs::hide("var_2N")
      shinyjs::hide("var_2D")
    } else {
      shinyjs::show("var_2N")
      shinyjs::show("var_2D")
    }  
  })
  
  
  observe({
    if(!input$mva_checkbox_out_ptm){
     shinyjs::hide("mva_report_grep") 
    }else{
     shinyjs::show("mva_report_grep")
    }
    })
    
  observe({
    if(!input$mva_checkbox_report_accession){
      shinyjs::hide("mva_report_accession") 
    }else{
      shinyjs::show("mva_report_accession")
    }
  })
  
  observe({
    if(input$check_mva==0){
      shinyjs::disable("start_mva") 
    }else{
      shinyjs::enable("start_mva")
    }
  })
  
  
}



