
#----------------------------------------------------------------------
#----------------------------------------------------------------------
 hide_enable <- function(session, input){
   
   
   observeEvent(input$design_file, {
      shinyjs::enable("action_load_design")
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
    if (!input$checkbox_n11){
      shinyjs::hide("protein_norm_list")
    }else {
      shinyjs::show("protein_norm_list")
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
      if (is.null(input$dpmsr_set_file) || input$dpmsr_set_file == "") {
        shinyjs::disable("load_dpmsr_set")
      } else {
        shinyjs::enable("load_dpmsr_set")
      }
    })
    
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
  
}



