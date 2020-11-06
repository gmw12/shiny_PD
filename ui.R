
library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinyWidgets)
library(rhandsontable)
library(rgl)
library(DT)
library(shinyalert)

source("Shiny_Startup_v1.R")


shinyUI(
  tagList(
  fluidPage(
  useShinyjs(),
  useShinyalert(),
  setBackgroundColor("DarkGray", shinydashboard = TRUE),
  setBackgroundColor("LightGray", shinydashboard = FALSE),
  titlePanel("Proteome Discoverer Data Processing"),
  
  navlistPanel(widths=c(1,11), id = "nlp1",
    tabPanel("Load Page", value = "load", align="center",
             br(),
             br(),
             br(),
             br(),
             br(),
             tags$h1("Loading app...") 
    ),
    tabPanel("Load Design", value = "tp_load_design", align="center",
                        hr(),
                        tags$h1("Choose and Load the study design file..."),
                        br(),
                        shinyFilesButton('design_file', label='Choose Design File', title='Please select excel design file', multiple=FALSE,
                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        br(),                            
                        br(),
                        actionButton("action_load_design", label = "Load Design File", width = 200, 
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        br(),
                        hr(),
                        br(),
                        tags$h3("Or load saved dpmsr_set file..."),
                       shinyFilesButton('dpmsr_set_file', label='Choose dpmsr_set File', title='Choose dpmsr_set File', multiple=FALSE,
                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                       br(),
                       br(),
                       actionButton("action_load_dpmsr_set", label = "Load dpmsr_set", width = 200, 
                                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        br(), 
                        br(),
                        hr(),                         
                        br(), 
                        tags$h3("If needed clear data and functions..."),
                        br(),
                       actionButton("action_clear", label = "Clear Data", width = 200, 
                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
           ), #end tab panel
  
    tabPanel("Load Data", value = "tp_load_data", align="center",
             tags$h1("Choose type and load data file"),
             hr(),
             fluidRow(align = "left",
                       column(width=4, offset =2,
                         radioButtons("radio_input", label = h3("Select Input Data Type"),
                                      choices = list("Protein_Peptide" = 1, "Protein" = 2, "Peptide" = 3),
                                      selected = 3)
                         ),
                        column(width=5, offset =0,
                             radioButtons("radio_output", label = h3("Select Output Data Type"),
                                          choices = list("Protein" = 1, "Peptide" = 2),
                                          selected = 2)
                        )
                 ),
             hr(),
             fluidRow(align="left",
               column(width=3, offset =1,
                      checkboxInput("checkbox_tmt", label = "SPQC Normalized TMT sets"),
                      checkboxInput("checkbox_isoform", label = "Use Peptide Isoform?"),
               ),
               column(width=3, offset =1,   
                      checkboxInput("checkbox_norm_ptm", label = "Normalize on PTM?"),
                      textInput("peptide_norm_grep", label="Normalize PTM grep", value = "Enter value here")
               ),
              column(width=3, offset = 1,       
                      checkboxInput("checkbox_impute_ptm", label = "Impute Distribution based on PTM?"),
                      textInput("peptide_impute_grep", label="Impute PTM grep", value = "Enter value here")
               )
             ),
             fluidRow(        
                         selectInput("razor", label = h5("Peptides to Use"), 
                                     choices = list("Razor", "Unique", "Shared"), 
                                     selected = "Razor")
                    ),
              hr(),
            fluidRow(
                      textInput("fileprefix", label="Enter file prefix for data output", value = "Enter value here"),
 
                      shinyFilesButton('raw_files', label='Select PD Text Export Files', title='Please select raw data files', multiple=TRUE,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                      br(),
                      br(),
                      
                      actionButton("load_data", label = "Load Data File(s)",
                                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                     )
        ), #end tab panel
    
    tabPanel("Overview", value = "tp_overview",
          hr(),
          fluidRow(  
            column(width=3, offset =1,
               imageOutput("mass_accuracy"), 
               imageOutput("inj_summary"), 
               imageOutput("peptide_MZ"),
               imageOutput("faims_psm_cv"),
               imageOutput("peptide_PI"), 
               imageOutput("peptide_cruc") 
            ),
            column(width=3, offset =1,
                imageOutput("peptide_RT"),
                imageOutput("feature_width"),
                imageOutput("peptide_Charge"),
                imageOutput("peptide_aainfo"), 
                imageOutput("peptide_ai") 
               ),
            column(width=3, offset =1,
                   rHandsontableOutput("project_overview")   
            )
          )
    ),  


      tabPanel("Filters", value = "tp_filters", align="center",
               tags$h1("Apply Raw Peptide Filters"),
               hr(),
               fluidRow(
                 column(width=3, offset = 1,
                        fluidRow(align = "left",
                          textOutput("text1"),
                          tags$head(tags$style("#text1{color: blue; font-size: 16px; font-style: bold;}")),
                          textOutput("text2"),
                          tags$head(tags$style("#text2{color: blue; font-size: 16px; font-style: bold;}")),
                          textOutput("text3"),
                          tags$head(tags$style("#text3{color: blue; font-size: 16px; font-style: bold;}")),
                          textOutput("text4"),
                          tags$head(tags$style("#text4{color: blue; font-size: 16px; font-style: bold;}"))
                           )
                        ),
                  column(width=6, offset = 1,
                    fluidRow(align = "left",
                          numericInput("minimum_measured_all", label="Enter minimum # measured values (all samples)", value = 2),
                          hr(),
                          checkboxInput("checkbox_require_x", label = "Require X% measured values in at least one group"),
                          numericInput("require_x_cutoff", label="Enter X% measured values (decimal)", value = 0.8),
                          hr(),
                          checkboxInput("checkbox_filtercv", label = "Filter on Specific Group CV"),
                          selectInput("text_filtercvgroup", label="Enter group for CV filter", 
                                      choices = list("SPQC"), selected = "SPQC"),
                          numericInput("number_filtercvcutoff", label="Enter CV% for cutoff", value = "Enter value here"),
                          hr(),
                          br(),
                          actionButton("filter_apply", label = "Apply Filters", width = 300,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                       )
                      )
                 )
               ), #end tab panel


      tabPanel("Normalize", value = "tp_normalize",
               fluidRow(
                 column(width=4, offset =1,
                     br(),
                     hr(),
                     textOutput("text_n1"),
                     tags$head(tags$style("#text_n1{color: blue; font-size: 20px; font-style: bold;}")),
                     checkboxInput("checkbox_n1", label = "Sample Loading - Total", value = TRUE),
                     checkboxInput("checkbox_n2", label = "Trimmed Mean - 10%"),
                     checkboxInput("checkbox_n3", label = "SL TMM"),
                     checkboxInput("checkbox_n4", label = "Quantile"),
                     checkboxInput("checkbox_n5", label = "Linear Regression"),
                     checkboxInput("checkbox_n6", label = "LOESS"),
                     checkboxInput("checkbox_n7", label = "VSN"),
                     checkboxInput("checkbox_n8", label = "Median of Total Intensity"),
                     checkboxInput("checkbox_n9", label = "Median Intensity"),
                     checkboxInput("checkbox_n10", label = "Average Intensity"),
                     checkboxInput("checkbox_n11", label = "Protein"),
                     br(),
                     textInput("protein_norm_list", label="Protein Norm List", value = "Enter value here"),
                     hr(),
                     br(),
                     actionButton("norm1", label = "Apply Normalization", width = 300,
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                     ),
               column(width=7, offset =0,
                      "Raw Data Total Intensity",
                      br(),
                      fluidRow(
                        imageOutput("raw_bar"),
                      ),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      fluidRow(
                        imageOutput("raw_ptm_bar")
                      )
               )
          )
       ),


    tabPanel("Impute", value = "tp_impute",
             fluidRow(
               column(width=4, offset =1,
                      textOutput("text_i1"),
                      tags$head(tags$style("#text_i1{color: blue; font-size: 20px; font-style: bold;}")),
                      br(),
                      radioButtons("radio_impute", label=NULL,
                                   choices = list("Duke/BottomX" = 1, "Floor" = 2, "Minimum" = 3, "Average/Group" = 4,
                                                  "Average/Global" = 9, "KNN"= 5, "LocalLeastSquares" = 6, "MLE" = 7, "BottomX" = 8
                                   ),
                                   selected = 1),
                      dropdownButton(
                        numericInput("impute_floor", label="Impute Floor Intensity", value = "Enter value here"),
                        numericInput("missing_cutoff", label="%minimum  measured values in group to allow missing values to be imputed in measured range", value = 50, width = '100%'),
                        checkboxInput("checkbox_misaligned", label = "Misaligned Filter"),
                        numericInput("misaligned_cutoff", label="%missing values to be considered for misalignment if average > intensity cutoff", value = 50, width='100%'),
                        numericInput("intensity_cutoff_mean_sd", label="#standard deviations (+/-) from mean for intensity cuttof", value = -0.5, width = '100%'),
                        textOutput("text_i2"),
                        br(),
                        actionButton("adjust_intensity_cutoff", label = "Calc Intensity Cutoff", width = 300,
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        textOutput("text_impute_ptm"),
                        br(),
                        numericInput("bottom_x", label="Bottom X%", value = "5"),
                        circle = TRUE, status = "info", icon = icon("gear"), width = "500px", size = "sm",
                        tooltip = tooltipOptions(title = "More Imputation Options")
                        ),
                      hr(),
                      actionButton("impute1", label = "Apply Imputation", width = 300,
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                    ),
               column(width=7, offset =0,
                         "Total Intensity Histogram",
                      br(),
                      imageOutput("histogram"),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      imageOutput("ptm_histogram")

                      )
             )
    ), #end tab panel

    
    tabPanel("TMT SPQC", value = "tp_tmt",
             fluidRow(
               column(width=2, offset =0,
               numericInput("tmt_sets", label="TMT Sets for IRS Normalization", value = "Enter value here"),
               numericInput("tmt_channels", label="TMT Channels", value = "Enter value here"),
               checkboxInput("checkbox_tmt_filter", label = "Filter peptides by stdev of average %CV"),
               numericInput("tmt_filter_sd", label="Stdev for filter", value = "Enter value here"),
               br(),
               numericInput("TMT_SPQC_bottom_x", label="Bottom X%", value = "5"),
               hr(),
               actionButton("tmt_irs_go", label = "Apply IRS Norm", width = 200,
                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
               hr(),
               actionButton("tmt_irs_qc", label = "Start Data QC", width = 200,
                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
               ),
               
               column(width=5, offset =0,
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      "Total Intensity Histogram",
                      br(),
                      imageOutput("histogram_tmt")   
               ),        
               column(width=5, offset =0,
                      selectInput("tmt_step", label = h5("TMT IRS Step Barplot"), 
                                  choices = list("Step0", "Step1", "Step2", "Step3", "Step4", "Step5", "Step6", "Step7", "Final"), 
                                  selected = "Step0"),
                      br(),
                      imageOutput("tmt_sl")   
               )
               
               
             )
    ), #end tab panel
    
    
    
    tabPanel("QC", value="tp_qc",
             navbarPage("QC:", id ="np5",
                        tabPanel("CV",
                                 rHandsontableOutput("data_CV"),
                                 br(),
                                 imageOutput("cv_plot")
                        ),  
                        tabPanel("ADH",
                                 imageOutput("adh_plot"),
                                 br(),
                                 imageOutput("adh_boxplot")
                                 ),
                        tabPanel("Protein Spike",
                                 fluidRow(
                                   column(width=6, offset =0,
                                    rHandsontableOutput("data_proteinQC")
                                   ),
                                   column(width=3, offset =0,
                                    textInput("protein_spike_list", label="Protein Spike Accession(s)")  
                                   ),
                                   column(width=3, offset =0,
                                          actionButton("calc_protein_spike", label = "Recalculate Protein Spike", width = 200,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")  
                                   )
                                 ),
                                 hr(),
                                 fluidRow(
                                   column(width=2, offset =0,
                                          selectInput("norm_type", label = h5("Normalization Type"), 
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)),
                                   column(width=10, offset =0,
                                          imageOutput("qc_spike_plot"))
                                 ),
                                 hr(),
                                 fluidRow(
                                   rHandsontableOutput("protein_qc_spike_levels"),
                                   hr(),
                                   actionButton("update_qc_spike_levels", label = "Save Table", width = 100,
                                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                 )
                                 ),
                        tabPanel("Norm Comparison",
                                 fluidRow(
                                   column(width=2, offset =0,
                                          selectInput("plot_select", label = NULL, 
                                                      choices = list("barplot", "boxplot"), 
                                                      selected = "barplot"),
                                          checkboxInput("checkbox_nc1", label = "Sample Loading - Total"),
                                          checkboxInput("checkbox_nc2", label = "Trimmed Mean - 10%"),
                                          checkboxInput("checkbox_nc3", label = "SL TMM"),
                                          checkboxInput("checkbox_nc4", label = "Quantile"),
                                          checkboxInput("checkbox_nc5", label = "Linear Regression"),
                                          checkboxInput("checkbox_nc6", label = "LOESS"),
                                          checkboxInput("checkbox_nc7", label = "VSN"),
                                          checkboxInput("checkbox_nc8", label = "Median of Total Intensity"),
                                          checkboxInput("checkbox_nc9", label = "Median Intensity"),
                                          checkboxInput("checkbox_nc10", label = "Average Intensity"),
                                          checkboxInput("checkbox_nc11", label = "Protein")
                                   ),
                                   column(width=4, offset =0,
                                          imageOutput("impute_plot"),
                                          br(),
                                          imageOutput("sl_plot"),
                                          br(),
                                          imageOutput("quantile_plot"),
                                          br(),
                                          imageOutput("lr_plot"),
                                          br(),
                                          imageOutput("ai_plot"),
                                          br(),
                                          imageOutput("mi_plot")
                                   ),
                                   column(width=4, offset =1,
                                          imageOutput("sltmm_plot"),
                                          br(),
                                          imageOutput("tmm_plot"),
                                          br(),
                                          imageOutput("vsn_plot"),
                                          br(),
                                          imageOutput("loess_plot"),
                                          br(),
                                          imageOutput("ti_plot"),
                                          br(),
                                          imageOutput("protein_plot")
                                   )
                                 )),
                        
                          tabPanel("Selected Proteins",
                                 fluidRow(
                                   column(width=2, offset =0,
                                          div(style = "font-size:12px;",
                                              selectInput("protein_select", label = "Select barplot type", 
                                                          choices = list("ADH", "BiraA", "Bait", "Carbox", "Avidin", 
                                                                         "Protein1", "Protein2", "Protein3", "Protein4"), 
                                                          selected = "ADH"),
                                              
                                              textInput("adh_list", label="ADH Accession"),
                                              textInput("bait_list", label="Bait Accession"),
                                              textInput("avidin_list", label="Avidin Accession"),
                                              textInput("carbox_list", label="Carbox Accession"),
                                              textInput("bira_list", label="BirA Accession"),
                                              textInput("protein1_list", label="Protein1 Accession"),
                                              textInput("protein2_list", label="Protein2 Accession"),
                                              textInput("protein3_list", label="Protein3 Accession"),
                                              textInput("protein4_list", label="Protein4 Accession")
                                          ),
                                          actionButton("protein_select_plots", label = "Create/Display Plots", width = 300,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                          
                                   ),
                                   column(width=4, offset =0,
                                          imageOutput("impute_plot_select"),
                                          br(),
                                          imageOutput("sl_plot_select"),
                                          br(),
                                          imageOutput("quantile_plot_select"),
                                          br(),
                                          imageOutput("lr_plot_select"),
                                          br(),
                                          imageOutput("ai_plot_select"),
                                          br(),
                                          imageOutput("mi_plot_select")
                                   ),
                                   column(width=4, offset =1,
                                          imageOutput("sltmm_plot_select"),
                                          br(),
                                          imageOutput("tmm_plot_select"),
                                          br(),
                                          imageOutput("vsn_plot_select"),
                                          br(),
                                          imageOutput("loess_plot_select"),
                                          br(),
                                          imageOutput("ti_plot_select"),
                                          br(),
                                          imageOutput("protein_plot_select")
                                   )
                                 )),
                        
                        tabPanel("One Protein",
                                 fluidRow( 
                                   column(width=2, offset =0,
                                          selectInput("select_oneprotein_norm", label = "Normalization", width = 150,
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=1, offset =0,
                                          textInput("oneprotein_accession", label="Accession", value = "0", width = 100)
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("oneprotein_show", label = "Show Graph", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   )
                                 ),
                                 fluidRow(
                                   hr(),
                                   column(width=4, offset =1,
                                          imageOutput("oneprotein_plot_select")
                                   ), 
                                   column(width=4, offset =1,
                                          rHandsontableOutput("oneprotein_stats")
                                   )
                                 )),
                        
                        tabPanel("One Peptide",
                                 fluidRow( 
                                   column(width=2, offset =0,
                                          selectInput("select_onepeptide_norm", label = "Normalization", width = 150,
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=1, offset =0,
                                          textInput("onepeptide_sequence", label="Sequence", value = "0", width = 100)
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("onepeptide_show", label = "Show Graph", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   )
                                 ),
                                 fluidRow(
                                   hr(),
                                   column(width=4, offset =1,
                                          imageOutput("onepeptide_plot_select")
                                   ), 
                                   column(width=4, offset =1,
                                          rHandsontableOutput("onepeptide_stats")
                                   )
                                 )),
                        
                        tabPanel("Data", value = "tp10",
                                 fluidRow( 
                                   
                                   column(width=2, offset =0,
                                          selectInput("select_final_data", label = "Normalization", width = 150,
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=1, offset =0,  
                                          numericInput("data_topn", label="TopN", value = 0, width = 100)
                                   ),
                                   column(width=1, offset =0,
                                          textInput("data_accession", label="Accession", value = "0", width = 100)
                                   ),
                                   column(width=1, offset =0,
                                          textInput("data_description", label="Description", value = "0", width = 100)
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("data_show", label = "Filter Data", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   )
                                 ),
                                 fluidRow(
                                   hr(),
                                   tags$head(tags$style("#data_final{color: blue;
                                 font-size: 12px;
                                 }"
                                   )
                                   ),
                                   rHandsontableOutput("data_final")
                                 )
                        )
             ) #end of qc navbarpage
          ), #end of qc tabpanel
    
   
    
    tabPanel("Load Data", value = "tp_customer_load", align="center",
             hr(),
             tags$h1("Welcome to DPMSR Proteome Discoverer Data Visualization Tool"),
             br(),
             br(),
             br(),
             br(),
             tags$h2("Choose and Load the study design file..."),
             br(),
             fileInput("customer_dpmsr_set", "Choose dpmsr_set file",
                       multiple = FALSE,
                       accept = c(".dpmsr_set", ".txt")),
            br(),
            br(),
            br(),
            actionButton("load_customer_dpmsr_set", label = "Load dpmsr_set", width = 200, 
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    ),
    
    
    
    tabPanel("Stats", value = "tp_stats",
             navbarPage("Stats:", id ="nbp_stats",
                        tabPanel("Setup", id="tp_stats_setup", 
                                 fluidRow(
                                   column(width=1, offset =0,
                                          selectInput("comp_number", label = "Comp #", width = 150,
                                                      choices = list(1,2,3,4,5,6,7,8,9,10,11,12), 
                                                      selected = 1)
                                   ),
                                   column(width=1, offset =0,
                                          selectInput("select_final_data_stats", label = "Normalization", width = 150,
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=1, offset =0,
                                          numericInput("pvalue_cutoff", label="pvalue cutoff", value = .05)
                                   ),
                                   column(width=1, offset =0,
                                          numericInput("foldchange_cutoff", label="FC cutoff", value = 2)
                                   ),
                                   column(width=2, offset =0,
                                          numericInput("missing_factor", label="Measured % (decimal)", value = 0.6)
                                   ),
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_spqc", label = "SPQC Group?",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                                  )
                                        ),
                                   column(width=1, offset =0,
                                     dropdownButton(
                                           textOutput("stats_gear1"),
                                           checkboxInput("peptide_missing_filter", label = "Refilter peptides by requiring X% measured values in one group?"),
                                           numericInput("peptide_missing_factor", label="Peptide X% measured cutoff (decimal)", value = 0.8),
                                           checkboxInput("peptide_cv_filter", label = "Refilter peptides by requiring X %CV one group?"),
                                           numericInput("peptide_cv_factor", label="Peptide X CV% cutoff", value = 100),
                                           textOutput("stats_gear2"),
                                           checkboxInput("stats_spqc_cv_filter", label = "Filter by SPQC CV"),
                                            numericInput("stats_spqc_cv_filter_factor", label="SPQC %CV Cutoff", value = 50),
                                            checkboxInput("stats_comp_cv_filter", label = "Require one group CV below cuttoff"),
                                            numericInput("stats_comp_cv_filter_factor", label="Comp %CV Cutoff", value = 50),
                                            checkboxInput("pair_comp", label = "Pairwise Comparisons"),
                                            checkboxInput("stats_peptide_minimum", label = "Require minimum # of peptides per protein", value = 0),
                                            numericInput("stats_peptide_minimum_factor", label="Peptide Minimum", value = 1),
                                            h5('Extra Stats'),
                                            checkboxInput("checkbox_adjpval", label = "Include bonferroni adjusted pvalue?"),
                                            selectInput("padjust_options", label = "p.adjust method", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), 
                                                       selected = "bonferroni"),
                                            checkboxInput("checkbox_cohensd", label = "Include Cohen's D?"),
                                            checkboxInput("checkbox_cohensd_hedges", label = "Use Hedge's Correction (low N)?"),
                                            checkboxInput("checkbox_limmapvalue", label = "Include Limma Pvalue?"),
                                            h5('Final Excel Options'),
                                            checkboxInput("checkbox_report_ptm", label = "Report Only PTM?"),
                                            textInput("peptide_report_grep", label="Report PTM grep", value = "Enter value here"),
                                            checkboxInput("checkbox_report_accession", label = "Report Specific Accession(s) Only"),
                                            textInput("report_accession", label="Protein Accessions for Final Report", value = "Enter value"),
                                            circle = TRUE, status = "info", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see more options!")
                                     )
                                   ),
                                   column(width=2, offset =0,
                                          textInput("final_stats_name", label="Final Stats Excel Name", 
                                                    value = str_c("Final_CompHere_stats.xlsx"))
                                   )
                                 ),
                                 hr(),
                                 fluidRow(
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_1N", label = "Comp1_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_1D", label = "Comp1_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )                
                                   ),
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_2N", label = "Comp2_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_2D", label = "Comp2_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )
                                   ),
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_3N", label = "Comp3_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_3D", label = "Comp3_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )
                                   ),
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_4N", label = "Comp4_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_4D", label = "Comp4_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )
                                   ) ,            
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_5N", label = "Comp5_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_5D", label = "Comp5_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )
                                   ),            
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_6N", label = "Comp6_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_6D", label = "Comp6_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ) 
                                   )
                                 ),
                                 fluidRow(
                                   hr()
                                 ),
                                 fluidRow(
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_7N", label = "Comp7_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_7D", label = "Comp7_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )                     
                                   ),
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_8N", label = "Comp8_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_8D", label = "Comp8_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )  
                                   ),
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_9N", label = "Comp9_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_9D", label = "Comp9_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )  
                                   ),
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_10N", label = "Comp10_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_10D", label = "Comp10_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )  
                                   ) ,            
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_11N", label = "Comp11_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_11D", label = "Comp11_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )  
                                   ),            
                                   column(width=2, offset =0,
                                          pickerInput(inputId = "comp_12N", label = "Comp12_N",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          ),
                                          pickerInput(inputId = "comp_12D", label = "Comp12_D",  choices = "None", 
                                                      options = list(`actions-box` = TRUE,size = 10,
                                                                     `selected-text-format` = "count > 3"),  multiple = TRUE
                                          )  
                                   )
                                 ),
                                 fluidRow(align="center",
                                          hr(),
                                          #rHandsontableOutput("stat2_table"),
                                          tags$head(tags$style("#stat2_N_1{color: blue; font-size: 16px; font-style: bold;}")),
                                          
                                          actionButton("check_stats", label = "Set Comparisons", width = 300,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          br(),
                                          br(),
                                          actionButton("start_stats", label = "Start Anlaysis", width = 300,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          br(),
                                          br(),
                                          actionButton("save_stats", label = "Create Excel Output File", width = 300,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          br(),
                                          br(),
                                          downloadButton('download_stats_excel')
                                 )
                        ),
                        
                        tabPanel("Graphs", id="tp_stats_graph1", 
                                 fluidRow(
                                   column(width=3, offset =0,
                                          pickerInput(inputId = "stats_plot_comp", label = "Comparison(s) to plot",  choices = "None", 
                                                      options = list(`actions-box` = TRUE, size = 100,
                                                                     `selected-text-format` = "count > 5"),  multiple = TRUE)
                                   ),
                                   column(width=2, offset =0,
                                          checkboxInput("stats_plot_spqc", label = "Add SPQC?")
                                   ),
                                   column(width=4, offset =0,
                                          actionButton("create_stats_plots", label = "Create Plots", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   )
                                 ),
                                 fluidRow(
                                   column(width=6, offset =0,
                                          dropdownButton(
                                            textInput("stats_barplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                            textInput("stats_barplot_title", label="plot title", value = "Total Summed Intensity", width = 200),
                                            sliderInput("stats_barplot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 11),
                                            sliderInput("stats_barplot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("stats_barplot", width = 600, height = 400)
                                          ),
                                          downloadButton('download_stats_barplot')
                                   ),  
                                   
                                   column(width=6, offset =0,
                                          dropdownButton(
                                            textInput("stats_boxplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                            textInput("stats_boxplot_title", label="plot title", value = "Total Summed Intensity", width = 200),
                                            sliderInput("stats_boxplot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 11),
                                            sliderInput("stats_boxplot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("stats_boxplot", width = 600, height = 400)
                                          ),
                                          downloadButton('download_stats_boxplot')
                                   )  
                                 ),
                                 hr(),
                                 fluidRow(
                                   column(width=6, offset =0,
                                          dropdownButton(
                                            selectInput("stats_pca2d_x", label = "pca xaxis", choices = list("PC1", "PC2", "PC3", "PC4", "PC5"), 
                                                        selected = "PC1"),
                                            selectInput("stats_pca2d_y", label = "pca yaxis", choices = list("PC1", "PC2", "PC3", "PC4", "PC5"), 
                                                        selected = "PC2"),
                                            textInput("stats_pca2d_title", label="plot title", value = "pca2d", width = 200),
                                            sliderInput("stats_pca2d_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 11),
                                            sliderInput("stats_pca2d_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20),
                                            sliderInput("stats_pca2d_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 20, value = 4),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("stats_pca2d", width = 600, height = 400,
                                                       hover = hoverOpts("plot_pca2d_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("hover_pca2d_info")
                                          ),
                                          downloadButton('download_stats_pca2d')
                                   ),  
                                   
                                   column(width=6, offset =0,
                                          dropdownButton(
                                            textInput("stats_pca3d_title", label="plot title", value = "pca3d", width = 200),
                                            sliderInput("stats_pca3d_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 11),
                                            sliderInput("stats_pca3d_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20),
                                            sliderInput("stats_pca3d_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            rglwidgetOutput("stats_pca3d", width = 600, height = 400)
                                          ),
                                          downloadButton('download_stats_pca3d')
                                   )  
                                 ),
                                 hr(),
                                 fluidRow(
                                   column(width=6, offset =0,
                                          dropdownButton(
                                            textInput("stats_cluster_title", label="plot title", value = "cluster", width = 200),
                                            sliderInput("stats_cluster_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 11),
                                            sliderInput("stats_cluster_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20),
                                            colourpicker::colourInput("cluster_high_color", "Select High Color", "#FF3366"),
                                            colourpicker::colourInput("cluster_low_color", "Select Low Color", "#009933"),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("stats_cluster", width = 600, height = 400)
                                          ),
                                          downloadButton('download_stats_cluster')
                                   ),  
                                   
                                   column(width=6, offset =0,
                                          dropdownButton(
                                            textInput("stats_heatmap_title", label="plot title", value = "heatmap", width = 200),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("stats_heatmap", width = 600, height = 400)
                                          ),
                                          downloadButton('download_stats_heatmap')
                                   )  
                                 )  
                        ),   # end graphs
                        
                        tabPanel("Volcano", id="tp_stats_volcano", 
                                 fluidRow(
                                   column(width=3, offset =0,
                                          actionButton("create_stats_volcano", label = "Create Volcano Plots", width = 300,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   ), 
                                   column(width=2, offset = 0,
                                          checkboxInput("stats_volcano_fixed_axis", label = "Fix x and y axis for all plots?")
                                          ),
                                   column(width=1, offset = 0,
                                           numericInput("stats_volcano_y_axis", label = "y axis", value=10)
                                   ),                                   
                                   column(width=1, offset = 0,
                                           numericInput("stats_volcano_x_axis", label = "x axis", value=5)
                                   )
                                 ), 
                                 fluidRow(
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano1_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano1_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano1_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano1_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano1_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano1_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano1_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano1_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano1_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano1_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano1')
                                   ), # end col
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano2_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano2_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano2_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano2_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano2_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano2_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano2_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano2_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano2_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano2_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano2')
                                   ) # end col
                                 ), #end fr
                                 
                                 fluidRow(
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano3_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano3_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano3_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano3_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano3_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano3_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano3_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano3_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano3_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano3_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano3')
                                   ), # end col
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano4_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano4_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano4_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano4_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano4_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano4_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano4_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano4_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano4_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano4_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano4')
                                   ) # end col
                                 ), #end fr   
                                 
                                 fluidRow(
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano5_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano5_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano5_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano5_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano5_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano5_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano5_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano5_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano5_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano5_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano5')
                                   ), # end col
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano6_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano6_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano6_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano6_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano6_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano6_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano6_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano6_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano6_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano6_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano6')
                                   ) # end col
                                 ), #end fr   
                                 
                                 fluidRow(
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano7_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano7_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano7_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano7_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano7_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano7_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano7_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano7_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano7_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano7_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano7')
                                   ), # end col
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano8_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano8_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano8_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano8_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano8_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano8_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano8_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano8_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano8_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano8_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano8')
                                   ) # end col
                                 ), #end fr
                                 
                                 fluidRow(
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano9_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano9_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano9_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano9_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano9_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano9_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano9_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano9_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano9_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano9_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano9')
                                   ), # end col
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano10_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano10_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano10_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano10_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano10_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano10_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano10_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano10_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano10_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano10_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano10')
                                   ) # end col
                                 ), #end fr   
                                 
                                 fluidRow(
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano11_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano11_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano11_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano11_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano11_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano11_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano11_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano11_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano11_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano11_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano11')
                                   ), # end col
                                   column(width=5, offset =0,
                                          dropdownButton(
                                            textInput("volcano12_stats_plot_title", label="plot title", 
                                                      value = "Volcano", width = 200),
                                            textInput("volcano12_stats_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                            textInput("volcano12_stats_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                            colourpicker::colourInput("volcano12_stats_dot_color", "Select Color", "blue"),
                                            sliderInput("volcano12_stats_plot_dot_size", label = h5("Point Size"), min = 1, 
                                                        max = 10, value = 2),
                                            sliderInput("volcano12_stats_plot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 20),
                                            sliderInput("volcano12_stats_plot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20), 
                                            
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("volcano12_stats_plot", width = 600, height = 600,
                                                       hover = hoverOpts("volcano12_stats_hover", delay = 100, delayType = "debounce")),
                                            uiOutput("volcano12_stats_hover_info")
                                          ),
                                          downloadButton('download_stats_volcano12')
                                   ) # end col
                                 ) #end fr   
                        ), # end volcano
                        

                        
                        tabPanel("Data", value = "tp_stats_data",
                                 fluidRow( 
                                   column(width=1, offset =0,  
                                          numericInput("stats_data_topn", label="TopN", value = 0, width = 100)
                                   ),
                                   column(width=1, offset =0,
                                          textInput("stats_data_accession", label="Accession", value = "0", width = 100)
                                   ),
                                   column(width=1, offset =0,
                                          textInput("stats_data_description", label="Description", value = "0", width = 200)
                                   ),
                                   column(width=3, offset =0,
                                          selectInput("stats_select_data_comp", label = "comparison", 
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=2, offset =0,
                                          checkboxInput("stats_add_filters", label="Apply stat filters (from setup)?", value = 0)
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("stats_data_show", label = "Filter Data", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   ),
                                   column(width=2, offset =0,
                                          textInput("stats_data_filename", label="File Name", value = "my_data.xlsx", width = 250)
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("stats_data_save", label = "Save Data", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   )
                                 ),
                                 
                                 fluidRow(
                                  verbatimTextOutput('stats_data_final_protein')
                                 ),
                                 
                                 fluidRow(
                                   column(width=12, offset =0,
                                   hr(),
                                   tags$head(tags$style("#stats_data_final{color: blue;
                                                           font-size: 12px;
                                                           }"
                                   )
                                   ),
                                   DT::dataTableOutput("stats_data_final", width='100%')
                                 )
                              ),
                              br(),
                              
                              fluidRow(
                                     actionButton("stats_data_update", label = "Remove from Stats", width = 200,
                                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                              )
                              
                              
                              
                        ), #end data
               
                        
                        tabPanel("Protein Plots", value = "tp_stats_oneprotein",
                                 fluidRow(
                                   column(width=1, offset =0,
                                          textInput("stats_oneprotein_accession", label="Accession", value = "0", width = 100)
                                   ),
                                   column(width=3, offset =0,
                                          selectInput("stats_oneprotein_plot_comp", label = "comparison", 
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   
                                   column(width=1, offset =0,
                                          checkboxInput("stats_oneprotein_plot_spqc", label = "Add SPQC?")
                                   ),
                                   column(width=1, offset =0,
                                          checkboxInput("stats_use_zscore", label = "Use zscore?")
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("create_stats_oneprotein_plots", label = "Create Plots", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   ),
                                   column(width=2, offset =1,
                                          textInput("stats_oneprotein_data_filename", label="File Name", value = "my_protein_data.xlsx", width = 250)
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("stats_oneprotein_data_save", label = "Save Data", width = 100,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   )
                                 ),
                                 fluidRow(
                                   column(width=12, offset =0,
                                          dropdownButton(
                                            textInput("stats_oneprotein_barplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                            textInput("stats_oneprotein_barplot_title", label="plot title", value = "Total Summed Intensity", width = 200),
                                            sliderInput("stats_oneprotein_barplot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 11),
                                            sliderInput("stats_oneprotein_barplot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("stats_oneprotein_barplot", width = 1200, height = 400)
                                          ),
                                          downloadButton('download_stats_oneprotein_barplot')
                                   )
                                   ),
                                  fluidRow(
                                   column(width=12, offset =0,
                                          dropdownButton(
                                            textInput("stats_oneprotein_grouped_barplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                            textInput("stats_oneprotein_grouped_barplot_x_axis_label", label="x axis label", value = "Sequence", width = 200),
                                            textInput("stats_oneprotein_grouped_barplot_title", label="plot title", value = "Total Summed Intensity", width = 200),
                                            sliderInput("stats_oneprotein_grouped_barplot_label_size", label = h5("Label Size"), min = 1, 
                                                        max = 50, value = 11),
                                            sliderInput("stats_oneprotein_grouped_barplot_title_size", label = h5("Title Size"), min = 10, 
                                                        max = 50, value = 20),
                                            circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                            tooltip = tooltipOptions(title = "Click to see inputs !")
                                          ),
                                          div(
                                            style = "position:relative",
                                            plotOutput("stats_oneprotein_grouped_barplot", width = 1200, height = 400)
                                          ),
                                          downloadButton('download_stats_oneprotein_grouped_barplot')
                                   )
                                   
                                   
                                   ),
                                 
                                 fluidRow(
                                   column(width=12, offset =0,
                                          hr(),
                                          tags$head(tags$style("#oneprotein_peptide_table{color: blue;
                                                           font-size: 12px;
                                                           }"
                                          )
                                          ),
                                          DT::dataTableOutput("oneprotein_peptide_table", width='100%')
                                   )
                                 )
                                 # 
                                 # 
                                 # fluidRow(
                                 #  br(),
                                 #  hr(),
                                 #  rHandsontableOutput("oneprotein_peptide_table")
                                 # )
              ), #end tab panel
              
              tabPanel("Peptide Plots", value = "tp_stats_onepeptide",
                       fluidRow(
                         column(width=1, offset =0,
                                textInput("stats_onepeptide_accession", label="Accession", value = "0", width = 100)
                         ),
                         column(width=2, offset =0,
                                textInput("stats_onepeptide_sequence", label="Peptide Sequence", width = 250)
                         ),
                         column(width=2, offset =0,
                                textInput("stats_onepeptide_modification", label="Peptide Modification", width = 250)
                         ),
                         column(width=3, offset =0,
                                selectInput("stats_onepeptide_plot_comp", label = "comparison", 
                                            choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                            selected = 1)
                         ),
                         
                         column(width=1, offset =0,
                                checkboxInput("stats_onepeptide_plot_spqc", label = "Add SPQC?")
                         ),
                         column(width=1, offset =0,
                                checkboxInput("stats_onepeptide_use_zscore", label = "Use zscore?")
                         ),
                         column(width=1, offset =0,
                              numericInput("stats_onepeptide_residue", label="Residue", value = 0, width = 100)
                         ),
                         column(width=1, offset =0,
                                actionButton("create_stats_onepeptide_plots", label = "Create Plots", width = 100,
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                         )
                       ),
                       fluidRow(
                         column(width=12, offset =0,
                                dropdownButton(
                                  textInput("stats_onepeptide_barplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                  textInput("stats_onepeptide_barplot_title", label="plot title", value = "Total Summed Intensity", width = 200),
                                  sliderInput("stats_onepeptide_barplot_label_size", label = h5("Label Size"), min = 1, 
                                              max = 50, value = 11),
                                  sliderInput("stats_onepeptide_barplot_title_size", label = h5("Title Size"), min = 10, 
                                              max = 50, value = 20),
                                  circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                  tooltip = tooltipOptions(title = "Click to see inputs !")
                                ),
                                div(
                                  style = "position:relative",
                                  plotOutput("stats_onepeptide_barplot", width = 1200, height = 400)
                                ),
                                downloadButton('download_stats_onepeptide_barplot')
                         )
                       ),
                       fluidRow(
                         column(width=12, offset =0,
                                dropdownButton(
                                  textInput("stats_onepeptide_grouped_barplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                  textInput("stats_onepeptide_grouped_barplot_x_axis_label", label="x axis label", value = "Sequence", width = 200),
                                  textInput("stats_onepeptide_grouped_barplot_title", label="plot title", value = "Total Summed Intensity", width = 200),
                                  sliderInput("stats_onepeptide_grouped_barplot_label_size", label = h5("Label Size"), min = 1, 
                                              max = 50, value = 11),
                                  sliderInput("stats_onepeptide_grouped_barplot_title_size", label = h5("Title Size"), min = 10, 
                                              max = 50, value = 20),
                                  sliderInput("stats_onepeptide_grouped_barplot_axis_angle", label = h5("X-axis angle"), min = 0, 
                                              max = 90, value = 45),
                                  numericInput("stats_onepeptide_grouped_barplot_axis_vjust", label="X-axis vjust", value = 0.5),
                                  circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                  tooltip = tooltipOptions(title = "Click to see inputs !")
                                ),
                                div(
                                  style = "position:relative",
                                  plotOutput("stats_onepeptide_grouped_barplot", width = 1200, height = 400)
                                ),
                                downloadButton('download_stats_onepeptide_grouped_barplot')
                         )
                       ),
                       
                       fluidRow(
                         column(width=12, offset =0,
                                hr(),
                                tags$head(tags$style("#onepeptide_peptide_table{color: blue;
                                                           font-size: 12px;
                                                           }"
                                )
                                ),
                                DT::dataTableOutput("onepeptide_peptide_table", width='100%')
                         ),
                         column(width=2, offset =0,
                                textInput("stats_onepeptide_data_filename", label="File Name", value = "my_peptide_data.xlsx", width = 250)
                         ),
                         column(width=1, offset =0,
                                actionButton("stats_onepeptide_data_save", label = "Save Data", width = 100,
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                         )
                       )

              ) #end tab panel              
              
              
              
              
              
             ) # end of navbar panel
    ),  #end of tab panel
    
    
    
    
    
    
    tabPanel("Pathway", value = "pathway",
      navbarPage("Pathway:", id ="path",
        tabPanel("Organism", value = "tp_save", align="center",
                 tags$h1("Select organism for pathway analysis/enrichment..."),
                 hr(),
                 selectInput("select_organism", label = "organism", 
                             choices = list("Human", "Mouse", "Rat", "Danio", "Arabidopsis", "Ecoli"), 
                             selected = "Human"),
                 br(),
                 actionButton("set_pathway", label = "Set Pathway", width = 300, 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        ),
        
        tabPanel("WikiPathways", id="wiki",
                 fluidRow( 
                   column(width=3, offset =0,
                          selectInput("select_data_comp_wiki", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=3, offset =0,  
                          radioButtons("wiki_direction", label="Fold Change Direction", choices = list("Up", "Down", "UpDown"),  selected = "Up", width = 200)
                   ),
                   column(width=1, offset =1,
                          actionButton("wiki_show", label = "Find WikiPathway", width = 150,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                   )
                 ),
                 
                 fluidRow(
                   hr(),
                  tags$head(tags$style("#data_final{color: blue;
                                 font-size: 12px;
                                  }"
                     )
                   ),
                   rHandsontableOutput("wiki_table")
                   ),
                  br(),
                 downloadButton('download_wiki_table')
          ),
        
        tabPanel("Go Profile", id="go_profile",
                 fluidRow( 
                   column(width=2, offset =0,
                          selectInput("select_data_comp_profile", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=3, offset =0,  
                          radioButtons("profile_direction", label="Fold Change Direction", choices = list("Up", "Down", "UpDown"),  selected = "Up", width = 200)
                   ),
                   column(width=2, offset =0,
                          selectInput("select_ont_profile", label = "ontology", 
                                      choices = list("CC", "BP", "MF"), 
                                      selected = "BP")
                   ),
                   column(width=2, offset =0,
                          selectInput("select_level_profile", label = "ontology level", 
                                      choices = list("1", "2", "3", "4", "5", "6"), 
                                      selected = "3")
                   ),
                   column(width=1, offset =0,
                          actionButton("profile_show", label = "Find Go Profile", width = 150,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                   )
                 ),  
                 fluidRow(
                   hr(),
                   column(width=6, offset =0,
                          tags$head(tags$style("#data_final{color: blue;
                                     font-size: 12px;
                                     }"
                          )
                          ),
                          rHandsontableOutput("go_profile_table"),
                          br(),
                          downloadButton('download_go_profile_table')
                   ),
                   column(width=6, offset =0,
                          plotOutput("go_profile_plot"),
                          br(),
                          downloadButton('download_go_profile_plot')
                   )
                 )
        ),
        
        tabPanel("Go Analysis", id="go",
                 fluidRow( 
                   column(width=3, offset =0,
                          selectInput("select_data_comp_go", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=3, offset =0,  
                          radioButtons("go_direction", label="Fold Change Direction", choices = list("Up", "Down", "UpDown"),  selected = "Up", width = 200)
                   ),
                   column(width=1, offset =0,
                          selectInput("select_ont_go", label = "ontology", 
                                      choices = list("CC", "BP", "MF"), 
                                      selected = "BP")
                   ),
                   # column(width=1, offset =0,
                   #        selectInput("select_algorithm", label = "algorithm", 
                   #                    choices = list("classic", "elim", "weight", "weight01"), 
                   #                    selected = "classic")
                   # ),
                   column(width=1, offset =0,
                          actionButton("go_show", label = "Go Analysis", width = 100,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                   )
                 ),
                 
                 fluidRow(
                   hr(),
                     tags$head(tags$style("#data_final{color: blue;
                                 font-size: 12px;
                                  }"
                          )
                          ),
                          rHandsontableOutput("go_table"),
                           downloadButton('download_go_table')
                 )
           ), #end of go analysis
        
        tabPanel("Go Volcano", id="go_volcano",
                 fluidRow(
                   column(width=2, offset =0,
                          textInput("go_volcano_id", label="GO ID", value = "", width = 200),
                          textInput("plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                          textInput("plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                          textInput("plot_title", label="plot title", value = "Go Volcano", width = 200),
                          colourpicker::colourInput("volcano_dot_color", "Select Color", "blue"),
                          actionButton("go_volcano", label = "Volcano", width = 100,
                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                          ),  
                   column(width=2, offset =0,
                          sliderInput("plot_dot_size", label = h5("Point Size"), min = 1, 
                                      max = 10, value = 2),
                          sliderInput("plot_label_size", label = h5("Label Size"), min = 1, 
                                      max = 50, value = 20),
                          sliderInput("plot_title_size", label = h5("Title Size"), min = 10, 
                                      max = 50, value = 20),   
                          downloadButton('download_go_volcano')
                   ),
                   column(width=8, offset =0,  
                     div(
                       style = "position:relative",
                       plotOutput("volcano_go_plot", width = 800, height = 600,
                                  hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")),
                       uiOutput("hover_info")
                     )
                   )
                 ),
                 fluidRow(
                   column(width=12, offset =0,
                          hr(),
                          tags$head(tags$style("#volcano_data_final{color: blue;
                                                           font-size: 12px;
                                                           }"
                          )
                          ),
                          DT::dataTableOutput("volcano_data_final", width='100%'),
                          downloadButton('download_go_volcano_table')
                   )
                 )
        ), #end of go volcano 
        
        
        tabPanel("StringDB", id="string",
                 fluidRow( 
                   column(width=2, offset =0,
                          actionButton("setup_string", label = "String Setup", width = 150,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                   ),
                   column(width=3, offset =0,
                          selectInput("select_data_comp_string", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=3, offset =0,  
                          radioButtons("string_direction", label="Fold Change Direction", choices = list("Up", "Down", "UpDown"),  selected = "Up", width = 200)
                   ),
                   column(width=1, offset =0,
                          selectInput("protein_number", label = "Max #Proteins", 
                                      choices = list(10, 50, 100, 150, 200, 250, 300, 350, 400), 
                                      selected = 200)
                   ),                   
                   column(width=2, offset =0,
                          actionButton("go_string", label = "String Analysis", width = 150,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                   )
                 ),
                 fluidRow(
                   hr(),
                   # string link depreciated - save incase comes back
                   # tags$head(tags$style("#string_link{color: blue;
                   #               font-size: 12px;
                   #                }"
                   # )
                   # ),
                   # textOutput("string_link"),
                   downloadButton('download_string_plot'),
                   plotOutput("string_plot")  
                   #rHandsontableOutput("string_table")
                 )
        ), #end of string analysis
        
        tabPanel("StringDB_enrich", id="string_enrich",
                 fluidRow( 
                   column(width=3, offset =0,
                          selectInput("select_data_comp_string_enrich", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=3, offset =0,  
                          radioButtons("string_enrich_direction", label="Fold Change Direction", choices = list("Up", "Down", "UpDown"),  selected = "Up", width = 200)
                   ),
                   column(width=1, offset =0,
                          selectInput("select_string_enrich", label = "Enrichment", 
                                      choices = list("Process", "Component", "Function", "KEGG", "Pfam", "InterPro"), 
                                      selected = "Process")
                   ),
                   column(width=2, offset =0,
                          actionButton("go_string_enrich", label = "String Enrichment", width = 150,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                   )
                 ),
                 
                 fluidRow(
                   hr(),
                   tags$head(tags$style("#data_final{color: blue;
                                 font-size: 12px;
                                  }"
                   )
                   ),
                   rHandsontableOutput("string_table"),
                   downloadButton('download_string_enrich_table')
                 )
        ) #end of string analysis        
        
        
      ) # end of navbarpage
    ), # end of tab panel Pathway
    
    
    
    tabPanel("Phos", value = "tp_phos",
             navbarPage("Phosphorylation:", id ="np_phos",
                        tabPanel("Format Fasta", id="motif",
                                 fluidRow(
                                   shinyFilesButton('motif_fasta', label='Select Motif-X FASTA', title='Please select motif-x formated text file', multiple=FALSE,
                                                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                   br(),
                                   textOutput("fasta"),
                                   tags$head(tags$style("#fasta{color: blue; font-size: 16px; font-style: bold;}")),
                                   br(),
                                   selectInput("accession_split", label="Split character after accession:", 
                                               choices = list("Space", "Bar"),
                                               selected = "Bar"),
                                   br(),
                                   actionButton("parse_fasta", label = "Format fasta", width = 150,
                                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),  
                                   br(),
                                   br(),
                                   
                                   rHandsontableOutput("fasta_example")
                                 )
                        ),
                        tabPanel("MotifX", id="motif",
                                 fluidRow(
                                   column(width=2, offset =0,
                                          selectInput("select_data_comp_motif", label = "comparison", 
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=2, offset =0,  
                                          numericInput("pval_motif", label="MotifX pval (x1e-5)", value =1, width = 150)
                                   ),
                                   column(width=3, offset =0,  
                                          textInput("protein_filter", label="Specific Protein Filter Accession", value ="", width = 250)
                                   ),
                                   column(width=1, offset =0,
                                          actionButton("motif_show", label = "Send to MotifX", width = 150,
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                   )
                                 ),
                                 
                                 fluidRow(
                                   hr(),
                                   tags$head(tags$style("#data_final{color: blue;
                                 font-size: 12px;
                                 }"
                                   )
                                   ),
                                   rHandsontableOutput("motif_table"),
                                   br(),
                                   br(),
                                   downloadButton('download_motifx_excel')
                                 )
                        )           
                        
             )
    ), # end of tab panel Phos  
   
   
    
    
    tabPanel("Save", value = "tp_save", align="center",
               tags$h1("Save dpmsr_set file..."),
               hr(),
               br(),
               textInput("dpmsr_set_name", label="File Name", value ="dpmsr_set_filename", width = 250),
               br(),
               actionButton("save_dpmsr_set", label = "Create new dpmsr_set", width = 300,
                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
              br(),
              br(),
              br(),
              downloadButton('download_new_dpmsr_set'),
             hr(),
             br(),
             textInput("dpmsr_set_name_customer", label="File Name", value ="dpmsr_set_customer_filename", width = 250),
             br(),
             actionButton("save_customer_dpmsr_set", label = "Create customer dpmsr_set", width = 300,
                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
             br(),
             br(),
             br(),
             
             # shinyDirButton('new_save_directory', label='Choose Directory', title='Please select directory for files',
             #            style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
             # actionButton("test_new_save", label = "test new save", width = 300,
             #              style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
             
    ),# end of tab panel Save   
    
    
    tabPanel("Report", value = "tp_report", align="center",
             tags$h1("Create report from templates..."),
             hr(),
             column(width=3, offset =0,
                    selectInput("select_final_data_report", label = "Normalization", width = 150,
                                choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                selected = 1),
                    br(),
                    br(),
                    actionButton("create_report", label = "Create Report", width = 200,
                                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                    ),
             column(width=9, offset =0, align="left",
                    textOutput("report1_output")
             )
    ) # end of tab panel Save      
    
    
    
    
    
     
 ) #end of navlistpanel
))) #end of ui
