library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinyWidgets)
library(rhandsontable)
library(rgl)

shinyUI(fluidPage(
  useShinyjs(),
  setBackgroundColor("DarkGray", shinydashboard = TRUE),
  setBackgroundColor("LightGray", shinydashboard = FALSE),
  titlePanel("Proteome Discoverer Data Processing"),
  
  navlistPanel(widths=c(1,11), id = "nlp1",
          
    tabPanel("Load Design", value = "tp1", align="center",
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
                        fluidRow(
                         radioButtons("radio_input", label = h3("Select Input Data Type"),
                                      choices = list("Protein_Peptide" = 1, "Protein" = 2, "Peptide" = 3),
                                      selected = 1),
                         br(),
                         hr(),
                         checkboxInput("checkbox_isoform", label = "Use Peptide Isoform?"),
                         checkboxInput("checkbox_tmt", label = "SPQC Normalized TMT sets"),
                         selectInput("razor", label = h5("Peptides to Use"), 
                                     choices = list("Razor", "Unique", "Shared"), 
                                     selected = "Razor"),
                    )),
                     column(width=5, offset =0,
                            radioButtons("radio_output", label = h3("Select Output Data Type"),
                                         choices = list("Protein" = 1, "Peptide" = 2),
                                         selected = 1),
                            br(),
                            br(),
                            hr()
                           )
                       ),

            fluidRow(
                     hr(),
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
               imageOutput("peptide_PI"), 
               imageOutput("peptide_cruc"), 
               imageOutput("faims_psm_cv")
            ),
            column(width=3, offset =1,
                imageOutput("peptide_RT"),
                imageOutput("peptide_MZ"),
                imageOutput("peptide_Charge"),
                imageOutput("peptide_aainfo"), 
                imageOutput("peptide_ai"), 
               ),
            column(width=3, offset =1,
                   rHandsontableOutput("project_overview")   
            )
          )
    ),  


      tabPanel("Filters", value = "tp4", align="center",
               tags$h1("Apply Filters"),
               hr(),
               fluidRow(
                 column(width=5, offset = 1,
                        fluidRow(align = "left",
                          textOutput("text1"),
                          tags$head(tags$style("#text1{color: blue; font-size: 16px; font-style: bold;}")),
                          textOutput("text2"),
                          tags$head(tags$style("#text2{color: blue; font-size: 16px; font-style: bold;}")),
                          textOutput("text3"),
                          tags$head(tags$style("#text3{color: blue; font-size: 16px; font-style: bold;}")),
                          textOutput("text4"),
                          tags$head(tags$style("#text4{color: blue; font-size: 16px; font-style: bold;}")),
                           )
                        ),
                  column(width=5, offset = 1,
                    fluidRow(align = "left",
                          "Fixed:  entry must have at least 2 data points",
                          hr(),
                          checkboxInput("checkbox_require2", label = "Require two data points in one group"),
                          hr(),
                          checkboxInput("checkbox_filtercv", label = "Filter on Specific Group CV"),
                          textInput("text_filtercvgroup", label="Enter group for CV filter", value = "Enter value here"),
                          numericInput("number_filtercvcutoff", label="Enter CV% for cutoff", value = "Enter value here"),
                          hr(),
                          br(),
                          actionButton("filter_apply", label = "Apply Filters", width = 300,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                       )
                      )
                 )
               ), #end tab panel


      tabPanel("Normalize", value = "tp5",
               fluidRow(
                 column(width=5, offset =1,
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
                     checkboxInput("checkbox_norm_ptm", label = "Normalized on Specific Modification Only"),
                     textInput("filter_grep", label="Filter(grep) for Modification", value = "Enter value here"),
                     br(),
                     actionButton("norm1", label = "Apply Normalization", width = 300,
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                     ),
               column(width=6, offset =0,
                      "Raw Data Total Intensity",
                      br(),
                      imageOutput("raw_bar")
                      
               )
          )
       ),


    tabPanel("Impute", value = "tp6",
             fluidRow(
               column(width=5, offset =1,
                      textOutput("text_i1"),
                      tags$head(tags$style("#text_i1{color: blue; font-size: 20px; font-style: bold;}")),
                      br(),
                      radioButtons("radio_impute", label=NULL,
                                   choices = list("Duke" = 1, "Floor" = 2, "Minimum" = 3,"Average" = 4,
                                                  "KNN"= 5, "LocalLeastSquares" = 6, "MLE" = 7, "BottomX" = 8
                                   ),
                                   selected = 1),
                      numericInput("bottom_x", label="Bottom X%", value = "5"),
                      hr(),
                      checkboxInput("checkbox_impute_ptm", label = "Impute with Data from Specific Modification Only"),
                      hr(),
                      numericInput("impute_floor", label="Impute Floor Intensity", value = "Enter value here"),
                      checkboxInput("checkbox_missing_50", label = "If >50% values missing from a group then apply intensity cutoff"),
                      textOutput("text_i2"),
                      hr(),
                      actionButton("impute1", label = "Apply Imputation", width = 300,
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                    ),
               column(width=6, offset =0,
                         "Total Intensity Histogram",
                      br(),
                      imageOutput("histogram")

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
               hr(),
               actionButton("tmt_irs_go", label = "Apply IRS Norm", width = 200,
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
    
    
    tabPanel("Stats", value = "tp7",
             fluidRow(
               column(width=3, offset =0,
                 sliderInput("comp_number", label = h5("Enter Number of Comparisons"), min = 1, 
                             max = 6, value = 2)
               ),
               column(width=3, offset =0,
                  numericInput("pvalue_cutoff", label="Enter pvalue cutoff", value = "Enter value here"),
                  numericInput("foldchange_cutoff", label="Enter folchange cutoff", value = "Enter value here")
               ),
               column(width=3, offset =0,
                 checkboxInput("pair_comp", label = "Use Pairwise Comparisons"),
                 checkboxInput("checkbox_out_ptm", label = "Report Specific Modification Only"),
                 textInput("report_grep", label="Filter(grep) for Modification", value = "Enter value here"),
               ),
               column(width=3, offset =0,
                      checkboxInput("checkbox_report_accession", label = "Report Specific Accession(s) Only"),
                      textInput("report_accession", label="Protein Accessions for Final Report", value = "Enter value here"),
               ) 
               
               
             ),
             hr(),
             fluidRow(
               column(width=4, offset =0,
                 selectInput("comp_1N", label = h5("Comp1 Numerator/Denominator"), 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_1D", label = NULL, 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_4N", label = h5("Comp4 Numerator/Denominator"), 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_4D", label = NULL, 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1)
                  ),
               column(width=4, offset =0,
                 selectInput("comp_2N", label = h5("Comp2 Numerator/Denominator"), 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_2D", label = NULL, 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_5N", label = h5("Comp5 Numerator/Denominator"), 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_5D", label = NULL, 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1)
               ),
               column(width=4, offset =0,
                 selectInput("comp_3N", label = h5("Comp3 Numerator/Denominator"), 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_3D", label = NULL, 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_6N", label = h5("Comp6 Numerator/Denominator"), 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1),
                 selectInput("comp_6D", label = NULL, 
                             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                             selected = 1)
               ),
             ),
             fluidRow(align="center",
               hr(),
               actionButton("start_stats", label = "Apply Plots/Stats", width = 300,
                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                     )
             ),

    
    tabPanel("QC", value="tp8",
             navbarPage("QC:", id ="np5",
                        tabPanel("CV",
                                 rHandsontableOutput("data_CV"),
                                 br(),
                                 imageOutput("cv_plot")
                                 ),  
                        tabPanel("ADH",
                                 imageOutput("adh_plot")),
                        tabPanel("Protein Spike",
                                 rHandsontableOutput("data_proteinQC"),
                                 hr(),
                                 fluidRow(
                                  column(width=2, offset =0,
                                     selectInput("norm_type", label = h5("Normalization Type"), 
                                                 choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                 selected = 1)),
                                  column(width=10, offset =0,
                                     imageOutput("qc_spike_plot"))
                                 ))
             
             )),

    tabPanel("Plots", value = "tp9",
             navbarPage("Plots:", id ="np6",
                      tabPanel("Norm Comparison",
                                 fluidRow(
                                   column(width=2, offset =0,
                                     selectInput("plot_select", label = NULL, 
                                                 choices = list("barplot", "boxplot", "cluster", "pca2d", "heatmap",
                                                                "MDS", "density", "heatmap"), 
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
                                   ),
                                 )),
                        
      
                        
                        tabPanel("Volcano",
                                 fluidRow(
                                   selectInput("volcano_select", label = NULL, 
                                                 choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                 selected = 1),
                                   ),
                                 fluidRow(
                                   column(width=4, offset =0,
                                     imageOutput("volcano_plot1"),
                                     imageOutput("volcano_plot4")
                                   ),
                                   column(width=4, offset =0,
                                      imageOutput("volcano_plot2"),
                                      imageOutput("volcano_plot5")
                                   ),
                                     column(width=4, offset =0,  
                                     imageOutput("volcano_plot3"),
                                     imageOutput("volcano_plot6")
                             
                                 )),
                         ),
                        
                        
                        tabPanel("Selected Proteins",
                                 fluidRow(
                                   column(width=2, offset =0,
                                          div(style = "font-size:12px;",
                                          selectInput("protein_select", label = "Select barplot type", 
                                                      choices = list("ADH", "BiraA", "Bait", "Carbox", "Avidin", 
                                                                     "Protein1", "Protein2", "Protein3", "Protein4"), 
                                                      selected = "ADH"),
                                          
                                          textInput("adh_list", label="ADH Accession", value = "Enter value here"),
                                          textInput("bait_list", label="Bait Accession", value = "Enter value here"),
                                          textInput("avidin_list", label="Avidin Accession", value = "Enter value here"),
                                          textInput("carbox_list", label="Carbox Accession", value = "Enter value here"),
                                          textInput("bira_list", label="BirA Accession", value = "Enter value here"),
                                          textInput("protein1_list", label="Protein1 Accession", value = "Enter value here"),
                                          textInput("protein2_list", label="Protein2 Accession", value = "Enter value here"),
                                          textInput("protein3_list", label="Protein3 Accession", value = "Enter value here"),
                                          textInput("protein4_list", label="Protein4 Accession", value = "Enter value here"),
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
                                   ),
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
                                        imageOutput("oneprotein_plot_select"),
                                  ), 
                                 column(width=4, offset =1,
                                        rHandsontableOutput("oneprotein_stats"),
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
                                        imageOutput("onepeptide_plot_select"),
                                 ), 
                                 column(width=4, offset =1,
                                        rHandsontableOutput("onepeptide_stats"),
                                 )
                               ))
                        
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
                numericInput("data_pvalue", label="pvalue", value = 0, width = 100)
            ),
            column(width=1, offset =0,  
                numericInput("data_foldchange1", label="foldchange up", value = 0, width = 100)
            ),
            column(width=1, offset =0,  
                   numericInput("data_foldchange2", label="foldchange dn", value = 0, width = 100)
            ),
            column(width=3, offset =0,
                  selectInput("select_data_comp", label = "comparison", 
                           choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                           selected = 1)
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
          ),
    
    
    
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
                   
                   column(width=2, offset =0,
                          selectInput("select_final_data_wiki", label = "Normalization", width = 150,
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=1, offset =0,  
                          numericInput("pval_filter_wiki", label="  pvalue", value = 0.05, width = 100)
                   ),
                   column(width=1, offset =0,  
                          numericInput("fc_filter_wiki", label="    FC", value = 2, width = 100)
                   ),
                   
                   column(width=2, offset =0,
                          selectInput("select_data_comp_wiki", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=1, offset =0,
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
                   )
          ),
        
        tabPanel("Go Profile", id="go_profile",
                 fluidRow( 
                   
                   column(width=2, offset =0,
                          selectInput("select_final_data_profile", label = "Normalization", width = 150,
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=1, offset =0,  
                          numericInput("pval_filter_profile", label="  pvalue", value = 0.05, width = 100)
                   ),
                   column(width=1, offset =0,  
                          numericInput("fc_filter_profile", label="    FC", value = 2, width = 100)
                   ),
                   
                   column(width=2, offset =0,
                          selectInput("select_data_comp_profile", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
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
                          rHandsontableOutput("go_profile_table")
                   ),
                   column(width=6, offset =0,
                          plotOutput("go_profile_plot")               
                   )
                 )
        ),
        
        tabPanel("Go Analysis", id="go",
                 fluidRow( 
                   
                   column(width=1, offset =0,
                          selectInput("select_final_data_go", label = "Normalization", width = 150,
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=1, offset =0,  
                          numericInput("pval_filter_go", label="  pvalue", value = 0.05, width = 200)
                   ),
                   column(width=1, offset =0,  
                          numericInput("fc_filter_go", label="    FC", value = 2, width = 200)
                   ),
                   
                   column(width=3, offset =0,
                          selectInput("select_data_comp_go", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
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
                          rHandsontableOutput("go_table")
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
                          downloadButton('Download')
                   ),
                   column(width=8, offset =0,  
                     div(
                       style = "position:relative",
                       plotOutput("volcano_go_plot", width = 800, height = 600,
                                  hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")),
                       uiOutput("hover_info")
                     )
                   )
                 )
        ), #end of go volcano 
        
        
        tabPanel("StringDB", id="string",
                 fluidRow( 
                   column(width=1, offset =0,
                          selectInput("select_final_data_string", label = "Normalization", width = 150,
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=2, offset =0,
                          actionButton("setup_string", label = "String Setup", width = 150,
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                   ),
                   column(width=1, offset =0,  
                          numericInput("pval_filter_string", label="  pvalue", value = 0.05, width = 200)
                   ),
                   column(width=1, offset =0,  
                          numericInput("fc_filter_string", label="    FC", value = 2, width = 200)
                   ),
                   column(width=3, offset =0,
                          selectInput("select_data_comp_string", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=1, offset =0,
                          selectInput("protein_number", label = "#Proteins", 
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
                   tags$head(tags$style("#string_link{color: blue;
                                 font-size: 12px;
                                  }"
                   )
                   ),
                   textOutput("string_link"),
                   plotOutput("string_plot")  
                   #rHandsontableOutput("string_table")
                 )
        ), #end of string analysis
        
        tabPanel("StringDB_enrich", id="string_enrich",
                 fluidRow( 
                   column(width=1, offset =0,  
                          numericInput("pval_filter_string_enrich", label="  pvalue", value = 0.05, width = 200)
                   ),
                   column(width=1, offset =0,  
                          numericInput("fc_filter_string_enrich", label="    FC", value = 2, width = 200)
                   ),
                   
                   column(width=3, offset =0,
                          selectInput("select_data_comp_string_enrich", label = "comparison", 
                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                      selected = 1)
                   ),
                   column(width=1, offset =0,
                          selectInput("select_string_enrich", label = "Enrichment", 
                                      choices = list("Process", "Component", "Function", "KEGG", "Pfam", "InterPro"), 
                                      selected = "Process")
                   ),
                   column(width=1, offset =0,
                          selectInput("select_methodMT", label = "methodMT", 
                                      choices = list("fdr", "bonferroni"), 
                                      selected = "fdr")
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
                   rHandsontableOutput("string_table") 
                 )
        ) #end of string analysis        
        
        
      ) # end of navbarpage
    ), # end of tab panel Pathway
    
    
    
    tabPanel("Phos", value = "tp_phos",
             navbarPage("Phosphorylation:", id ="np_phos",
                        tabPanel("MotifX", id="motif",
                                 fluidRow( 
                                   
                                   column(width=1, offset =0,
                                          selectInput("select_final_data_motif", label = "Normalization", width = 150,
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=1, offset =0,  
                                          numericInput("pval_filter", label="  pvalue", value = 0.05, width = 100)
                                   ),
                                   column(width=1, offset =0,  
                                          numericInput("fc_filter", label="    FC", value = 2, width = 100)
                                   ),
                                   
                                   column(width=2, offset =0,
                                          selectInput("select_data_comp_motif", label = "comparison", 
                                                      choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                                      selected = 1)
                                   ),
                                   column(width=1, offset =0,  
                                          textInput("protein_filter", label="Specific Protein", value ="", width = 150)
                                   ),
                                   column(width=2, offset =0,
                                          shinyFilesButton('motif_fasta', label='Select Motif-X FASTA', title='Please select motif-x formated text file', multiple=FALSE,
                                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          textOutput("fasta"),
                                          tags$head(tags$style("#fasta{color: blue; font-size: 16px; font-style: bold;}"))
                                   ),
                                   column(width=2, offset =0,  
                                          numericInput("pval_motif", label="MotifX pval (x1e-5)", value =1, width = 150)
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
                                   rHandsontableOutput("motif_table")
                                 )
                        )           
                        
             )
    ), # end of tab panel Extra   
   
    tabPanel("MVA", value = "tp_mva",
        navbarPage("MVA:", id ="tp_,mva2",
            tabPanel("Setup", id="tp_mva_setup", 
             fluidRow(
               column(width=1, offset =0,
                      selectInput("mva_comp", label = "Comp #", width = 150,
                                  choices = list(1,2,3,4,5,6,7,8,9,10,11,12), 
                                  selected = 1)
               ),
               column(width=1, offset =0,
                      selectInput("select_final_data_mva", label = "Normalization", width = 150,
                                  choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                  selected = 1)
               ),
               column(width=1, offset =0,
                      numericInput("mva_pvalue_cutoff", label="pvalue cutoff", value = .05)
               ),
               column(width=1, offset =0,
                      numericInput("mva_foldchange_cutoff", label="FC cutoff", value = 2)
               ),
               column(width=2, offset =0,
                      checkboxInput("mva_pair_comp", label = "Use Pairwise Comparisons"),
                      checkboxInput("mva_checkbox_out_ptm", label = "Report Specific Modification Only")
               ),
               column(width=2, offset =0,
                      textInput("mva_report_grep", label="Filter(grep) for Modification", value = "Enter value")
               ),
               column(width=2, offset =0,
                      checkboxInput("mva_checkbox_report_accession", label = "Report Specific Accession(s) Only")
               ), 
               column(width=2, offset =0,
                      textInput("mva_report_accession", label="Protein Accessions for Final Report", value = "Enter value")
               ) 
             ),
             hr(),
             fluidRow(
               column(width=2, offset =0,
                      pickerInput(inputId = "var_1N", label = "Comp1_N",  choices = "None", 
                        options = list(`actions-box` = TRUE,size = 10,
                          `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_1D", label = "Comp1_D",  choices = "None", 
                        options = list(`actions-box` = TRUE,size = 10,
                                       `selected-text-format` = "count > 3"),  multiple = TRUE
                      )                
                     ),
               column(width=2, offset =0,
                      pickerInput(inputId = "var_2N", label = "Comp2_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_2D", label = "Comp2_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )
                     ),
               column(width=2, offset =0,
                      pickerInput(inputId = "var_3N", label = "Comp3_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_3D", label = "Comp3_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )
                  ),
               column(width=2, offset =0,
                      pickerInput(inputId = "var_4N", label = "Comp4_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_4D", label = "Comp4_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )
               ) ,            
               column(width=2, offset =0,
                      pickerInput(inputId = "var_5N", label = "Comp5_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_5D", label = "Comp5_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )
               ),            
               column(width=2, offset =0,
                      pickerInput(inputId = "var_6N", label = "Comp6_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_6D", label = "Comp6_D",  choices = "None", 
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
                      pickerInput(inputId = "var_7N", label = "Comp7_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_7D", label = "Comp7_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )                     
               ),
               column(width=2, offset =0,
                      pickerInput(inputId = "var_8N", label = "Comp8_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_8D", label = "Comp8_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )  
               ),
               column(width=2, offset =0,
                      pickerInput(inputId = "var_9N", label = "Comp9_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_9D", label = "Comp9_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )  
               ),
               column(width=2, offset =0,
                      pickerInput(inputId = "var_10N", label = "Comp10_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_10D", label = "Comp10_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )  
               ) ,            
               column(width=2, offset =0,
                      pickerInput(inputId = "var_11N", label = "Comp11_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_11D", label = "Comp11_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )  
               ),            
               column(width=2, offset =0,
                      pickerInput(inputId = "var_12N", label = "Comp12_N",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      ),
                      pickerInput(inputId = "var_12D", label = "Comp12_D",  choices = "None", 
                                  options = list(`actions-box` = TRUE,size = 10,
                                                 `selected-text-format` = "count > 3"),  multiple = TRUE
                      )  
               )
             ),
             fluidRow(align="center",
                      hr(),
                      #rHandsontableOutput("stat2_table"),
                      tags$head(tags$style("#stat2_N_1{color: blue; font-size: 16px; font-style: bold;}")),

                      actionButton("check_mva", label = "Check Comparisons", width = 300,
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                      br(),
                      br(),
                      actionButton("start_mva", label = "Start Anlaysis", width = 300,
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
             )
            ),
            
            tabPanel("Graphs", id="tp_mva_graph1", 
                     fluidRow(
                       column(width=8, offset =0,
                            pickerInput(inputId = "mva_plot_comp", label = "Comp #",  choices = "None", 
                                          options = list(`actions-box` = TRUE, size = 100,
                                                         `selected-text-format` = "count > 5"),  multiple = TRUE)
                            ),
                            # sliderInput("mva_plot_comp", label = "Comp #", min = 1, max = 12,
                            #        value = 1)
                            #  ),
                       column(width=4, offset =0,
                              actionButton("create_mva_plots", label = "Create Plots", width = 100,
                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                       ),
                            ),
                     fluidRow(
                       column(width=6, offset =0,
                          dropdownButton(
                                textInput("mva_barplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                textInput("mva_barplot_title", label="plot title", value = "Barplot", width = 200),
                                sliderInput("mva_barplot_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 11),
                                sliderInput("mva_barplot_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20),
                                 circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                 tooltip = tooltipOptions(title = "Click to see inputs !")
                                       ),
                            div(
                              style = "position:relative",
                              plotOutput("mva_barplot", width = 600, height = 400)
                               ),
                              downloadButton('download_mva_barplot')
                            ),  

                       column(width=6, offset =0,
                              dropdownButton(
                                textInput("mva_boxplot_y_axis_label", label="y axis label", value = "Intensity", width = 200),
                                textInput("mva_boxplot_title", label="plot title", value = "Boxplot", width = 200),
                                sliderInput("mva_boxplot_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 11),
                                sliderInput("mva_boxplot_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20),
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("mva_boxplot", width = 600, height = 400)
                              ),
                              downloadButton('download_mva_boxplot')
                            ),  
                          ),
                     hr(),
                     fluidRow(
                       column(width=6, offset =0,
                              dropdownButton(
                                selectInput("mva_pca2d_x", label = "pca xaxis", choices = list("PC1", "PC2", "PC3", "PC4", "PC5"), 
                                            selected = "PC1"),
                                selectInput("mva_pca2d_y", label = "pca yaxis", choices = list("PC1", "PC2", "PC3", "PC4", "PC5"), 
                                            selected = "PC2"),
                                textInput("mva_pca2d_title", label="plot title", value = "pca2d", width = 200),
                                sliderInput("mva_pca2d_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 11),
                                sliderInput("mva_pca2d_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20),
                                sliderInput("mva_pca2d_dot_size", label = h5("Point Size"), min = 1, 
                                            max = 20, value = 4),
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("mva_pca2d", width = 600, height = 400,
                                           hover = hoverOpts("plot_pca2d_hover", delay = 100, delayType = "debounce")),
                                uiOutput("hover_pca2d_info")
                              ),
                              downloadButton('download_mva_pca2d')
                       ),  
                       
                       column(width=6, offset =0,
                              dropdownButton(
                                textInput("mva_pca3d_title", label="plot title", value = "pca3d", width = 200),
                                sliderInput("mva_pca3d_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 11),
                                sliderInput("mva_pca3d_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20),
                                sliderInput("mva_pca3d_dot_size", label = h5("Point Size"), min = 1, 
                                            max = 10, value = 2),
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                rglwidgetOutput("mva_pca3d", width = 600, height = 400)
                              ),
                              downloadButton('download_mva_pca3d')
                       ),  
                     ),
                     hr(),
                     fluidRow(
                       column(width=6, offset =0,
                              dropdownButton(
                                textInput("mva_cluster_title", label="plot title", value = "cluster", width = 200),
                                sliderInput("mva_cluster_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 11),
                                sliderInput("mva_cluster_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20),
                                colourpicker::colourInput("cluster_high_color", "Select High Color", "#FF3366"),
                                colourpicker::colourInput("cluster_low_color", "Select Low Color", "#009933"),
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("mva_cluster", width = 600, height = 400)
                              ),
                              downloadButton('download_mva_cluster')
                       ),  
                       
                       column(width=6, offset =0,
                              dropdownButton(
                                textInput("mva_heatmap_title", label="plot title", value = "heatmap", width = 200),
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("mva_heatmap", width = 600, height = 400)
                              ),
                              downloadButton('download_mva_heatmap')
                       ),  
                     )  
                    ),   # end graphs
              
            tabPanel("Volcano", id="tp_mva_volcano", 
                     fluidRow(
                       column(width=1, offset =0,
                              actionButton("create_mva_volcano", label = "Create Volcano Plots", width = 300,
                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), 
                       )
                     ), 
                     fluidRow(
                       column(width=5, offset =0,
                              dropdownButton(
                                textInput("volcano1_mva_plot_title", label="plot title", 
                                          value = "Volcano", width = 200),
                                textInput("volcano1_mva_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                textInput("volcano1_mva_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                colourpicker::colourInput("volcano1_mva_dot_color", "Select Color", "blue"),
                                sliderInput("volcano1_mva_plot_dot_size", label = h5("Point Size"), min = 1, 
                                            max = 10, value = 2),
                                sliderInput("volcano1_mva_plot_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 20),
                                sliderInput("volcano1_mva_plot_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20), 
                                
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("volcano1_mva_plot", width = 600, height = 600,
                                           hover = hoverOpts("volcano1_mva_hover", delay = 100, delayType = "debounce")),
                                uiOutput("volcano1_mva_hover_info")
                              ),
                              downloadButton('download_mva_volcano1')
                       ), # end col
                       column(width=5, offset =0,
                              dropdownButton(
                                textInput("volcano2_mva_plot_title", label="plot title", 
                                          value = "Volcano", width = 200),
                                textInput("volcano2_mva_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                textInput("volcano2_mva_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                colourpicker::colourInput("volcano2_mva_dot_color", "Select Color", "blue"),
                                sliderInput("volcano2_mva_plot_dot_size", label = h5("Point Size"), min = 1, 
                                            max = 10, value = 2),
                                sliderInput("volcano2_mva_plot_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 20),
                                sliderInput("volcano2_mva_plot_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20), 
                                
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("volcano2_mva_plot", width = 600, height = 600,
                                           hover = hoverOpts("volcano2_mva_hover", delay = 100, delayType = "debounce")),
                                uiOutput("volcano2_mva_hover_info")
                              ),
                              downloadButton('download_mva_volcano2')
                       ) # end col
                       ), #end fr
                       
                     fluidRow(
                       column(width=5, offset =0,
                              dropdownButton(
                                textInput("volcano3_mva_plot_title", label="plot title", 
                                          value = "Volcano", width = 200),
                                textInput("volcano3_mva_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                textInput("volcano3_mva_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                colourpicker::colourInput("volcano3_mva_dot_color", "Select Color", "blue"),
                                sliderInput("volcano3_mva_plot_dot_size", label = h5("Point Size"), min = 1, 
                                            max = 10, value = 2),
                                sliderInput("volcano3_mva_plot_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 20),
                                sliderInput("volcano3_mva_plot_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20), 
                                
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("volcano3_mva_plot", width = 600, height = 600,
                                           hover = hoverOpts("volcano3_mva_hover", delay = 100, delayType = "debounce")),
                                uiOutput("volcano3_mva_hover_info")
                              ),
                              downloadButton('download_mva_volcano3')
                       ), # end col
                       column(width=5, offset =0,
                              dropdownButton(
                                textInput("volcano4_mva_plot_title", label="plot title", 
                                          value = "Volcano", width = 200),
                                textInput("volcano4_mva_plot_y_axis_label", label="y axis label", value = "-log_pvalue", width = 200),
                                textInput("volcano4_mva_plot_x_axis_label", label="x axis label", value = "log_FC", width = 200),
                                colourpicker::colourInput("volcano4_mva_dot_color", "Select Color", "blue"),
                                sliderInput("volcano4_mva_plot_dot_size", label = h5("Point Size"), min = 1, 
                                            max = 10, value = 2),
                                sliderInput("volcano4_mva_plot_label_size", label = h5("Label Size"), min = 1, 
                                            max = 50, value = 20),
                                sliderInput("volcano4_mva_plot_title_size", label = h5("Title Size"), min = 10, 
                                            max = 50, value = 20), 
                                
                                circle = TRUE, status = "danger", icon = icon("gear"), width = "300px", size = "sm",
                                tooltip = tooltipOptions(title = "Click to see inputs !")
                              ),
                              div(
                                style = "position:relative",
                                plotOutput("volcano4_mva_plot", width = 600, height = 600,
                                           hover = hoverOpts("volcano4_mva_hover", delay = 100, delayType = "debounce")),
                                uiOutput("volcano4_mva_hover_info")
                              ),
                              downloadButton('download_mva_volcano4')
                       ) # end col
                     ) #end fr    
            ), # end volcano
            
            tabPanel("Data", value = "tp_mva_data",
                     fluidRow( 
                       column(width=1, offset =0,  
                              numericInput("mva_data_topn", label="TopN", value = 0, width = 100)
                       ),
                       column(width=1, offset =0,
                              textInput("mva_data_accession", label="Accession", value = "0", width = 100)
                       ),
                       column(width=1, offset =0,
                              textInput("mva_data_description", label="Description", value = "0", width = 100)
                       ),
                       column(width=1, offset =0,  
                              numericInput("mva_data_pvalue", label="pvalue", value = 0, width = 100)
                       ),
                       column(width=1, offset =0,  
                              numericInput("mva_data_foldchange1", label="foldchange up", value = 0, width = 100)
                       ),
                       column(width=1, offset =0,  
                              numericInput("mva_data_foldchange2", label="foldchange dn", value = 0, width = 100)
                       ),
                       column(width=3, offset =0,
                              selectInput("mva_select_data_comp", label = "comparison", 
                                          choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                                          selected = 1)
                       ),
                       column(width=1, offset =0,
                              actionButton("mva_data_show", label = "Filter Data", width = 100,
                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                       )
                     ),
                     fluidRow(
                       hr(),
                       tags$head(tags$style("#mva_data_final{color: blue;
                                 font-size: 12px;
                                 }"
                       )
                       ),
                       rHandsontableOutput("mva_data_final")
                     )
            ) #end data
            
          ) # end of navbar panel
    ),  #end of tab panel
    
    
    
    
    tabPanel("Save", value = "tp_save", align="center",
               tags$h1("Save dpmsr_set file..."),
               hr(),
               br(),
               br(),
               actionButton("save_dpmsr_set", label = "Save dpmsr_set", width = 300,
                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
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
)) #end of ui
