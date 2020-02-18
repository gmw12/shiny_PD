library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinyWidgets)
library(rhandsontable)

ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor("DarkGray", shinydashboard = TRUE),
  setBackgroundColor("LightGray", shinydashboard = FALSE),
  titlePanel("Proteome Discoverer Data Processing"),

  
  navlistPanel(widths=c(1,11), id = "nlp1",
          
    tabPanel("Load Design", value = "tp1", align="center",
                        tags$h1("Choose and Load the study design file..."),
                        hr(),
                        shinyFilesButton('design_file', label='Choose Design File', title='Please select excel design file', multiple=FALSE,
                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        br(),                            
                        hr(),
                        actionButton("action_load_design", label = "Load Design File",
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
             column(width=3, offset =1,
                    fluidRow(
                      rHandsontableOutput("project_overview"))
             ),
             
            column(width=6, offset =0,
                    imageOutput("mass_accuracy"),
                    br(),
                    imageOutput("inj_summary")
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
                     checkboxInput("checkbox_n1", label = "Sample Loading - Total"),
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
                                                  "KNN"= 5, "LocalLeastSquares" = 6, "MLE" = 7, "Bottom5" = 8
                                   ),
                                   selected = 1),
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
                                                 choices = list("barplot", "boxplot", "cluster", "pca2d", "pca3d",
                                                                "MDS", "density", "heatmap3"), 
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
    
    tabPanel("Extra", value = "tp11",
      navbarPage("Extra:", id ="np_extra",
        tabPanel("save_load", value = "tp_save", align="center",
                 tags$h1("Save/Load dpmsr_set files..."),
                 hr(),
                 br(),
                 br(),
                 actionButton("save_dpmsr_set", label = "Save dpmsr_set", width = 300,
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 hr(),
                 shinyFilesButton('dpmsr_set_file', label='Choose Design File', title='Choose dpmsr_set File', multiple=FALSE,
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 
                 # fileInput("dpmsr_set_file2", "Choose dpmsr_set File", multiple = FALSE, accept = "dpmsr",
                 #           width = NULL, buttonLabel = "Browse...",
                 #           placeholder = "No file selected"),                               
                 br(),
                 br(),
                 br(),
                 actionButton("load_dpmsr_set", label = "Load dpmsr_set", width = 300, 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        ),
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
                    column(width=2, offset =0,
                           shinyFilesButton('motif_fasta', label='Select Motif-X FASTA', title='Please select motif-x formated text file', multiple=FALSE,
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                           textOutput("fasta"),
                           tags$head(tags$style("#fasta{color: blue; font-size: 16px; font-style: bold;}")),
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
    )
    
    
 ) #end of navlistpanel
) #end of ui