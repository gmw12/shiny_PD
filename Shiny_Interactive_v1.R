interactive_go_volcano <- function(session, input, output)
{
    volcano_data <- create_go_volcano(input, output, session)
    
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
      # left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
      # top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
      left_px <- left_pct * (hover$range$right - hover$range$left)
      top_px <- top_pct * (hover$range$bottom - hover$range$top)
      
      #cat(file=stderr(), str_c("hoverrr=", hover$range$right, "   hoverrl=", hover$range$left, "  left_pct=", left_pct), "\n")
      #cat(file=stderr(), str_c("hoverrb=", hover$range$bottom, "   hoverrt=", hover$range$top, "  top_pct=", top_pct), "\n")
      #cat(file=stderr(), str_c("leftpx=", left_px, "   toppx=", top_px), "\n", "\n"  )
      
      # create style property fot tooltip
      # background color is set so tooltip is a bit transparent
      # z-index is set so we are sure are tooltip will be on top
      style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                      "left:", 10, "px; top:", 10, "px;")
      # actual tooltip created as wellPanel
      wellPanel(
        style = style,
        p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                      "<b> FC: </b>", point$foldchange, "<br/>",
                      "<b> pvalue: </b>", point$pvalue, "<br/>")))
      )
    })
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#match(namesdf, names(df))

interactive_barplot <- function(session, input, output, df, namex, color_list)
{
  datay <- colSums(df, na.rm = TRUE)
  df2 <- data.frame(namex)
  df2$Total_Intensity <- datay
  colnames(df2) <- c("Sample", "Total_Intensity")
  df2$Sample <- factor(df2$Sample, levels = df2$Sample)
  ymax <- max(datay)
  
  create_mva_barplot <- reactive({
    ggplot(data=df2, aes(x=Sample, y=Total_Intensity)) +
      geom_bar(stat="identity", fill=color_list)+ theme_classic() + 
      ggtitle(input$mva_barplot_title) + 
      ylab(input$mva_barplot_y_axis_label) +
      xlab(NULL) +
      #scale_y_discrete(labels = NULL) +
      coord_cartesian(ylim=NULL, expand = TRUE) +
      theme(plot.title = element_text(hjust = 0.5, size=input$mva_barplot_title_size), 
            axis.title = element_text(size=input$mva_barplot_label_size, color="black"),
            axis.text.x = element_text(size=input$mva_barplot_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$mva_barplot_label_size,  color="black"),
      ) 
  })
  
  output$mva_barplot <- renderPlot({
    req(create_mva_barplot())
    create_mva_barplot()
  })
  
  output$download_mva_barplot <- downloadHandler(
    filename = function(){
      str_c("MVA_Barplot_", dpmsr_set$y$mva$groups$comp_name[as.numeric(input$mva_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_mva_barplot())
      ggsave(file, plot = create_mva_barplot(), device = 'png')
    }
  )
  
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_boxplot <- function(session, input, output, df, namex, color_list)
{
  colnames(df) <- namex
  df3 <- log2(df) %>% gather(Sample, Intensity, colnames(df))
  df3$Sample <- factor(df3$Sample, levels = rev(namex))
  
  create_mva_boxplot <- reactive({
    ggplot(data=df3, aes(x=Sample, y=Intensity)) +
      geom_boxplot(notch = TRUE, outlier.colour="red", outlier.shape=1,
                   outlier.size=1, fill=rev(color_list)) + theme_classic() + 
      coord_flip()+
      xlab(input$mva_boxplot_x_axis_label) +
      ggtitle(input$mva_boxplot_title) + 
      theme(plot.title = element_text(hjust = 0.5, size=input$mva_boxplot_title_size), 
            axis.title = element_text(size=input$mva_boxplot_label_size, color="black"),
            axis.text.x = element_text(size=input$mva_boxplot_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$mva_boxplot_label_size,  color="black"),
      ) 
  })
  
  output$mva_boxplot <- renderPlot({
    req(create_mva_boxplot())
    create_mva_boxplot()
  })
  
  output$download_mva_boxplot <- downloadHandler(
    filename = function(){
      str_c("MVA_Boxplot_", dpmsr_set$y$mva$groups$comp_name[as.numeric(input$mva_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_mva_boxplot())
      ggsave(file, plot = create_mva_boxplot(), device = 'png')
    }
  )
  
}
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_pca2d <- function(session, input, output, df, namex, color_list, groupx)
{
  x_transpose <- t(df)
  x_transpose <-data.frame(x_transpose)
  row.names(x_transpose) <- NULL
  x_transpose <-cbind(groupx, x_transpose)
  x_pca <- prcomp(x_transpose[,-1], scale=TRUE)
  test_this <-x_transpose[,1]
  x_gr <- factor(unlist(test_this))
  summary(x_gr)
  df_out <- as.data.frame(x_pca$x)
  df_out_test <<- df_out
  df_xgr <- data.frame(x_gr)
  df_xgr_test <<- df_xgr
  #df_xgr$x_gr <- as.character(df_xgr$x_gr)
  
  hover_data <- data.frame(cbind(namex, df_out[[input$mva_pca2d_x]], df_out[[input$mva_pca2d_y]]), stringsAsFactors = FALSE  )
  colnames(hover_data) <- c("Sample", "get(input$mva_pca2d_x)", "get(input$mva_pca2d_y)")
  hover_data$`get(input$mva_pca2d_x)` <- as.numeric(hover_data$`get(input$mva_pca2d_x)`)
  hover_data$`get(input$mva_pca2d_y)` <- as.numeric(hover_data$`get(input$mva_pca2d_y)`)
  
  hover_data_test <<- hover_data
  
  create_mva_pca2d <- reactive({
    ggplot(df_out, aes(x=get(input$mva_pca2d_x), y=get(input$mva_pca2d_y), color=x_gr )) +
      geom_point(alpha=0.8, size=input$mva_pca2d_dot_size) +
      theme(legend.title=element_blank()) +
      ggtitle(input$mva_pca2d_title) + 
      ylab(input$mva_pca2d_y) +
      xlab(input$mva_pca2d_x) +
      scale_color_manual(values = rev(unique(color_list))) +
      theme(plot.title = element_text(hjust = 0.5, size=input$mva_pca2d_title_size), 
            axis.title = element_text(size=input$mva_pca2d_label_size, color="black"),
            axis.text.x = element_text(size=input$mva_pca2d_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$mva_pca2d_label_size,  color="black"),
    ) 
  })
  
  output$mva_pca2d <- renderPlot({
    req(create_mva_pca2d())
    create_mva_pca2d()
  })
  
  output$download_mva_pca2d <- downloadHandler(
    filename = function(){
      str_c("MVA_pca2d_", dpmsr_set$y$mva$groups$comp_name[as.numeric(input$mva_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_mva_pca2d())
      ggsave(file, plot = create_mva_pca2d(), device = 'png')
    }
  )
  
  output$hover_pca2d_info <- renderUI({
    hover <- input$plot_pca2d_hover
    point <- nearPoints(hover_data, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)

    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", 10, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Sample: </b>", point$Sample, "<br/>")))
    )
  })
  
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_pca3d <- function(session, input, output, df, namex, color_list, groupx)
{
  x_transpose <- t(df)
  x_transpose <-data.frame(x_transpose)
  row.names(x_transpose) <- NULL
  x_transpose <-cbind(groupx, x_transpose)
  x_pca <- prcomp(x_transpose[,-1], scale=TRUE)
  test_this <-x_transpose[,1]
  x_gr <- factor(unlist(test_this))
  summary(x_gr)

  create_mva_pca3d <- reactive({
    pca3d(x_pca, 
          group=x_gr,
          new=FALSE,
          legend = "right",
          palette = rev(unique(color_list)), 
          radius = input$mva_pca3d_dot_size,
          title = input$mva_pca3d_title)
  })
  
  output$mva_pca3d <- renderRglwidget ({
    try(rgl.close())
    req(create_mva_pca3d())
    create_mva_pca3d()
    rglwidget()
  })
  
  output$download_mva_pca3d <- downloadHandler(
    filename = function(){
      str_c("MVA_pca3d_", dpmsr_set$y$mva$groups$comp_name[as.numeric(input$mva_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_mva_pca3d())
      snapshotPCA3d(file)
      #ggsave(file, plot = create_mva_pca3d(), device = 'png')
    }
  )
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_cluster <- function(session, input, output, df, namex)
{
  colnames(df) <- namex
  
  df <- t(df)
  df <-data.frame(df)
  #row.names(df) <- NULL
  #df <- na.omit(df)
  df[] <- lapply(df, as.numeric)
  df <- scale(df)
  
  
  create_mva_cluster <- reactive({
    distance <- get_dist(df, method="euclidean")
    fviz_dist(distance,  show_labels = TRUE, gradient = list(low = input$cluster_low_color, mid = "white", high = input$cluster_high_color)) +
      ggtitle(input$mva_cluster_title) +
      theme(plot.title = element_text(hjust = 0.5, size=input$mva_cluster_title_size), 
            axis.text.x = element_text(size=input$mva_cluster_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$mva_cluster_label_size,  color="black"))
      
  })
  
  output$mva_cluster <- renderPlot({
    req(create_mva_cluster())
    create_mva_cluster()
  })
  
  output$download_mva_cluster <- downloadHandler(
    filename = function(){
      str_c("MVA_cluster_", dpmsr_set$y$mva$groups$comp_name[as.numeric(input$mva_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_mva_cluster())
      ggsave(file, plot = create_mva_cluster(), device = 'png')
    }
  )
}


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_heatmap <- function(session, input, output, df, namex, groupx)
{
  colnames(df) <- namex
  
  df <- log2(df)
  df <- data.matrix(df)
  
  create_mva_heatmap <- reactive({
    ## Row- and column-wise clustering 
    hr <- hclust(as.dist(1-cor(t(df), method="pearson")), method="complete")
    hc <- hclust(as.dist(1-cor(df, method="spearman")), method="complete") 
    ## Tree cutting
    mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
    ## Plot heatmap 
    mycol <- redgreen(75)
    #png(filename="erasemyheatmap.png", units="px", width = 1776, height = 1146)  
    heatmap.2(df, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, labCol=groupx, 
              scale="row", density.info="none", trace="none", RowSideColors=mycolhc, main = input$mva_heatmap_title) 
  })
  
  output$mva_heatmap <- renderPlot({
    req(create_mva_heatmap())
    create_mva_heatmap()
  })
  
  output$download_mva_heatmap <- downloadHandler(
    filename = function(){
      str_c("MVA_heatmap_", dpmsr_set$y$mva$groups$comp_name[as.numeric(input$mva_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_mva_heatmap())
      png(filename=file, units="px", width = 1776, height = 1146)  
      create_mva_heatmap()
      dev.off()
    }
  )
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_mva_volcano1 <- function(session, input, output)
{
  df <- dpmsr_set$data$mva$final
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[1]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[1]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[1]))
  
  df <- cbind(df$Accession, df_fc, df_pval)
  colnames(df) <- c("Accession", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano1_mva_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano1_mva_plot_dot_size, color = input$volcano1_mva_dot_color) +
      xlab(input$volcano1_mvaplot_x_axis_label) + ylab(input$volcano1_mva_plot_y_axis_label) +
      scale_colour_gradient(low = input$volcano1_mva_dot_color, high = input$volcano1_mva_dot_color) +
      ggtitle(input$volcano1_mva_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano1_mva_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano1_mva_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano1_mva_plot<- renderPlot({
    req(volcano1_mva_plot())
    volcano1_mva_plot()
  })
  
  output$download_mva_volcano1 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$mva$groups$comp_name[1], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano1_mva_plot())
      ggsave(file, plot = volcano1_mva_plot(), device = 'png')
    }
  )
  
  output$volcano1_mva_hover_info <- renderUI({
    hover <- input$volcano1_mva_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)

    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", 10, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_mva_volcano2 <- function(session, input, output)
{
  df <- dpmsr_set$data$mva$final
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[2]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[2]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[1]))
  
  df <- cbind(df$Accession, df_fc, df_pval)
  colnames(df) <- c("Accession", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano2_mva_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano2_mva_plot_dot_size, color = input$volcano2_mva_dot_color) +
      xlab(input$volcano2_mvaplot_x_axis_label) + ylab(input$volcano2_mva_plot_y_axis_label) +
      scale_colour_gradient(low = input$volcano2_mva_dot_color, high = input$volcano2_mva_dot_color) +
      ggtitle(input$volcano2_mva_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano2_mva_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano2_mva_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano2_mva_plot<- renderPlot({
    req(volcano2_mva_plot())
    volcano2_mva_plot()
  })
  
  output$download_mva_volcano2 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$mva$groups$comp_name[2], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano2_mva_plot())
      ggsave(file, plot = volcano2_mva_plot(), device = 'png')
    }
  )
  
  output$volcano2_mva_hover_info <- renderUI({
    hover <- input$volcano2_mva_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", 10, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_mva_volcano3 <- function(session, input, output)
{
  df <- dpmsr_set$data$mva$final
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[3]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[3]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[1]))
  
  df <- cbind(df$Accession, df_fc, df_pval)
  colnames(df) <- c("Accession", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano3_mva_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano3_mva_plot_dot_size, color = input$volcano3_mva_dot_color) +
      xlab(input$volcano3_mvaplot_x_axis_label) + ylab(input$volcano3_mva_plot_y_axis_label) +
      scale_colour_gradient(low = input$volcano3_mva_dot_color, high = input$volcano3_mva_dot_color) +
      ggtitle(input$volcano3_mva_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano3_mva_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano3_mva_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano3_mva_plot<- renderPlot({
    req(volcano3_mva_plot())
    volcano3_mva_plot()
  })
  
  output$download_mva_volcano3 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$mva$groups$comp_name[3], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano3_mva_plot())
      ggsave(file, plot = volcano3_mva_plot(), device = 'png')
    }
  )
  
  output$volcano3_mva_hover_info <- renderUI({
    hover <- input$volcano3_mva_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", 10, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_mva_volcano4 <- function(session, input, output)
{
  df <- dpmsr_set$data$mva$final
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[4]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[4]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$mva$groups$pval[1]))
  
  df <- cbind(df$Accession, df_fc, df_pval)
  colnames(df) <- c("Accession", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano4_mva_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano4_mva_plot_dot_size, color = input$volcano4_mva_dot_color) +
      xlab(input$volcano4_mvaplot_x_axis_label) + ylab(input$volcano4_mvap_lot_y_axis_label) +
      scale_colour_gradient(low = input$volcano4_mva_dot_color, high = input$volcano4_mva_dot_color) +
      ggtitle(input$volcano4_mva_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano4_mva_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano4_mva_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano4_mva_plot<- renderPlot({
    req(volcano4_mva_plot())
    volcano4_mva_plot()
  })
  
  output$download_mva_volcano4 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$mva$groups$comp_name[4], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano4_mva_plot())
      ggsave(file, plot = volcano4_mva_plot(), device = 'png')
    }
  )
  
  output$volcano4_mva_hover_info <- renderUI({
    hover <- input$volcano4_mva_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", 10, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}

