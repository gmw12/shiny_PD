interactive_go_volcano <- function(session, input, output)
{
    volcano_data <- create_go_volcano(input, output, session)
    
    testdf <- data.frame(cbind(dpmsr_set$data$stats$final$Accession, dpmsr_set$data$stats$final$Description))
    colnames(testdf) <- c("Accession", "Description")
    volcano_data <- merge(x=volcano_data, y=testdf, by.x="Accession", by.y="Accession")
    
    
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
                      "<b> Description: </b>", point$Description, "<br/>",
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
  
  create_stats_barplot <- reactive({
    ggplot(data=df2, aes(x=Sample, y=Total_Intensity)) +
      geom_bar(stat="identity", fill=color_list)+ theme_classic() + 
      ggtitle(input$stats_barplot_title) + 
      ylab(input$stats_barplot_y_axis_label) +
      xlab(NULL) +
      #scale_y_discrete(labels = NULL) +
      coord_cartesian(ylim=NULL, expand = TRUE) +
      theme(plot.title = element_text(hjust = 0.5, size=input$stats_barplot_title_size), 
            axis.title = element_text(size=input$stats_barplot_label_size, color="black"),
            axis.text.x = element_text(size=input$stats_barplot_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$stats_barplot_label_size,  color="black"),
      ) 
  })
  
  output$stats_barplot <- renderPlot({
    req(create_stats_barplot())
    create_stats_barplot()
  })
  
  output$download_stats_barplot <- downloadHandler(
    filename = function(){
      str_c("stats_Barplot_", dpmsr_set$y$stats$groups$comp_name[as.numeric(input$stats_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_stats_barplot())
      ggsave(file, plot = create_stats_barplot(), device = 'png')
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
  
  create_stats_boxplot <- reactive({
    ggplot(data=df3, aes(x=Sample, y=Intensity)) +
      geom_boxplot(notch = TRUE, outlier.colour="red", outlier.shape=1,
                   outlier.size=1, fill=rev(color_list)) + theme_classic() + 
      coord_flip()+
      xlab(input$stats_boxplot_x_axis_label) +
      ggtitle(input$stats_boxplot_title) + 
      theme(plot.title = element_text(hjust = 0.5, size=input$stats_boxplot_title_size), 
            axis.title = element_text(size=input$stats_boxplot_label_size, color="black"),
            axis.text.x = element_text(size=input$stats_boxplot_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$stats_boxplot_label_size,  color="black"),
      ) 
  })
  
  output$stats_boxplot <- renderPlot({
    req(create_stats_boxplot())
    create_stats_boxplot()
  })
  
  output$download_stats_boxplot <- downloadHandler(
    filename = function(){
      str_c("stats_Boxplot_", dpmsr_set$y$stats$groups$comp_name[as.numeric(input$stats_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_stats_boxplot())
      ggsave(file, plot = create_stats_boxplot(), device = 'png')
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
  
  hover_data <- data.frame(cbind(namex, df_out[[input$stats_pca2d_x]], df_out[[input$stats_pca2d_y]]), stringsAsFactors = FALSE  )
  colnames(hover_data) <- c("Sample", "get(input$stats_pca2d_x)", "get(input$stats_pca2d_y)")
  hover_data$`get(input$stats_pca2d_x)` <- as.numeric(hover_data$`get(input$stats_pca2d_x)`)
  hover_data$`get(input$stats_pca2d_y)` <- as.numeric(hover_data$`get(input$stats_pca2d_y)`)
  
  hover_data_test <<- hover_data
  
  create_stats_pca2d <- reactive({
    ggplot(df_out, aes(x=get(input$stats_pca2d_x), y=get(input$stats_pca2d_y), color=x_gr )) +
      geom_point(alpha=0.8, size=input$stats_pca2d_dot_size) +
      theme(legend.title=element_blank()) +
      ggtitle(input$stats_pca2d_title) + 
      ylab(input$stats_pca2d_y) +
      xlab(input$stats_pca2d_x) +
      scale_color_manual(values = rev(unique(color_list))) +
      theme(plot.title = element_text(hjust = 0.5, size=input$stats_pca2d_title_size), 
            axis.title = element_text(size=input$stats_pca2d_label_size, color="black"),
            axis.text.x = element_text(size=input$stats_pca2d_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$stats_pca2d_label_size,  color="black"),
    ) 
  })
  
  output$stats_pca2d <- renderPlot({
    req(create_stats_pca2d())
    create_stats_pca2d()
  })
  
  output$download_stats_pca2d <- downloadHandler(
    filename = function(){
      str_c("stats_pca2d_", dpmsr_set$y$stats$groups$comp_name[as.numeric(input$stats_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_stats_pca2d())
      ggsave(file, plot = create_stats_pca2d(), device = 'png')
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

  create_stats_pca3d <- reactive({
    pca3d(x_pca, 
          group=x_gr,
          new=FALSE,
          legend = "right",
          palette = rev(unique(color_list)), 
          radius = input$stats_pca3d_dot_size,
          title = input$stats_pca3d_title)
  })
  
  output$stats_pca3d <- renderRglwidget ({
    try(rgl.close())
    req(create_stats_pca3d())
    create_stats_pca3d()
    rglwidget()
  })
  
  output$download_stats_pca3d <- downloadHandler(
    filename = function(){
      str_c("stats_pca3d_", dpmsr_set$y$stats$groups$comp_name[as.numeric(input$stats_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_stats_pca3d())
      snapshotPCA3d(file)
      #ggsave(file, plot = create_stats_pca3d(), device = 'png')
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
  
  
  create_stats_cluster <- reactive({
    distance <- get_dist(df, method="euclidean")
    fviz_dist(distance,  show_labels = TRUE, gradient = list(low = input$cluster_low_color, mid = "white", high = input$cluster_high_color)) +
      ggtitle(input$stats_cluster_title) +
      theme(plot.title = element_text(hjust = 0.5, size=input$stats_cluster_title_size), 
            axis.text.x = element_text(size=input$stats_cluster_label_size, angle = 90,  color="black"),
            axis.text.y = element_text(size=input$stats_cluster_label_size,  color="black"))
      
  })
  
  output$stats_cluster <- renderPlot({
    req(create_stats_cluster())
    create_stats_cluster()
  })
  
  output$download_stats_cluster <- downloadHandler(
    filename = function(){
      str_c("stats_cluster_", dpmsr_set$y$stats$groups$comp_name[as.numeric(input$stats_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_stats_cluster())
      ggsave(file, plot = create_stats_cluster(), device = 'png')
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
  
  create_stats_heatmap <- reactive({
    ## Row- and column-wise clustering 
    hr <- hclust(as.dist(1-cor(t(df), method="pearson")), method="complete")
    hc <- hclust(as.dist(1-cor(df, method="spearman")), method="complete") 
    ## Tree cutting
    mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
    ## Plot heatmap 
    mycol <- redgreen(75)
    #png(filename="erasemyheatmap.png", units="px", width = 1776, height = 1146)  
    heatmap.2(df, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, labCol=groupx, 
              scale="row", density.info="none", trace="none", RowSideColors=mycolhc, main = input$stats_heatmap_title) 
  })
  
  output$stats_heatmap <- renderPlot({
    req(create_stats_heatmap())
    create_stats_heatmap()
  })
  
  output$download_stats_heatmap <- downloadHandler(
    filename = function(){
      str_c("stats_heatmap_", dpmsr_set$y$stats$groups$comp_name[as.numeric(input$stats_plot_comp)],
            ".png", collapse = " ")
    },
    content = function(file){
      req(create_stats_heatmap())
      png(filename=file, units="px", width = 1776, height = 1146)  
      create_stats_heatmap()
      dev.off()
    }
  )
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_stats_volcano <- function(session, input, output, i)
{
  df <- dpmsr_set$data$stats$final
  df <- subset(df, df[ , dpmsr_set$y$stats$groups$mf[i]] >= input$missing_factor )
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[i]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[i]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  df <- cbind(df$Accession, df$Description, df_fc, df_pval)
  colnames(df) <- c("Accession", "Description", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$Description <- as.character(df$Description)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano_stats_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input[[str_c("volcano",i,"_stats_plot_dot_size")]], color = input[[str_c("volcano",i,"_stats_dot_color")]] ) +
      xlab(input[[str_c("volcano",i,"_stats_plot_x_axis_label")]]) + 
      ylab(input[[str_c("volcano",i,"_stats_plot_y_axis_label")]]) +
      scale_colour_gradient(low = input[[str_c("volcano", i, "_stats_dot_color")]], high = input[[str_c("volcano", i, "_stats_dot_color")]] ) +
      ggtitle(input[[str_c("volcano",i,"_stats_plot_title")]])+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input[[str_c("volcano",i,"_stats_plot_title_size")]], hjust = 0.5),
            axis.title = element_text(size=input[[str_c("volcano",i,"_stats_plot_label_size")]], color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  plot_name <- str_c("volcano", i, "_stats_plot")
  download_name <- str_c("download_stats_volcano", i)
  hover_name <- str_c("volcano", i, "_stats_hover_info")
  hover_stats_name <- str_c("volcano", i, "_stats_hover")
  
  output[[plot_name]]<- renderPlot({
    req(volcano_stats_plot())
    volcano_stats_plot()
  })
  
  output[[download_name]] <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$stats$groups$comp_name[i], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano_stats_plot())
      ggsave(file, plot = volcano_stats_plot(), device = 'png')
    }
  )
  
  output[[hover_name]] <- renderUI({
    hover <- input[[hover_stats_name]]
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    #cat(file=stderr(), str_c("top_pct = ", top_pct), "\n")
    #cat(file=stderr(), str_c("top_px = ", top_px), "\n")
    
    if(top_pct > 0.3){
      top_custom <- 10
    }else{
      top_custom <- 200
    }
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", top_custom, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> Description: </b>", point$Description, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}



#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_stats_volcano1 <- function(session, input, output)
{
  df <- dpmsr_set$data$stats$final
  df <- subset(df, df[ , dpmsr_set$y$stats$groups$mf[1]] >= input$missing_factor )
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  df <- cbind(df$Accession, df$Description, df_fc, df_pval)
  colnames(df) <- c("Accession", "Description", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$Description <- as.character(df$Description)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano1_stats_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano1_stats_plot_dot_size, color = input$volcano1_stats_dot_color) +
      xlab(input$volcano1_statsplot_x_axis_label) + ylab(input$volcano1_stats_plot_y_axis_label) +
      scale_colour_gradient(low = input$volcano1_stats_dot_color, high = input$volcano1_stats_dot_color) +
      ggtitle(input$volcano1_stats_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano1_stats_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano1_stats_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano1_stats_plot<- renderPlot({
    req(volcano1_stats_plot())
    volcano1_stats_plot()
  })
  
  output$download_stats_volcano1 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$stats$groups$comp_name[1], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano1_stats_plot())
      ggsave(file, plot = volcano1_stats_plot(), device = 'png')
    }
  )
  
  output$volcano1_stats_hover_info <- renderUI({
    hover <- input$volcano1_stats_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)

    cat(file=stderr(), str_c("top_pct = ", top_pct), "\n")
    cat(file=stderr(), str_c("top_px = ", top_px), "\n")
    
    if(top_pct > 0.3){
      top_custom <- 10
    }else{
      top_custom <- 200
    }
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", top_custom, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> Description: </b>", point$Description, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_stats_volcano2 <- function(session, input, output)
{
  df <- dpmsr_set$data$stats$final
  df <- subset(df, df[ , dpmsr_set$y$stats$groups$mf[2]] >= input$missing_factor )
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[2]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[2]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  df <- cbind(df$Accession, df$Description, df_fc, df_pval)
  colnames(df) <- c("Accession", "Description", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$Description <- as.character(df$Description)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano2_stats_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano2_stats_plot_dot_size, color = input$volcano2_stats_dot_color) +
      xlab(input$volcano2_statsplot_x_axis_label) + ylab(input$volcano2_stats_plot_y_axis_label) +
      scale_colour_gradient(low = input$volcano2_stats_dot_color, high = input$volcano2_stats_dot_color) +
      ggtitle(input$volcano2_stats_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano2_stats_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano2_stats_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano2_stats_plot<- renderPlot({
    req(volcano2_stats_plot())
    volcano2_stats_plot()
  })
  
  output$download_stats_volcano2 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$stats$groups$comp_name[2], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano2_stats_plot())
      ggsave(file, plot = volcano2_stats_plot(), device = 'png')
    }
  )
  
  output$volcano2_stats_hover_info <- renderUI({
    hover <- input$volcano2_stats_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    if(top_pct > 0.3){
      top_custom <- 10
    }else{
      top_custom <- 200
    }
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", top_custom, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> Description: </b>", point$Description, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_stats_volcano3 <- function(session, input, output)
{
  df <- dpmsr_set$data$stats$final
  df <- subset(df, df[ , dpmsr_set$y$stats$groups$mf[3]] >= input$missing_factor )
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[3]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[3]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  df <- cbind(df$Accession, df$Description, df_fc, df_pval)
  colnames(df) <- c("Accession", "Description", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$Description <- as.character(df$Description)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano3_stats_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano3_stats_plot_dot_size, color = input$volcano3_stats_dot_color) +
      xlab(input$volcano3_statsplot_x_axis_label) + ylab(input$volcano3_stats_plot_y_axis_label) +
      scale_colour_gradient(low = input$volcano3_stats_dot_color, high = input$volcano3_stats_dot_color) +
      ggtitle(input$volcano3_stats_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano3_stats_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano3_stats_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano3_stats_plot<- renderPlot({
    req(volcano3_stats_plot())
    volcano3_stats_plot()
  })
  
  output$download_stats_volcano3 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$stats$groups$comp_name[3], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano3_stats_plot())
      ggsave(file, plot = volcano3_stats_plot(), device = 'png')
    }
  )
  
  output$volcano3_stats_hover_info <- renderUI({
    hover <- input$volcano3_stats_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    if(top_pct > 0.3){
      top_custom <- 10
    }else{
      top_custom <- 200
    }
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", top_custom, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> Description: </b>", point$Description, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_stats_volcano4 <- function(session, input, output)
{
  df <- dpmsr_set$data$stats$final
  df <- subset(df, df[ , dpmsr_set$y$stats$groups$mf[4]] >= input$missing_factor )
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[4]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[4]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  df <- cbind(df$Accession, df$Description, df_fc, df_pval)
  colnames(df) <- c("Accession", "Description", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$Description <- as.character(df$Description)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano4_stats_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano4_stats_plot_dot_size, color = input$volcano4_stats_dot_color) +
      xlab(input$volcano4_statsplot_x_axis_label) + ylab(input$volcano4_statsp_lot_y_axis_label) +
      scale_colour_gradient(low = input$volcano4_stats_dot_color, high = input$volcano4_stats_dot_color) +
      ggtitle(input$volcano4_stats_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano4_stats_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano4_stats_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano4_stats_plot<- renderPlot({
    req(volcano4_stats_plot())
    volcano4_stats_plot()
  })
  
  output$download_stats_volcano4 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$stats$groups$comp_name[4], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano4_stats_plot())
      ggsave(file, plot = volcano4_stats_plot(), device = 'png')
    }
  )
  
  output$volcano4_stats_hover_info <- renderUI({
    hover <- input$volcano4_stats_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    if(top_pct > 0.3){
      top_custom <- 10
    }else{
      top_custom <- 200
    }
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", top_custom, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> Description: </b>", point$Description, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_stats_volcano5 <- function(session, input, output)
{
  df <- dpmsr_set$data$stats$final
  df <- subset(df, df[ , dpmsr_set$y$stats$groups$mf[5]] >= input$missing_factor )
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[5]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[5]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  df <- cbind(df$Accession, df$Description, df_fc, df_pval)
  colnames(df) <- c("Accession", "Description", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$Description <- as.character(df$Description)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano5_stats_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano5_stats_plot_dot_size, color = input$volcano5_stats_dot_color) +
      xlab(input$volcano5_statsplot_x_axis_label) + ylab(input$volcano5_statsp_lot_y_axis_label) +
      scale_colour_gradient(low = input$volcano5_stats_dot_color, high = input$volcano5_stats_dot_color) +
      ggtitle(input$volcano5_stats_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano5_stats_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano5_stats_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano5_stats_plot<- renderPlot({
    req(volcano5_stats_plot())
    volcano5_stats_plot()
  })
  
  output$download_stats_volcano5 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$stats$groups$comp_name[5], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano5_stats_plot())
      ggsave(file, plot = volcano5_stats_plot(), device = 'png')
    }
  )
  
  output$volcano5_stats_hover_info <- renderUI({
    hover <- input$volcano5_stats_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    if(top_pct > 0.3){
      top_custom <- 10
    }else{
      top_custom <- 200
    }
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", top_custom, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> Description: </b>", point$Description, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

interactive_stats_volcano6 <- function(session, input, output)
{
  df <- dpmsr_set$data$stats$final
  df <- subset(df, df[ , dpmsr_set$y$stats$groups$mf[6]] >= input$missing_factor )
  df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[6]))
  df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[6]))
  
  #df_fc <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$fc[1]))
  #df_pval <- df %>% dplyr::select(contains(dpmsr_set$y$stats$groups$pval[1]))
  
  df <- cbind(df$Accession, df$Description, df_fc, df_pval)
  colnames(df) <- c("Accession", "Description", "fc", "fc2", "pval")
  df$Accession <- as.character(df$Accession)
  df$Description <- as.character(df$Description)
  df$log_pvalue <- -log(as.numeric(df$pval), 10)
  df$log_fc <- log(as.numeric(df$fc2), 2)
  
  volcano6_stats_plot <- reactive({
    ggplot(df, aes(x = log_fc, y = log_pvalue)) +
      theme_minimal() +
      geom_point(alpha=0.4, size=input$volcano6_stats_plot_dot_size, color = input$volcano6_stats_dot_color) +
      xlab(input$volcano6_statsplot_x_axis_label) + ylab(input$volcano6_statsp_lot_y_axis_label) +
      scale_colour_gradient(low = input$volcano6_stats_dot_color, high = input$volcano6_stats_dot_color) +
      ggtitle(input$volcano6_stats_plot_title)+    
      xlim(-max(df$log_fc), max(df$log_fc)) +
      theme(plot.title = element_text(size=input$volcano6_stats_plot_title_size, hjust = 0.5),
            axis.title = element_text(size=input$volcano6_stats_plot_label_size, color="black"),
            axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10,  color="black"),
            legend.position = "none")
  })
  
  output$volcano6_stats_plot<- renderPlot({
    req(volcano6_stats_plot())
    volcano6_stats_plot()
  })
  
  output$download_stats_volcano6 <- downloadHandler(
    filename = function(){
      str_c("Volcano_", dpmsr_set$y$stats$groups$comp_name[6], ".png", collapse = " ")
    },
    content = function(file){
      req(volcano6_stats_plot())
      ggsave(file, plot = volcano6_stats_plot(), device = 'png')
    }
  )
  
  output$volcano6_stats_hover_info <- renderUI({
    hover <- input$volcano6_stats_hover
    point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    left_px <- left_pct * (hover$range$right - hover$range$left)
    top_px <- top_pct * (hover$range$bottom - hover$range$top)
    
    if(top_pct > 0.3){
      top_custom <- 10
    }else{
      top_custom <- 200
    }
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", 10, "px; top:", top_custom, "px;")
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Accession: </b>", point$Accession, "<br/>",
                    "<b> Description: </b>", point$Description, "<br/>",
                    "<b> FC: </b>", point$fc, "<br/>",
                    "<b> pvalue: </b>", point$pval, "<br/>")))
    )
  })
}