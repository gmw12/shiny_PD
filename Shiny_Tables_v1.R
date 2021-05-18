   
protein_table <- function(session, input, output, filter_df){

  pval_cols <- colnames(filter_df %>% dplyr::select(contains("pval") ) )
  pval_col_numbers <- list(match(pval_cols, names(filter_df)))
  pval_col_numbers <- unlist(pval_col_numbers)
  sample_cols <- c(colnames(filter_df %>% dplyr::select(contains("Normalized"))),
                   colnames(filter_df %>% dplyr::select(contains("Imputed"))) )
  sample_col_numbers <- list(match(sample_cols, names(filter_df)))
  sample_col_numbers <- unlist(sample_col_numbers)
  cv_cols <- colnames(filter_df %>% dplyr::select(contains("CV") ) )
  mf_cols <- colnames(filter_df %>% dplyr::select(contains("MF") ) )
  stat_col <- ncol(filter_df) -1
  
 stats_DT <-  DT::datatable(filter_df,
                               rownames = FALSE,
                               extensions = c("FixedColumns"), #, "Buttons"),
                               #editable = list(target='cell', disable = list(columns=c(0:(stat_col-1))) ) ,
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
    
    stats_DT <- stats_DT %>%  formatRound(columns=c(sample_col_numbers), digits=0)
    #stats_DT <- stats_DT %>%  formatRound(columns=c(pval_col_numbers), digits=2)
    
 return(stats_DT)   
    
}




#--------------------------------




peptide_table <- function(session, input, output, filter_df){
  
  pval_cols <- colnames(filter_df %>% dplyr::select(contains("pval") ) )
  pval_col_numbers <- list(match(pval_cols, names(filter_df)))
  pval_col_numbers <- unlist(pval_col_numbers)
  sample_cols <- c(colnames(filter_df %>% dplyr::select(contains("Normalized"))),
                   colnames(filter_df %>% dplyr::select(contains("Imputed"))) )
  sample_col_numbers <- list(match(sample_cols, names(filter_df)))
  sample_col_numbers <- unlist(sample_col_numbers)
  cv_cols <- colnames(filter_df %>% dplyr::select(contains("CV") ) )
  mf_cols <- colnames(filter_df %>% dplyr::select(contains("MF") ) )
  stat_col <- ncol(filter_df) -1
  
  
  stats_DT <-  DT::datatable(filter_df,
                             rownames = FALSE,
                             extensions = c("FixedColumns"), #, "Buttons"),
                             #editable = list(target='cell', disable = list(columns=c(0:(stat_col-1))) ) ,
                             options=list(
                               #dom = 'Bfrtipl',
                               autoWidth = TRUE,
                               scrollX = TRUE,
                               scrollY=500,
                               scrollCollapse=TRUE,
                               columnDefs = list(list(targets = c(0), visibile = TRUE, "width"='30', className = 'dt-center'),
                                                 list(targets = c(1), visible = TRUE, "width"='20', className = 'dt-center'),
                                                 list(
                                                   targets = c(2),
                                                   width = '250',
                                                   render = JS(
                                                     "function(data, type, row, meta) {",
                                                     "return type === 'display' && data.length > 35 ?",
                                                     "'<span title=\"' + data + '\">' + data.substr(0, 35) + '...</span>' : data;",
                                                     "}")
                                                 ),
                                                 list(
                                                   targets = c(4),
                                                   width = '200',
                                                   render = JS(
                                                     "function(data, type, row, meta) {",
                                                     "return type === 'display' && data.length > 20 ?",
                                                     "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                     "}")
                                                 ),
                                                 list(
                                                   targets = c(8),
                                                   width = '150',
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
  
  stats_DT <- stats_DT %>%  formatRound(columns=c(sample_col_numbers), digits=0)
  #stats_DT <- stats_DT %>%  formatRound(columns=c(pval_col_numbers), digits=2)
  
  return(stats_DT)   
  
}
