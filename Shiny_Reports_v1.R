report_templates <- function(session, input){
  
 digest_strap_mini <- readLines("Shiny_Report_STrap_Mini.txt", warn = FALSE)
 
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_01", as.character(dpmsr_set$y$sample_number-1))
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_02", as.character(max(dpmsr_set$design$Replicate)))
 samples_list <- str_c(dpmsr_set$y$sample_groups$Group, collapse = ", ")
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_03", samples_list)
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_04", as.character(unique(dpmsr_set$design$`QC Spike Level`)[1]))
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_05", as.character(unique(dpmsr_set$design$`QC Spike Level`)[2]))
 volume <- trunc(as.numeric(dpmsr_set$x$strap_binding_buffer),0)
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_06", as.character(volume))
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_07", as.character(dpmsr_set$x$final_volume))
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_08", as.character(dpmsr_set$x$ADH_recon_conc  ))
 digest_strap_mini <- str_replace_all(digest_strap_mini, "input_09", as.character(dpmsr_set$x$SPQC_aliquot ))
 }
  
dpmsr_set$y$sample_groups
sample_groups <- dpmsr_set$y$sample_groups[dpmsr_set$y$sample_groups$Group != "SPQC",]


samples_list <- str_c(sample_groups$Group, collapse = ", ")
