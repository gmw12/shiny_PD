set_pathway <- function(input, output, session){
  
  cat(file=stderr(), "Set Pathway..." , "\n")
  tax_choice <- input$select_organism

  if(!is.null(dpmsr_set$tax_choice)){
    if (tax_choice == dpmsr_set$tax_choice){
      cat(file=stderr(), "Tax choice same as previous..." , "\n")
    }else{
      cat(file=stderr(), "Tax has changed..." , "\n")
      dpmsr_set$x$pathway_set <<- 0
    }
  }
  
  cat(file=stderr(), str_c("Pathway tax choice...", tax_choice), "\n")
  
  string_dir <- str_c(getwd(), "/")
  
  load(file=str_c(string_dir,"Pathway_wp2gene_",tax_choice), envir = .GlobalEnv)
  load(file=str_c(string_dir,"Pathway_myGENE2GO_",tax_choice), envir = .GlobalEnv)
  
  if (tax_choice == "Human"){
    cat(file=stderr(), str_c("Load tax library...", tax_choice ), "\n")
    library(org.Hs.eg.db)
    tax_db <<- org.Hs.eg.db}
  
  if (tax_choice == "Mouse"){
    cat(file=stderr(), str_c("Load tax library...", tax_choice ), "\n")
    library(org.Mm.eg.db)
    tax_db <<- org.Mm.eg.db}
  
  if (tax_choice == "Rat"){
    cat(file=stderr(), str_c("Load tax library...", tax_choice ), "\n")
    library(org.Rn.eg.db)
    tax_db <<- org.Rn.eg.db} 
  

  
  if (dpmsr_set$x$pathway_set == 1) {
    cat(file=stderr(), "Pathway previously set, skipping string download..." , "\n")
  }else{
    cat(file=stderr(), "Set String..." , "\n")
    setup_string(session, input, output)
    cat(file=stderr(), "String setup complete..." , "\n")
  }
  
  dpmsr_set$x$pathway_set <<- 1
  dpmsr_set$tax_choice <<- tax_choice
  
  gc()
  cat(file=stderr(), "Pathway setup complete..." , "\n")
}


#-----------------------------------------------------------------------------------------------------------------------------------

set_pathway_old <- function(input, output, session){
  cat(file=stderr(), "Set Pathway...1" , "\n")
  tax_choice <- input$select_organism
  cat(file=stderr(), str_c("Pathway tax choice...", tax_choice), "\n")
  
  #---Wiki Setup----------------------------
  db_get <- function(tax_choice)
  {
    if (tax_choice == "Human"){
      library(org.Hs.eg.db)
      tax_db = org.Hs.eg.db}
    
    if (tax_choice == "Mouse"){
      library(org.Mm.eg.db)
      tax_db = org.Mm.eg.db}
    
    if (tax_choice == "Rat"){
      library(org.Rn.eg.db)
      tax_db = org.Rn.eg.db}   
    
    return(tax_db)
  }
  
  # rat org.Rn.eg.db
  #zebra fish org.Dr.eg.db
  #arabidopsis org.At.tair.db
  #yeast org.Sc.sgd.db
  
  
  #listOrganisms()
  
  if (!dir_exists(dpmsr_set$file$string)) {
    dpmsr_set$file$string <<- create_dir(str_c(dpmsr_set$file$data_dir,"/String"))
  }
  
  
  gmt_get <- function(tax_choice){
    if (tax_choice == "Human"){
      wp.gmt <- rWikiPathways::downloadPathwayArchive(date = "20220110", organism="Homo sapiens", format = "gmt", destpath = dpmsr_set$file$string)
      wp2gene <- clusterProfiler::read.gmt(str_c(dpmsr_set$file$string,wp.gmt))
    }
    if (tax_choice == "Mouse"){
      wp.gmt <- rWikiPathways::downloadPathwayArchive(date = "20220110", organism="Mus musculus" , format = "gmt", destpath = dpmsr_set$file$string)
      wp2gene <- clusterProfiler::read.gmt(str_c(dpmsr_set$file$string,wp.gmt))
    }
    if (tax_choice == "Rat"){
      wp.gmt <- rWikiPathways::downloadPathwayArchive(date = "20220110", organism="Rattus norvegicus" , format = "gmt", destpath = dpmsr_set$file$string)
      wp2gene <- clusterProfiler::read.gmt(str_c(dpmsr_set$file$string,wp.gmt))
    } 
    
    return(wp2gene)
  }
  
  cat(file=stderr(), "Set Pathway...2" , "\n")
  dpmsr_set$pathway$tax_db <<- db_get(tax_choice)
  
  dpmsr_set$pathway$wp2gene <<- gmt_get(tax_choice)
  
  if (version$major < 4){
    dpmsr_set$pathway$wp2gene <<- dpmsr_set$pathway$wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  }else {
    dpmsr_set$pathway$wp2gene <<- try(dpmsr_set$pathway$wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%") )
  }
  
  
  #if (class(dpmsr_set$pathway$wp2gene) == "try-error") {
  #  dpmsr_set$pathway$wp2gene <<- dpmsr_set$pathway$wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  #}
  
  #retired "ont"?
  #dpmsr_set$pathway$wp2gene <<- dpmsr_set$pathway$wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  
  
  
  
  cat(file=stderr(), "Set Pathway...3" , "\n")
  #---ViseaGo/topGo Setup----------------------------
  dpmsr_set$pathway$Uniprot <<- ViSEAGO::Uniprot2GO()
  
  cat(file=stderr(), "Set Pathway...4" , "\n")
  
  if (tax_choice == "Human"){
    dpmsr_set$pathway$myGENE2GO <<- ViSEAGO::annotate(
      "human", dpmsr_set$pathway$Uniprot)
  }
  if (tax_choice == "Mouse"){
    dpmsr_set$pathway$myGENE2GO <<- ViSEAGO::annotate(
      "mouse", dpmsr_set$pathway$Uniprot)
  }
  if (tax_choice == "Rat"){
    dpmsr_set$pathway$myGENE2GO <<- ViSEAGO::annotate(
      "rat", dpmsr_set$pathway$Uniprot)
  }
  
  #myGENE2GO<-ViSEAGO::annotate("human", Uniprot)
  
  
  # Display table of available organisms with Uniprot
  # connect to Uniprot-GOA
  # Uniprot<-ViSEAGO::Uniprot2GO()
  # ViSEAGO::available_organisms(Uniprot)
  
  
  cat(file=stderr(), str_c("Uniprot/ViSEAGO download has ", nrow(dpmsr_set$pathway$myGENE2GO@MF), " entries"), "\n")
  cat(file=stderr(), "Set Pathway...complete" , "\n")
  
  ###------------------------------------------------------
  
  cat(file=stderr(), "Setup String..." , "\n")
  
  
  dpmsr_set$string$string_db <<- NULL
  tax_choice <- input$select_organism
  cat(file=stderr(), str_c("organism...", tax_choice), "\n")
  
  #string_species <- get_STRING_species(version=10)
  # look up tax ID https://www.ncbi.nlm.nih.gov/taxonomy
  
  if (version$major < 4){
    string_version <- "10"
  }else{
    string_version <- "11.5"
  }
  
  cat(file=stderr(), str_c("string version...", string_version), "\n")
  
  if(input$select_organism=="Human"){
    dpmsr_set$string$string_db <<- STRINGdb$new(version=string_version, species=9606,
                                                score_threshold=0, input_directory=dpmsr_set$file$string)
  }
  
  if(input$select_organism=="Mouse"){
    dpmsr_set$string$string_db <<- STRINGdb$new( version=string_version, species=10090,
                                                 score_threshold=0, input_directory=dpmsr_set$file$string)
  } 
  
  if(input$select_organism=="Rat"){
    dpmsr_set$string$string_db <<- STRINGdb$new( version=string_version, species=10116,
                                                 score_threshold=0, input_directory=dpmsr_set$file$string)
  } 
  
  
  cat(file=stderr(), str_c("stringdb object created"), "\n")
  
  
  
}
