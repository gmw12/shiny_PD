
set_pathway <- function(input, output, session){
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
      wp.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt", destpath = dpmsr_set$file$string)
      wp2gene <- clusterProfiler::read.gmt(str_c(dpmsr_set$file$string,wp.gmt))
      }
    if (tax_choice == "Mouse"){
      wp.gmt <- rWikiPathways::downloadPathwayArchive(date = "20191010", organism="Mus musculus" , format = "gmt", destpath = dpmsr_set$file$string)
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
  

  #myGENE2GO<-ViSEAGO::annotate("human", Uniprot)
  
  cat(file=stderr(), str_c("Uniprot/ViSEAGO download has ", nrow(dpmsr_set$pathway$myGENE2GO@MF), " entries"), "\n")
  cat(file=stderr(), "Set Pathway...complete" , "\n")
  
  ###------------------------------------------------------
  
  cat(file=stderr(), "Setup String..." , "\n")
  
  
  dpmsr_set$string$string_db <<- NULL
  tax_choice <- input$select_organism
  cat(file=stderr(), str_c("organism...", tax_choice), "\n")
  #string_species <- get_STRING_species(version=10)
  
  if (version$major < 4){
    string_version <- "10"
  }else{
    string_version <- "11.0"
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
  
  cat(file=stderr(), str_c("stringdb object created"), "\n")
  
  
  
}
  
  
  