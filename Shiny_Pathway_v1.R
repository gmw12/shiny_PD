
set_pathway <- function(input, output, session){
  cat(file=stderr(), "Set Pathway...1" , "\n")
  tax_choice <- input$select_organism
  
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
  
  gmt_get <- function(tax_choice){
    if (tax_choice == "Human"){
      wp.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt", destpath = dpmsr_set$file$string)
      wp2gene <- clusterProfiler::read.gmt(str_c(dpmsr_set$file$string,wp.gmt))
      }
    if (tax_choice == "Mouse"){
      wp.gmt <- rWikiPathways::downloadPathwayArchive(organism="Mus musculus" , format = "gmt", destpath = dpmsr_set$file$string)
      wp2gene <- clusterProfiler::read.gmt(str_c(dpmsr_set$file$string,wp.gmt))
    }
    return(wp2gene)
  }
  
  cat(file=stderr(), "Set Pathway...2" , "\n")
  dpmsr_set$pathway$tax_db <<- db_get(tax_choice)
  
  dpmsr_set$pathway$wp2gene <<- gmt_get(tax_choice)
  dpmsr_set$pathway$wp2gene <<- dpmsr_set$pathway$wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")

  cat(file=stderr(), "Set Pathway...3" , "\n")
  #---ViseaGo/topGo Setup----------------------------
  dpmsr_set$pathway$Uniprot <<- ViSEAGO::Uniprot2GO()

    if (tax_choice == "Human"){
      dpmsr_set$pathway$myGENE2GO <<- ViSEAGO::annotate(
        "human", dpmsr_set$pathway$Uniprot)
    }
    if (tax_choice == "Mouse"){
      dpmsr_set$pathway$myGENE2GO <<- ViSEAGO::annotate(
        "mouse", dpmsr_set$pathway$Uniprot)
    }
  

  #myGENE2GO<-ViSEAGO::annotate("human", Uniprot)
  
  
  cat(file=stderr(), "Set Pathway...4" , "\n")
  
  
  
  
  
}
  
  
  