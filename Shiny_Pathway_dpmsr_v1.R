#listOrganisms()
#string_species <- get_STRING_species(version=10)
# look up tax ID https://www.ncbi.nlm.nih.gov/taxonomy

#i=1

library(org.Hs.eg.db) 
library(org.Mm.eg.db)
library(org.Rn.eg.db)
  
tax_list <- c("Human", "Mouse", "Rat")
tax_db_list <- c(org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db)
organism_list <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus")
go_list <- c("human", "mouse", "rat")
string_dir <- '/home/dpmsr/shiny_PD/'
  
for (i in (1:length(tax_list))){
  
  tax_choice <- tax_list[i]
  tax_db <- tax_db_list[i]
  organism_choice <- organism_list[i]
  go_choice <- go_list[i]
  
  cat(file=stderr(), str_c(tax_choice, " ", organism_choice, " ", go_choice), "\n" , "\n")
  
  wp.gmt <- rWikiPathways::downloadPathwayArchive(date = "20220110", organism=organism_choice, format = "gmt", destpath = string_dir)
  wp2gene <- clusterProfiler::read.gmt(str_c(string_dir, wp.gmt))

  if (version$major < 4){
    wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  }else {
    wp2gene <- try(wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%") )
  }
    
 
  Uniprot <- ViSEAGO::Uniprot2GO()
  myGENE2GO <- ViSEAGO::annotate(go_choice, Uniprot)
  
  save(wp2gene, file=str_c(string_dir,"Pathway_wp2gene_",tax_choice))
  save(myGENE2GO, file=str_c(string_dir,"Pathway_myGENE2GO_",tax_choice))
  
}
