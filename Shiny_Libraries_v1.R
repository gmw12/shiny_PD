cat(file=stderr(), "load libraries", "\n")

#shiny
library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinyWidgets)
library(shinyalert)

#parallel and background processing
library(callr)
library(doParallel)
library(parallel)
library(imp4p)
library(preprocessCore)

library(promises)
library(future)
plan(multisession)

#basic
library(tidyverse)
library(data.table)
library(fs)
library(miscTools)
 
#graphics
library(gplots)
library(rgl)
library(colourpicker)
library(randomcoloR)


#not updated on CRAN anymore
test_pca3d <- require(pca3d)
if (!test_pca3d) {
  cat(file = stderr(), "Package pca3d not found, adding from project directory", "\n")
  if(dir.exists('code/Packages/pca3d')){
    cat(file = stderr(), "Loading pca3d from primary", "\n")
    load_pca3d <- try(library("pca3d", lib.loc = str_c(getwd(), '/code/Packages/')), silent = TRUE)
  }else {
    cat(file = stderr(), "Loading pca3d from alternate", "\n")
    library("pca3d", lib.loc = str_c(getwd(), '/Packages/') )
  }
}
 
 
#read write
library(readxl)
library(openxlsx)
#library(tcltk)
 
# #tables
library(rhandsontable)
library(DT)




















#library(cluster)    # clustering algorithms
#library(grid)

#library(effsize)

#library(ggpubr)
#library(pca3d)

##library(limma)
##library(edgeR)

# library(gridExtra)
# library(MASS)
# library(pcaMethods)
# #library(vsn)
# #library(robustbase)
# library(factoextra) # clustering algorithms & visualization

#library(ViSEAGO)
#library(topGO)
#library(clusterProfiler)
#library(GSEABase)
#library(rWikiPathways)
#library(STRINGdb)
#library(igraph)

#library(dplyr)
#library(stringr)
#library(tidyr)
#library(tibble)
