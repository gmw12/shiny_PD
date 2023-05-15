cat(file=stderr(), "Shiny_install.R.... starting....", "\n")

package_list <- c('devtools', 'tidyr', 'httr', 'png', 'tidyverse', 'dplyr', 'fs', 'effsize',
                  'colourpicker', 'tibble', 'stringr', 'readxl', 'randomcoloR', 'gplots',
                  'ggpubr', 'rgl', 'pca3d', 'robustbase', 'cluster', 'factoextra',
                  'igraph', 'shiny', 'shinyWidgets', 'shinyFiles', 'rhandsontable', 
                  'shinyjs', 'shinyalert', 'DT', 'ggraph', 'imp4p', 'Peptides',
                  'flexdashboard', 'openxlsx', 'stringi', 'jsonlite', 'remotes', 
                  'BiocManager', 'rAmCharts', 'future', 'promises')


biocmanager_list = c('impute', 'ViSEAGO', 'topGO', 'clusterProfiler', 'GSEABase', 'rWikiPathways', 
                     'STRINGdb', 'limma', 'edgeR', 'pcaMethods', 'gridExtra', 'MASS', 'vsn',
                     'preprocessCore', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db')


cat(file=stderr(), "Shiny_install.R.... package_list", "\n")
# loop to install require packages
for (pack in package_list){
  print(pack)
  if(pack %in% rownames(installed.packages())) {   
    print("not installing")
  }else{
    print("installing")
    install.packages(pack, dependencies = TRUE) 
  }
}


cat(file=stderr(), "Shiny_install.R.... bioconductor_list", "\n")
#loop to install required BioConductor packages
for (pack in biocmanager_list){
  print(pack)
  if(pack %in% rownames(installed.packages())) {   
    print("not installing")
  }else{
    print("installing")
    BiocManager::install(pack, dependencies = TRUE) 
  }
}        


cat(file=stderr(), "Shiny_install.R.... github", "\n")

#library(devtools)

#devtools::install_github('omarwagih/rmotifx', dependencies = TRUE) 

#devtools::install_github("jmwozniak/PTMphinder", dependencies = TRUE)



 
cat(file=stderr(), "Shiny_install.R.... END", "\n")

                              