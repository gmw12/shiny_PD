
install.packages('tidyr', dependencies = TRUE)


install.packages('tidyverse', dependencies = TRUE)
install.packages('dplyr', dependencies = TRUE)
install.packages('fs', dependencies = TRUE)
install.packages('effsize', dependencies = TRUE)
install.packages('colourpicker', dependencies = TRUE)

install.packages('tibble', dependencies = TRUE)
install.packages('stringr', dependencies = TRUE)
install.packages('readxl', dependencies = TRUE)
install.packages('randomcoloR', dependencies = TRUE)
install.packages('gplots', dependencies = TRUE) 
install.packages('ggpubr', dependencies = TRUE)

install.packages('rgl', dependencies = TRUE)
install.packages('pca3d', dependencies = TRUE)
install.packages('robustbase', dependencies = TRUE)
install.packages('cluster', dependencies = TRUE)    # clustering algorithms
install.packages('factoextra', dependencies = TRUE) # clustering algorithms & visualization
install.packages('igraph', dependencies = TRUE)

install.packages('shiny', dependencies = TRUE)
install.packages('shinyWidgets', dependencies = TRUE)
install.packages('shinyFiles', dependencies = TRUE)
install.packages('rhandsontable', dependencies = TRUE)
install.packages('shinyjs', dependencies = TRUE)
install.packages('shinyalert', dependencies = TRUE)
install.packages('DT', dependencies = TRUE)
install.packages('ggraph', dependencies = TRUE)

install.packages('imp4p', dependencies = TRUE)
install.packages('Peptides', dependencies = TRUE)

install.packages('flexdashboard', dependencies = TRUE)

install.packages('devtools')
require(devtools)
install_github('omarwagih/rmotifx')

install.packages("remotes")
remotes::install_github("jmwozniak/PTMphinder")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('ViSEAGO', dependencies = TRUE)
BiocManager::install('topGO', dependencies = TRUE)
BiocManager::install('clusterProfiler', dependencies = TRUE)
BiocManager::install('GSEABase', dependencies = TRUE)
BiocManager::install('rWikiPathways', dependencies = TRUE)
BiocManager::install('STRINGdb', dependencies = TRUE)
BiocManager::install('limma', dependencies = TRUE)
BiocManager::install('edgeR', dependencies = TRUE)
BiocManager::install('pcaMethods', dependencies = TRUE)
BiocManager::install('gridExtra', dependencies = TRUE)
BiocManager::install('MASS', dependencies = TRUE)
BiocManager::install('vsn', dependencies = TRUE)
BiocManager::install('preprocessCore', dependencies = TRUE)
