### R code from vignette source 'STRINGdb.Rnw'

###################################################
### code chunk number 1: initialization
###################################################
library(STRINGdb)
string_db <- STRINGdb$new( version="11", species=9606, 
                           score_threshold=200, input_directory="")


###################################################
### code chunk number 2: help
###################################################
STRINGdb$methods()              # To list all the methods available.
STRINGdb$help("get_graph")      # To visualize their documentation.


###################################################
### code chunk number 3: load_data
###################################################
data(diff_exp_example1)
head(diff_exp_example1)


###################################################
### code chunk number 4: map
###################################################
example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )


###################################################
### code chunk number 5: STRINGdb.Rnw:111-113
###################################################
options(SweaveHooks=list(fig=function()
  par(mar=c(2.1, 0.1, 4.1, 2.1))))


###################################################
### code chunk number 6: get_hits
###################################################
hits <- example1_mapped$STRING_id[1:200]  


###################################################
### code chunk number 7: plot_network
###################################################
getOption("SweaveHooks")[["fig"]]()
string_db$plot_network( hits )


###################################################
### code chunk number 8: add_diff_exp_color
###################################################
# filter by p-value and add a color column 
# (i.e. green down-regulated gened and red for up-regulated genes)
example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, pvalue<0.05), 
                                                        logFcColStr="logFC" )    


###################################################
### code chunk number 9: post_payload
###################################################
# post payload information to the STRING server
payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id, 
                                      colors=example1_mapped_pval05$color )


###################################################
### code chunk number 10: plot_halo_network
###################################################
getOption("SweaveHooks")[["fig"]]()
# display a STRING network png with the "halo"
string_db$plot_network( hits, payload_id=payload_id )


###################################################
### code chunk number 11: enrichment
###################################################
enrichment <- string_db$get_enrichment( hits )
head(enrichment, n=20)


###################################################
### code chunk number 12: background (eval = FALSE)
###################################################
## backgroundV <- example1_mapped$STRING_id[1:2000]   # as an example, we use the first 2000 genes                                                    
## string_db$set_background(backgroundV)


###################################################
### code chunk number 13: new_background_inst (eval = FALSE)
###################################################
## string_db <- STRINGdb$new( score_threshold=200, backgroundV = backgroundV )


###################################################
### code chunk number 14: enrichment
###################################################
annotations <- string_db$get_annotations( hits )
head(annotations, n=20)


###################################################
### code chunk number 15: clustering1
###################################################
# get clusters
clustersList <- string_db$get_clusters(example1_mapped$STRING_id[1:600])


###################################################
### code chunk number 16: STRINGdb.Rnw:230-232
###################################################
options(SweaveHooks=list(fig=function()
  par(mar=c(2.1, 0.1, 4.1, 2.1))))


###################################################
### code chunk number 17: clustering2
###################################################
getOption("SweaveHooks")[["fig"]]()
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}


###################################################
### code chunk number 18: proteins
###################################################
string_proteins <- string_db$get_proteins()


###################################################
### code chunk number 19: atmtp
###################################################
tp53 = string_db$mp( "tp53" )
atm = string_db$mp( "atm" )


###################################################
### code chunk number 20: neighbors (eval = FALSE)
###################################################
## string_db$get_neighbors( c(tp53, atm) )


###################################################
### code chunk number 21: interactions
###################################################
string_db$get_interactions( c(tp53, atm) )


###################################################
### code chunk number 22: paralogs (eval = FALSE)
###################################################
## # Get all homologs of TP53 in human.
## string_db$get_paralogs(tp53)


###################################################
### code chunk number 23: Closest homologs from other species (eval = FALSE)
###################################################
## # get the best hits of the following protein in all the STRING species
## string_db$get_homologs_besthits(tp53)


###################################################
### code chunk number 24: homologs_besthits in target species (eval = FALSE)
###################################################
## # get the homologs of the following two proteins in the mouse (i.e. species_id=10090)
## string_db$get_homologs_besthits(c(tp53, atm), target_species_id=10090, bitscore_threshold=60)

