df<- dpmsr_set$data$final$sltmm

testallgene <- AnnotationDbi::select(org.Hs.eg.db, keys=selection, columns = c("ENTREZID","GO", "GENENAME"), keytype = "UNIPROT")


testallgene <- testallgene[testallgene$ONTOLOGY=="BP",]
golookup <- data.frame(cbind(testallgene$UNIPROT, testallgene$GO))
goolookup <- golookup %>% distinct()
goolookup <- na.omit(goolookup)

v_results_df <- v_results_df[order(v_results_df$condition.pvalue  )]
firstgo <- v_results_df$GO.ID[1]
goproteins <- goolookup[goolookup$X2==firstgo,]


test1 <- goolookup %>% group_by(X1) %>% mutate(Y1 = paste0(X2, collapse = ",")) 
test1$X2 <- NULL
testtme <- unique(test1)
testtme$X1 <-as.character(testtme$X1)


mergedf <- merge(x=df, y=testtme, by.x="Accession", by.y="X1")


sub_df <- mergedf[grep("GO:0001732", mergedf$Y1),]
