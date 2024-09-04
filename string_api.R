cat(file=stderr(), "run string step 1", "\n")
input_fc_up <- log(1.2, 2)
input_fc_down <- log(1/1.2, 2)

input_pval <- 0.05
input_comp <- "KO_IH_v_KO_Nx"

cat(file=stderr(), "run string step 2", "\n")

df <- dpmsr_set$string[[input_comp]]
df <- subset(df, pvalue < input_pval)

cat(file=stderr(), str_c("length of dataframe...", nrow(df)), "\n")

cat(file=stderr(), "run string step 3", "\n")


  df <- subset(df, logFC >= input_fc_up | logFC <= input_fc_down )


df <- df[order(-df$logFC),]

cat(file=stderr(), str_c("length of dataframe...", nrow(df)), "\n")
cat(file=stderr(), "run string step 4", "\n")

if (nrow(df) > 200){
  hits <- df$STRING_id[1:input$protein_number]
}else{
  hits <- df$STRING_id
}

cat(file=stderr(), str_c("number of hits searched...", length(hits)), "\n")

cat(file=stderr(), "run string step 5", "\n")


cat(file=stderr(), "run string step 6", "\n")
string_file_name <- str_c(dpmsr_set$file$string, input_comp, ".png")


cat(file=stderr(), str_c("string file name... ", string_file_name ), "\n")


cat(file=stderr(), "run string step 7", "\n")
#dpmsr_set$string$string_db$get_png(hits, required_score=NULL, network_flavor="evidence", file=string_file_name, payload_id=NULL)

#save string png
cat(file=stderr(), "run string step 8", "\n")


library(httr)
library(jsonlite)

hits <- df$STRING_id

hit_list <- hits[1]
for(i in 2:length(hits)){
  hit_list <- str_c(hit_list,"%0d", hits[i])
}

network_api <- str_c("https://string-db.org/api/image/network?identifiers=", hit_list,"&species=10090","&caller_identity=DukeProteomics" )
res = GET(network_api)
#res = GET("https://string-db.org/api/image/network?identifiers=PTCH1")
res

res2 <- readPNG(res$content)
res3 <- rasterImage(res2, 100, 100, 100, 100)
writePNG(res2, target="testpng2.png")

network2_api <- str_c("https://string-db.org/api/tsv-no-header/get_link?identifiers=", hit_list,"&species=10090","&caller_identity=DukeProteomics" )
res2 = GET(network2_api)
res2
link_network <- rawToChar(res2$content)
link_network <- substr(link_network, 1, nchar(link_network)-2)
link_network



hit_list <- hits[1]
for(i in 2:length(hits)){
  hit_list <- str_c(hit_list,"%0d", hits[i])
}

bg <- dpmsr_set$string$string_db$backgroundV
bg_list <- bg[1]
for(i in 2:length(bg)){
  bg_list <- str_c(bg_list,"%0d", bg[i])
}

network3_api <- str_c("https://string-db.org/api/json/enrichment?identifiers=", hit_list,
                      "&background_string_indentifiers=", "",
                      "&species=10090",
                      "&caller_identity=DukeProteomics" )
enrich = GET(network3_api)

enrich_convert <- rawToChar(enrich$content)
test <- jsonlite::fromJSON(content(enrich, "text"), simplifyVector = FALSE)

enrichment <- data.frame(matrix(ncol=10, nrow = length(test)))
enrich_names <- c("category", "term", "genes", "genes_background", "ncbiTaxonId", "inputGenes", "preferredNames", "p_value", "fdr", "description")
colnames(enrichment) <- enrich_names

for (i in 1:length(test)) {
   enrichment[i,1] <- test[[i]]$category
   enrichment[i,2] <- test[[i]]$term
   enrichment[i,3] <- test[[i]]$number_of_genes
   enrichment[i,4] <- test[[i]]$number_of_genes_in_background
   enrichment[i,5] <- test[[i]]$ncbiTaxonId
   enrichment[i,6] <- str_c(unlist(test[[i]]$inputGenes), collapse = ",")
   enrichment[i,7] <- str_c(unlist(test[[i]]$preferredNames), collapse = ",")
   enrichment[i,8] <- test[[i]]$p_value
   enrichment[i,9] <- test[[i]]$fdr
   enrichment[i,10] <- test[[i]]$description
}


