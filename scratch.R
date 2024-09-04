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
#dpmsr_set$string$string_db$plot_network(hits)

cat(file=stderr(), "run string step 6", "\n")
string_file_name <- str_c(dpmsr_set$file$string, input_comp, ".png")
cat(file=stderr(), str_c("string file name... ", string_file_name ), "\n")


cat(file=stderr(), "run string step 7", "\n")
#dpmsr_set$string$string_db$get_png(hits, required_score=NULL, network_flavor="evidence", file=string_file_name, payload_id=NULL)

#save string png
cat(file=stderr(), "run string step 8", "\n")

# string_plot <- try(dpmsr_set$string$string_db$plot_network(hits, add_link = TRUE, add_summary = TRUE), silent = TRUE)
# cat(file=stderr(), "run string step 9 - plot object created", "\n")
# 
# if ( string_plot != "try-error"){
#   png(filename=string_file_name, units="px", width = 1200, height = 1200)
#   string_plot
#   dev.off()
#   cat(file=stderr(), "run string step 20 - saved plot", "\n")
# }else{
#   shinyalert("Oops!", "StringDB server failed to return data...", type = "error")
# }

png(filename=string_file_name, units="px", width = 1200, height = 1200)
dpmsr_set$string$string_db$plot_network(hits, add_link = FALSE, add_summary = TRUE)
dev.off()
cat(file=stderr(), "run string step 9 - saved plot", "\n")



graph <- dpmsr_set$string$string_db$get_subnetwork(hits)
full.graph <- dpmsr_set$string$string_db$get_graph()
# see how many proteins do you have    
vcount(full.graph)

# find top 200 proteins with the highest degree
top.degree.verticies <- names(tail(sort(degree(full.graph)), 200))

# extract the relevant subgraph
top.subgraph <- induced_subgraph(full.graph, top.degree.verticies)

# count the number of proteins in it
vcount(top.subgraph)










library(httr)
library(jsonlite)
res = GET("https://string-db.org/api/image/network?identifiers=PTCH1")
res

res2 <- readPNG(res$content)
res3 <- rasterImage(res2, 100, 100, 100, 100)
writePNG(res2, target="testpng.png")
