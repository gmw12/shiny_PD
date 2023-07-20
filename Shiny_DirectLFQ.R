library(reticulate)

use_python("/home/dpmsr/miniconda3/envs/directlfq/bin/python3")

directlfq <- import("directlfq.lfq_manager")

df <- dpmsr_set$data$normalized$impute

df_data <- df %>% dplyr::select(contains(c("Accession", "Sequence", "Modifications")))
df_data$index <- 1:nrow(df_data)

df_data$protein <- df_data$Accession
df_data$ion <- str_c(df_data$Sequence, "_", df_data$Modifications, "_", df_data$index)
df_data$Accession <- NULL
df_data$Modifications <- NULL
df_data$Sequence <- NULL
df_data$index <- NULL

df_data <- cbind(df_data, df[,14:ncol(df)])

directlfq_file <- str_c(dpmsr_set$file$data_dir, "directlfq.aq_reformat.tsv")
write.table(df_data, file = directlfq_file, sep = "\t", row.names = FALSE)

directlfq$run_lfq(directlfq_file)

testme <- directlfq$lfqnorm
