'if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDb", force = TRUE)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("ComplexHeatmap", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("hgu95av2.db", force = TRUE)
BiocManager::install("biomaRt", force = TRUE)
BiocManager::install("plotly", force = TRUE)
BiocManager::install("tidyverse", force = TRUE)
BiocManager::install("vctrs", force = TRUE)
BiocManager::install("matrixStats", force = TRUE)
BiocManager::install("cluster", force = TRUE)
BiocManager::install("factoextra", force = TRUE)
BiocManager::install("fpc", force = TRUE)
BiocManager::install("dendextend", force = TRUE)
BiocManager::install("cluster.stats", force = TRUE)
BiocManager::install("sva", force = TRUE)'

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
'library("org.Hs.eg.db")
library(tidyverse)
require(hgu95av2.db)'
library(stringr)
require(biomaRt)
library(matrixStats)
library(cluster)
library(fpc)
library(dendextend)
library(sva)


# setwd <- setwd("~/PycharmProjects/RNASeqProject/Results/OptiDonor_1/")
# setwd <- setwd("~/PycharmProjects/RNASeqProject/Results/OptiDonor_1/")

setwd <- setwd("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_STAR_GENCODE_TRUE_analysis/")
setwd <- setwd("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_STAR_GENCODE_TRUE_analysis/")

heatmap_plotter <- function(input_dir, output_dir, count_data, no_samples, no_percent, 
                            cluster_no, cv_max_threshold, cv_min_threshold, 
                            sum_threshold, removed_batch=NULL, by_method=NULL, normalize_to_ref=FALSE, normalize_gene="ENSG00000089157") {
  #'Debugging'
# input_dir <- setwd
# output_dir <- sprintf("%s/R_outputs", input_dir)
# count_data <- read.csv(sprintf("%s/final_countdown_count.csv", input_dir))
# removed_batch = c("AD371", "AD374" ,"BFII.110", "BFII.215", "BFII.216", "MRB003.23",	"MRB005.23")
# if (!is.null(removed_batch)) {
#   count_data <- count_data %>% dplyr::select(-c(removed_batch))
# }
# no_samples <- 16
# no_percent <- 0.1
# cluster_no <- 4
# cv_max_threshold <- 100
# cv_min_threshold <- 20
# sum_threshold <- no_samples*10
# by_method = "ensembl_gene" #"hgnc_symbol" "ensembl_transcript"
# normalize_to_ref <- TRUE
# normalize_gene <- "ENSG00000089157"
# 
# if (!is.null(removed_batch)) {
#   count_data <- count_data %>% dplyr::select(-c(removed_batch))
# }
# rownames(count_data) <- count_data$ensembl_gene
# count_data <- count_data %>% dplyr::select(-c(ensembl_gene))
  
count_data <- count_data[which(rowSums(count_data) > sum_threshold), ]
count_data <- na.omit(count_data)
  if (by_method == "hgnc_symbol") {
  count_data$ensemblt <- gsub("\\..*","", rownames(count_data))
  
  ensembl95 <- useEnsembl(biomart = 'genes', 
                          dataset = 'hsapiens_gene_ensembl', version = 95)
  attri = c("ensembl_transcript_id","ensembl_gene_id", "hgnc_symbol")
  
  annotLookup <- getBM(attributes=attri, filter="ensembl_transcript_id", 
                       values=count_data$ensemblt, mart=ensembl95, uniqueRows=TRUE)
  annotLookup$ensemblt <- annotLookup$ensembl_transcript_id
  annotLookup <- annotLookup %>% dplyr::select(-c(ensembl_transcript_id))
  count_data$ensembl_transcript <- rownames(count_data)
  count_data <- merge(count_data, annotLookup, by="ensemblt")
  count_data <- count_data[!duplicated(count_data[c("hgnc_symbol")]), ]
  range_samples <- (1:no_samples+1)
  } else {
  range_samples <- (1:no_samples)
}
if (normalize_to_ref == TRUE) {
  normalize_gene_df <- count_data[normalize_gene,]
  normalize_gene_df <- normalize_gene_df[, 1:no_samples]
  count_data <- count_data[,1:no_samples]/normalize_gene_df[col(count_data[,1:no_samples])]
  df <- as.data.frame(lapply(count_data, log10))
  rownames(df) <- rownames(count_data)
  clean_df <- df %>%
    filter_all(all_vars(is.finite(.)))
  count_data <- clean_df

}
i <- c(range_samples)
count_data[ , i] <- apply(count_data[ , i], 2,            # Specify own function within apply
                          function(x) as.numeric(as.character(x)))
count_data$var <- apply(count_data[, range_samples], 1, var)
count_data$sd <- apply(count_data[, range_samples], 1, sd)
count_data$mean <- apply(count_data[, range_samples], 1, mean)
count_data$cv <- count_data$sd / count_data$mean
count_data$sum <- apply(count_data[, range_samples], 1, sum)
sum_totals <- colSums(count_data[, range_samples])

#for gene linearity
if (by_method == "hgnc_symbol") {
  count_data <- na.omit(count_data)
  rownames(count_data) <- count_data$hgnc_symbol
  file1 <- "final_data_svf_with_genes"
  write.csv(count_data, file = sprintf("%s/%s.csv", output_dir, file1))
  
  count_data_var <- count_data[order(-count_data$var),]

  count_data_var <- count_data_var %>% dplyr::select(-c("ensemblt", "ensembl_transcript", "var", "sd", "mean", "cv", "sum", "hgnc_symbol")) # 
  file2 <- "final_data_svf_transcription_annotation_var"
  write.csv(count_data_var, file = sprintf("%s/%s.csv", output_dir, file2))
  count_data_var <- count_data_var %>% dplyr::select(-c("ensembl_gene_id"))
  file3 <- "final_data_svf_transcription_annotation"
  write.csv(count_data_var, file = sprintf("%s/%s.csv", output_dir, file3))
  
  count_data_cv <- count_data[order(-count_data$cv),]
  count_data_cv <- count_data_cv %>% dplyr::select(-c("ensemblt", "ensembl_transcript", "var", "sd", "mean", "cv", "sum", "hgnc_symbol", "ensembl_gene_id")) # 
  file_cv <- "final_data_svf_transcription_annotation_cv"
  write.csv(count_data_var, file = sprintf("%s/%s.csv", output_dir, file_cv))
}
if (by_method == "ensembl_transcript" | by_method == "ensembl_gene") {
  count_data <- na.omit(count_data)
  count_data_var <- count_data[order(-count_data$var),]
  # count_data_var <- count_data_var %>% dplyr::select(-c("var", "sd", "mean", "cv", "sum")) #
  file4 <- "final_data_svf_transcription"
  write.csv(count_data_var, file = sprintf("%s/%s.csv", output_dir, file4))
}


#for Heatmap

count_data_var_filtered <- count_data[which(count_data$cv < cv_max_threshold), ]
count_data_var_filtered <- count_data_var_filtered[which(count_data_var_filtered$cv > cv_min_threshold), ]
# count_data_var_filtered <- count_data_var_filtered[which(count_data_var_filtered$sum > sum_threshold), ]
count_data_var_filtered <- count_data_var_filtered[order(-count_data_var_filtered$var),]

if (by_method == "hgnc_symbol") {
  count_data_var_filtered <- count_data_var_filtered %>% dplyr::select(-c("ensemblt", "ensembl_transcript", "var", "sd", "mean", "cv", "sum", "hgnc_symbol")) # 
  file5 <- "final_data_svf_transcription_annotation_var_filtered"
  write.csv(count_data_var_filtered, file = sprintf("%s/%s.csv", output_dir, file5))
  count_data_var_filtered <- count_data_var_filtered %>% dplyr::select(-c("ensembl_gene_id"))
}
if (by_method == "ensembl_transcript" | by_method == "ensembl_gene")  {
  count_data_var_filtered <- count_data_var_filtered %>% dplyr::select(-c("var", "sd", "mean", "cv", "sum")) # 
  file6 <- "final_data_svf_transcription_var_filtered"
  write.csv(count_data_var_filtered, file = sprintf("%s/%s.csv", output_dir, file6))
}

mat.z_var <- t(apply(count_data_var_filtered, 1, scale))
colnames(mat.z_var) <- colnames(count_data_var_filtered)
mat.z_var <- na.omit(mat.z_var)
file7 <- "final_data_svf_mat_z"
write.csv(mat.z_var, file = sprintf("%s/%s.csv", output_dir, file7))
# top percent
shown_genes <- ceiling(nrow(mat.z_var)*no_percent)
mat.z_var <- head(mat.z_var, shown_genes)
count_data_var_filtered_top10 <- head(count_data_var_filtered, shown_genes)

file8 <-  sprintf("final_data_svf_mat_z_var_%s", shown_genes)
write.csv(mat.z_var, file = sprintf("%s/%s.csv", output_dir, file8))

# #Choosing the number of clusters to be used
# #"elbow method"
# wcss <- vector()
# for (i in 1:10) wcss[i] <- sum(kmeans(mat.z_var, i)$withinss)
# 
# # Plot the WCSS against the number of clusters
# plot(1:10, wcss, type = "b", xlab = "Number of clusters (k)", ylab = "Within-cluster sum of squares")
# 
# #"silhouette method"
# avg_sil_width <- vector()
# for (i in 2:10) {
#   clusters <- kmeans(mat.z_var, i)
#   avg_sil_width[i] <- mean(silhouette(clusters$cluster, dist(mat.z_var)))
# }
# 
# # Plot the average silhouette width against the number of clusters
# 
# plot(1:10, avg_sil_width, type = "b", xlab = "Number of clusters (k)", ylab = "Average silhouette width")
# 
# # ROWS - batches
# wcss <- apply(mat.z_var, 1, function(x) {
#   wcss_row <- vector()
#   for (i in 1:10) {
#     if (length(unique(x)) >= i) {
#       wcss_row[i] <- sum(kmeans(x, i)$withinss)
#     } else {
#       wcss_row[i] <- NA
#     }
#   }
#   return(wcss_row)
# })
# 
# # Calculate the average within-cluster sum of squares for each row
# avg_wcss <- apply(wcss, 1, function(x) mean(x, na.rm = TRUE))
# 
# # Plot the average within-cluster sum of squares against the number of clusters
# plot(1:10, avg_wcss, type = "b", xlab = "Number of clusters (k)", ylab = "Average within-cluster sum of squares")
# 
# if (!is.null(removed_batch)) {
#   file9 <- sprintf("HeatMap_final_data_svf_transcription_annotation_var_%s_%dclusters_wo_%s_%d_batches_by_%s"
#                    , shown_genes, cluster_no, removed_batch, no_samples, by_method)
#   png(sprintf("%s/%s.png", output_dir, file9), height=4600, width=4000, res=500)
#   } else {
file10 <- sprintf("HeatMap_final_data_svf_transcription_annotation_var_%s_%dclusters_%d_batches_by_%s"
                     , shown_genes, cluster_no, no_samples, by_method)
png(sprintf("%s/%s.png", output_dir, file10), height=4600, width=4000, res=500)
  # }

g1 <- Heatmap(mat.z_var, cluster_rows = TRUE, cluster_columns = TRUE,
              column_labels = colnames(mat.z_var), name = "Z-score",
              row_labels = rownames(mat.z_var), show_row_names = FALSE,
              row_km = 8, column_km = cluster_no,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 1),
              col = RColorBrewer::brewer.pal(9, "RdBu"))

ht <- draw(g1)
clusters <- row_order(ht)

num_clusters <- length(clusters)

gene_lists <- list()

for (i in 1:num_clusters) {
  cluster_indices <- clusters[[as.character(4)]] #it is the index in the dat frame not the value
  cluster_genes <- rownames(mat.z_var)[cluster_indices]
  gene_lists[[i]] <- cluster_genes
  
  # Write cluster genes to a CSV file
  file11 <- sprintf("cluster%d_genes_%s_var.csv", i, shown_genes)
  cluster_file <- sprintf("%s/%s", output_dir, file11)
  cluster_data <- mat.z_var[cluster_indices, ]
  cluster_data <- merge(count_data_var_filtered, cluster_data, by.x = 0, by.y = 0)
  write.csv(cluster_data, file = cluster_file)
}

print(g1)
dev.off()
}


input_dir <- setwd
output_dir <- sprintf("%s/20240704R_outputs_with_groups", input_dir)
count_data <- read.csv(sprintf("%s/final_countdown_count.csv", input_dir))
removed_batch <- NULL
# removed_batch = c("AD371", "AD374", "BFII.110", "BFII.215", "BFII.216", "MRB003.23",	"MRB005.23") #"AD371",  ,
if (!is.null(removed_batch)) {
  count_data <- count_data %>% dplyr::select(-c(removed_batch))
}


order_columns <- c("ensembl_gene", "AD369", "AD370", "AD371", "AD372", "AD373", 
                   "AD374", "AD376", "AD377", "BFII.109", "AD379", "AD380", 
                   "AD382", "AD383", "AD384", "BFII.110", "BFII.112", 
                   "BFII.215", "AD385", "AD386", "AD387", "BFII.216", 
                   "MRB003.23", "MRB005.23", "AD371.2", "AD376.2", "AD384.2", 
                   "AD387.2", "AD388", "AD389", "AD390", "AD392", "AD393", 
                   "BFII.702", "MRB006.23")
count_data <- count_data[, order_columns]
rownames(count_data) <- count_data$ensembl_gene
count_data <- count_data %>% dplyr::select(-c(ensembl_gene))
count_data_mtx = as.matrix(count_data)
batch <- c(rep(1, 9), rep(2, 8), rep(3, 6), rep(4, 11))
group<- c(0,0,1,0,0,0,2,0,0,0,0,0,0,3,0,0,0,0,0,4,0,0,0,1,2,3,4,0,0,0,0,0,0,0)

adjusted_count_data_mtx <- ComBat_seq(count_data_mtx, batch=batch, group=group)
adjusted_count_data <- as.data.frame(adjusted_count_data_mtx)
file_combat <- "final_data_svf_combat_seq"
write.csv(adjusted_count_data, file = sprintf("%s/%s.csv", output_dir, file_combat))


if (!is.null(removed_batch)) {
  no_samples <- ncol(adjusted_count_data)
} else {no_samples <- ncol(count_data)}
# no_samples <- 16
no_percent <- 1
cluster_no <- 4
cv_max_threshold <- 10
cv_min_threshold <- 0
sum_threshold <- no_samples*5
by_method = "ensembl_transcript"#"ensembl_gene" #"hgnc_symbol" "ensembl_transcript"  



heatmap_plotter(input_dir = input_dir, output_dir = output_dir, count_data = adjusted_count_data, 
                no_samples = no_samples, no_percent = no_percent, cluster_no = cluster_no, 
                cv_max_threshold = cv_max_threshold, cv_min_threshold = cv_min_threshold, sum_threshold = sum_threshold, 
                removed_batch = removed_batch, by = by_method, normalize_to_ref = FALSE, normalize_gene = "ENSG00000089157")


ensembl_geneID_to_hgnc <- function(data) {
  ensembl95 <- useEnsembl(biomart = 'genes', 
                          dataset = 'hsapiens_gene_ensembl', version = 95)
  attri = c("ensembl_gene_id", "hgnc_symbol")
  
  annotLookup <- getBM(attributes=attri, filter="ensembl_gene_id", 
                       values=data$Gene, mart=ensembl95, uniqueRows=TRUE)
  data$ensembl_gene_id <- data$Gene
  data <- data %>% dplyr::select(-c(Gene))
  data <- merge(data, annotLookup, by="ensembl_gene_id")
  return(data)
}

data <- read.csv("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_STAR_analysis/R_outputs/20240108_Linear_reg_Autophagy_log(IDO activity)_IC50_FGF-7/significant_features.csv")
hgnc_data <- ensembl_geneID_to_hgnc(data = data)
write.csv(hgnc_data, file = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_STAR_analysis/R_outputs/20240108_Linear_reg_Autophagy_log(IDO activity)_IC50_FGF-7/significant_features_hgnc.csv", row.names = FALSE)


