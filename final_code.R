library(biomaRt)
library(tidyverse)
library(tibble)
library(readr)
library(stringr)
library(MatrixGenerics)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(plotly)
library(ClassDiscovery)
# Load Data
load_expression <- function(filepath) {
  return(readr::read_csv(filepath)) #changing the source file's header seemed like the least work-intensive solution, please use the file provided in the repo to see the change
}
# first filter function 4.1
filter_20 <- function(df){
  tibble <- as_tibble(df)
  newtib <- dplyr::mutate(tibble, keep_rows1 = (rowSums(dplyr::select(tibble, starts_with("GSM")) > log2(15)) / (ncol(tibble)-1)) >= 0.2)
  return(dplyr::filter(newtib, keep_rows1 == TRUE))
}
# second filter function 4.2
chisqfilter <- function(tib, N){
  chisc_lower <- qchisq(0.01/2 ,df = N-1)
  chisc_upper <- qchisq(0.01/2, df = N-1, lower.tail = FALSE)
  new_tib <- dplyr::mutate(tib, keep_row_upper = population_variance > (((N-1)*row_var)/chisc_lower))%>%
    dplyr::mutate(tib, keep_row_lower = population_variance < (((N-1)*row_var)/chisc_upper))%>%
    dplyr::select(-contains("row_var"))%>%
    dplyr::select(-contains("population_variance"))
  return(dplyr::filter(new_tib, keep_row_upper == TRUE | keep_row_lower == TRUE))
}
# 4.3 function
coeff_filter <- function(tibble){
  newtib <- dplyr::mutate(tibble, row_coeff = (rowSds(as.matrix(dplyr::select(tibble, starts_with("GSM"))))/rowMeans(dplyr::select(tibble, starts_with("GSM"))) >= 0.186))
  return(dplyr::filter(newtib, row_coeff == TRUE))
}
# implementing our functions
data <- load_expression("/projectnb/bf528/users/group_7/project_1/proj1_intensity_data.csv")
#(data) we filtered out outlier #GSM972097
filtered <- filter_20(dplyr::select(data, -contains("GSM972097")))
data_frame <- as.data.frame(filtered)
data_frame$keep_rows <- NULL
data_frame$row_var = apply(data_frame[,-1], 1, var)
row.names(data_frame) <- data_frame[,1]
data_frame[,1] <- NULL
df1 <- as.data.frame(filtered)
df1 <- df1[,-1]
df1$keep_rows <- NULL
mat1 <- as.matrix(df1)
gene_list <- c(mat1)
median_row_variance <- median(data_frame$row_var)
N <- ncol(data_frame)-1
data_frame$population_variance <- median_row_variance
probeids <- as.vector(filtered[,1])
data_frame$keep_rows1 <- NULL
data_frame2 <- cbind(data_frame, probeids)
data_tib <- as_tibble(data_frame2)
data_tib_ordered <- dplyr::relocate(data_tib, any_of(c("probeids","row_var", "population_variance")))
filtered2 <- chisqfilter(data_tib_ordered, N)
#4.3
final_filter_tib <- coeff_filter(filtered2)
final_filter <- as.data.frame(final_filter_tib)
row.names(final_filter) <- final_filter[,1]
final_filter[,1] <- NULL
final_filter[,134:136] <- NULL
df_4 <- final_filter
colnames(df_4) <- as.character(substring(colnames(final_filter),1,9))
final_filter <- df_4
#write to file
write.table(final_filter, file = 'expression_data_after_filtering.csv', sep = ',')
#5
cor_matrix <- cor(final_filter, method = 'pearson')
dist_matrix <- as.dist(1 - cor_matrix)
cluster_tree <- hclust(dist_matrix, method = 'ward.D2')
cluster_assignments <- cutree(cluster_tree, k = 2)
p<- read.csv('/projectnb/bf528/users/group_7/project_1/code/proj_metadata.csv', header = TRUE, row.names = 1)
subtype <- p[,c('geo_accession','cit.coloncancermolecularsubtype')]
matching_rows <- subtype[subtype$geo_accession %in% colnames(final_filter), ]
annotations <- matching_rows$cit.coloncancermolecularsubtype
colors <- c('blue', 'red')
sample_colors <- ifelse(annotations == 'C3', colors[2], colors[1])
heatmap(data.matrix(final_filter), ColSideColors = sample_colors, Colv = as.dendrogram(cluster_tree))
legend(x='right', legend=c('C3(57-1=56)', 'C4(76+1=77)'),fill=c('red','blue'), title='Tumor Subtype',cex = 0.6)
#finding the clusters.
cluster_1 <- final_filter[, cluster_assignments == 1]
cluster_2 <- final_filter[, cluster_assignments == 2]
#carry out the welch t-test (FDR)
test_results <- apply(final_filter, 1, function(x) t.test(x ~ cluster_assignments))
p_values <- sapply(test_results, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = 'fdr')
results <- data.frame(probeset_ID = rownames(final_filter),t_statistic = sapply(test_results, function(x) x$statistic),
                      p_value = p_values,
                      adjusted_p_value = adjusted_p_values
)
# write out the data frame to file
write.csv(results, 'differential_expression_results.csv', row.names = FALSE)
# number of differentially expressed genes
differentially_expressed_genes <- sum(results$adjusted_p_value < 0.05)
cat('Number of differentially expressed genes:', differentially_expressed_genes, '\n')
# Find the top differentially expressed genes
topmost <- head(results[order(results$adjusted_p_value), ], 5)$probeset_ID
print(topmost)
