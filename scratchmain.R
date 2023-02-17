library(biomaRt)
library(tidyverse)
library(matrixStats)
library(ggplot2)
library(ggdendro)
library(plotly)
library(ClassDiscovery)

load_expression <- function(filepath) {
  return(readr::read_csv(filepath)) #changing the source file's header seemed like the least work-intensive solution, please use the file provided in the repo to see the change
}
filter_20 <- function(df){
  tibble <- as_tibble(df)
  newtib <- dplyr::mutate(tibble, keep_rows1 = (rowSums(dplyr::select(tibble, starts_with("GSM")) > log2(15)) / (ncol(tibble)-1)) >= 0.05)
  return(dplyr::filter(newtib, keep_rows1 == TRUE))
}

chisqfilter <- function(tib, N){
  chisc_lower <- qchisq(0.01/2 ,df = N-1)
  chisc_upper <- qchisq(0.01/2, df = N-1, lower.tail = FALSE)
  new_tib <- dplyr::mutate(tib, keep_row_upper = population_variance > (((N-1)*row_var)/chisc_lower))%>%
              dplyr::mutate(tib, keep_row_lower = population_variance < (((N-1)*row_var)/chisc_upper))%>%
              dplyr::select(-contains("row_var"))%>%
              dplyr::select(-contains("population_variance"))
  return(dplyr::filter(new_tib, keep_row_upper == TRUE | keep_row_lower == TRUE))
}

coeff_filter <- function(tibble){
  newtib <- dplyr::mutate(tibble, row_coeff = (rowSds(as.matrix(dplyr::select(tibble, starts_with("GSM"))))/rowMeans(dplyr::select(tibble, starts_with("GSM"))) >= 0.186))
  return(dplyr::filter(newtib, row_coeff == TRUE))
}

write_to_file <- function(tibble, destination){
  write(tibble, file=destination)
}