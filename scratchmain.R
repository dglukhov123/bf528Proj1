library(biomaRt)
library(tidyverse)
library(MatrixGenerics)

load_expression <- function(filepath) {
  return(readr::read_csv(filepath)) #changing the source file's header seemed like the least work-intensive solution, please use the file provided in the repo to see the change
}
filter_20 <- function(df){
  tibble <- as_tibble(df)
  newtib <- dplyr::mutate(tibble, keep_rows1 = (rowSums(dplyr::select(tibble, starts_with("GSM")) > log2(15)) / (ncol(tibble)-1)) >= 0.2)
  return(dplyr::filter(newtib, keep_rows1 == TRUE))
}

calc_var <- function(tibble){
  newtib <- dplyr::mutate(tibble, row_var = (rowVars(dplyr::select(tibble, starts_with("GSM")), rows)))
  return(newtib)
}

remove_upper <- function(tibble){
  newtib <- dplyr::mutate(tibble, keep_rows2 = (rowVars(dplyr::select(tibble, starts_with("GSM")) < .995) / (ncol(tibble)-1)) >= 0.2)
  return(dplyr::filter(newtib, keep_rows2 == TRUE))
}
remove_lower <- function(tibble){
  newtib <- dplyr::mutate(tibble, keep_rows3 = (rowVars(dplyr::select(tibble, starts_with("GSM")) > .005) / (ncol(tibble)-1)) >= 0.2)
  return(dplyr::filter(newtib, keep_rows2 == TRUE))
}