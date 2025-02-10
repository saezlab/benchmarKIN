# libraries
library(dplyr)
library(httr)
library(readxl)
library(usethis)
library(stringr)

## code to prepare `tyrosine` dataset ------
# URL of the file to download
tyrosine <- "https://zenodo.org/records/14824013/files/benchmark_data_tyrosine.csv?download=1"

# Perform the GET request to fetch the file content
response <- httr::GET(tyrosine)
# Read the content directly into R
file_content <- httr::content(response, as = "raw")
# Create a temporary file connection
temp_file <- base::tempfile()
# Write the content to the temporary file
base::writeBin(file_content, temp_file)

# Benchmark data
# load logFC and annotations
tyrosineData <- utils::read.csv(temp_file)
tyrosineData <- tibble::as_tibble(tyrosineData)
colnames(tyrosineData) <- stringr::str_remove(colnames(tyrosineData), "X")

## code to prepare `tyrosine` metadata ------
# URL of the Excel file to download
tyrosine_meta <- "https://zenodo.org/records/14824013/files/benchmark_metadata_tyrosine.csv?download=1"

# Perform the GET request to fetch the file content
response <- httr::GET(tyrosine_meta)
# Read the content directly into R
file_content <- httr::content(response, as = "raw")
# Create a temporary file connection
temp_file <- base::tempfile()
# Write the content to the temporary file
base::writeBin(file_content, temp_file)

# Benchmark data
# Meta data
metaData <- utils::read.csv(temp_file)
tyrosineMeta <- tibble::as_tibble(metaData) %>%
  dplyr::mutate(treatment = dplyr::case_when(
    stringr::str_detect(id, "PMID_24362263|Dasatinib") ~ "Dasatinib",
    stringr::str_detect(id, "QSAV|EGF") ~ "EGF",
    stringr::str_detect(id, "HRG") ~ "HRG",
    stringr::str_detect(id, "sec") ~ "TCRstimulation",
    stringr::str_detect(id, "Imatinib") ~ "Imatinib",
    stringr::str_detect(id, "Bosutinib") ~ "Bosutinib",
    stringr::str_detect(id, "Nilotinib") ~ "Nilotinib"
  )) %>%
  dplyr::mutate(PMID = dplyr::case_when(
    stringr::str_detect(id, "PMID_24362263") ~ "24362263",
    stringr::str_detect(id, "QSAV") ~ "17389395",
    stringr::str_detect(id, "24H") ~ "17016520",
    stringr::str_detect(id, "sec") ~ "25147952",
    stringr::str_detect(id, "Log2") ~ "24804581",
    stringr::str_detect(id, "min") ~ "30257219"
  ))

## add kinase class
kinase_class <- "https://zenodo.org/records/14824013/files/kinase_class.csv?download=1"

# Perform the GET request to fetch the file content
response <- httr::GET(kinase_class)
# Read the content directly into R
file_content <- httr::content(response, as = "raw")
# Create a temporary file connection
temp_file <- base::tempfile()
# Write the content to the temporary file
base::writeBin(file_content, temp_file)

kin_class <- utils::read.csv(temp_file)
tyrosineMeta <- tyrosineMeta %>%
  dplyr::left_join(kin_class, by = c("target" = "source"))

## Merge and save `perturbData` dataset ------
usethis::use_data(tyrosineData, overwrite = TRUE)
usethis::use_data(tyrosineMeta, overwrite = TRUE)
