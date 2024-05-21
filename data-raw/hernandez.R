# libraries
library(dplyr)
library(httr)
library(readxl)
library(usethis)

## code to prepare `hernandez` dataset ------
# URL of the Excel file to download
hernandez <- "https://zenodo.org/records/5645208/files/annotations.xlsx?download=1"

# Perform the GET request to fetch the file content
response <- httr::GET(hernandez)
# Read the content directly into R
file_content <- httr::content(response, as = "raw")
# Create a temporary file connection
temp_file <- base::tempfile()
# Write the content to the temporary file
base::writeBin(file_content, temp_file)

# Benchmark data
# load logFC and annotations
phosphosites <- readxl::read_excel(temp_file, sheet = "Phosphosites")
logFC <- readxl::read_excel(temp_file, sheet = "Foldchanges", na = "NA")

perturbData <- dplyr::tibble(cbind(ID = paste0(phosphosites$gene_name, "|",
                                           phosphosites$residues, phosphosites$positions, "|",
                                           phosphosites$ensg, "|",
                                           phosphosites$ensp),
                               logFC))

## code to prepare `hernandez` metadata ------
# URL of the Excel file to download
hernandez_meta <- "https://zenodo.org/records/5645208/files/benchmark_data.xlsx?download=1"

# Perform the GET request to fetch the file content
response <- httr::GET(hernandez_meta)
# Read the content directly into R
file_content <- httr::content(response, as = "raw")
# Create a temporary file connection
temp_file <- base::tempfile()
# Write the content to the temporary file
base::writeBin(file_content, temp_file)

# Benchmark data
# Meta data
metaData <- readxl::read_excel(temp_file, sheet = "KinaseConditionPairs")
metaData <- tibble::as_tibble(metaData) %>%
  dplyr::rename(id = Condition, target = Kinase, sign = Regulation) %>%
  dplyr::filter(id %in% colnames(perturbData))
metaData <- metaData %>%
  dplyr::mutate(sign = base::replace(sign, sign == "up", 1)) %>%
  dplyr::mutate(sign = base::replace(sign, sign == "down", -1))
metaData$sign <- base::as.double(metaData$sign)

#manually change perturbation for 702_225 as this should be a knock in according to publication
metaData$sign[metaData$id == "702_225"] <- 1

## rename meta_data kinases
metaData <- metaData %>%
  dplyr::mutate(target = dplyr::recode(target,
                                       "ABL" = "ABL1"))

## Remove experiments with consistent bad performance
exp_to_remove <- c("705_225", "1298_272", "1288_272", "1291_272", "387_117", "1289_272", "1290_272", "1308_272", "699_225")
hernandezMeta <- metaData %>%
  dplyr::filter(!id %in% exp_to_remove)

hernandezData <- perturbData[colnames(perturbData) %in% c("ID", metaData$id)]

## Merge and save `perturbData` dataset ------
usethis::use_data(hernandezData, overwrite = TRUE)
usethis::use_data(hernandezMeta, overwrite = TRUE)
