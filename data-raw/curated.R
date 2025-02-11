# libraries
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(usethis)

## code to prepare `curated` dataset, combination from phosphositeplus, PTMsigDB and the gold standard set from GPS 5.0 ------
# URL of the file to download
curated_file <- "https://zenodo.org/records/14824013/files/curated_library.csv?download=1"

# Perform the GET request to fetch the file content
response <- httr::GET(curated_file)
# Read the content directly into R
file_content <- httr::content(response, as = "raw")
# Create a temporary file connection
temp_file <- base::tempfile()
# Write the content to the temporary file
base::writeBin(file_content, temp_file)

# load curated targets
curated <- read.csv(curated_file)

## Prepare data ---------------------------
## Include ENSEMBL gene name
source_keys <- unique(curated$target_protein)
translate_df <- AnnotationDbi::select(org.Hs.eg.db, keys=toupper(source_keys),
                                      columns=c("ENSEMBL","SYMBOL"), keytype="SYMBOL")

## Manual renaming for kinases where the alias was not found
translate_df <- translate_df %>%
  dplyr::distinct(SYMBOL, .keep_all = T) %>%
  dplyr::rename("target_protein" = SYMBOL)


## Rename kinases in prior knowledge
curated <- dplyr::left_join(curated, translate_df,
                                  by = "target_protein", relationship = "many-to-many") %>%
  dplyr::distinct()

## Save processed phosphositeplus
usethis::use_data(curated, overwrite = TRUE)
