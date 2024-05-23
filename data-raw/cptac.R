library(dplyr)
library(purrr)
library(usethis)

# load cptac processed phosphoproteomics data
file_paths <- list.files("../kinase_benchmark/data/datasets/CPTAC_phospho/final", full.names = T, pattern = "original") #currently only locale
cptacData <- purrr::map(file_paths, readRDS)
names(cptacData) <- sub(".*/([^_]+)_.*", "\\1", file_paths) %>% toupper()

# turn to long format to decrease file size
cptac_long <- purrr::map_dfr(names(cptacData), function(cancer_type){
  cptacData[[cancer_type]] %>%
    tibble::rownames_to_column("id") %>%
    tidyr::pivot_longer(!id, names_to = "sample", values_to = "statistic") %>%
    tibble::add_column(cancer = cancer_type) %>%
    dplyr::filter(!is.na(statistic))
})

cptacData <- cptac_long

## Save processed cptacData
usethis::use_data(cptacData, compress = "xz", overwrite = TRUE)

# load goldS Standard set based on proteomics
cptacGS <- readRDS("../kinase_benchmark/data/misc/benchmarking_gold_standard_top5per.Rds")

## Save processed cptacData
usethis::use_data(cptacGS, overwrite = TRUE)
