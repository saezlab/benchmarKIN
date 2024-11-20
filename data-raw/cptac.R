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
cptacGS_2pt5 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/protein_2pt5percent.Rds")
cptacGS_5 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/protein_5percent.Rds")
cptacGS_10 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/protein_10percent.Rds")
cptacGS_15 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/protein_15percent.Rds")

## Save processed cptacData
usethis::use_data(cptacGS_2pt5, overwrite = TRUE)
usethis::use_data(cptacGS_5, overwrite = TRUE)
usethis::use_data(cptacGS_10, overwrite = TRUE)
usethis::use_data(cptacGS_15, overwrite = TRUE)

# load goldS Standard set based on activating sites
cptacGS_act_2pt5 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/actsite_2pt5percent.Rds")
cptacGS_act_5 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/actsite_5percent.Rds")
cptacGS_act_10 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/actsite_10percent.Rds")
cptacGS_act_15 <- readRDS("../kinase_benchmark/data/tumor_benchmark/GSsets/actsite_15percent.Rds")

## Save processed cptacData
usethis::use_data(cptacGS_act_2pt5, overwrite = TRUE)
usethis::use_data(cptacGS_act_5, overwrite = TRUE)
usethis::use_data(cptacGS_act_10, overwrite = TRUE)
usethis::use_data(cptacGS_act_15, overwrite = TRUE)

