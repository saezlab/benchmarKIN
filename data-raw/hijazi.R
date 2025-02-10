# libraries
library(tidyverse)
library(usethis)

## code to prepare `hijazi` dataset ------
# currently only locale
raw_HL60 <- "../kinase_benchmark/data/datasets/hijazi/HL60_fc.tsv"
p_HL60 <- "../kinase_benchmark/data/datasets/hijazi/HL60_pval.tsv"
raw_MCF7 <- "../kinase_benchmark/data/datasets/hijazi/MCF7_fc.tsv"
p_MCF7 <- "../kinase_benchmark/data/datasets/hijazi/MCF7_pval.tsv"

targets_file <- "../kinase_benchmark/data/datasets/hijazi/targets_hijazi.tsv"
inhibitory_file <- "../kinase_benchmark/data/datasets/hijazi/inhibitor_selectivity.tsv"

## Prepare data ---------------------------
# LogFC of benchmark data HL60
HL60_df <- readr::read_tsv(raw_HL60, col_types = cols())
HL60_pval_df <- readr::read_tsv(p_HL60, col_types = cols())

HL60p_long <- HL60_pval_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "pVal") %>%
  mutate(experiment = str_remove(experiment, "\\.p\\.value"))

HL60FC_long <-HL60_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "logFC") %>%
  mutate(experiment = str_remove(experiment, "\\.fold"))

HL60_mono <- full_join(HL60FC_long, HL60p_long, by = c("sh.index.sites", "experiment")) %>%
  tidyr::separate_rows(sh.index.sites, sep = ";") %>%
  filter(!sh.index.sites == "")

HL60_mono_df <- HL60_mono %>%
  dplyr::group_by(sh.index.sites, experiment) %>%
  filter(pVal == min(pVal)) %>%
  ungroup() %>%
  dplyr::select(-pVal) %>%
  distinct(sh.index.sites, experiment, .keep_all = T) %>%
  pivot_wider(names_from = "experiment", values_from = "logFC") %>%
  mutate(protein = map_chr(str_split(sh.index.sites, "\\("), 1), .after = sh.index.sites) %>%
  mutate(aa = map_chr(str_split(sh.index.sites, "\\("), 2) %>% str_remove("\\)"), .after = protein) %>%
  mutate(id = paste0(protein, "_", aa, "|", protein, "|", aa), .before = sh.index.sites) %>%
  dplyr::select(-protein, -aa, -sh.index.sites) %>%
  # Fix gene names which were converted to dates
  dplyr::mutate(id = stringr::str_replace_all(id, "09/09/18 00:00:00", "SEPT9")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/07/18 00:00:00", "SEPT7")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/06/18 00:00:00", "SEPT6")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/05/18 00:00:00", "SEPT5")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/02/18 00:00:00", "SEPT2"))


# LogFC of benchmark data MCF7
MCF7_df <- read_tsv(raw_MCF7, col_types = cols())
MCF7_pval_df <- read_tsv(p_MCF7, col_types = cols())

MCF7p_long <- MCF7_pval_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "pVal") %>%
  mutate(experiment = str_remove(experiment, "\\.p\\.value"))

MCF7FC_long <-MCF7_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "logFC") %>%
  mutate(experiment = str_remove(experiment, "\\.fold"))

MCF7_mono <- full_join(MCF7FC_long, MCF7p_long, by = c("sh.index.sites", "experiment")) %>%
  tidyr::separate_rows(sh.index.sites, sep = ";") %>%
  filter(!sh.index.sites == "")

MCF7_mono_df <- MCF7_mono %>%
  dplyr::group_by(sh.index.sites, experiment) %>%
  filter(pVal == min(pVal)) %>%
  ungroup() %>%
  dplyr::select(-pVal) %>%
  distinct(sh.index.sites, experiment, .keep_all = T) %>%
  pivot_wider(names_from = "experiment", values_from = "logFC") %>%
  mutate(protein = map_chr(str_split(sh.index.sites, "\\("), 1), .after = sh.index.sites) %>%
  mutate(aa = map_chr(str_split(sh.index.sites, "\\("), 2) %>% str_remove("\\)"), .after = protein) %>%
  mutate(id = paste0(protein, "_", aa, "|", protein, "|", aa), .before = sh.index.sites) %>%
  dplyr::select(-protein, -aa, -sh.index.sites) %>%
  # Fix gene names which were converted to dates
  dplyr::mutate(id = stringr::str_replace_all(id, "09/09/18 00:00:00", "SEPT9")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/07/18 00:00:00", "SEPT7")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/06/18 00:00:00", "SEPT6")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/05/18 00:00:00", "SEPT5")) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "09/02/18 00:00:00", "SEPT2"))

phospho_df <- full_join(HL60_mono_df, MCF7_mono_df, by = "id") %>%
  mutate(id = str_remove_all(id, ";")) %>%
  column_to_rownames("id")

phospho_df[phospho_df == 0] <- NA
hijaziData <- phospho_df[!rowSums(is.na(phospho_df)) == ncol(phospho_df),]
hijaziData <- hijaziData %>%
  tibble::rownames_to_column("ID")

# MetaData
metaData <- data.frame(id = colnames(phospho_df)) %>%
  mutate(cell_line = map_chr(str_split(id, "\\."), 1)) %>%
  mutate(drug = map_chr(str_split(id, "\\."), 2)) %>%
  add_column(sign = -1)

targets_manual <- read_tsv(targets_file, col_types = cols()) %>%
  mutate(target = case_when(
    is.na(Target_TTD) & !is.na(Manual) ~ Manual,
    !is.na(Target_TTD) & is.na(Manual) ~ Target_TTD,
    is.na(Target_TTD) & is.na(Manual) ~ NA
  )) %>%
  filter(!is.na(target))

hijaziMeta <- left_join(metaData, targets_manual %>% dplyr::select(drug, target), by = "drug") %>%
  filter(!is.na(target))

## Select targets based on drug selectivity data
drug_selectivity <- read_tsv(inhibitory_file, col_types = cols())

drug_selectivity_df <- drug_selectivity %>%
  pivot_longer(!kinase, names_to = "drug_id", values_to = "selectivity") %>%
  mutate(dataset = map_chr(str_split(drug_id, "\\.\\."), 1)) %>%
  mutate(drug = map_chr(str_split(drug_id, "\\.\\."), 2)) %>%
  mutate(drug = str_remove_all(drug, "-")) %>%
  mutate(drug = recode(drug,
                       "Silmitasertib" = "CX4945",
                       "Pictilisib" = "GDC0941",
                       "Ribociclib" = "LEE011",
                       "Abemaciclib" = "LY2835219",
                       "Amuvatinib" = "MP470"))

targets_discoverX <- drug_selectivity_df %>%
  filter(dataset == "DiscoverX") %>%
  filter(selectivity < 50) %>%
  distinct(kinase, drug) %>%
  group_by(drug) %>%
  summarise(target_discoverX = paste(kinase, collapse = ";"))

hijaziMeta <- left_join(hijaziMeta, targets_discoverX, by = "drug") %>%
  filter(!is.na(target))  %>%
  tibble::add_column(PMID = 31959955)

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
hijaziMeta_separated <- hijaziMeta %>%
  tidyr::separate_rows(target, sep = ";") %>%
  dplyr::left_join(kin_class, by = c("target" = "source"))

# Collapse target and source back together
hijaziMeta_collapsed <- hijaziMeta_separated %>%
  dplyr::group_by(dplyr::across(-c(target, class))) %>%  # Group by all other columns
  dplyr::summarize(
    target = paste(target, collapse = ";"),
    class = paste(class, collapse = ";"),
    .groups = "drop"
  )

## do the same for discoverX targets
hijaziMeta_separated_disc <- hijaziMeta_collapsed %>%
  tidyr::separate_rows(target_discoverX, sep = ";") %>%
  dplyr::left_join(kin_class %>% rename("class_discoverX" = class), by = c("target_discoverX" = "source"))

# Collapse target and source back together
hijaziMeta_collapsed_disc <- hijaziMeta_separated_disc %>%
  dplyr::group_by(dplyr::across(-c(target_discoverX, class_discoverX))) %>%  # Group by all other columns
  dplyr::summarize(
    target_discoverX = paste(target_discoverX, collapse = ";"),
    class_discoverX = paste(class_discoverX, collapse = ";"),
    .groups = "drop"
  )

hijaziMeta <- hijaziMeta_collapsed_disc %>%
  dplyr::rename("treatment" = drug)

hijaziData <- hijaziData[colnames(hijaziData) %in% c("ID", unique(hijaziMeta$id))]

## Merge and save `perturbData` dataset ------
usethis::use_data(hijaziData, overwrite = TRUE)
usethis::use_data(hijaziMeta, overwrite = TRUE)
