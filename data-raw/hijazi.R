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
  dplyr::select(-protein, -aa, -sh.index.sites)


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
  dplyr::select(-protein, -aa, -sh.index.sites)

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
hijaziMeta <- tidyr::separate_rows(hijaziMeta, target, sep = ";") %>%
  tibble::add_column(PMID = 31959955)

## Merge and save `perturbData` dataset ------
usethis::use_data(hijaziData, overwrite = TRUE)
usethis::use_data(hijaziMeta, overwrite = TRUE)
