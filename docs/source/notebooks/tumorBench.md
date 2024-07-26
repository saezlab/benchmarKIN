# Tumor-based benchmark

## Introduction

This benchmark is based on the cptac dataset which contains genomics,
proteomics and phosphoproteomics information of multiple patients from
different cancer types. Here, the proteomics data was used to identify a
gold standard set of deregulated kinases per patient which are expected
to have an increased or decreased activity.

The phosphoproteomic data is then used to infer kinase activities and to
evaluate which method and kinase-substrate library is able to
recapitulate the gold standard set of kinases.

For this, we provide the normalized abundance of phosphorylation sites
for each patient across different cancer types.

## Getting set up

We first load the required packages to run the benchmark.

    library(benchmarKIN)
    library(dplyr)
    #> 
    #> Attaching package: 'dplyr'
    #> The following objects are masked from 'package:stats':
    #> 
    #>     filter, lag
    #> The following objects are masked from 'package:base':
    #> 
    #>     intersect, setdiff, setequal, union

## CPTAC data

We then take a look at the phosphoproteomics data from CPTAC. This data
is stored as a long data frame in cptacData. For our analysis we will
transform it to a list of wide data frames for each cancer type which
can be used for activity inference.

    cptac_list <- purrr::map(unique(cptacData$cancer), function(cancer_type){
      cptacData %>%
        dplyr::filter(cancer == cancer_type) %>%
        tidyr::pivot_wider(names_from = "sample", values_from = "statistic", values_fill = NA) %>%
        dplyr::select(-cancer) %>%
        tibble::column_to_rownames("id")
    })

    names(cptac_list) <- unique(cptacData$cancer)

## Activity estimation

Here, we will test the z-score (as implemented by RoKAI) in combination
with the PhosphoSitePlus kinase-substrate library as in example.

In theory any other method can also be tested with that data.

### Kinase-substrate library

For that we have already processed a version of PhosphoSitePlus
(accessed: 19/04/2023) that can be mapped to our phosphorylation site
ids.

    head(phosphositeplus)
    #> # A tibble: 6 × 8
    #>   source  target      target_protein position   mor sequence ENSEMBL ENSEMBLPROT
    #>   <chr>   <chr>       <chr>          <chr>    <dbl> <chr>    <chr>   <chr>      
    #> 1 EIF2AK1 EIF2S1_S52  EIF2S1         S52          1 MILLSEL… ENSG00… <NA>       
    #> 2 EIF2AK1 EIF2S1_S49  EIF2S1         S49          1 IEGMILL… ENSG00… <NA>       
    #> 3 PRKCD   HDAC5_S259  HDAC5          S259         1 FPLRKTA… ENSG00… <NA>       
    #> 4 PRKCD   PTPRA_S204  PTPRA          S204         1 PLLARSP… ENSG00… <NA>       
    #> 5 PRKCD   BCL2_S70    BCL2           S70          1 RDPVART… ENSG00… <NA>       
    #> 6 PRKCD   HNRNPK_S302 HNRNPK         S302         1 GRGGRGG… ENSG00… <NA>

After mapping the phosphorylation sites we can bring it into the right
format to run the run\_zscore function.

    pps_id <- cptacData$id %>% 
      unique()
    map_pps <- data.frame(pps_id = pps_id) %>%
      dplyr::mutate(ENSEMBL = purrr::map_chr(stringr::str_split(pps_id, "\\|"), 1)) %>%
      dplyr::mutate(ENSEMBL = purrr::map_chr(stringr::str_split(ENSEMBL, "\\."), 1))  %>%
      dplyr::mutate(sequence = purrr::map_chr(stringr::str_split(pps_id, "\\|"), 4))
    mapped_phosphositeplus <- phosphositeplus %>%
      dplyr::left_join(map_pps, by = c("ENSEMBL", "sequence"), relationship = "many-to-many") %>%
      dplyr::mutate(target = pps_id) %>%
      dplyr::filter(!is.na(pps_id)) %>%
      dplyr::select(source, target, mor) %>%
      dplyr::distinct()

### Activity inference using the z-score

The z-score calculates a score for each kinase by aggregating the change
in abundance of the direct targets in relation to changes in the
non-targets.

    act_scores <- purrr::map(cptac_list, ~list(zscore_ppsp = run_zscore(mat = .x, network = mapped_phosphositeplus)))

## CPTAC benchmark

### Gold standard set of kinases

For each cancer type, a gold standard set of kinases was developed based
on the protein levels for each kinase. For the gold standard set kinases
in the top and bottom 5% relative to the normal distribution of the
protein levels were selected. Here we can look at the kinases which were
located in the top 5% and are supposed to have a high activity for
patients which BRCA. Additionally we also provide the gold standard sets
using the top and bottom 2.5%, 10% and 15% (cptacGS\_2pt5, cptacGS\_10
and cptacGS\_15) of kinases relative to the normal distribution.

    data(cptacGS_5)
    str(cptacGS_5$BRCA$GS_pos_pairs)
    #> List of 61
    #>  $ AKT1    : chr [1:5] "X18BR017" "X09BR005" "X11BR058" "X20BR006" ...
    #>  $ BLVRA   : chr [1:6] "X11BR050" "X01BR025" "X05BR004" "X11BR074" ...
    #>  $ CDK7    : chr [1:3] "X11BR015" "X11BR013" "X11BR011"
    #>  $ EEF2K   : chr [1:2] "X11BR047" "X05BR009"
    #>  $ NADK    : chr [1:8] "X11BR024" "X01BR030" "X01BR043" "X01BR042" ...
    #>  $ NEK9    : chr [1:5] "X20BR008" "X01BR033" "X11BR072" "X18BR019" ...
    #>  $ PAK1    : chr [1:6] "X05BR043" "X05BR026" "X11BR024" "X03BR011" ...
    #>  $ PKN1    : chr [1:4] "X14BR008" "X11BR003" "X11BR047" "X09BR001"
    #>  $ PRKCA   : chr [1:7] "X05BR045" "X05BR009" "X05BR029" "X01BR044" ...
    #>  $ PRKCD   : chr [1:4] "X18BR010" "X11BR006" "X11BR028" "X11BR075"
    #>  $ PRKCQ   : chr [1:8] "X06BR003" "X01BR043" "X05BR043" "X01BR042" ...
    #>  $ PRKD2   : chr [1:4] "X05BR016" "X11BR012" "X11BR027" "X01BR040"
    #>  $ RPS6KA3 : chr [1:4] "X01BR009" "X21BR002" "X01BR026" "X21BR010"
    #>  $ RPS6KA4 : chr [1:6] "X11BR018" "X11BR032" "X11BR036" "X09BR005" ...
    #>  $ RPS6KB1 : chr [1:6] "X11BR028" "X05BR026" "X20BR008" "X11BR022" ...
    #>  $ RPS6KB2 : chr [1:3] "X01BR033" "X03BR010" "X11BR017"
    #>  $ ATR     : chr [1:4] "X05BR042" "X11BR016" "X604" "X01BR042"
    #>  $ BAZ1B   : chr [1:3] "X11BR023" "X01BR031" "X11BR004"
    #>  $ CDK1    : chr [1:7] "X05BR045" "X01BR018" "X05BR042" "X11BR036" ...
    #>  $ CDK9    : chr [1:3] "X11BR049" "X11BR058" "X01BR023"
    #>  $ CHEK2   : chr [1:8] "X01BR027" "X20BR005" "X20BR002" "X05BR042" ...
    #>  $ CSNK1D  : chr [1:5] "X11BR043" "X11BR049" "X01BR032" "X01BR030" ...
    #>  $ DAPK1   : chr [1:4] "X11BR023" "X18BR002" "X11BR075" "X21BR002"
    #>  $ DAPK2   : chr [1:7] "X20BR002" "X01BR031" "X06BR005" "X11BR036" ...
    #>  $ DCK     : chr [1:3] "X05BR009" "X11BR075" "X21BR001"
    #>  $ EIF2AK2 : chr [1:9] "X11BR023" "X20BR002" "X01BR033" "X03BR010" ...
    #>  $ GRK2    : chr [1:2] "X03BR005" "X11BR054"
    #>  $ GSK3B   : chr [1:6] "X09BR005" "X05BR045" "X11BR032" "X11BR010" ...
    #>  $ LMTK2   : chr [1:2] "X01BR023" "X01BR020"
    #>  $ MAP2K2  : chr [1:6] "X20BR001" "X20BR008" "X05BR009" "X05BR029" ...
    #>  $ MAP2K4  : chr [1:5] "X20BR001" "X06BR014" "X18BR016" "X21BR002" ...
    #>  $ MAP2K5  : chr [1:7] "X06BR003" "X03BR010" "X03BR002" "X21BR002" ...
    #>  $ MAPK13  : chr [1:2] "X01BR017" "X09BR001"
    #>  $ MAPK1   : chr "X05BR005"
    #>  $ MAPK3   : chr [1:2] "X11BR022" "X05BR005"
    #>  $ MAPK6   : chr [1:5] "X05BR042" "X01BR043" "X05BR029" "X05BR043" ...
    #>  $ MAPKAPK5: chr [1:3] "X01BR043" "X05BR026" "X01BR040"
    #>  $ MTOR    : chr [1:2] "X11BR024" "X01BR044"
    #>  $ NEK6    : chr [1:6] "X01BR017" "X21BR010" "X20BR008" "X01BR033" ...
    #>  $ NME1    : chr [1:2] "X06BR003" "X11BR055"
    #>  $ OXSR1   : chr [1:2] "X11BR047" "X11BR054"
    #>  $ PDPK1   : chr [1:8] "X18BR004" "X20BR001" "X03BR010" "X01BR025" ...
    #>  $ PFKL    : chr [1:7] "X20BR001" "X05BR001" "X03BR005" "X18BR003" ...
    #>  $ PKN2    : chr [1:3] "X06BR006" "X11BR030" "X01BR043"
    #>  $ PKN3    : chr "X01BR040"
    #>  $ PRKACA  : chr [1:5] "X01BR025" "X01BR040" "X03BR002" "X01BR042" ...
    #>  $ PRKCE   : chr [1:3] "X11BR027" "X11BR059" "X14BR014"
    #>  $ PRKCI   : chr [1:2] "X01BR044" "X01BR023"
    #>  $ PRKG1   : chr "X01BR033"
    #>  $ RIPK2   : chr [1:8] "X01BR001" "X05BR045" "X11BR053" "X20BR005" ...
    #>  $ RPS6KA6 : chr "X11BR054"
    #>  $ SIK2    : chr [1:2] "X20BR002" "X11BR060"
    #>  $ SLK     : chr [1:8] "X01BR025" "X05BR003" "X11BR072" "X01BR042" ...
    #>  $ SRC     : chr [1:5] "X11BR053" "X11BR030" "X01BR026" "X03BR006" ...
    #>  $ SRPK1   : chr [1:5] "X11BR003" "X20BR007" "X05BR042" "X05BR029" ...
    #>  $ STK39   : chr [1:6] "X11BR025" "X01BR032" "X03BR010" "X05BR009" ...
    #>  $ SYK     : chr [1:5] "X11BR049" "X01BR017" "X11BR050" "X11BR038" ...
    #>  $ TBK1    : chr [1:2] "X11BR047" "X11BR038"
    #>  $ TEC     : chr [1:2] "X11BR006" "X01BR026"
    #>  $ TK1     : chr [1:5] "X01BR017" "X01BR027" "X03BR004" "X11BR022" ...
    #>  $ ULK1    : chr [1:4] "X11BR025" "X11BR054" "X05BR029" "X13BR009"

Additionally we developed an alternative gold standard set of kinases
using the activating sites on kinases. Hereby, the same thresholds were
applied to the kinase phosphosite levels. This gold standard set can
also be accessed and used in the benchmark (cptacGS\_act\_2pt5,
cptacGS\_act\_5, cptacGS\_act\_10 and cptacGS\_act\_15).

### Area under the receiver operator curve

The inferred kinase activities can then be evaluated using the
benchmarkROC() function. Hereby, it is important to note that the kinase
symbol needs to match with the gold standard set (here gene symbols).

    performance <- benchmarkROC(score_lists = act_scores, GS_list = cptacGS_5)

We can then plot the areas under the receiver operator curve (AUROCs)
and can compare this with other methods or prior knowledge resources.

    auroc_p <- performance %>%
      ggplot2::ggplot(ggplot2::aes(x=method, y=auroc)) +
      ggplot2::geom_boxplot(outlier.size=1, lwd=0.6) +
      ggplot2::geom_hline(yintercept = 0.5)
    auroc_p

![](/private/var/folders/th/nbdnn8l96tx88tt8nm212dpw0000gn/T/RtmpHj2LOB/preview-3b944d36a746.dir/tumorBench_files/figure-markdown_strict/plot-1.png)
