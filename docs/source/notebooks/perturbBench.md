# Peturbation-based benchmark

## Introduction

This benchmark is based on experiments where specific kinases are
expected to have an increased or decreased activity after being
perturbed in various cell line experiments.

We provide the log fold changes of phosphorylation sites after the
perturbation which can be used to infer changes in kinase activities
using different methods or kinase-substrate libraries.

The methods are then evaluated based on whether they are able to
recapitulate the peturbed kinases based on their inferred activities.

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

Since we call the get\_performances function from the decoupler-py
package we also need to set up reticulate.

    library(reticulate)

    # Use reticulate to install the Python package
    py_install("git+https://github.com/saezlab/decoupler-py.git", pip = TRUE)

    # If you are using a conda environment make sure to use the correct path for your
    # python version and install the decoupler-py package directly into your environment

    # use_python("path-to-python", required = TRUE)  #you can find the correct path by activating your environment and typing "which python" into the terminal
    # py_install("git+https://github.com/saezlab/decoupler-py.git", pip = TRUE, envname = "path-to-environment")

## Perturbation data

The perturbation data can be extracted using the load\_perturbData()
function. This data frame contains the logFC for the perturbation
collection previously published by Hernandez-Armenta et al. and the
perturbation experiments from Hijazi et al.

    mat <- load_perturbData()
    mat[1:5, 20:23]
    #>                            273_75      276_75     279_75     282_75
    #> ARF5_S103|ARF5|S103            NA  0.04380698         NA         NA
    #> M6PR_S267|M6PR|S267   -0.05630667 -0.01278283 0.01527238 -0.3346162
    #> FKBP4_S118|FKBP4|S118          NA  0.45777031         NA  0.1405713
    #> FKBP4_S15|FKBP4|S15            NA          NA         NA         NA
    #> FKBP4_S258|FKBP4|S258          NA          NA         NA         NA

Additional information about each experiment can be found in the meta
data. This contains information about the targeted kinase(s), the PMID,
etc.

    meta <- load_meta()
    head(meta)
    #>        id                     Description              Control Time_min
    #> 1 102_117                          RG7356              Vehicle       30
    #> 2 105_117                          RG7356              Vehicle       90
    #> 3  111_45 BI 4834 (on Mitosis Nocodazole) Mitosis (Nocodazole)       15
    #> 4  114_48                         Starved              Insulin       NA
    #> 5  117_48                          Torin1              Insulin       NA
    #> 6  120_48                       Rapamycin              Insulin       NA
    #>       PMID target sign            class cell_line treatment
    #> 1 22777824   AKT1    1 Serine/Threonine      <NA>      <NA>
    #> 2 22777824   AKT1    1 Serine/Threonine      <NA>      <NA>
    #> 3 21857030   PLK1   -1 Serine/Threonine      <NA>      <NA>
    #> 4 21659604   MTOR   -1 Serine/Threonine      <NA>      <NA>
    #> 5 21659604   MTOR   -1 Serine/Threonine      <NA>      <NA>
    #> 6 21659604   MTOR   -1 Serine/Threonine      <NA>      <NA>

## Kinase activity inference

This data can be used to test any method for kinase activity inference.

We will use the z-score (as implemented by RoKAI) and the curated
kinase- substrate library as in example.

### Kinase-substrate library

For that we have already processed a version of a combination of curated
libraries, namely PhosphoSitePlus, PTMsigDB (excluding iKiP-DB) and the
gold standard set of GPS 5.0, that can be mapped to our phosphorylation
site ids.

    head(curated)
    #>    source      target target_protein position mor        sequence
    #> 1 EIF2AK1  EIF2S1_S52         EIF2S1      S52   1 MILLSELSRRRIRSI
    #> 2 EIF2AK1  EIF2S1_S49         EIF2S1      S49   1 IEGMILLSELSRRRI
    #> 3   PRKCD  HDAC5_S259          HDAC5     S259   1 FPLRKTASEPNLKVR
    #> 4   PRKCD  PTPRA_S204          PTPRA     S204   1 PLLARSPSTNRKYPP
    #> 5   PRKCD    BCL2_S70           BCL2      S70   1 RDPVARTSPLQTPAA
    #> 6   PRKCD HNRNPK_S302         HNRNPK     S302   1 GRGGRGGSRARNLPL
    #>           ENSEMBL
    #> 1 ENSG00000134001
    #> 2 ENSG00000134001
    #> 3 ENSG00000108840
    #> 4 ENSG00000132670
    #> 5 ENSG00000171791
    #> 6 ENSG00000165119

After mapping the phosphorylation sites we can bring it into the right
format to run the run\_zscore function.

    curated$target <- paste(curated$target, curated$target_protein, curated$position, sep = "|")

    curatedLib <- curated %>%
      dplyr::select(source, target, mor) %>%
      dplyr::distinct()

    head(curatedLib)
    #>    source                  target mor
    #> 1 EIF2AK1   EIF2S1_S52|EIF2S1|S52   1
    #> 2 EIF2AK1   EIF2S1_S49|EIF2S1|S49   1
    #> 3   PRKCD   HDAC5_S259|HDAC5|S259   1
    #> 4   PRKCD   PTPRA_S204|PTPRA|S204   1
    #> 5   PRKCD       BCL2_S70|BCL2|S70   1
    #> 6   PRKCD HNRNPK_S302|HNRNPK|S302   1

### Activity inference using the z-score

The z-score calculates a score for each kinase by aggregating the change
in abundance of the direct targets in relation to changes in the
non-targets.

    act_scores <- run_zscore(mat = mat, network = curatedLib)

    act_scores[1:5, 1:5]
    #>               15_3       18_3        21_3       24_3      96_36
    #> PRKCD    3.1802463  2.3143338  3.42871525  1.2961965 -4.2716360
    #> CAMK2A   0.6594898  0.7286016  0.07022767  1.9630066         NA
    #> CSNK2A1  1.2340623  1.2400087  1.38745758  2.1008802 -0.1351479
    #> LATS1   -1.6980139 -2.0950547 -0.86222145  0.3599908         NA
    #> CDK7     0.1237668 -0.6063455 -0.92708849 -2.0542065  0.4655647

## Perturbation benchmark

### Area under the receiver operator curve

The inferred kinase activities can now be evaluated using the
run\_perturbBench() function. Hereby, it is important to note that the
kinase symbol needs to match with the meta data (here gene symbols).

    performance <- run_perturbBench(act = act_scores, meta = meta, method_id = "zscore_ppsp")

We can then plot the areas under the receiver operator curve (AUROCs)
and can compare this with other methods or prior knowledge resources.

    auroc_p <- performance %>%
      ggplot2::ggplot(ggplot2::aes(x=method, y=auroc)) +
      ggplot2::geom_boxplot(outlier.size=1, lwd=0.6) +
      ggplot2::geom_hline(yintercept = 0.5)

    auroc_p

![](/private/var/folders/th/nbdnn8l96tx88tt8nm212dpw0000gn/T/RtmpFl3WvE/preview-3e9d2799757d.dir/perturbBench_files/figure-markdown_strict/plot-1.png)

### Area under the receiver operator curve

Additionally we can calculate the rank and scaled rank of the perturbed
kinase(s) in their respective experiments. A lower rank/scaled rank
indicates a better performance of the method

    rank <- run_rank(act = act_scores, meta = meta)
    median(rank$scaled_rank)
    #> [1] 0.2354167
