#' RoKAI_zscore
#'
#' @description
#' Calculates Z-score for kinase activity  described in RoKAI.
#'
#' @details
#' RoKAI_zscore infers regulator activities
#'
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `source`: Source nodes of `network`.
#'  2. `condition`: Condition representing each column of `mat`.
#'  3. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_zscore_RoKAI(mat, net, minsize=0)
run_zscore <- function(mat,
                       network,
                       minsize = 5
) {
  network_filtered <- network %>%
    dplyr::filter(target %in% rownames(mat)) %>%
    dplyr::group_by(source) %>%
    dplyr::filter(dplyr::n() >= minsize)

  # Analysis ----------------------------------------------------------------
  kin_sub <- network_filtered %>%
    tidyr::pivot_wider(values_from = mor, names_from = target, values_fill = 0) %>%
    tibble::column_to_rownames("source")

  scores <- purrr::map_dfr(1:ncol(mat), function(i_mat){
    V <- mat[i_mat] %>%
      tidyr::drop_na()

    S <- stats::sd(V[,1])

    valid_pps <- base::intersect(rownames(V), colnames(kin_sub))

    kin_sub_f <- kin_sub[, valid_pps]
    V_f <- V %>%
      dplyr::filter(rownames(V) %in% valid_pps)
    V_f <- V_f[valid_pps,] %>%
      base::as.matrix()

    kinaseScores <- (base::as.matrix(kin_sub_f) %*% V_f) / (S * base::sqrt(base::abs(base::as.matrix(kin_sub_f)) %*% base::rep(1, base::length(V_f))))
    kinaseScores <- kinaseScores[!is.na(kinaseScores),]

    data.frame(source = names(kinaseScores), condition = colnames(V), score = base::as.numeric(base::as.vector(kinaseScores)))
  })

  act <- scores %>%
    tidyr::pivot_wider(names_from = "condition", values_from = "score") %>%
    tibble::column_to_rownames("source") %>%
    base::as.data.frame()

  return(act)
}
