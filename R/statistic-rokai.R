#' run_zscore
#'
#' @description
#' Calculates kinase activities from phosphoproteomics data
#' using the Z-score as implemented in RoKAI.
#'
#' @details
#' RoKAI zscore algorithm re-implementation in R
#'
#' @param mat Matrix with phosphorylation sites as rows and experiments as columns.
#' @param network Data frame containing kinase-substrate interactions.
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @return A long format tibble of the activity scores for each source (kinase)
#'  across the samples.
#' @export
#' @import dplyr purrr tibble tidyr
#' @importFrom stats sd
#' @examples
#' # Create random network and matrix
#' set.seed(123)
#' net <- data.frame(source = rep(c("A", "B", "C"), each = 5),
#'                   target = sample(rep(c("A", "B", "C", "D", "E"), each = 3)),
#'                   mor = 1)
#' net <- dplyr::distinct(net)
#'
#' mat <- data.frame(exp1 = runif(5, min = -2, max = 2))
#' rownames(mat) <- c("A", "B", "C", "D", "E")
#'
#' # Activity estimation
#' run_zscore(mat = mat, network = net, minsize = 2)

run_zscore <- function(mat,
                       network,
                       minsize = 5
) {
  # Analysis ----------------------------------------------------------------
  scores <- purrr::map_dfr(1:ncol(mat), function(i_mat){
    V <- mat[i_mat] %>%
      tidyr::drop_na()

    network_filtered <- network %>%
      dplyr::filter(target %in% rownames(V)) %>%
      dplyr::group_by(source) %>%
      dplyr::filter(dplyr::n() >= minsize)

    kin_sub <- network_filtered %>%
      tidyr::pivot_wider(values_from = mor, names_from = target, values_fill = 0) %>%
      tibble::column_to_rownames("source")

    S <- stats::sd(V[,1])

    valid_pps <- base::intersect(rownames(V), colnames(kin_sub))

    #Prevent coercion to 1-dimensional object in case there is only one result
    kin_sub_f <- kin_sub[, valid_pps, drop=FALSE]
    V_f <- V %>%
      dplyr::filter(rownames(V) %in% valid_pps)
    V_f <- V_f[valid_pps,] %>%
      base::as.matrix()

    kinaseScores <- (base::as.matrix(kin_sub_f) %*% V_f) / (S * base::sqrt(base::abs(base::as.matrix(kin_sub_f)) %*% base::rep(1, base::length(V_f))))
    kinaseScores <- kinaseScores[!is.na(kinaseScores),]

    #Prevent error if a single experiment has no overlap with the network
    if (length(kinaseScores) > 0) {
      data.frame(source = names(kinaseScores), condition = colnames(V), score = base::as.numeric(base::as.vector(kinaseScores)))
    }
  })
  #Prevent error if NO experiment had any overlap with the network
  if (length(scores) > 0) {
    act <- scores %>%
      tidyr::pivot_wider(names_from = "condition", values_from = "score") %>%
      tibble::column_to_rownames("source") %>%
      base::as.data.frame()

    return(act)
  }
}
