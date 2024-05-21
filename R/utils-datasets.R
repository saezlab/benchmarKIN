#' load_perturbData
#'
#' This function returns a matrix of logFC of phosphorylation sites
#' from perturbation experiments. These can be used to infer kinase
#' activities.
#'
#' @return data.frame with phosphorylation sites as rownames and perturbation
#' experiments as columns.
#' @export
#'
load_perturbData <- function(){
  hernandez <- hernandezData
  hijazi <- hijaziData

  mat <- dplyr::full_join(hernandezData, hijaziData, by = "ID")

  rownames(mat) <- mat$ID
  mat <- mat[, -which(names(mat) == "ID")]

  return(mat)
}


#' load_meta
#'
#' This function  extract regulatory phosphites from the phonemes network, e.i.
#' psites that are differentially regulated on proteins that are found in
#' the phonemes network

#' @return list with two elements: the regulatory psites with predicted mode of
#' regulations and the kinases that catalyse their phosphorilation
#' @export
#'
load_meta <- function(){
  hernandez <- base::as.data.frame(hernandezMeta)
  hijazi <- base::as.data.frame(hijaziMeta)

  meta <- dplyr::bind_rows(hernandez, hijazi)

  return(meta)
}
