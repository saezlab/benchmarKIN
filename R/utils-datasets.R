#' load_perturbData
#'
#' @description
#' Loads data from the perturbation collection
#' published by Hernandez-Armenta with the Hijazi datasets.
#'
#' @return data frame containing logFC of different
#' perturbation experiments
#' @export
#' @importFrom dplyr full_join
#' @examples
#' load_perturbData()
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
#' @description
#' Loads meta data from the perturbation collection
#' published by Hernandez-Armenta with the Hijazi datasets.
#'
#' @return data frame containing additional information
#' about the perturbation experiments
#' @export
#' @importFrom dplyr bind_rows
#' @examples
#' load_meta()
load_meta <- function(){
  hernandez <- base::as.data.frame(hernandezMeta)
  hijazi <- base::as.data.frame(hijaziMeta)

  meta <- dplyr::bind_rows(hernandez, hijazi)

  return(meta)
}
