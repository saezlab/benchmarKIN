#' load_perturbData
#'
#' @description
#' Loads data from the perturbation collection
#' published by Hernandez-Armenta with the Hijazi datasets.
#'
#' @return data frame containing logFC of different
#' perturbation experiments
#' @param dataset whether to load the "hijazi", "hernandez" or "both" datasets
#' @export
#' @importFrom dplyr full_join
#' @examples
#' load_perturbData()
load_perturbData <- function(dataset = "both"){
  hernandez <- hernandezData
  hijazi <- hijaziData

  mat <- dplyr::full_join(hernandezData, hijaziData, by = "ID")

  rownames(mat) <- mat$ID
  mat <- mat[, -which(names(mat) == "ID")]

  if (dataset == "hernandez"){
    rownames(hernandez) <- hernandez$ID
    hernandez <- hernandez[, -which(names(hernandez) == "ID")]
    return(hernandez)
  } else if (dataset == "hijazi"){
    rownames(hijazi) <- hijazi$ID
    hijazi <- hijazi[, -which(names(hijazi) == "ID")]
    return(hijazi)
  } else if (dataset == "both"){
    return(mat)
  }
}

#' load_meta
#'
#' @description
#' Loads meta data from the perturbation collection
#' published by Hernandez-Armenta with the Hijazi datasets.
#'
#' @return data frame containing additional information
#' about the perturbation experiments
#' @param dataset whether to load the "hijazi", "hernandez" or "both" datasets
#' @param targets_hijazi whether to use a "manual" list of targets for the
#' inhibitors or whether to consider all kinases targets based on the "selectivity"
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom tidyr separate_rows
#' @examples
#' load_meta()
load_meta <- function(dataset = "both", targets_hijazi = "manual"){
  hernandez <- base::as.data.frame(hernandezMeta)
  hijazi <- base::as.data.frame(hijaziMeta)

  if (targets_hijazi == "manual"){
    hijazi <- hijazi[!colnames(hijazi) == "target_discoverX"]
  } else if (targets_hijazi == "selectivity"){
    hijazi <- hijazi[!colnames(hijazi) == "target"]
    colnames(hijazi) <- c("id", "cell_line", "drug", "sign", "target", "PMID")
  }

  hijazi <- tidyr::separate_rows(hijazi, target, sep = ";")
  hijazi <- hijazi[!is.na(hijazi$target),]

  meta <- dplyr::bind_rows(hernandez, hijazi)

  if (dataset == "hernandez"){
    return(hernandez)
  } else if (dataset == "hijazi"){
    return(hijazi)
  } else if (dataset == "both"){
    return(meta)
  }
}
