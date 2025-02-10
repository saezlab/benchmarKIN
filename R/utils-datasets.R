#' load_perturbData
#'
#' @description
#' Loads data from the perturbation collection
#' published by Hernandez-Armenta with the Hijazi datasets.
#'
#' @return data frame containing logFC of different
#' perturbation experiments
#' @param dataset whether to load the "hijazi", "hernandez", "tyrosine" or "all" datasets
#' @export
#' @importFrom dplyr full_join
#' @examples
#' load_perturbData()
load_perturbData <- function(dataset = "all"){
  hernandez <- hernandezData
  hijazi <- hijaziData
  tyrosine <- tyrosineData

  mat <- dplyr::full_join(hernandezData, hijaziData, by = "ID")
  mat <- dplyr::full_join(mat, tyrosineData, by = "ID")

  rownames(mat) <- mat$ID
  mat <- mat[, -which(names(mat) == "ID")]
  colnames(mat) <- sub("^X", "", colnames(mat))

  if (dataset == "hernandez"){
    rownames(hernandez) <- hernandez$ID
    hernandez <- hernandez[, -which(names(hernandez) == "ID")]
    return(hernandez)
  } else if (dataset == "hijazi"){
    rownames(hijazi) <- hijazi$ID
    hijazi <- hijazi[, -which(names(hijazi) == "ID")]
    return(hijazi)
  } else if (dataset == "tyrosine"){
    rownames(tyrosine) <- tyrosine$ID
    tyrosine <- tyrosine[, -which(names(tyrosine) == "ID")]
    return(hijazi)
  } else if (dataset == "all"){
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
#' @param dataset whether to load the "hijazi", "hernandez", "tyrosine" or "all" datasets
#' @param targets_hijazi whether to use a "manual" list of targets for the
#' inhibitors or whether to consider all kinases targets based on the "selectivity"
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom tidyr separate_rows
#' @examples
#' load_meta()
load_meta <- function(dataset = "all", targets_hijazi = "manual"){
  hernandez <- base::as.data.frame(hernandezMeta)
  hijazi <- base::as.data.frame(hijaziMeta)
  tyrosine <- base::as.data.frame(tyrosineMeta)
  tyrosine$PMID <- as.numeric(tyrosine$PMID)

  if (targets_hijazi == "manual"){
    hijazi <- hijazi[!colnames(hijazi) %in% c("target_discoverX", "class_discoverX")]
  } else if (targets_hijazi == "selectivity"){
    hijazi <- hijazi[!colnames(hijazi) %in% c("target", "class")]
    colnames(hijazi) <- c("id", "cell_line", "treatment", "PMID", "sign", "target", "class")
  }

  hijazi <- hijazi[c("id", "cell_line", "treatment", "PMID", "target",  "sign", "class")]
  hijazi <- tidyr::separate_rows(hijazi, target, class, sep = ";")
  hijazi <- hijazi[!is.na(hijazi$target),]

  meta <- dplyr::bind_rows(hernandez, hijazi, tyrosine)

  if (dataset == "hernandez"){
    return(hernandez)
  } else if (dataset == "hijazi"){
    return(hijazi)
  } else if (dataset == "tyrosine"){
    return(tyrosine)
  } else if (dataset == "all"){
    return(meta)
  }
}
