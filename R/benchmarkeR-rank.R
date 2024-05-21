#' run_rank
#'
#' This function returns a matrix of logFC of phosphorylation sites
#' from perturbation experiments. These can be used to infer kinase
#' activities.
#'
#' @return data.frame with phosphorylation sites as rownames and perturbation
#' experiments as columns.
#' @export
#'
run_rank <- function(act, meta){
  ## Get rank ---------------------------
  method_act_long <- act %>%
    tibble::rownames_to_column("kinase") %>%
    tidyr::pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
    dplyr::filter(!is.na(score))

  # Rename columns in meta and remove experiments with unknown targets (e.g. several members)
  obs <- meta %>%
    dplyr::rename("perturb" = target) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(targets = paste(perturb, collapse = ";"), sign = unique(sign)) %>%
    dplyr::rename("perturb" = targets) %>%
    dplyr::filter(!perturb == "" | !is.na(perturb)) %>%
    dplyr::rename("sample" = id)

  method_act_long <- method_act_long %>%
    dplyr::filter(sample %in% obs$sample)

  ## Get rank
  rank_df <- purrr::map_dfr(unique(method_act_long$sample), function(exp){
    direction <- obs %>%
      dplyr::filter(sample == exp) %>%
      dplyr::pull(sign) %>%
      base::unique()

    act_df <- method_act_long %>%
      dplyr::filter(sample == exp)

    if (direction == 1){
      act_df <- act_df %>%
        dplyr::arrange(desc(score))
    } else if (direction == -1){
      act_df <- act_df %>%
        dplyr::arrange(score)
    }


    if (exp %in% obs$sample){
      targets <- obs %>%
        dplyr::filter(sample == exp) %>%
        dplyr::pull(perturb) %>%
        stringr::str_split(";") %>%
        base::unlist()

      rank <- purrr::map_dbl(targets, function(target){
        position <- base::which(act_df$kinase %in% target)
        if (base::length(position) == 0){
          position <- NA
        }
        position
      })

      data.frame(sample = exp,
                 targets = targets,
                 rank = rank,
                 kinases_act = nrow(act_df),
                 all_kinases_act = paste(act_df %>%
                                           dplyr::pull(kinase), collapse = ";")) %>%
        dplyr::mutate(scaled_rank = rank/kinases_act)
    }
  }) %>%
    dplyr::filter(!is.na(rank))

  return(rank_df)
}
