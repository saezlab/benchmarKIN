#' run_rank
#'
#' @description
#' Calculates the rank of the perturbed kinases in the respective experiment
#'
#' @param act Activity matrix with kinases as rows and experiments as columns.
#' @param meta Data frame containing sample (experiment) information and perturb (target) in each experiment.
#'
#' @return Data frame containing the rank of each perturbed kinase based on its activity.
#' @export
#' @import tibble tidyr dplyr stringr
#'
#' @examples

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
