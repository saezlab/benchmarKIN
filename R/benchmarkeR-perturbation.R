#' run_perturbBench
#'
#' This function returns a matrix of logFC of phosphorylation sites
#' from perturbation experiments. These can be used to infer kinase
#' activities.
#'
#' @return data.frame with phosphorylation sites as rownames and perturbation
#' experiments as columns.
#' @export
#'
run_perturbBench <- function(act, meta, scale_data = T, rm_experiments = F, n_iter = 1000, ...){
  if(scale_data){
    act <- scale_scores(act)
  }

  res <- prepareBench(act = act, meta = meta, rm_experiments = rm_experiments)
  act <- res$act
  obs <- res$obs

  act_py <- reticulate::r_to_py(act)
  obs_py <- reticulate::r_to_py(obs)

  # Call the Python function
  performances <- dc$get_performances(
    res = act_py,
    obs = obs_py,
    groupby = NULL,
    by = 'experiment',
    metrics = c('auroc', 'auprc', 'mcauroc', 'mcauprc'),
    n_iter = n_iter,
    min_exp = 1
  )

  # Convert the result back to R if needed
  performances_r <- reticulate::py_to_r(performances)

  return(performances_r)
}

#' prepareBench
#'
#' This function returns a matrix of logFC of phosphorylation sites
#' from perturbation experiments. These can be used to infer kinase
#' activities.
#'
#' @return data.frame with phosphorylation sites as rownames and perturbation
#' experiments as columns.
#' @export
#'
prepareBench <- function(act, meta, rm_experiments = F, method_id = "method"){
  # Rename columns in meta and remove experiments with unknown targets (e.g. several members)
  target_df <- meta %>%
    dplyr::rename("perturb" = target) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(targets = paste(perturb, collapse = ";"), sign = unique(sign)) %>%
    dplyr::rename("perturb" = targets) %>%
    dplyr::filter(!perturb == "" | !is.na(perturb)) %>%
    dplyr::rename("sample" = id)

  # change direction according to perturbation
  df <- purrr::map_dfr(target_df$sample, function(experiment){
    if (experiment %in% colnames(act)){
      mat <- t(act[,experiment]) * target_df$sign[target_df$sample == experiment]
      data.frame(mat) %>%
        tibble::add_column(experiment = experiment, .before = 1)
    } else {
      df <- data.frame(matrix(NA, ncol = nrow(mat) + 1))
      colnames(df) <- c("experiment", rownames(mat))
      df
    }
  }) %>% dplyr::filter(!is.na(experiment))

  colnames(df) <- c("experiment", rownames(act))
  df <- df[df$experiment %in% target_df$sample,]

  # filter out experiments where no activity was inferred for perturbed kinase
  # or should they be kept as background
  if (rm_experiments){
    df_filtered <- purrr::map_dfr(1:nrow(df), function(i){
      tmp <- df[i,]
      targets <- target_df %>%
        dplyr::filter(sample %in% tmp$experiment) %>%
        dplyr::pull(perturb) %>%
        stringr::str_split(";") %>%
        base::unlist()

      sum_act <- base::sum(tmp[,colnames(tmp) %in% targets], na.rm = T)

      if(!is.na(sum_act) & !sum_act == 0){
        df[i,]
      } else {
        df[i,] %>% dplyr::mutate(experiment = "remove")
      }
    })

    df <- df_filtered %>%
      dplyr::filter(!experiment == "remove")
  }

  df <- df %>%
    tibble::column_to_rownames("experiment")

  # filter out meta to fit to experiments in matrix
  target_df <- target_df %>%
    dplyr::filter(sample %in% rownames(df)) %>%
    dplyr::mutate(sign = 1)

  # Replace 0 with NA
  df[df == 0] <- NA

  # Create a list with the method name and scores for the input
  res <- list(df)
  names(res) <- method_id

  # Split the 'perturb' column by ";"
  target_df <- target_df %>%
    dplyr::mutate(perturb = stringr::str_split(perturb, ";")) %>%
    tibble::column_to_rownames("sample")

  return(list(act = res, obs = target_df))
}

#' scale_scores
#'
#' This function returns a matrix of logFC of phosphorylation sites
#' from perturbation experiments. These can be used to infer kinase
#' activities.
#'
#' @return data.frame with phosphorylation sites as rownames and perturbation
#' experiments as columns.
#' @export
#'
scale_scores <- function(mat, scaling = "sd"){
  if (scaling == "max"){
    map_dfc(colnames(mat), function(col_i){
      mat[col_i]/max(abs(mat[,col_i]), na.rm = T)
    })
  } else if (scaling == "sd") {
    scale(mat, center = FALSE, scale = TRUE)[,]
  } else {
    mat
  }
}




