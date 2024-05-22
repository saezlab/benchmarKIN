#' run_perturbBench
#'
#' @description
#' Runs perturbation benchmark (based on decoupler).
#'
#' @param act Activity matrix with kinases as rows and experiments as columns.
#' @param meta Data frame containing sample (experiment) information and perturb (target) in each experiment.
#' @param scale_data Whether scores should be scaled per experiment.
#' @param rm_bg Whether to remove experiments where no activity was inferred for the target kinase from the background set.
#' @param n_iter Number of iterations for the subsampling
#' @param method_id ID for the method currently being tested.
#' @param metric Performance metric to be returned (either auroc or auprc)
#'
#' @import reticulate dplyr
#'
#' @return Data frame containing the performance metric for the inferred kinase activities.
#' @export
#' @examples
#' # Create random meta and matrix
#' set.seed(123)
#'
#' mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
#'                   exp2 = runif(5, min = -3, max = 2))
#' rownames(mat) <- c("A", "B", "C", "D", "E")
#'
#' meta <- data.frame(id = c("exp1", "exp2"),
#'                    target = c("E", "A"),
#'                    sign = c(1, -1))
#'
#' # run benchmark
#' res <- run_perturbBench(act = mat, meta = meta, method_id = "test")

run_perturbBench <- function(act, meta, scale_data = T, rm_bg = F, n_iter = 1000, method_id = "method", metric = "auroc", ...){
  if(scale_data){
    act <- scale_scores(act)
  }

  res <- prepareBench(act = act, meta = meta, rm_bg = rm_bg, method_id = method_id)
  act <- res$act
  obs <- res$obs

  act_py <- reticulate::r_to_py(act)
  obs_py <- reticulate::r_to_py(obs)

  # Import the decoupler package
  dc <- reticulate::import("decoupler")

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
  performances_r <- reticulate::py_to_r(performances) %>%
    dplyr::select(method, metric, score)

  if (metric == "auroc"){
    performances_r <- performances_r %>%
      dplyr::filter(metric == "mcauroc") %>%
      dplyr::select(method, score) %>%
      dplyr::rename("auroc" = score)
  } else if (metric == "auprc"){
    performances_r <- performances_r %>%
      dplyr::filter(metric == "mcauprc") %>%
      dplyr::select(method, score) %>%
      dplyr::rename("auprc" = score)
  }

  return(performances_r)
}

#' prepareBench
#'
#' @description
#' Prepares the input for the decoupler benchmark function.
#'
#' @param act Activity matrix with kinases as rows and experiments as columns.
#' @param meta Data frame containing sample (experiment) information and perturb (target) in each experiment.
#' @param rm_bg Whether to remove experiments where no activity was inferred for the target kinase from the background set.
#' @param method_id ID for the method currently being tested.
#'
#' @import dplyr purrr tibble
#'
#' @return list of processed activity scores and meta data which can be used as direct input to the decoupler benchmark function
#' @export
#' @examples
#' set.seed(123)
#'
#' mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
#'                   exp2 = runif(5, min = -3, max = 2))
#' rownames(mat) <- c("A", "B", "C", "D", "E")
#'
#' meta <- data.frame(id = c("exp1", "exp2"),
#'                    target = c("A", "C"),
#'                    sign = c(1, -1))
#'
#' out_list <- prepareBench(act = mat, meta = meta, method_id = "test")

prepareBench <- function(act, meta, rm_bg = F, method_id = "method"){
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
  # or whether they should be kept as background
  if (rm_bg){
    df <- remove_bg(mat = df, meta = target_df)
  }

  rownames(df) <- NULL
  df <- df %>%
    tibble::column_to_rownames("experiment")

  # Replace 0 with NA
  df[df == 0] <- NA

  # filter out meta to fit to experiments in matrix
  target_df <- target_df %>%
    dplyr::filter(sample %in% rownames(df)) %>%
    dplyr::mutate(sign = 1) %>%
    dplyr::mutate(perturb = stringr::str_split(perturb, ";")) %>%
    tibble::column_to_rownames("sample")

  # Create a list with the method name and scores for the input
  res <- list(df)
  names(res) <- method_id

  return(list(act = res, obs = target_df))
}

#' scale_scores
#'
#' @description
#' Scales inferred kinase activity scores per experiment. Standardizes standard
#' deviation or maximum score across experiments.
#'
#' @param mat Matrix with kinases as rows and experiments as columns.
#' @param scaling Scaling based on max or standard deviation.
#'
#' @return Scaled matrix with kinases as rows and experiments as columns.
#' @export
#' @examples
#' # Create random network and matrix
#' set.seed(123)
#'
#' mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
#'                   exp2 = runif(5, min = -3, max = 2))
#' rownames(mat) <- c("A", "B", "C", "D", "E")
#'
#' # Scaling
#' scale_scores(mat = mat)

scale_scores <- function(mat, scaling = "sd"){
  if (scaling == "max"){
    map_dfc(colnames(mat), function(col_i){
      mat[col_i]/max(abs(mat[,col_i]), na.rm = T)
    })
  } else if (scaling == "sd") {
    base::scale(mat, center = FALSE, scale = TRUE)[,]
  } else {
    mat
  }
}

#' remove_bg
#'
#' @description
#' Identifies experiments where no activity could be inferred for any of the
#' target kinases and removes the experiment.
#'
#' @param mat Matrix with experiments as rows and kinases as columns.
#' @param meta Data frame containing sample (experiment) information and perturb (target) in each experiment
#'
#' @import dplyr purrr
#' @importFrom stringr str_split
#'
#' @return Matrix after filtering out experiments where no activity for the target was inferred.
#' @export
#' @examples
#'  Create random network and matrix
#'  set.seed(123)
#'
#'  mat <- data.frame(experiment = c("exp1", "exp2"),
#'                    A = runif(2, min = -2, max = 2),
#'                    B = runif(2, min = -3, max = 2),
#'                    C = runif(2, min = -2, max = 2),
#'                    D = runif(2, min = -3, max = 2))
#'
#'  meta <- data.frame(sample = c("exp1", "exp2"),
#'                     perturb = c("A", "F"))
#'
#' # Remove experiment from background
#' mat_filtered <- remove_bg(mat = mat, meta = meta)

remove_bg <- function(mat, meta){
  # mark experiments where no activity was inferred for target
  df_filtered <- purrr::map_dfr(1:nrow(mat), function(i){
    tmp <- mat[i,]
    targets <- meta %>%
      dplyr::filter(sample %in% tmp$experiment) %>%
      dplyr::pull(perturb) %>%
      stringr::str_split(";") %>%
      base::unlist()

    sum_act <- base::sum(tmp[,colnames(tmp) %in% targets], na.rm = T)

    if(!is.na(sum_act) & !sum_act == 0){
      mat[i,]
    } else {
      mat[i,] %>% dplyr::mutate(experiment = "remove")
    }
  })

  # Filter out the complete experiment
  df <- df_filtered %>%
    dplyr::filter(!experiment == "remove")

  return(df)
}
