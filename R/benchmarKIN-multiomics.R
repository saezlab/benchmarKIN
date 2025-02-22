#' benchmarkROC
#'
#' @description
#' Evaluates inferred kinase activities from the cptac data. Gold standard sets have
#' been defined based on the proteomics data of cptac.
#'
#' @param score_lists A list of lists containing kinase scores for each method for each of cancer types.
#' @param GS_list A list of data frames containing gold standard (GS) positive and negative pairs of kinase for each cancer type.
#' @param cts A character vector specifying the cancer types or experimental conditions to be considered. Default is "all", indicating all cancer types in score_lists.
#' @param min_num An integer specifying the minimum number of measurements required for each kinase. Default is 30.
#' @param kins A list of potential kinases to consider, which should match the kinase names in the GS lists. Default is "all", indicating all kinases present in the GS lists.
#' @param sub_samp_frac A numeric value specifying the fraction of samples to be used for benchmarking. Default is 0.8.
#' @param ROC_obj A logical value indicating whether to return ROC objects along with AUROC values. Default is FALSE.
#' @param metric A character value specifying the evaluation metric.
#' @param ... Additional arguments to be passed to the internal functions, such as the `.ROC_sampler` function.
#'
#' @return Data frame containing the performance metric for the inferred kinase activities.
#'
#' @export
#'
benchmarkROC <- function(score_lists, GS_list, cts="all", min_num=30, kins="all", sub_samp_frac=0.8, ROC_obj=F, metric = "auroc", ...){
  #filter to cancer types of (cts) interest if specified; otherwise use all cts in score_lists
  if(length(cts)==1){
    if(cts=="all"){
      cts <- names(score_lists)
    }
  }
  score_lists <- score_lists[cts]
  GS_list <- GS_list[cts]
  #set list of potential kinases to consider (need to be in Gold Standard (GS), filter to user list if provided)
  if(length(unlist(kins))==1){
    if(kins=="all"){
      kins <- list()
      for(i in 1:length(cts)){
        kins[[cts[i]]] <- names(GS_list[[cts[i]]]$GS_pos_pairs)
      }
    }
  } else {
    for(i in 1:length(cts)){
      kins[[cts[i]]] <- intersect(kins[[cts[i]]], names(GS_list[[cts[i]]]$GS_pos_pairs))
    }
  }
  kins <- kins[cts]
  #convert subset to scores for kinases of relevance, filter out kinase scores without enough measurements (min_num), convert to Z-scores
  score_listsZ <- score_lists
  for(i in 1:length(score_listsZ)){
    for(j in 1:length(score_listsZ[[i]])){
      score_listsZ[[i]][[j]] <- score_listsZ[[i]][[j]][rowSums(!is.na(score_listsZ[[i]][[j]])) >= min_num, ]
      kins[[names(score_listsZ)[i]]] <- intersect(kins[[names(score_listsZ)[i]]], rownames(score_listsZ[[i]][[j]]))
    }
    for(j in 1:length(score_listsZ[[i]])){
      score_listsZ[[i]][[j]] <- score_listsZ[[i]][[j]][kins[[names(score_listsZ)[i]]], ]
    }
    score_listsZ[[i]] <- lapply(score_listsZ[[i]], function(x) as.data.frame(t(scale(t(x)))))
  }
  #get Z-scores for pairs in GS sets for each ct
  GSvals <- list()
  for(i in 1:length(cts)){
    GSvals[[cts[i]]] <- lapply(score_listsZ[[cts[i]]], .getGSvals, GS_list[[cts[i]]]$GS_pos_pairs, GS_list[[cts[i]]]$GS_neg_pairs)
    #make sure that kinase has GS scores for at least GS positive pair and one GS negative pair for all methods >> update kinase list and use to subset later
    for(j in 1:(length(GSvals[[cts[i]]]))){
      kins[[cts[i]]] <- intersect(kins[[cts[i]]], names(GSvals[[cts[i]]][[j]]$GSpos_values))
    }
  }
  all_score_list_pos <- list()
  all_score_list_neg <- list()
  samps <- list()
  for(i in 1:(length(GSvals[[1]]))){
    for(j in 1:length(GSvals)){
      all_score_list_pos[[names(GSvals[[1]])[i]]] <- c(all_score_list_pos[[names(GSvals[[1]])[i]]], unlist(GSvals[[j]][[i]]$GSpos_values[kins[[j]]]))
      all_score_list_neg[[names(GSvals[[1]])[i]]] <- c(all_score_list_neg[[names(GSvals[[1]])[i]]], unlist(GSvals[[j]][[i]]$GSneg_values[kins[[j]]]))
    }
    #unlist(c(brca_scores_GSvals[[i]]$GSpos_values[brca_kins], ccrcc_scores_GSvals[[i]]$GSpos_values[ccrcc_kins], gbm_scores_GSvals[[i]]$GSpos_values[gbm_kins], hnscc_scores_GSvals[[i]]$GSpos_values[hnscc_kins], lscc_scores_GSvals[[i]]$GSpos_values[lscc_kins], luad_scores_GSvals[[i]]$GSpos_values[luad_kins], ucec_scores_GSvals[[i]]$GSpos_values[ucec_kins]))
    #all_score_list_neg[[names(brca_scores_GSvals)[i]]] <- unlist(c(brca_scores_GSvals[[i]]$GSneg_values[brca_kins], ccrcc_scores_GSvals[[i]]$GSneg_values[ccrcc_kins], gbm_scores_GSvals[[i]]$GSneg_values[gbm_kins], hnscc_scores_GSvals[[i]]$GSneg_values[hnscc_kins], lscc_scores_GSvals[[i]]$GSneg_values[lscc_kins], luad_scores_GSvals[[i]]$GSneg_values[luad_kins], ucec_scores_GSvals[[i]]$GSneg_values[ucec_kins]))
    samps[[names(GSvals[[1]][i])]] <- sub_samp_frac*c(length(all_score_list_pos[[i]]),length(all_score_list_neg[[i]]))
  }
  #samps <- 0.8*c(257,240)
  all_roc_list <- list()
  all_roc <- numeric()
  all_roc_sd <- numeric()
  for(i in 1:length(all_score_list_pos)){
    suppressMessages({
      all_roc_list[[names(all_score_list_pos)[i]]] <- .ROC_sampler(all_score_list_pos[[i]], all_score_list_neg[[i]], set_size = samps[[i]], ...)
    })
    all_roc[i] <- mean(all_roc_list[[names(all_score_list_pos)[i]]]$sample_AUROCs)
    names(all_roc)[i] <- names(all_score_list_pos)[i]
    all_roc_sd[i] <- sd(all_roc_list[[names(all_score_list_pos)[i]]]$sample_AUROCs)
    names(all_roc_sd)[i] <- names(all_score_list_pos)[i]
  }
  all_roc_tab <- cbind(all_roc, all_roc_sd)
  op <- list(kins, GSvals, all_score_list_pos, all_score_list_neg, all_roc_list, all_roc_tab)
  names(op) <- c("evaluation_kinases", "list_of_standard_set_values", "all_positive_pair_values", "all_negative_pair_values", "ROC_results", "ROC_summary_table")

  if (metric == "auroc"){
    op <- purrr::map_dfr(names(op$ROC_results), function(method_id){
      data.frame(method = method_id,
                 auroc = op$ROC_results[[method_id]]$sample_AUROCs)

    })
  }

  return(op)
}

#' .getGSvals
#'
#' @description
#' This function processes a dataframe of kinase scores, optionally scales them,
#' and extracts values based on specified positive and negative gene sets.
#'
#' @param score_df A dataframe of scores with kinases as row names and samples/conditions as column names.
#' @param GSpos_list A named list of positive gene sets for each kinase.
#' @param GSneg_list A named list of negative gene sets for each kinase.
#' @param scale A logical indicating whether to scale the scores by z-scores across each row. Default is FALSE.
#'
#' @return A list with the following elements:
#'   \item{Z-scores}{Dataframe of (optionally scaled) scores.}
#'   \item{GSpos_values}{List of extracted positive gene set values for each kinase.}
#'   \item{GSneg_values}{List of extracted negative gene set values for each kinase.}
#'   \item{number_GSpos_per_kinase}{Numeric vector of counts of positive gene sets per kinase.}
#'   \item{number_GSneg_per_kinase}{Numeric vector of counts of negative gene sets per kinase.}
#'
#' @export
#'
.getGSvals <- function(score_df, GSpos_list, GSneg_list, scale=F){
  if(scale==T){
    scoreZ_df <- t(scale(t(score_df)))
  } else {
    scoreZ_df <- score_df
  }
  pos_vals <- list()
  neg_vals <- list()
  kins <- intersect(rownames(scoreZ_df), names(GSpos_list))
  kins_pos <- numeric(length(kins))
  names(kins_pos) <- kins
  kins_neg <- kins_pos
  for(i in 1:length(kins)){
    to_add <- as.numeric(scoreZ_df[kins[i], GSpos_list[[kins[i]]]])
    if(sum(!is.na(to_add)) > 0){
      names(to_add) <- paste0(kins[i], "_", GSpos_list[[kins[i]]])
    }
    to_add <- to_add[!is.na(to_add)]
    kins_pos[i] <- length(to_add)
    to_add_neg <- as.numeric(scoreZ_df[kins[i], GSneg_list[[kins[i]]]])
    if(sum(!is.na(to_add_neg)) > 0){
      names(to_add_neg) <- paste0(kins[i], "_", GSneg_list[[kins[i]]])
    }
    to_add_neg <- to_add_neg[!is.na(to_add_neg)]
    kins_neg[i] <- length(to_add_neg)
    if(kins_pos[i] > 0 & kins_neg[i] > 0){
      pos_vals[[kins[i]]] <- to_add
      neg_vals[[kins[i]]] <- to_add_neg
    }
  }
  out_list <- list(scoreZ_df, pos_vals, neg_vals, kins_pos, kins_neg)
  names(out_list) <- c("Z-scores", "GSpos_values", "GSneg_values", "number_GSpos_per_kinase", "number_GSneg_per_kinase")
  return(out_list)
}

#' .ROC_sampler
#'
#' Computes AUROC (Area Under the Receiver Operating Characteristic curve) for given positive and negative gene sets.
#'
#' @param GSpos A vector of positive gene sets.
#' @param GSneg A vector of negative gene sets.
#' @param Npert An integer specifying the number of permutations/samples to generate. Default is 1000.
#' @param set_size Either a single value or a vector of two values specifying the sizes of the sampled positive and negative sets. Default is "GSpos_size". If "all", uses all elements in GSpos and GSneg.
#' @param ROC_obj A logical value indicating whether to return ROC objects along with AUROC values. Default is FALSE.
#' @param seed An integer used to set the random seed for reproducibility. Default is 1.
#' @param ... Additional arguments passed to `pROC::roc` function.
#'
#' @return A list containing:
#' \item{sample_AUROCs}{A numeric vector of AUROC values for each sample.}
#' \item{sample_ROC_objects}{(Optional) A list of ROC objects for each sample. Returned if `ROC_obj` is TRUE.}
#' \item{GS_pos_samples}{A list of sampled positive gene sets.}
#' \item{GS_neg_samples}{A list of sampled negative gene sets.}
#'
#' @export
#'
.ROC_sampler <- function(GSpos, GSneg, Npert=1000, set_size= "GSpos_size", ROC_obj=F, seed=1, ...){
  set.seed(seed)
  GS_pos_samples <- list()
  GS_neg_samples <- list()
  sample_ROC <- list()
  auroc <- numeric()
  if(length(set_size) == 1) { if(set_size=="all"){
    GS_pos_samples <- GSpos
    GS_neg_samples <- GSneg
    sample_ROC <- pROC::roc(controls=GS_neg_samples, cases=GS_pos_samples, ...)
    auroc <- as.numeric(sample_ROC$auc)
  } else {
    if((length(set_size) == 1) & (set_size=="GSpos_size")){
      GS_pos_samples <- GSpos
      for(i in 1:Npert){
        GS_neg_samples[[i]] <- sample(GSneg, length(GSpos), replace = F)
        sample_ROC[[i]] <- pROC::roc(controls=GS_neg_samples[[i]], cases=GS_pos_samples, ...)
        auroc[i] <- as.numeric(sample_ROC[[i]]$auc)
      }
    }}} else {
      for(i in 1:Npert){
        GS_pos_samples[[i]] <- sample(GSpos, set_size[1], replace = F)
        GS_neg_samples[[i]] <- sample(GSneg, set_size[2], replace = F)
        sample_ROC[[i]] <- pROC::roc(controls=GS_neg_samples[[i]], cases=GS_pos_samples[[i]], ...)
        auroc[[i]] <- as.numeric(sample_ROC[[i]]$auc)
      }

    }
  if(ROC_obj==T){
    out_list <- list(auroc, sample_ROC, GS_pos_samples, GS_neg_samples)
    names(out_list) <- c("sample_AUROCs", "sample_ROC_objects", "GS_pos_samples", "GS_neg_samples")
  } else {
    out_list <- list(auroc, GS_pos_samples, GS_neg_samples)
    names(out_list) <- c("sample_AUROCs", "GS_pos_samples", "GS_neg_samples")
  }
  return(out_list)
}
