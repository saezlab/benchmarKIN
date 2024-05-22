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
    all_roc_list[[names(all_score_list_pos)[i]]] <- .ROC_sampler(all_score_list_pos[[i]], all_score_list_neg[[i]], set_size = samps[[i]], ...)
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
.ROC_sampler <- function(GSpos, GSneg, Npert=1000, set_size="GSpos_size", ROC_obj=F, seed=1, ...){
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
