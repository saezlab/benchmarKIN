% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmarkeR-multiomics.R
\name{.ROC_sampler}
\alias{.ROC_sampler}
\title{.ROC_sampler}
\usage{
.ROC_sampler(
  GSpos,
  GSneg,
  Npert = 1000,
  set_size = "GSpos_size",
  ROC_obj = F,
  seed = 1,
  ...
)
}
\arguments{
\item{GSpos}{A vector of positive gene sets.}

\item{GSneg}{A vector of negative gene sets.}

\item{Npert}{An integer specifying the number of permutations/samples to generate. Default is 1000.}

\item{set_size}{Either a single value or a vector of two values specifying the sizes of the sampled positive and negative sets. Default is "GSpos_size". If "all", uses all elements in GSpos and GSneg.}

\item{ROC_obj}{A logical value indicating whether to return ROC objects along with AUROC values. Default is FALSE.}

\item{seed}{An integer used to set the random seed for reproducibility. Default is 1.}

\item{...}{Additional arguments passed to \code{pROC::roc} function.}
}
\value{
A list containing:
\item{sample_AUROCs}{A numeric vector of AUROC values for each sample.}
\item{sample_ROC_objects}{(Optional) A list of ROC objects for each sample. Returned if \code{ROC_obj} is TRUE.}
\item{GS_pos_samples}{A list of sampled positive gene sets.}
\item{GS_neg_samples}{A list of sampled negative gene sets.}
}
\description{
Computes AUROC (Area Under the Receiver Operating Characteristic curve) for given positive and negative gene sets.
}
