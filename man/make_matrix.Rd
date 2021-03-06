% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dash_utils.R
\name{make_matrix}
\alias{make_matrix}
\title{Constructs a matrix containing PSI values of the given ASE events}
\usage{
make_matrix(
  se,
  event_list,
  sample_list = colnames(se),
  method = c("PSI", "logit", "Z-score"),
  depth_threshold = 10,
  logit_max = 5,
  na.percent.max = 0.1
)
}
\arguments{
\item{se}{A NxtIRF SummarizedExperiment}

\item{event_list}{A character vector containing the row names of ASE events (as given by the
\code{EventName} column of differential ASE results table using \code{limma_ASE()} or \code{DESeq_ASE()})}

\item{sample_list}{(default = \code{colnames(se)}) A list of sample names referring to the subset
of samples in the given experiment to be included in the returned matrix}

\item{method}{The values to be returned (default = "PSI"). It can alternately be "logit"
which returns logit-transformed PSI values, or "Z-score" which returns Z-score-transformed PSI values}

\item{depth_threshold}{(default = 10) If any PSI is derived from raw values with both isoforms represented
by reads below this value, it is given as NA as the uncertainty of PSI would be deemed too highlight}

\item{logit_max}{The max or min logit values to be capped at. Because logit(0) == -Inf and logit(1) = Inf,
this function caps logit values using logit_max. (Only relevant if \code{method = "logit"})}

\item{na.percent.max}{(default = 0.1) The maximum number of NA values in an event for the PSI
values for each event to be returned. Most heatmap functions will spring an error if there are too
many NA values in any given row. This option caps the number of NA values to avoid returning this error.}
}
\value{
A matrix of PSI (or alternate) values, with columns as samples and rows as ASE events.
}
\description{
This function takes an input SummarizedExperiment \code{se}, a list of alternative
splicing events \code{event_list}, and (optionally) a list of sample names \code{sample_list}.
It returns a matrix containing PSI values with columns as samples and rows as ASE events.
}
