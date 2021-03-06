% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuildRef.R
\name{GenerateMappabilityReads}
\alias{GenerateMappabilityReads}
\title{Generates Mappability reads from a genome FASTA file}
\usage{
GenerateMappabilityReads(
  fasta = "genome.fa",
  ah_genome = "",
  reference_path,
  read_len = 70,
  read_stride = 10,
  error_pos = 35,
  verbose = FALSE,
  localHub = FALSE,
  ah = AnnotationHub(localHub = localHub)
)
}
\arguments{
\item{fasta}{The path to the user-supplied genome fasta file}

\item{ah_genome}{The name of the AnnotationHub record containing the genome 2bit file.
Leave blank to use user-supplied \code{fasta} file.}

\item{reference_path}{The directory to store the reference files}

\item{read_len}{The nucleotide length of the generated reads}

\item{read_stride}{The nucleotide distance between sequentially generated
reads}

\item{error_pos}{The position of the generated error within the reads}

\item{verbose}{Whether additional status messages are shown}

\item{localHub}{See \code{?AnnotationHub::AnnotationHub()}. Setting \code{TRUE} will run \code{AnnotationHub()}
in offline mode}

\item{ah}{An AnnotationHub object containing the records \code{ah_genome} and/or \code{ah_transcriptome}
records to be used.}
}
\value{
None. Output reads are stored in \code{MappabilityReads.genome.fa}
in the given reference path. If from an AnnotationHub object, \code{genome}
is replaced by the name of the AnnotationHub record
}
\description{
This function generates sequential reads from a supplied FASTA file. These
reads of length \code{read_len = 70} are generated from a starting point spaced
\code{read_stride = 10} nucleotides apart, alternating between forward and reverse
strand, and contains a single nucleotide error \code{error_pos = 35} nucleotides
from the start of the read.\cr\cr
This function replicates the functionality of generateReadsError.pl in
the original IRFinder
}
\seealso{
\url{https://github.com/williamritchie/IRFinder/blob/master/bin/util/generateReadsError.pl}
}
