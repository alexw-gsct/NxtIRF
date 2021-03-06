% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CollateData.R
\name{CollateData}
\alias{CollateData}
\title{Processes data from IRFinder output}
\usage{
CollateData(
  Experiment,
  reference_path,
  output_path,
  IRMode = c("SpliceOverMax", "SpliceMax"),
  localHub = FALSE,
  ah = AnnotationHub(localHub = localHub),
  low_memory_mode = FALSE,
  samples_per_block = 16,
  n_threads = 1
)
}
\arguments{
\item{Experiment}{A 2-column data frame (generated by \verb{FindSamples()), with the first column designating the sample names, and the 2nd column containing the primary IRFinder output file (of type }.txt.gz`). A third optional column
can contain the coverage files of the corresponding samples.
NB: all other columns are ignored.\cr
We recommend using the output of \code{\link[=Find_IRFinder_Output]{Find_IRFinder_Output()}} here.}

\item{reference_path}{THe path to the reference generated by BuildReference()}

\item{output_path}{The path for the output files to be generated by this function.}

\item{IRMode}{The algorithm to calculate 'splice abundance' in IR quantification.
The original algorithm by Middleton et al (2017) proposes \code{SpliceMax}, which
calculates the number of mapped splice events that share the boundary coordinate
of either the left or right flanking exon (SpliceLeft, SpliceRight) and defines
splice abundance as the larger of the two values. NxtIRF proposes a new algorithm,
\code{SpliceOverMax}, to account for the possibility that the major isoform shares neither
boundary, but arises from either of the flanking "exon islands". Exon islands are
contiguous regions covered by exons from any transcript (except those designated
as 'retained_intron' or 'sense_intronic'), and are separated by obligate
intronic regions (genomic regions that are introns for all transcripts). Note for
introns that are internal to a single exon island (i.e. akin to "known-exon"
introns from IRFinder), \code{SpliceOverMax} uses \code{findOverlaps()} to summate competing
mapped splice reads.}

\item{localHub}{See \code{?AnnotationHub::AnnotationHub()}. Setting \code{TRUE} will run \code{AnnotationHub()}
in offline mode}

\item{ah}{An AnnotationHub object containing the records \code{ah_genome} and/or \code{ah_transcriptome}
records to be used.}

\item{low_memory_mode}{Use this mode in memory-limited systems with many samples (> 16).
CollateData will write to file for every N samples as defined by \code{samples_per_block = N}.
Memory usage is often a problem for large datasets or using multiple cores. If you
experience crashes due to running out of memory, set this to true and make sure
\code{n_threads = 1}}

\item{samples_per_block}{How many samples to process per thread. Use in conjunction
with low_memory_mode to lower memory requirements}

\item{n_threads}{The number of threads to use.
records to be used.}
}
\value{
None. \code{CollateData()} writes to the directory given by \code{output_path}
}
\description{
CollateData unifies a list of IRFinder output files belonging to an experiment.
It is assumed every sample is analysed using the same IRFinder reference.
The combination of junction counts and IR quantification from IRFinder is used
to calculate percentage spliced in (PSI) of alternative splice events, and percent
intron retention (PIR) of retained introns. Also, QC information is extracted,
and data is collated into fst files for fast downstream access such as \code{MakeSE()}
}
\examples{
\donttest{
# Collate experiment from IRFinder output contained in "./IRFinder_output"
# Using laptop (low memory, low CPU cores)
CollateData(
    Experiment = Find_IRFinder_Output("./IRFinder_output"),
    reference_path = "./Reference",
    output_path = "./NxtIRF_output",
    low_memory_mode = TRUE,
    n_threads = 1
)
# Run the same, using a server with lots of memory and cores
CollateData(
    Experiment = Find_IRFinder_Output("./IRFinder_output"),
    reference_path = "./Reference",
    output_path = "./NxtIRF_output",
    low_memory_mode = FALSE,
    n_threads = 8
)
}
}
\seealso{
"NxtIRF: 2 - Building an Experiment"
(accessed via \code{browseVignettes("NxtIRF")})
}
