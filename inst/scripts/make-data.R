#' Making Mappability Exclusion Region BED files from FASTA file
#'
#' Making a Mappability Region Exclusion BED file requires 3 steps.
#'   1) Generate mappability reads using `run_IRFinder_GenerateMapReads()`, 
#'      using the primary assembly genome FASTA as `genome.fa`
#'   2) Align `out.fa` to the corresponding genome using a genome splice-aware aligner such as STAR
#'   3) Process the aligned BAM file using `run_IRFinder_MapExclusionRegions()`
#' @examples
#' \dontrun{
#' FASTA = "/path/to/genome.fa"
#' run_IRFinder_GenerateMapReads(
#'      genome.fa = FASTA,
#'      out.fa = "/path/to/mappability_reads.fa",
#'      read_len = 70,
#'      read_stride = 10,
#'      error_pos = 35
#'  )
#'
#' # Now run STAR in the command line, e.g.:
#' # STAR \
#' # --genomeDir /path/to/Reference \
#' # --genomeLoad NoSharedMemory \
#' # --runThreadN 4 --outStd SAM --outSAMmode NoQS \
#' # --outSAMattributes None \
#' # --outFilterMultimapNmax 1 \
#' # --readFilesIn /path/to/mappability_reads.fa \
#' # > /path/to/genome_fragments.sam
#'
#' run_IRFinder_MapExclusionRegions(
#'      bamfile = "/path/to/genome_fragments.sam",
#'      output_file = "/path/to/MappabilityExclusionBED.txt",
#'      threshold = 4,
#'      includeCov = FALSE
#'  )
#' }
NULL
