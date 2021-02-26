#' Searches a specified path for files of a particular file pattern
#'
#' This convenience function identifies files with the specified suffix.
#' The sample name is assumed to be the name of the file minus its suffix
#' (`use_subdir = FALSE`), but can also be changed to be the names 
#' of the parent directory (`use_subdir = TRUE`). See example below.
#'
#' @param sample_path The path in which to recursively search for files
#'   that match the given `suffix`
#' @param suffix A vector of or or more strings that specifies the file suffix 
#'   (e.g. '.bam' denotes BAM files, whereas ".txt.gz" denotes gzipped txt files).
#' @param suffix_type A vector of string that determines the column
#'   names of the files retrieved by `suffix`. Must be the same length as `suffix`
#' @param use_subdir Whether to assume the directory name containing
#'   the found files denote the sample name. If `FALSE`, the base name
#'   of the file is assumed to be the sample name. See below example.
#' @return A 2-column data frame with the first column containing
#'   the sample name, and the second column being the file path.
#' @examples
#'\dontrun{
#'  # If your sample files are contained like this:
#'  #   /path/to/project/samples/Apple.BAM and /path/to/project/samples/Banana.BAM
#'  df = FindSamples(sample_path = '/path/to/project/samples', suffix = ".BAM",
#'   use_subdir = FALSE)
#'  # returns df with 
#'  #   df$sample = c("Apple", "Banana")
#'  #   df$path = c("/path/to/project/samples/Apple.BAM",
#'  #     "/path/to/project/samples/Banana.BAM")
#'
#'  # Alternately, if your sample files are contained like this:
#'  #   /path/to/project/samples/Apple/Aligned.BAM and
#'  #   /path/to/project/samples/Banana/Aligned.BAM
#'  df = FindSamples(sample_path = '/path/to/project/samples', suffix = ".BAM",
#'   use_subdir = TRUE)
#'  # returns df with 
#'  #   df$sample = c("Apple", "Banana"), 
#'  #   df$path = c("/path/to/project/samples/Apple/Aligned.BAM",
#'  #     "/path/to/project/samples/Banana/Aligned.BAM")
#' }
#' @md
#' @export
FindSamples <- function(sample_path, suffix = ".txt.gz", 
            suffix_type = "path", use_subdir = FALSE) {
    
    assert_that(dir.exists(sample_path),
        msg = "Given path does not exist")
    assert_that(length(suffix) > 0 && length(suffix) == length(suffix_type),
        msg = "suffix must be of length greater than zero and same length as suffix_type")

    df.list = list()
    for(i in seq_len(length(suffix))) {
        files_found = list.files(pattern = paste0("\\", suffix[i], "$"),
        path = normalizePath(sample_path), full.names = TRUE, recursive = TRUE)    
        if(length(files_found) > 0) {
            df = data.frame(sample = "", path = files_found)
            if(use_subdir) {
                df$sample = basename(dirname(df$path))
            } else {
                df$sample = sub(suffix[i],"",basename(df$path))
            }
            colnames(df)[2] = suffix_type[i]
            df.list[[i]] = df
        } else {
            df.list[[i]] = NULL
        }
    }
    if(length(df.list) == 0) return(as.data.frame(df.list))
    final = df.list[[1]]
    if(length(suffix) > 1) {
        for(i in seq(2, length(suffix))) {
            if(!is.null(df.list[[i]])) {
                final = suppressMessages(left_join(final, df.list[[i]]))
            }
        }
    }
    return(final)
}

#' Convenience Function to find BAM files in a certain folder
#' 
#' Runs FindSamples to find BAM files. Assumes file names are names of samples.
#' 
#' @param sample_path The path in which to recursively search for BAM files
#' @param ... Additional parameters to pass into `FindSamples()`
#' @return A 2-column data frame with the first column containing
#'   the sample name, and the second column being the BAM file path.
#' @seealso [FindSamples()]
#' @md
#' @export
Find_Bams <- function(sample_path, ...) {
    return(FindSamples(sample_path, ".bam", "BAM", ...))
}

#' Convenience Function to find IRFinder output and COV files
#' 
#' Runs FindSamples to find IRFinder .txt.gz and .cov files. 
#' Assumes file names are names of samples.
#' 
#' @param sample_path The path in which to recursively search for BAM files
#' @param ... Additional parameters to pass into `FindSamples()`
#' @return A 3-column data frame with the first column containing
#'   the sample name, and the second column being the IRFinder main output path,
#'   and the third column being the COV file path.
#' @seealso [FindSamples()]
#' @md
#' @export
Find_IRFinder_Output <- function(sample_path, ...) {
    return(FindSamples(sample_path, c(".txt.gz", ".cov"), c("irf_file", "cov_file"), ...))
}

#' A wrapper function to call NxtIRF/IRFinder
#'
#' This function calls IRFinder on one or more BAM files.
#' @param bamfiles The file names of 1 or more BAM files
#' @param sample_names The sample names of the given BAM files. Must
#'   be a vector of the same length as `bamfiles`
#' @param reference_path The directory of the NxtIRF reference
#' @param output_path The directory where NxtIRF/IRFinder output
#'   should be stored
#' @param n_threads The number of threads to use. On Linux / Windows, this will
#'   use OpenMP from within the C++ subroutine. On Macs, BiocParallel
#'   MulticoreParam will be used on single-threaded NxtIRF/IRFinder
#' @param run_featureCounts Whether this function will run 
#'   `Rsubread::featureCounts()` on the BAM files. If so, the output will be
#'   saved to "main.FC.Rds" in the output directory as a list object
#' @param localHub Set as TRUE to disable AnnotationHub online mode
#' @param ah An AnnotationHub object.
#' @return None. `IRFinder()` will save output to `output_path`. \cr\cr
#'   sample.txt.gz: The main IRFinder output file containing the quantitation
#'   of IR and splice junctions, as well as QC information\cr\cr
#'   sample.cov: Contains coverage information in compressed binary. This
#'   format is 5-10X faster than BigWig format (see [GetCoverage()])\cr\cr
#'   main.FC.Rds: A single file containing gene counts for the whole dataset
#'   (only if `run_featureCounts == TRUE`)
#' @md
#' @export
IRFinder <- function(
        bamfiles = "Unsorted.bam", 
        sample_names = "sample1",
        reference_path = "./Reference",
        output_path = "./IRFinder_Output",
        n_threads = 1,
        run_featureCounts = FALSE,
        localHub = FALSE,
        ah = AnnotationHub(localHub = localHub)
        ) {

    assert_that(length(bamfiles) == length(sample_names),
        msg = "Number of BAM files and sample names must be the same")

    assert_that(dir.exists(dirname(output_path)),
        msg = paste(dirname(output_path), " - path does not exist"))

    if(!dir.exists(output_path)) dir.create(output_path)

    s_output = file.path(normalizePath(output_path), sample_names)
    
    run_IRFinder_multithreaded(
        reference_path = reference_path,
        bamfiles = bamfiles,
        output_files = s_output,
        max_threads = n_threads,
        run_featureCounts = run_featureCounts,
        localHub = localHub,
        ah = ah
    )
}

#' Processes data from IRFinder output
#'
#' CollateData unifies a list of IRFinder output files belonging to an experiment.
#' It is assumed every sample is analysed using the same IRFinder reference.
#' The combination of junction counts and IR quantification from IRFinder is used
#' to calculate percentage spliced in (PSI) of alternative splice events, and percent
#' intron retention (PIR) of retained introns. Also, QC information is extracted,
#' and data is collated into fst files for fast downstream access such as `MakeSE()`
#'
#' @param Experiment A 2-column data frame (generated by `FindSamples()), with
#'   the first column designating the sample names, and the 2nd column containing
#'   the primary IRFinder output file (of type `.txt.gz`). A third optional column
#'   can contain the coverage files of the corresponding samples.
#' @param reference_path THe path to the reference generated by BuildReference()
#' @param output_path The path for the output files to be generated by this function.
#' @param IRMode The algorithm to calculate 'splice abundance' in IR quantification.
#'   The original algorithm by Middleton et al (2017) proposes `SpliceMax`, which
#'   calculates the number of mapped splice events that share the boundary coordinate
#'   of either the left or right flanking exon (SpliceLeft, SpliceRight) and defines 
#'   splice abundance as the larger of the two values. NxtIRF proposes a new algorithm,
#'   `SpliceOverMax`, to account for the possibility that the major isoform shares neither 
#'   boundary, but arises from either of the flanking "exon islands". Exon islands are 
#'   contiguous regions covered by exons from any transcript (except those designated
#'   as 'retained_intron' or 'sense_intronic'), and are separated by obligate
#'   intronic regions (genomic regions that are introns for all transcripts). Note for
#'   introns that are internal to a single exon island (i.e. akin to "known-exon"
#'   introns from IRFinder), `SpliceOverMax` uses `findOverlaps()` to summate competing
#'   mapped splice reads.
#' @param localHub See `?AnnotationHub::AnnotationHub()`. Setting `TRUE` will run `AnnotationHub()`
#'   in offline mode
#' @param ah An AnnotationHub object containing the records `ah_genome` and/or `ah_transcriptome`
#'   records to be used.
#' @param low_memory_mode Use this mode in memory-limited systems with many samples (> 16).
#'   CollateData will write to file for every N samples as defined by `samples_per_block = N`.
#'   Memory usage is often a problem for large datasets or using multiple cores. If you
#'   experience crashes due to running out of memory, set this to true and make sure
#'   `n_threads = 1`
#' @param samples_per_block How many samples to process per thread. Use in conjunction
#'   with low_memory_mode to lower memory requirements
#' @param n_threads The number of threads to use.
#'   records to be used.
#' @return None. `CollateData()` writes to the directory given by `output_path`
#' @md
#' @export
CollateData <- function(Experiment, reference_path, output_path,
        IRMode = c("SpliceOverMax", "SpliceMax"), 
        localHub = FALSE, ah = AnnotationHub(localHub = localHub), 
        low_memory_mode = FALSE, samples_per_block = 16, n_threads = 1) {
################################################################################
        
    assert_that(is(Experiment, "data.frame"),
        msg = "Experiment object needs to be a data frame")
    assert_that(ncol(Experiment) >= 2,
        msg = paste("Experiment needs to contain at least two columns,",
        "with the first 2 columns containing",
        "(1) sample name and (2) IRFinder output"))
    assert_that(file.exists(file.path(reference_path, "settings.Rds")),
        msg = paste(file.path(reference_path, "settings.Rds"),
            "does not exist"))
    
    coverage_files = ""
    
    Experiment.df = as.data.frame(Experiment)
    
    if(ncol(Experiment) > 2 && all(file.exists(Experiment.df[,3]))) {
        coverage_files = Experiment.df[,3]
    }
    if(!IsCOV(coverage_files)){
        message("Some coverage files do not exist or are corrupted")
        coverage_files = ""
    }
    if(length(coverage_files) != nrow(Experiment)) {
        message("The number of coverage files must equal the number of samples")
        coverage_files = ""
    }
    
    IRMode = match.arg(IRMode)
    assert_that(IRMode != "",
        msg = "IRMode must be either 'SpliceOverMax' (default) or 'SpliceMax'")

    BPPARAM = BiocParallel::bpparam()
    n_threads_to_use = as.numeric(n_threads)
    assert_that(!is.na(n_threads_to_use),
        msg = "n_threads must be a numeric value")
    assert_that(n_threads_to_use <= (parallel::detectCores() - 2) | n_threads_to_use == 1,
        msg = paste("NxtIRF does not support more threads than",
        "parallel::detectCores() - 2") )
    
    if(Sys.info()["sysname"] == "Windows") {
        BPPARAM_mod = BiocParallel::SnowParam(n_threads_to_use)
        message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
    } else {
        BPPARAM_mod = BiocParallel::MulticoreParam(n_threads_to_use)
        message(paste("Using MulticoreParam", BPPARAM_mod$workers, "threads"))
    }

    settings = readRDS(file.path(reference_path, "settings.Rds"))
    if(settings$ah_genome != "") {
        genome = .fetch_AH(settings$ah_genome, localHub = localHub)
    } else {
        genome = rtracklayer::TwoBitFile(file.path(reference_path, 
            "resource", "genome.2bit"))
    }
    
    colnames(Experiment)[c(1,2)] = c("sample", "path")
    assert_that(all(vapply(Experiment$path, file.exists, logical(1))),
        msg = "Some files in Experiment do not exist")
    assert_that(dir.exists(dirname(output_path)),
        msg = paste("Parent directory of output path:", 
        dirname(output_path), "needs to exist"))
        
    # Create a subdirectory for each sample within output_path
    base_output_path = normalizePath(dirname(output_path)) 
    norm_output_path = file.path(base_output_path, basename(output_path))
    if(!dir.exists(norm_output_path)) {
        dir.create(norm_output_path)
    }
    temp_output_path = file.path(norm_output_path, "temp")
    if(!dir.exists(temp_output_path)) {
        dir.create(temp_output_path)
    }        
    if(!dir.exists(file.path(norm_output_path, "samples"))) {
        dir.create(file.path(norm_output_path, "samples"))
    }

    df.internal = as.data.table(Experiment[,c(1,2)])
    df.internal$paired = FALSE
    df.internal$strand = 0
    df.internal$depth = 0
    df.internal$mean_frag_size = 0
    df.internal$directionality_strength = 0
    df.internal$Intergenic_Fraction = 0
    df.internal$rRNA_Fraction = 0
    df.internal$NonPolyA_Fraction = 0
    df.internal$Mitochondrial_Fraction = 0
    df.internal$Unanno_Jn_Fraction = 0
    df.internal$NMD_Jn_Fraction = 0
    df.internal$Fraction_Splice_Reads = 0
    df.internal$Fraction_Span_Reads = 0

    df.internal$IRBurden_clean = 0
    df.internal$IRBurden_exitrons = 0
    df.internal$IRBurden_clean_unstranded = 0
    df.internal$IRBurden_exitrons_unstranded = 0
    df.internal$IRBurden_antisense = 0
    
    n_jobs = min(ceiling(nrow(df.internal) / samples_per_block), 
        BPPARAM_mod$workers)
    jobs = NxtIRF.SplitVector(seq_len(nrow(df.internal)), n_jobs)    
    n_jobs = length(jobs)
    
################################################################################
    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(0.01, message = "Compiling Sample Stats")
    }
    message("Compiling Sample Stats")
    df.internal = suppressWarnings(rbindlist(
        BiocParallel::bplapply(
            seq_len(n_jobs),
            function(x, jobs, df.internal) {
                suppressPackageStartupMessages({
                    requireNamespace("data.table")
                    requireNamespace("stats")
                })
                work = jobs[[x]]
                block = df.internal[work]
                for(i in seq_len(length(work))) {
                    data.list = get_multi_DT_from_gz(
                        normalizePath(block$path[i]), 
                        c("BAM", "Directionality", "QC")) 

                    stats = data.list$BAM
                    direct = data.list$Directionality
                    QC = data.list$QC
                    if(stats$Value[3] == 0 & stats$Value[4] > 0) {
                        block$paired[i] = TRUE
                        block$depth[i] = stats$Value[4]
                        block$mean_frag_size[i] = stats$Value[2] / 
                            stats$Value[4]
                    } else if(stats$Value[3] > 0 && 
                            stats$Value[4] / stats$Value[3] / 1000) {
                        block$paired[i] = TRUE
                        block$depth[i] = stats$Value[4]
                        block$mean_frag_size[i] = stats$Value[2] / 
                            stats$Value[4]
                    } else {
                        block$paired[i] = FALSE
                        block$depth[i] = stats$Value[3]
                        block$mean_frag_size[i] = stats$Value[2] / 
                            stats$Value[3]
                    }
                    block$strand[i] = direct$Value[9]

                    # QC
                    block$directionality_strength[i] = direct$Value[8]
                    block$Intergenic_Fraction[i] =
                        QC$Value[QC$QC == "Intergenic Reads"] / 
                            block$depth[i]
                    block$rRNA_Fraction[i] =    
                        QC$Value[QC$QC == "rRNA Reads"] / 
                            block$depth[i]
                    block$NonPolyA_Fraction[i] =
                        QC$Value[QC$QC == "NonPolyA Reads"] / 
                            block$depth[i]
                    block$Mitochondrial_Fraction[i] =
                        QC$Value[QC$QC == "Mitochondrial Reads"] / 
                            block$depth[i]
                    block$Unanno_Jn_Fraction[i] =
                        QC$Value[QC$QC == "Unannotated Junctions"] / 
                        (QC$Value[QC$QC == "Unannotated Junctions"] +
                        QC$Value[QC$QC == "Annotated Junctions"])
                    block$NMD_Jn_Fraction[i] =
                        QC$Value[QC$QC == "NMD Junctions"] / 
                        QC$Value[QC$QC == "Annotated Junctions"]
                    block$Fraction_Splice_Reads[i] =
                        QC$Value[QC$QC == "Annotated Junctions"] / 
                        block$depth[i]
                    block$Fraction_Span_Reads[i] =
                        QC$Value[QC$QC == "Spans Reads"] / 
                            block$depth[i]

# IRBurden calculations                            
                    block$IRBurden_clean_unstranded[i] =
                        QC$Value[QC$QC == "Non-Directional Clean IntronDepth Sum"] / 
                        (QC$Value[QC$QC == "Non-Directional Clean IntronDepth Sum"] +
                        QC$Value[QC$QC == "Annotated Junctions"])
                    block$IRBurden_exitrons_unstranded[i] =
                        QC$Value[QC$QC == "Non-Directional Known-Exon IntronDepth Sum"] / 
                        (QC$Value[QC$QC == "Non-Directional Known-Exon IntronDepth Sum"] +
                        QC$Value[QC$QC == "Annotated Junctions"])
                    block$IRBurden_antisense[i] =
                        QC$Value[QC$QC == "Non-Directional Anti-Sense IntronDepth Sum"] / 
                        (QC$Value[QC$QC == "Non-Directional Anti-Sense IntronDepth Sum"] +
                        QC$Value[QC$QC == "Annotated Junctions"])
                    if(block$strand[i] != 0) {
                        block$IRBurden_clean[i] =
                            QC$Value[QC$QC == "Directional Clean IntronDepth Sum"] / 
                            (QC$Value[QC$QC == "Directional Clean IntronDepth Sum"] +
                            QC$Value[QC$QC == "Annotated Junctions"])
                        block$IRBurden_exitrons[i] =
                            QC$Value[QC$QC == "Directional Known-Exon IntronDepth Sum"] / 
                            (QC$Value[QC$QC == "Directional Known-Exon IntronDepth Sum"] +
                            QC$Value[QC$QC == "Annotated Junctions"])
                    }
                }
                return(block)
            }, jobs = jobs, df.internal = df.internal, BPPARAM = BPPARAM_mod
        )
    ))
    
    if(any(df.internal$strand == 0)) {
        runStranded = FALSE
    } else {
        runStranded = TRUE
    }
    
    # TODO: Simply check version of reference used to generate the IRFinder files
    
    # Compile junctions and IR lists first, save to temp files
    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(0.15, message = "Compiling Junction List")
    }
    message("Compiling Junction List")       
    # Compile junc.common via merge
    junc.list = suppressWarnings(BiocParallel::bplapply(
        seq_len(n_jobs),
        function(x, jobs, df.internal, temp_output_path) {
            suppressPackageStartupMessages({
                requireNamespace("data.table")
                requireNamespace("stats")
            })
            work = jobs[[x]]
            block = df.internal[work]
            junc.segment = NULL
            for(i in seq_len(length(work))) {
            # junc = get_multi_DT_from_gz(block$path[i], c("JC_seqname"))
            # junc = junc$JC_seqname
                junc = suppressWarnings(data.table::as.data.table(
                data.table::fread(block$path[i], skip = "JC_seqname")))
                data.table::setnames(junc, "JC_seqname", "seqnames")
            if(is.null(junc.segment)) {
                junc.segment = junc[,seq_len(4), with = FALSE]
            } else {
                junc.segment = merge(junc.segment, 
                    junc[,seq_len(4), with = FALSE], 
                    all = TRUE)
            }
            # Write temp file
            fst::write.fst(as.data.frame(junc), 
                file.path(temp_output_path, paste(block$sample[i], 
                "junc.fst.tmp", sep=".")))        
            }
            return(junc.segment)
        }, jobs = jobs, df.internal = df.internal, temp_output_path = temp_output_path, BPPARAM = BPPARAM_mod
    ))
    junc.common = NULL
    for(i in seq_len(length(junc.list))) {
        if(is.null(junc.common)) {
            junc.common = junc.list[[i]]
        } else {
            junc.common = merge(junc.common, junc.list[[i]], 
                all = TRUE, by = colnames(junc.common))
        }
    }
    rm(junc.list)
    gc()
    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(0.15, message = "Compiling Intron Retention List")
    }
    
    message("Compiling Intron Retention List")
    if(!runStranded) {
        irf = suppressWarnings(as.data.table(
        fread(df.internal$path[i], skip = "Nondir_")))
        setnames(irf, c("Nondir_Chr", "Start", "End", "Strand"), 
            c("seqnames","start","end", "strand"))
    } else {
        irf = suppressWarnings(as.data.table(
        fread(df.internal$path[i], skip = "Dir_Chr")))
        setnames(irf, c("Dir_Chr", "Start", "End", "Strand"), 
            c("seqnames","start","end", "strand"))          
    }
    irf.common = irf[,seq_len(6), with = FALSE]
    rm(irf)
    irf.list = suppressWarnings(BiocParallel::bplapply(
        seq_len(n_jobs),
        function(x, jobs, df.internal, temp_output_path, runStranded) {
            suppressPackageStartupMessages({
                requireNamespace("data.table")
                requireNamespace("stats")
                requireNamespace("openssl")
            })
            work = jobs[[x]]
            block = df.internal[work]
            irf.md5 = NULL
            for(i in seq_len(length(work))) {

                if(!runStranded) {
                    irf = suppressWarnings(data.table::as.data.table(
                    fread(block$path[i], skip = "Nondir_")))
                    setnames(irf, c("Nondir_Chr", "Start", "End", "Strand"), 
                        c("seqnames","start","end", "strand"))
                } else {
                    irf = suppressWarnings(data.table::as.data.table(
                    fread(block$path[i], skip = "Dir_Chr")))
                    setnames(irf, c("Dir_Chr", "Start", "End", "Strand"), 
                        c("seqnames","start","end", "strand"))          
                }
                if(is.null(irf.md5)) {
                    irf.md5 = openssl::md5(paste(irf$Name, collapse=" "))
                } else {
                    irf.md5 = unique(c(irf.md5, 
                        openssl::md5(paste(irf$Name, collapse=" "))))
                }
                fst::write.fst(as.data.frame(irf), 
                file.path(temp_output_path, 
                    paste(block$sample[i], "irf.fst.tmp", sep=".")))
            }
            return(irf.md5)
        }, jobs = jobs, df.internal = df.internal, temp_output_path = temp_output_path, 
        runStranded = runStranded, BPPARAM = BPPARAM_mod
    ))
    irf.md5.check = NULL
    for(i in seq_len(length(irf.list))) {
        if(is.null(irf.md5.check)) {
            irf.md5.check = irf.list[[i]]
        } else {
            irf.md5.check = unique(irf.md5.check, irf.list[[i]])
        }
    }
    rm(irf.list)
    gc()

    assert_that(length(irf.md5.check) == 1,
        msg = paste(
            "MD5 check of IRFinder introns are not the same.",
            "Perhaps some samples were processed by a different reference"
        )
    )

    irf.common[, start := start + 1]
    junc.common[, start := start + 1]

# Reassign +/- based on junctions.fst annotation
    # Annotate junctions
    
    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(0.15, 
            message = "Tidying up splice junctions and intron retentions")
    }  
    message("Tidying up splice junctions and intron retentions")
    
    candidate.introns = as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst")))

    junc.strand = unique(candidate.introns[, c("seqnames", "start", "end", "strand")])

    junc.common[, c("strand") := NULL]
    junc.common = unique(junc.common)
    junc.common = merge(junc.common, junc.strand, 
        all = TRUE, by = c("seqnames", "start", "end"))
    junc.common[is.na(get("strand")), c("strand") := "*"]

    left.gr = GRanges(seqnames = junc.common$seqnames, 
        ranges = IRanges(start = junc.common$start, 
        end = junc.common$start + 1), strand = "+")
    right.gr = GRanges(seqnames = junc.common$seqnames, 
        ranges = IRanges(start = junc.common$end - 1, 
        end = junc.common$end), strand = "+")
        
    left.seq = getSeq(genome, left.gr)
    right.seq = getSeq(genome, right.gr)

    junc.common$motif_pos = paste0(as.character(left.seq), as.character(right.seq))
    junc.common$motif_infer_strand = "n"
    junc.common[ get("motif_pos") %in% c("GTAG", "GCAG", "ATAC", "ATAG"), 
        c("motif_infer_strand") := "+"]
    junc.common[ get("motif_pos") %in% c("CTAC", "CTGC", "GTAT", "CTAT"), 
        c("motif_infer_strand") := "-"]
    junc.common[ get("motif_pos") %in% c("GTAC"), c("motif_infer_strand") := "n"]       
    # Do not accept un-annotated GTACs - too confusing
    # Exclude non-splice motifs (that are also not annotated - i.e. strand == "*")
    junc.common = junc.common[get("motif_infer_strand") != "n" | 
        get("strand") != "*"]

    # Use motif_infer_strand (only for *)
    junc.common[get("strand") == "*", c("strand") := get("motif_infer_strand")]
    junc.common$motif_infer_strand = NULL
    
    # Should splicing across gene groups be allowed? Exclude
    Genes = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path, "fst", "Genes.fst"))
    )
   
    # Assign region names to junctions:
    junc.common[, c("Event") := paste0(get("seqnames"), ":", 
        get("start"), "-", get("end"), "/", get("strand"))]
    
    candidate.introns[, c("transcript_biotype_2") := get("transcript_biotype")]
    candidate.introns[
        !(get("transcript_biotype") %in% c("protein_coding", "processed_transcript",
        "lincRNA", "antisense", "nonsense_mediated_decay")), 
        c("transcript_biotype_2") := "other"]

    candidate.introns[, c("transcript_biotype_2") := 
        factor(get("transcript_biotype_2"), c("protein_coding", "processed_transcript",
        "lincRNA", "antisense", "other", "nonsense_mediated_decay"), 
        ordered = TRUE)]
        
    if("transcript_support_level" %in% colnames(candidate.introns)) {
        setorderv(candidate.introns, c("transcript_biotype_2", "transcript_support_level"))
    } else {
        setorder(candidate.introns, "transcript_biotype_2")    
    }
    introns.unique = unique(candidate.introns, 
        by = c("seqnames", "start", "end", "width", "strand"))
    setorderv(introns.unique, c("seqnames", "start", "end", "strand"))

    junc.annotation = introns.unique[junc.common, 
        c("seqnames", "start", "end", "strand", "transcript_id", 
        "intron_number", "gene_name", "gene_id", "transcript_biotype"),
        on = c("seqnames", "start", "end", "strand")]
    
    # Use Exon Groups file to designate exon groups to all junctions
    Exon.Groups = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path, "fst", "Exons.Group.fst")),
        keep.extra.columns = TRUE)
    
    # Always calculate stranded for junctions
    # if(!runStranded) {
        # Exon.Groups = Exon.Groups[strand(Exon.Groups) == "*"]
    # } else {
        # Exon.Groups = Exon.Groups[strand(Exon.Groups) != "*"]    
    # }
    Exon.Groups.S = Exon.Groups[strand(Exon.Groups) != "*"]    
    
    junc.common.left = copy(junc.common)
    junc.common.left[, c("start") := get("start") - 1]
    junc.common.left[, c("end") := get("start") + 1]
    OL = suppressWarnings(
        findOverlaps(
            .grDT(junc.common.left), 
            .grDT(Exon.Groups.S)
        )
    )
    junc.common[, c("gene_group_left") := NA]
    junc.common[, c("exon_group_left") := NA]
    junc.common$gene_group_left[from(OL)] = Exon.Groups.S$gene_group[to(OL)]
    junc.common$exon_group_left[from(OL)] = Exon.Groups.S$exon_group[to(OL)]

    junc.common.right = copy(junc.common)
    junc.common.right[, c("end") := get("end") + 1]
    junc.common.right[, c("start") := get("end") - 1]
    OL = suppressWarnings(
        findOverlaps(
            .grDT(junc.common.right), 
            .grDT(Exon.Groups.S)
        )
    )
    junc.common[, c("gene_group_right") := NA]
    junc.common[, c("exon_group_right") := NA]
    junc.common$gene_group_right[from(OL)] = Exon.Groups.S$gene_group[to(OL)]
    junc.common$exon_group_right[from(OL)] = Exon.Groups.S$exon_group[to(OL)]
    
    junc.common[, c("JG_up") := ""]
    junc.common[, c("JG_down") := ""]
    junc.common[get("strand") == "+" & 
        !is.na(get("gene_group_left")) & !is.na(get("exon_group_left")), 
        c("JG_up") := paste(get("gene_group_left"), 
            get("exon_group_left"), sep="_")]
    junc.common[get("strand") == "-" & 
        !is.na(get("gene_group_right")) & !is.na(get("exon_group_right")), 
        c("JG_up") := paste(get("gene_group_right"), get("exon_group_right"), sep="_")]
    junc.common[get("strand") == "+" & 
        !is.na(get("gene_group_right")) & !is.na(get("exon_group_right")), 
        c("JG_down") := paste(get("gene_group_right"), get("exon_group_right"), sep="_")]
    junc.common[get("strand") == "-" & 
        !is.na(get("gene_group_left")) & !is.na(get("exon_group_left")), 
        c("JG_down") := paste(get("gene_group_left"), get("exon_group_left"), sep="_")]

    junc.common$gene_group_left = NULL
    junc.common$gene_group_right = NULL
    junc.common$exon_group_left = NULL
    junc.common$exon_group_right = NULL
    
    irf.common.left = copy(irf.common)
    irf.common.left[, c("start") := get("start") - 1]
    irf.common.left[, c("end") := get("start") + 1]
    OL = suppressWarnings(
        findOverlaps(
            .grDT(irf.common.left), 
            .grDT(Exon.Groups.S)
        )
    )
    irf.common[, c("gene_group_left") := NA]
    irf.common[, c("exon_group_left") := NA]
    irf.common$gene_group_left[from(OL)] = Exon.Groups.S$gene_group[to(OL)]
    irf.common$exon_group_left[from(OL)] = Exon.Groups.S$exon_group[to(OL)]

    irf.common.right = copy(irf.common)
    irf.common.right[, c("end") := get("end") + 1]
    irf.common.right[, c("start") := get("end") - 1]
    OL = suppressWarnings(
        findOverlaps(
            .grDT(irf.common.right), 
            .grDT(Exon.Groups.S)
        )
    )
    irf.common[, c("gene_group_right") := NA]
    irf.common[, c("exon_group_right") := NA]
    irf.common$gene_group_right[from(OL)] = Exon.Groups.S$gene_group[to(OL)]
    irf.common$exon_group_right[from(OL)] = Exon.Groups.S$exon_group[to(OL)]
    
    irf.common[, c("JG_up") := ""]
    irf.common[, c("JG_down") := ""]
    irf.common[get("strand") == "+" & 
        !is.na(get("gene_group_left")) & !is.na(get("exon_group_left")), 
        c("JG_up") := paste(get("gene_group_left"), 
            get("exon_group_left"), sep="_")]
    irf.common[get("strand") == "-" & 
        !is.na(get("gene_group_right")) & !is.na(get("exon_group_right")), 
        c("JG_up") := paste(get("gene_group_right"), get("exon_group_right"), sep="_")]
    irf.common[get("strand") == "+" & 
        !is.na(get("gene_group_right")) & !is.na(get("exon_group_right")), 
        c("JG_down") := paste(get("gene_group_right"), get("exon_group_right"), sep="_")]
    irf.common[get("strand") == "-" & 
        !is.na(get("gene_group_left")) & !is.na(get("exon_group_left")), 
        c("JG_down") := paste(get("gene_group_left"), get("exon_group_left"), sep="_")]
    
    if(!runStranded) {
        Exon.Groups.S = Exon.Groups[strand(Exon.Groups) == "*"]
    } else {
        Exon.Groups.S = Exon.Groups[strand(Exon.Groups) != "*"]    
    }
    irf.common.left = copy(irf.common)
    irf.common.left[, c("start") := get("start") - 1]
    irf.common.left[, c("end") := get("start") + 1]
    OL = suppressWarnings(
        findOverlaps(
            .grDT(irf.common.left), 
            .grDT(Exon.Groups.S)
        )
    )
    irf.common[, c("gene_group_left") := NA]
    irf.common[, c("exon_group_left") := NA]
    irf.common$gene_group_left[from(OL)] = Exon.Groups.S$gene_group[to(OL)]
    irf.common$exon_group_left[from(OL)] = Exon.Groups.S$exon_group[to(OL)]

    irf.common.right = copy(irf.common)
    irf.common.right[, c("end") := get("end") + 1]
    irf.common.right[, c("start") := get("end") - 1]
    OL = suppressWarnings(
        findOverlaps(
            .grDT(irf.common.right), 
            .grDT(Exon.Groups.S)
        )
    )
    irf.common[, c("gene_group_right") := NA]
    irf.common[, c("exon_group_right") := NA]
    irf.common$gene_group_right[from(OL)] = Exon.Groups.S$gene_group[to(OL)]
    irf.common$exon_group_right[from(OL)] = Exon.Groups.S$exon_group[to(OL)]
    
    irf.common[, c("IRG_up") := ""]
    irf.common[, c("IRG_down") := ""]
    irf.common[get("strand") == "+" & 
        !is.na(get("gene_group_left")) & !is.na(get("exon_group_left")), 
        c("IRG_up") := paste(get("gene_group_left"), 
            get("exon_group_left"), sep="_")]
    irf.common[get("strand") == "-" & 
        !is.na(get("gene_group_right")) & !is.na(get("exon_group_right")), 
        c("IRG_up") := paste(get("gene_group_right"), get("exon_group_right"), sep="_")]
    irf.common[get("strand") == "+" & 
        !is.na(get("gene_group_right")) & !is.na(get("exon_group_right")), 
        c("IRG_down") := paste(get("gene_group_right"), get("exon_group_right"), sep="_")]
    irf.common[get("strand") == "-" & 
        !is.na(get("gene_group_left")) & !is.na(get("exon_group_left")), 
        c("IRG_down") := paste(get("gene_group_left"), get("exon_group_left"), sep="_")]
    irf.common$gene_group_left = NULL
    irf.common$gene_group_right = NULL
    irf.common$exon_group_left = NULL
    irf.common$exon_group_right = NULL

    irf.common[, c("EventRegion") := 
        paste0(get("seqnames"), ":", get("start"), "-", 
            get("end"), "/", get("strand"))]

    Splice.Anno = as.data.table(read.fst(file.path(reference_path, "fst", "Splice.fst")))
    candidate.introns[, c("Event1a") := get("Event")]
    candidate.introns[, c("Event2a") := get("Event")]
    Splice.Anno[candidate.introns, on = "Event1a", 
        c("up_1a") := paste(get("i.gene_group_stranded"), 
        get("i.exon_group_stranded_upstream"), sep="_")]
    Splice.Anno[candidate.introns, on = "Event1a", 
        c("down_1a") := paste(get("i.gene_group_stranded"), 
        get("i.exon_group_stranded_downstream"), sep="_")]
    Splice.Anno[candidate.introns, on = "Event2a", 
        c("down_2a") := paste(get("i.gene_group_stranded"), 
        get("i.exon_group_stranded_downstream"), sep="_")]
    
    Splice.Anno[get("EventType") %in% c("MXE", "SE", "ALE", "A3SS"),
        c("JG_up") := get("up_1a")]
    Splice.Anno[get("EventType") %in% c("SE", "AFE", "A5SS"),
        c("JG_down") := get("down_1a")]
    Splice.Anno[get("EventType") %in% c("MXE"),
        c("JG_down") := get("down_2a")]
    
    Splice.Anno$up_1a = NULL
    Splice.Anno$down_1a = NULL
    Splice.Anno$down_2a = NULL
    Splice.Anno[, c("strand") := 
        tstrsplit(get("Event1a"), split="/")[[2]]]

    # Save irf.common, Splice.Anno
  
    if(!dir.exists(file.path(norm_output_path, "annotation"))) {
        dir.create(file.path(norm_output_path, "annotation"))
    }          
  
    write.fst(as.data.frame(junc.common), 
        file.path(norm_output_path, "annotation", "Junc.fst"))
    write.fst(as.data.frame(irf.common), 
        file.path(norm_output_path, "annotation", "IR.fst"))
    write.fst(as.data.frame(Splice.Anno), 
        file.path(norm_output_path, "annotation", "Splice.fst"))

    # make rowEvent brief here
    irf.anno.brief = irf.common[, c("Name", "EventRegion")]
    setnames(irf.anno.brief, "Name", "EventName")
    irf.anno.brief[, c("EventType") := "IR"]
    irf.anno.brief = irf.anno.brief[, c("EventName", "EventType", "EventRegion")]
    
    splice.anno.brief = Splice.Anno[, c("EventName", "EventType", "EventRegion")]
    
    rowEvent = rbind(irf.anno.brief, splice.anno.brief)    
    
    write.fst(rowEvent, file.path(norm_output_path, "rowEvent.brief.fst"))
    
    item.todo = c("Included", "Excluded", "Depth", "Coverage", "minDepth", 
        "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc", "junc_PSI", "junc_counts")

    se_output_path = norm_output_path
    
    # Implement filters:
    # Is_Protein_Coding: at least one of the options 
    #   is a valid protein_coding transcript
    # Triggers_NMD: at least one isoform is pure NMD
    # TSL: maximum TSL of either isoform

    IR_NMD = as.data.table(read.fst(
        file.path(reference_path, "fst", "IR.NMD.fst")))
    Splice.Options = as.data.table(read.fst(
        file.path(reference_path, "fst", "Splice.options.fst")))
    Transcripts = as.data.table(read.fst(
        file.path(reference_path, "fst", "Transcripts.fst")))
    Splice.Options[Splice.Anno, on = "EventID", 
        c("EventName") := get("i.EventName")]
    Splice.Options[Transcripts, on = "transcript_id", 
        c("transcript_biotype") := get("i.transcript_biotype")]

    Splice.Options.Summary = copy(Splice.Options)
    Splice.Options.Summary[, 
        c("tsl_min") := min(get("transcript_support_level")), 
        by = c("EventID", "isoform")]
    Splice.Options.Summary[, 
        c("any_is_PC") := any(get("is_protein_coding")), 
        by = c("EventID", "isoform")]  
    Splice.Options.Summary[, 
        c("all_is_NMD") := all(grepl("decay", get("transcript_biotype"))), 
        by = c("EventID", "isoform")]  

    rowEvent.Extended = copy(rowEvent)

    rowEvent.Extended[get("EventType") == "IR", 
        c("intron_id") := tstrsplit(get("EventName"), split="/")[[2]]]
    rowEvent.Extended[, c("Inc_Is_Protein_Coding") := FALSE]
    rowEvent.Extended[, c("Exc_Is_Protein_Coding") := FALSE]
    rowEvent.Extended[IR_NMD, on = "intron_id", 
        c("Exc_Is_Protein_Coding") := TRUE]
    rowEvent.Extended[IR_NMD, on = "intron_id", 
        c("Inc_Is_Protein_Coding") := (get("i.intron_type") == "CDS")]

    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "A"], 
        on = "EventName", c("Inc_Is_Protein_Coding") := get("i.any_is_PC")]
    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "B"], 
        on = "EventName", c("Exc_Is_Protein_Coding") := get("i.any_is_PC")]

    rowEvent.Extended[, c("Inc_Is_NMD") := FALSE]
    rowEvent.Extended[, c("Exc_Is_NMD") := FALSE]
    rowEvent.Extended[IR_NMD[!is.na(get("splice_is_NMD"))], on = "intron_id", 
        c("Exc_Is_NMD") := get("i.splice_is_NMD")]
    rowEvent.Extended[IR_NMD, on = "intron_id", 
        c("Inc_Is_NMD") := get("i.IRT_is_NMD")]
    rowEvent.Extended[get("EventType") == "IR" & 
        get("Exc_Is_Protein_Coding") == FALSE, 
        c("Exc_Is_NMD") := NA]
    rowEvent.Extended[get("EventType") == "IR" & 
        get("Inc_Is_Protein_Coding") == FALSE, 
        c("Inc_Is_NMD") := NA]

    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "A"], 
        on = "EventName", c("Inc_Is_NMD") := get("i.all_is_NMD")]
    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "B"], 
        on = "EventName", c("Exc_Is_NMD") := get("i.all_is_NMD")]

    rowEvent.Extended[candidate.introns, on = "intron_id", 
        c("Inc_TSL") := get("i.transcript_support_level")]
    rowEvent.Extended[candidate.introns, on = "intron_id", 
        c("Exc_TSL") := get("i.transcript_support_level")]

    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "A"], 
        on = "EventName", c("Inc_TSL") := get("i.tsl_min")]
    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "B"], 
        on = "EventName", c("Exc_TSL") := get("i.tsl_min")]

    # define Event1 / Event2
    rowEvent.Extended[get("EventType") == "IR", 
        c("Event1a") := get("EventRegion")]
    rowEvent.Extended[Splice.Anno, on = "EventName",
        c("Event1a", "Event2a", "Event1b", "Event2b") := 
        list(get("i.Event1a"), get("i.Event2a"), 
            get("i.Event1b"), get("i.Event2b"))]

    write.fst(rowEvent.Extended, file.path(se_output_path, "rowEvent.fst"))

# Write junc_PSI index
    junc_PSI = junc.common[, c("seqnames", "start", "end", "strand")]
    write.fst(junc_PSI, file.path(se_output_path, "junc_PSI_index.fst"))

    rm(candidate.introns, introns.unique)
    gc()


    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(0.15, message = "Generating NxtIRF FST files")
    }  
    message("Generating NxtIRF FST files")
      
    agg.list <- suppressWarnings(BiocParallel::bplapply(seq_len(n_jobs),
        function(x, jobs, df.internal, norm_output_path) {
            suppressPackageStartupMessages({
                requireNamespace("data.table")
                requireNamespace("stats")
            })
            
            # Read this from fst file
            rowEvent = as.data.table(read.fst(
                file.path(norm_output_path, "rowEvent.brief.fst")))
            junc.common = as.data.table(read.fst(
                file.path(norm_output_path, "annotation", "Junc.fst")))
            irf.common = as.data.table(read.fst(
                file.path(norm_output_path, "annotation","IR.fst")))
            Splice.Anno = as.data.table(read.fst(
                file.path(norm_output_path, "annotation","Splice.fst")))
      
            work = jobs[[x]]
            block = df.internal[work]
            
            Included = copy(rowEvent)
            Excluded = copy(rowEvent)
            Depth = copy(rowEvent)
            Coverage = copy(rowEvent)
            minDepth = copy(rowEvent)
            
            Up_Inc = rowEvent[get("EventType") %in% c("IR", "MXE", "SE")]
            Down_Inc = rowEvent[get("EventType") %in% c("IR", "MXE", "SE")]
            Up_Exc = rowEvent[get("EventType") %in% c("MXE")]        
            Down_Exc = rowEvent[get("EventType") %in% c("MXE")]        
                # for IR and SE, this defaults to rowEvent.Excluded
                
            junc_PSI = as.data.table(read.fst(
                file.path(se_output_path, "junc_PSI_index.fst")
            ))
            junc_counts = copy(junc_PSI)
            
            for(i in seq_len(length(work))) {
                junc = as.data.table(
                    read.fst(file.path(norm_output_path, "temp", 
                        paste(block$sample[i], "junc.fst.tmp", sep=".")))
                )
                junc[, c("start") := get("start") + 1]
                junc$strand = NULL

                junc = junc[junc.common, on = colnames(junc.common)[c(1,2,3)]]
                if(block$strand[i] == 0) {
                    junc$count = junc$total    
                } else if(block$strand[i] == -1) {
                    junc$count = 0
                    junc[get("strand") == "+", c("count") := get("neg")]
                    junc[get("strand") == "-", c("count") := get("pos")]
                    junc[get("strand") == "*", c("count") := get("total")]
                } else {
                    junc$count = 0
                    junc[get("strand") == "+", c("count") := get("pos")]
                    junc[get("strand") == "-", c("count") := get("neg")]
                    junc[get("strand") == "*", c("count") := get("total")]    
                }
                junc[is.na(get("count")), c("count") := 0]
                junc = junc[,c("seqnames", "start", "end", "strand", "Event", "count")]
                junc = cbind(junc, junc.common[, c("JG_up", "JG_down")])
                junc[, c("SO_L") := 0]
                junc[, c("SO_R") := 0]
                junc[, c("SO_I") := 0]
        
                # SpliceLeft and SpliceRight calculations
                junc[, c("SL") := sum(get("count")), 
                    by = c("seqnames", "start", "strand")]
                junc[, c("SR") := sum(get("count")), 
                    by = c("seqnames", "end", "strand")]
        
                # first overlap any junction that has non-same-island junctions
                junc[get("JG_up") != get("JG_down") & 
                    get("JG_up") != "" & get("strand") == "+", 
                    c("SO_L") := sum(get("count")), by = "JG_up"]
                junc[get("JG_up") != get("JG_down") & 
                    get("JG_down") != "" & get("strand") == "+", 
                    c("SO_R") := sum(get("count")), by = "JG_down"]
                junc[get("JG_up") != get("JG_down") & 
                    get("JG_up") != "" & get("strand") == "-",
                    c("SO_R") := sum(get("count")), by = "JG_up"]
                junc[get("JG_up") != get("JG_down") & 
                    get("JG_down") != "" & get("strand") == "-", 
                    c("SO_L") := sum(get("count")), by = "JG_down"]

                # message("Calculating SpliceOver for annotated IR events")
                        
                # Then use a simple overlap method to account for the remainder
                junc.subset = junc[get("JG_up") == get("JG_down") & 
                    get("JG_up") != "" & get("JG_down") != ""]
                junc.from = makeGRangesFromDataFrame(as.data.frame(junc.subset))
                junc.to = makeGRangesFromDataFrame(as.data.frame(junc))
                
                OL = findOverlaps(junc.from, junc.to)

                splice.overlaps.DT = data.table(from = from(OL), to = to(OL))
                splice.overlaps.DT[, 
                    c("count") := junc$count[to(OL)]]
                splice.overlaps.DT[, 
                    c("count_sum") := sum(get("count")), by = "from"]
                splice.summa = unique(splice.overlaps.DT[, c("from", "count_sum")])        

                junc.subset[splice.summa$from, c("SO_I") := splice.summa$count_sum]

                junc[junc.subset, on = c("Event"), c("SO_I") := get("i.SO_I")]
        
                # For annotated junctions, take SpliceOver as max of SpliceLeft, SpliceRight, or SpliceOver
                
                junc[get("SO_L") < get("SO_I"), c("SO_L") := get("SO_I")]
                junc[get("SO_R") < get("SO_I"), c("SO_R") := get("SO_I")]

                # Finally, for extreme cases, make SO_L = SL if underestimates
                junc[get("SO_L") < get("SL"), c("SO_L") := get("SL")]
                junc[get("SO_R") < get("SR"), c("SO_R") := get("SR")]

                junc[, c("SO_I") := NULL]

                splice = copy(Splice.Anno)
                
                splice[, c("count_Event1a") := 0]
                splice[!is.na(get("Event1a")), 
                    c("count_Event1a") := junc$count[
                        match(get("Event1a"), junc$Event)]]
                splice[is.na(get("count_Event1a")), 
                    c("count_Event1a") := 0]
                splice[, c("count_Event2a") := 0]
                splice[!is.na(get("Event2a")), 
                    c("count_Event2a") := junc$count[
                        match(get("Event2a"), junc$Event)]]
                splice[is.na(get("count_Event2a")), 
                    c("count_Event2a") := 0]
                splice[, c("count_Event1b") := 0]
                splice[!is.na(get("Event1b")), 
                    c("count_Event1b") := junc$count[
                        match(get("Event1b"), junc$Event)]]
                splice[is.na(get("count_Event1b")), 
                    c("count_Event1b") := 0]
                splice[, c("count_Event2b") := 0]
                splice[!is.na(get("Event2b")), 
                    c("count_Event2b") := junc$count[
                        match(get("Event2b"), junc$Event)]]
                splice[is.na(get("count_Event2b")), 
                    c("count_Event2b") := 0]

                splice[, c("count_JG_up") := 0]
                splice[!is.na(get("JG_up")) & get("strand") == "+", 
                    c("count_JG_up") := junc$SO_L[match(get("JG_up"), junc$JG_up)]]
                splice[!is.na(get("JG_up")) & get("strand") == "-", 
                    c("count_JG_up") := junc$SO_R[match(get("JG_up"), junc$JG_up)]]
                splice[is.na(get("count_JG_up")), c("count_JG_up") := 0]
                splice[, c("count_JG_down") := 0]
                splice[!is.na(get("JG_down")) & get("strand") == "-", 
                    c("count_JG_down") := junc$SO_L[match(get("JG_down"), junc$JG_down)]]
                splice[!is.na(get("JG_down")) & get("strand") == "+", 
                    c("count_JG_down") := junc$SO_R[match(get("JG_down"), junc$JG_down)]]
                splice[is.na(get("count_JG_down")), c("count_JG_down") := 0]

                # Splice participation: sum of two events compared to JG_up / JG_down
                splice[, c("partic_up") := 0]
                splice[, c("partic_down") := 0]

                splice[get("EventType") %in% c("MXE", "SE", "ALE", "A3SS"), 
                    c("partic_up") := get("count_Event1a") + get("count_Event1b")]
                splice[get("EventType") %in% c("MXE"), 
                    c("partic_down") := get("count_Event2a") + get("count_Event2b")]
                splice[get("EventType") %in% c("SE"), 
                    c("partic_down") := get("count_Event2a") + get("count_Event1b")]
                splice[get("EventType") %in% c("AFE", "A5SS"), 
                    c("partic_down") := get("count_Event1a") + get("count_Event1b")]

                # Splice coverage = participation / max_JG
                
                splice[, c("cov_up") := 0]
                splice[get("count_JG_up") > 0, 
                    c("cov_up") := get("partic_up") / get("count_JG_up")]
                splice[, c("cov_down") := 0]
                splice[get("count_JG_down") > 0, 
                    c("cov_down") := get("partic_down") / get("count_JG_down")]
                splice[get("EventType") %in% c("MXE", "SE") & 
                    get("cov_up") < get("cov_down"), 
                    c("coverage") := get("cov_up")]
                splice[get("EventType") %in% c("MXE", "SE") & 
                    get("cov_up") >= get("cov_down"), 
                    c("coverage") := get("cov_down")]
                splice[get("EventType") %in% c("ALE", "A3SS"), 
                    c("coverage") := get("cov_up")]
                splice[get("EventType") %in% c("AFE", "A5SS"), 
                    c("coverage") := get("cov_down")]
                
                irf = as.data.table(
                        read.fst(file.path(norm_output_path, "temp", 
                            paste(block$sample[i], "irf.fst.tmp", sep=".")))
                )

                irf[, c("start") := get("start") + 1]
                irf = irf[irf.common, on = colnames(irf.common)[seq_len(6)], 
                    c("EventRegion") := get("i.EventRegion")]
                
                # Extra statistics:
                irf[, c("SpliceMax") := 0]
                irf[get("SpliceLeft") >= get("SpliceRight"), 
                    c("SpliceMax") := get("SpliceLeft")]
                irf[get("SpliceLeft") < get("SpliceRight"), 
                    c("SpliceMax") := get("SpliceRight")]

                irf[junc, on = c("seqnames", "start", "end", "strand"), 
                    c("SpliceOverLeft") := get("SO_L")]
                irf[junc, on = c("seqnames", "start", "end", "strand"), 
                    c("SpliceOverRight") := get("SO_R")]
                irf[get("SpliceOverLeft") >= get("SpliceOverRight"), 
                    c("SpliceOverMax") := get("SpliceOverLeft")]
                irf[get("SpliceOverLeft") < get("SpliceOverRight"), 
                    c("SpliceOverMax") := get("SpliceOverRight")]
                
                irf[, c("IROratio") := 0]
                irf[get("IntronDepth") < 1 & get("IntronDepth") > 0 & 
                    (get("Coverage") + get("SpliceOverMax")) > 0, 
                    c("IROratio") := get("Coverage") / (
                        get("Coverage") + get("SpliceOverMax"))]
                irf[get("IntronDepth") >= 1, 
                    c("IROratio") := get("IntronDepth") / 
                        (get("IntronDepth") + get("SpliceOverMax"))]

                irf[, c("TotalDepth") := get("IntronDepth") + get("SpliceOverMax")]

                splice.no_region = splice[!(get("EventRegion") %in% irf$EventRegion)]
                splice.no_region[, 
                    c("Depth1a") := irf$TotalDepth[match(get("Event1a"), irf$EventRegion)]]
                splice.no_region[, 
                    c("Depth2a") := irf$TotalDepth[match(get("Event2a"), irf$EventRegion)]]
                splice.no_region[, 
                    c("Depth1b") := irf$TotalDepth[match(get("Event1b"), irf$EventRegion)]]
                splice.no_region[, 
                    c("Depth2b") := irf$TotalDepth[match(get("Event2b"), irf$EventRegion)]]
                splice.no_region[, c("Depth") := 0]
                splice.no_region[get("count_JG_up") > get("count_JG_down"), 
                    c("Depth") := get("count_JG_up")]
                splice.no_region[get("count_JG_up") <= get("count_JG_down"), 
                    c("Depth") := get("count_JG_down")]
                splice.no_region[is.na(get("Depth1a")), c("Depth1a") := 0]
                splice.no_region[is.na(get("Depth1b")), c("Depth1b") := 0]
                splice.no_region[is.na(get("Depth2a")), c("Depth2a") := 0]
                splice.no_region[is.na(get("Depth2b")), c("Depth2b") := 0]
                splice.no_region[get("Depth1a") > get("Depth2a"), 
                    c("DepthA") := get("Depth1a")]
                splice.no_region[get("Depth1b") > get("Depth2b"), 
                    c("DepthB") := get("Depth1b")]
                splice.no_region[get("Depth1a") <= get("Depth2a"), 
                    c("DepthA") := get("Depth2a")]
                splice.no_region[get("Depth1b") <= get("Depth2b"), 
                    c("DepthB") := get("Depth2b")]
                splice.no_region[get("DepthA") > get("DepthB"), 
                    c("Depth") := get("DepthA")]
                splice.no_region[get("DepthA") <= get("DepthB"), 
                    c("Depth") := get("DepthB")]

                splice[, c("TotalDepth") := 0]
                splice[irf, on = "EventRegion", 
                    c("TotalDepth") := get("i.TotalDepth")]
                splice[splice.no_region, on = "EventName", 
                    c("TotalDepth") := get("i.Depth")]
                        
                file.remove(file.path(norm_output_path, "temp", 
                    paste(block$sample[i], "junc.fst.tmp", sep=".")))
                file.remove(file.path(norm_output_path, "temp", 
                    paste(block$sample[i], "irf.fst.tmp", sep=".")))
                
                # Do BuildSE here
                setnames(irf, "Name", "EventName")
                # Included
                Included[, c(block$sample[i]) := c(
                    irf$IntronDepth, 
                    0.5 * (splice$count_Event1a[splice$EventType %in% c("SE", "MXE")] + 
                        splice$count_Event2a[splice$EventType %in% c("SE", "MXE")]),
                    splice$count_Event1a[!splice$EventType %in% c("SE", "MXE")]
                )]
                
                if(IRMode == "SpliceOverMax") {
                    Excluded[, c(block$sample[i]) := c(
                        irf$SpliceOverMax,
                        0.5 * (splice$count_Event1b[splice$EventType %in% c("MXE")] + 
                            splice$count_Event2b[splice$EventType %in% c("MXE")]),
                        splice$count_Event1b[!splice$EventType %in% c("MXE")]
                    )]
                } else {
                    Excluded[, c(block$sample[i]) := c(
                        irf$SpliceMax,
                        0.5 * (splice$count_Event1b[splice$EventType %in% c("MXE")] + 
                            splice$count_Event2b[splice$EventType %in% c("MXE")]),
                        splice$count_Event1b[!splice$EventType %in% c("MXE")]
                    )]      
                }

                # Validity checking for IR, MXE, SE
                irf[get("strand") == "+", c("Up_Inc") := get("ExonToIntronReadsLeft")]
                irf[get("strand") == "-", c("Up_Inc") := get("ExonToIntronReadsRight")]
                irf[get("strand") == "+", c("Down_Inc") := get("ExonToIntronReadsRight")]
                irf[get("strand") == "-", c("Down_Inc") := get("ExonToIntronReadsLeft")]
                
                Up_Inc[, c(block$sample[i]) := 
                    c(irf$Up_Inc, 
                        splice$count_Event1a[splice$EventType %in% c("MXE", "SE")])]
                Down_Inc[, c(block$sample[i]) := 
                    c(irf$Down_Inc, 
                        splice$count_Event2a[splice$EventType %in% c("MXE", "SE")])]
                
                Up_Exc[, c(block$sample[i]) := 
                    splice$count_Event1b[splice$EventType %in% c("MXE")]]
                Down_Exc[, c(block$sample[i]) := 
                    splice$count_Event2b[splice$EventType %in% c("MXE")]]
                
                Depth[, c(block$sample[i]) := c(irf$TotalDepth, splice$TotalDepth)]
                Coverage[, c(block$sample[i]) := c(irf$Coverage, splice$coverage)]
                
                splice[get("EventType") %in% c("MXE", "SE") & 
                    get("cov_up") < get("cov_down"), 
                    c("minDepth") := get("count_JG_up")]
                splice[get("EventType") %in% c("MXE", "SE") & 
                    get("cov_up") >= get("cov_down"), 
                    c("minDepth") := get("count_JG_down")]
                splice[get("EventType") %in% c("ALE", "A3SS"), 
                    c("minDepth") := get("count_JG_up")]
                splice[get("EventType") %in% c("AFE", "A5SS"), 
                    c("minDepth") := get("count_JG_down")]
                
                minDepth[, c(block$sample[i]) := c(
                    irf$IntronDepth,
                    splice$minDepth)]

                junc[get("count") == 0, c("PSI") := 0]
                junc[get("SO_L") > get("SO_R"), 
                    c("PSI") := get("count") / get("SO_L")]
                junc[get("SO_R") >= get("SO_L") & get("SO_R") > 0, 
                    c("PSI") := get("count") / get("SO_R")]
                
                junc_PSI[junc, on = c("seqnames", "start", "end", "strand"),
                    c(block$sample[i]) := get("i.PSI")]
                junc_counts[junc, on = c("seqnames", "start", "end", "strand"),
                    c(block$sample[i]) := get("i.count")]
                    
            } # end FOR loop
            
            if(low_memory_mode == TRUE) {
                # Write BuildSE temp files
                value = t(as.matrix(Included[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Included", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(Excluded[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Excluded", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(Depth[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Depth", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(Coverage[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Coverage", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(minDepth[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("minDepth", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(Up_Inc[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Up_Inc", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(Down_Inc[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Down_Inc", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(Up_Exc[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Up_Exc", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(Down_Exc[, -c(1,2,3)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("Down_Exc", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)

                value = t(as.matrix(junc_PSI[, -c(1,2,3,4)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("junc_PSI", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
                value = t(as.matrix(junc_counts[, -c(1,2,3,4)]))
                fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
                    paste("junc_counts", as.character(x), "txt.gz", sep=".")), 
                    col.names = FALSE, row.names = FALSE)
            } else {
                final = list(
                    Included = Included[, -c(1,2,3)],
                    Excluded = Excluded[, -c(1,2,3)],
                    Depth = Depth[, -c(1,2,3)],
                    Coverage = Coverage[, -c(1,2,3)],
                    minDepth = minDepth[, -c(1,2,3)],
                    Up_Inc = Up_Inc[, -c(1,2,3)],
                    Down_Inc = Down_Inc[, -c(1,2,3)],
                    Up_Exc = Up_Exc[, -c(1,2,3)],
                    Down_Exc = Down_Exc[, -c(1,2,3)],
                    junc_PSI = junc_PSI[, -c(1,2,3,4)],
                    junc_counts = junc_counts[, -c(1,2,3,4)]
                )                
            }
            
        }, df.internal = df.internal, jobs = jobs, 
            norm_output_path = norm_output_path, BPPARAM = BPPARAM_mod
    ))
    gc()

    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(0.20, message = "Building Final SummarizedExperiment Object")
    }  
    message("Building Final SummarizedExperiment Object")
    
    if(low_memory_mode) {
        for(item in item.todo) {
            file.DT = data.table(file = list.files(pattern = item, 
                path = file.path(norm_output_path, "temp")))
            file.DT[, c("index") := 
                as.numeric(tstrsplit(file.DT, split=".", fixed = TRUE)[[2]])]
            setorder(file.DT, "index")
            mat = NULL
            for(x in seq_len(n_jobs)) {
                temp = t(fread(file.path(file.path(norm_output_path, "temp"), 
                    file.DT$file[x]), data.table = FALSE))
                colnames(temp) = df.internal$sample[jobs[[x]]]
                mat = cbind(mat, temp)
                file.remove(file.path(file.path(norm_output_path, "temp"), 
                    file.DT$file[x]))
            }
            outfile = file.path(se_output_path, paste(item, "fst", sep="."))
            write.fst(as.data.frame(mat), outfile)
        }
    } else {
        item.DTList = list()
        for(item in item.todo) {
            for(x in seq_len(n_jobs)) {
                if(x == 1) {
                    item.DTList[[item]] = agg.list[[x]][[item]]
                } else {
                    item.DTList[[item]] = cbind(item.DTList[[item]], agg.list[[x]][[item]])
                }
            }
            outfile = file.path(se_output_path, paste(item, "fst", sep="."))
            write.fst(as.data.frame(item.DTList[[item]]), outfile)
        }
    }
  
  # Rewrite junc_PSI and junc_count by adding in NxtIRF rownames (as first column)
    junc_index = fst::read.fst(file.path(
        se_output_path, "junc_PSI_index.fst"
    ))
    junc_PSI = fst::read.fst(file.path(
        se_output_path, "junc_PSI.fst"
    ))  
    junc_counts = fst::read.fst(file.path(
        se_output_path, "junc_counts.fst"
    ))
    junc_PSI$rownames = with(junc_index, 
        paste0(seqnames, ":", start, "-", end, "/", strand))
    junc_counts$rownames = junc_PSI$rownames
    fst::write.fst(cbind(junc_PSI[,ncol(junc_PSI),drop=FALSE],
        junc_PSI[,-ncol(junc_PSI)]),file.path(
        se_output_path, "junc_PSI.fst"
    ))
    fst::write.fst(cbind(junc_counts[,ncol(junc_counts),drop=FALSE],
        junc_counts[,-ncol(junc_counts)]),file.path(
        se_output_path, "junc_counts.fst"
    ))
    
    outfile = file.path(se_output_path, paste("stats", "fst", sep="."))
    write.fst(as.data.frame(df.internal), outfile)
    
    # Create barebones colData.Rds - save coverage files as well
    if(length(coverage_files) == nrow(df.internal) & IsCOV(coverage_files)) {
        df.files = data.table(
            sample = df.internal$sample,
            bam_file = "",
            irf_file = df.internal$path,
            cov_file = coverage_files
        )
    } else {
        df.files = data.table(
            sample = df.internal$sample,
            bam_file = "",
            irf_file = df.internal$path,
            cov_file = ""
        )
    }
    # 
    colData = list(
        df.files = df.files,
        df.anno = data.table(sample = df.internal$sample)
    )
    saveRDS(colData, file.path(se_output_path, "colData.Rds"))
    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(0.19, message = "NxtIRF Collation Finished")
    }  
    message("NxtIRF Collation Finished")
}


#' Constructs a SummarizedExperiment object from the collated data
#'
#' This function creates a SummarizedExperiment object from the data collated
#' from IRFinder output using CollateData().
#'
#' @param fst_path The output path given to CollateData() pointing to the
#'   collated data
#' @param colData A data frame containing the sample annotation information. Note
#'   that the first column must contain the sample names. If the names of only a subset
#'   of samples are given, then `MakeSE()` will construct the SE object based only on
#'   the samples given. Omit `colData` to generate an SE object based on the whole dataset.
#'   The colData can be set later using `SummarizedExperiment::colData()`
#' @param RemoveOverlapping (default = TRUE) Whether to filter out overlapping introns
#'   of IR events belonging to minor isoforms. MakeSE will try to identify which junctions
#'   belong to major isoforms, then select the junctions from non-overlapping minor isoforms
#'   in an iterative approach, until no non-overlapping introns remain. This is important
#'   to make sure IR events are not 'double-counted'
#'
#' @return A SummarizedExperiment object containing the following assays:
#'   "Included", "Excluded" : contains the IntronDepth and SpliceJunction counts of IR events.
#'   For alternative splicing events, they represent counts of Included and Excluded isoforms.
#'   For mutually exclusive events, Included represents counts for selection of the upstream
#'   cassette exon. For alternative 5' or 3' events and for alternate first / last exons, 
#'   Included represents splicing of the shorter splice junction.
#'   Other assays for advanced users include:
#'   "Depth": Sum of included and excluded counts
#'   "Coverage": Fraction of the intron covered by bases
#'   "MinDepth": Intron depth for IR events, or the count of the minor isoform for AS events
#'
#' @export
MakeSE = function(fst_path, colData, RemoveOverlapping = TRUE) {
    # Includes iterative filtering for IR events with highest mean PSI
        # To annotate IR events of major isoforms

    item.todo = c("rowEvent", "Included", "Excluded", "Depth", "Coverage", 
        "minDepth", "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")
    files.todo = file.path(normalizePath(fst_path), paste(item.todo, "fst", sep="."))
    assert_that(all(file.exists(files.todo)),
    msg = "FST File generation appears incomplete. Suggest run CollateData() again")

    assert_that(file.exists(file.path(fst_path, "colData.Rds")),
        msg = "colData.Rds does not exist in given path")
    colData.Rds = readRDS(file.path(fst_path, "colData.Rds"))
    assert_that("df.anno" %in% names(colData.Rds),
        msg = "colData.Rds must contain df.anno containing annotations")

    if(missing(colData)) {    
        colData = colData.Rds$df.anno
    } else {
        assert_that("sample" %in% colnames(colData),
            msg = "'sample' must be a column name in colData containing sample names")
        assert_that(all(colData$sample %in% colData.Rds$df.anno$sample),
            msg = "some samples in colData were not found")
    }
    colData = as.data.frame(colData)
    
    remove_na = NULL
    if(ncol(colData) > 1) {
        for(i in seq(2, ncol(colData))) {
            if(is(colData[,i], "character")) {
                colData[,i] = factor(unlist(colData[,i]))      
            } else if(is(colData[,i], "logical")) {
                colData[,i] <- factor(unlist(ifelse(colData[,i], "TRUE","FALSE")))                
            } else if(all(is.na(unlist(colData[,i])))) {
                remove_na = append(remove_na, i)
            }
        }
    }
    if(!is.null(remove_na)) {
        colData = colData[,-remove_na]
    }
  
    rowData = read.fst(files.todo[1])
    Included = as.matrix(read.fst(files.todo[2], columns = colData$sample))
    Excluded = as.matrix(read.fst(files.todo[3], columns = colData$sample))
    Depth = as.matrix(read.fst(files.todo[4], columns = colData$sample))
    Coverage = as.matrix(read.fst(files.todo[5], columns = colData$sample))
    minDepth = as.matrix(read.fst(files.todo[6], columns = colData$sample))
    Up_Inc = as.matrix(read.fst(files.todo[7], columns = colData$sample))
    Down_Inc = as.matrix(read.fst(files.todo[8], columns = colData$sample))
    Up_Exc = as.matrix(read.fst(files.todo[9], columns = colData$sample))
    Down_Exc = as.matrix(read.fst(files.todo[10], columns = colData$sample))

    rownames(Up_Inc) = rowData$EventName[rowData$EventType %in% c("IR", "MXE", "SE")]
    rownames(Down_Inc) = rowData$EventName[rowData$EventType %in% c("IR", "MXE", "SE")]
    rownames(Up_Exc) = rowData$EventName[rowData$EventType %in% c("MXE")]
    rownames(Down_Exc) = rowData$EventName[rowData$EventType %in% c("MXE")]
  
    # Annotate NMD direction
    rowData = as.data.table(rowData)
    rowData[, get("NMD_direction") := 0]
    rowData[get("Inc_Is_NMD") & !get("Exc_Is_NMD"), c("NMD_direction") := 1]
    rowData[!get("Inc_Is_NMD") & get("Exc_Is_NMD"), c("NMD_direction") := -1]
    rowData = as.data.frame(rowData)

    se = SummarizedExperiment(
        assays = SimpleList(
            Included = Included, Excluded = Excluded, 
            Depth = Depth, Coverage = Coverage, minDepth = minDepth
        ),
        rowData = rowData, colData = as.data.frame(colData[, -1, drop=FALSE], 
        row.names = colData$sample)
    )
    rownames(se) = rowData(se)$EventName

    metadata(se)$Up_Inc = Up_Inc
    metadata(se)$Down_Inc = Down_Inc
    metadata(se)$Up_Exc = Up_Exc
    metadata(se)$Down_Exc = Down_Exc

    if("df.files" %in% names(colData.Rds) &&
        "cov_file" %in% colnames(colData.Rds$df.files)) {
        metadata(se)$cov_file = colData.Rds$df.files$cov_file
        names(metadata(se)$cov_file) = colData.Rds$df.files$sample       
    }

    if(RemoveOverlapping == TRUE) {
        # Iterative filtering of IR
        tryCatch({
            junc_PSI = fst::read.fst(file.path(
                normalizePath(fst_path), "junc_PSI.fst"
            ))
            rownames(junc_PSI) = junc_PSI$rownames
            junc_PSI = junc_PSI[,-1,drop=FALSE]

            se.IR = se[rowData(se)$EventType == "IR",,drop = FALSE]
            se.IR = se.IR[rowData(se.IR)$EventRegion %in% rownames(junc_PSI),
                ,drop = FALSE]
            junc_PSI = junc_PSI[rowData(se.IR)$EventRegion,, drop = FALSE]

            if(nrow(se.IR) > 0) {
                message("Iterating through IR events to determine introns of main isoforms")

                se.IR.gr = NxtIRF.CoordToGR(rownames(junc_PSI))
                se.IR.gr.reduced = reduce(se.IR.gr)

                OL = findOverlaps(se.IR.gr, se.IR.gr.reduced)
                junc_PSI.group = as.data.table(junc_PSI)
                junc_PSI.group$group = to(OL)
                junc_PSI.group$means = rowMeans(junc_PSI)
                junc_PSI.group[, c("max_means") := max(get("means")), 
                    by = "group"]

                se.IR.final = se.IR[
                    junc_PSI.group$means == junc_PSI.group$max_means,,drop = FALSE]
                se.IR.excluded = se.IR[
                    junc_PSI.group$means != junc_PSI.group$max_means,,drop = FALSE]

                final.gr = NxtIRF.CoordToGR(rowData(se.IR.final)$EventRegion)
                excluded.gr = NxtIRF.CoordToGR(rowData(se.IR.excluded)$EventRegion)

                OL = findOverlaps(excluded.gr, final.gr)
                include = which(!(seq_len(length(excluded.gr))) %in% sort(unique(from(OL))))

                # Iteration to find events not overlapping with se.IR.final

                iteration = 0
                while(length(include) > 0 & length(final.gr) > 0) {
                    iteration = iteration + 1
                    message(paste("Iteration", iteration))
                    se.IR.excluded = se.IR.excluded[include,,drop = FALSE]
                    junc_PSI = junc_PSI[rowData(se.IR.excluded)$EventRegion,, drop = FALSE]
                    se.IR.gr = NxtIRF.CoordToGR(rownames(junc_PSI))
                    se.IR.gr.reduced = reduce(se.IR.gr)

                    OL = findOverlaps(se.IR.gr, se.IR.gr.reduced)
                    junc_PSI.group = as.data.table(junc_PSI)
                    junc_PSI.group$group = to(OL)
                    junc_PSI.group$means = rowMeans(junc_PSI)
                    junc_PSI.group[,  c("max_means") := max(get("means")), 
                        by = "group"]

                    if(length(which(junc_PSI.group$means == junc_PSI.group$max_means)) > 0) {
                        se.IR.final = rbind(se.IR.final, 
                            se.IR.excluded[junc_PSI.group$means == junc_PSI.group$max_means,
                                ,drop = FALSE])
                        se.IR.excluded = 
                            se.IR.excluded[junc_PSI.group$means != junc_PSI.group$max_means,
                                ,drop = FALSE]

                        final.gr = NxtIRF.CoordToGR(rowData(se.IR.final)$EventRegion)
                        excluded.gr = NxtIRF.CoordToGR(rowData(se.IR.excluded)$EventRegion)

                        OL = findOverlaps(excluded.gr, final.gr)
                        include = which(!(seq_len(length(excluded.gr))) %in% sort(unique(from(OL))))
                    } else {
                        final.gr = c()
                        include = c()
                    }
                }
            }

            se = rbind(se.IR[rownames(se.IR) %in% rownames(se.IR.final),,drop = FALSE],
                se[rowData(se)$EventType != "IR",,drop = FALSE])

        }, error = function(e) {
            message(paste("Iterative filtering of IR appears to have run into an error.",
                "Using RemoveOverlapping = FALSE"))
        })
    }
  
    # Encapsulate as NxtSE object
    se = as(se, "NxtSE")
    return(se)
}

#' Filtering for IR and Alternative Splicing Events
#'
#' This function implements filtering of IR or AS events based on customisable criteria
#' 
#' @details
#'   \strong{Annotation Filters}\cr\cr 
#'     \strong{Protein_Coding}: Filters for alternative splicing or IR events within protein reading
#'       frames. No additional parameters required.\cr \cr 
#'     \strong{NMD_Switching}: Filters for events in which one isoform is a predicted NMD substrate.\cr \cr 
#'     \strong{Transcript_Support_Level}: filters for events in which both isoforms have a TSL level
#'       below or equal to filterVars$minimum\cr \cr 
#'   \strong{Data Filters}\cr\cr 
#'     \strong{Depth}: Filters IR or alternative splicing events of transcripts that are "expressed"
#'       with adequate \code{Depth} as calculated by the sum of all splicing and IR reads spanning the event.
#'       Events with \code{Depth} below filterVars$minimum are excluded\cr\cr 
#'     \strong{Coverage}: Coverage means different things to IR and alternative splicing.\cr\cr 
#'       For \emph{IR}, Coverage refers to the percentage of the measured intron covered with
#'       reads. Introns of samples with an IntronDepth above \code{filterVars$minDepth} are
#'       assessed, with introns with coverage below \code{filterVars$minimum} are excluded.\cr\cr 
#'       For \emph{Alternative Splicing}, Coverage refers to the percentage of all splicing events
#'       observed across the genomic region that is compatible with either the included
#'       or excluded event. This prevents NxtIRF from doing differential analysis
#'       between two minor isoforms. Instead of IntronDepth, in AS events NxtIRF considers
#'       events where the spliced reads from both exonic regions exceed \code{filterVars$minDepth}.
#'       Then, events with a splicing coverage below \code{filterVars$minimum} are excluded.
#'       We recommend testing IR events for > 90% coverage and AS events for > 60% coverage
#'       as given in the default filters which can be accessed using \code{get_default_filters()}\cr\cr  
#'     \strong{Consistency}: Skipped exons (SE) and mutually exclusive exons (MXE) comprise reads of
#'       two contiguous splice junctions (for the included casette exon). Summating counts 
#'       from both junctions is misleading as there may be overlapping events (e.g. alternate
#'       first / last exons) that only rely on one splice event. To ensure the SE / MXE is the
#'       dominant event, we require both splice junctions to have comparable counts.\cr\cr
#'       Events are excluded if either of the upstream or downstream
#'       event is lower than total splicing events by a log-2 magnitude above filterVars$maximum.
#'       For example, if \code{filterVars$maximum = 2}, we require both upstream and downstream events
#'       to represent at least 1/(2^2) = 1/4 of the sum of upstream and downstream event.
#'       This is considered for each isoform of each event, as long as the total counts belonging
#'       to the considered isoform is above \code{filterVars$minDepth}.
#'
#'  We highly recommend using the default filters, which can be acquired using \code{get_default_filters()}
#' 
#' @param filterClass One of either \code{Annotation} or \code{Data}
#' @param filterType For filterClass \code{Annotation}, either one of 
#'   \code{Protein_Coding}, \code{NMD_Switching}, \code{Transcript_Support_Level}.
#'   For filterClass \code{Data}, either one of
#'   \code{Depth}, \code{Coverage}, \code{Consistency}.
#' @param filterVars A list of parameters, as explained below
#' @param filterObject the SummarizedExperiment to filter
#'
#' @return A vector of type \code{logical} designating which events to retain \code{TRUE} and which
#'   to remove \code{FALSE}.
#'
#' @export
runFilter <- function(filterClass, filterType, filterVars, filterObject) {
    # filterClass: can be one of 'Annotation', 'Data', 'Runtime'
    # filterType:
    # - Annotation:
    # - Data:
    # -     Depth: 1-minimum, 2-minCond, 3-pcTRUE
    # -     Coverage: 1-minimum, 1a-minDepth, 2-minCond, 3-pcTRUE
    # -     Consistency: 1-maximum, 1a-minDepth, 2-minCond, 3-pcTRUE
    # - for Consistency, maximum is the max(abs(log2_delta)) between comparison and calculated value

    filterResult = rep(TRUE, nrow(filterObject))
    usePC = ifelse("pcTRUE" %in% names(filterVars),
        filterVars$pcTRUE, 100)
    minDepth = ifelse("minDepth" %in% names(filterVars),
        filterVars$minDepth, 0)

    rowData = as.data.frame(rowData(filterObject))
    colData = as.data.frame(colData(filterObject))
    use_cond = ifelse(
            !is.null(names(filterVars)) && 
            "condition" %in% names(filterVars) && 
            filterVars$condition %in% colnames(colData), 
        TRUE, FALSE)
    if(filterClass == "Data") {
        if(filterType == "Depth") {
            message("Running Depth filter")
            minimum = ifelse("minimum" %in% names(filterVars),
                filterVars$minimum, 20)
            if(use_cond == TRUE) {
                cond_vec = unlist(colData[, which(colnames(colData) == filterVars$condition)])
                cond_vars = unique(cond_vec)
            }
            depth = as.matrix(assay(filterObject, "Depth"))
            sum_res = rep(0, nrow(filterObject))
            if(use_cond == TRUE) {
                for(cond in cond_vars) {
                    depth.subset = depth[, which(cond_vec == cond)]
                    sum = rowSums(depth.subset >= minimum)
                    sum_res = sum_res + 
                        ifelse(sum * 100 / ncol(depth.subset) >= usePC, 1, 0)
                }
                n_TRUE = ifelse(!is.na(suppressWarnings(as.numeric(filterVars$minCond))), 
                    as.numeric(filterVars$minCond), -1)
                if(n_TRUE == -1) n_TRUE = length(cond_vars)
                res = (sum_res >= n_TRUE)
            } else {
                sum = rowSums(depth >= minimum)
                res = ifelse(sum * 100 / ncol(depth) >= usePC, TRUE, FALSE)
            }
            if("EventTypes" %in% names(filterVars)) {
                res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
            }
            return(res)
        } else if(filterType == "Coverage") {
            message("Running Coverage filter")
            minimum = ifelse("minimum" %in% names(filterVars),
                filterVars$minimum, 20)
            if(use_cond == TRUE) {
                cond_vec = unlist(colData[, which(colnames(colData) == filterVars$condition)])
                cond_vars = unique(cond_vec)
            }
            cov = as.matrix(assay(filterObject, "Coverage"))
            depth = as.matrix(assay(filterObject, "minDepth"))
            
            # do not test if depth below threshold
            cov[depth < minDepth] = 1    

            sum_res = rep(0, nrow(filterObject))
            if(use_cond == TRUE) {
                for(cond in cond_vars) {
                    cov.subset = cov[, which(cond_vec == cond)]
                    sum = rowSums(cov.subset >= minimum / 100)
                    sum_res = sum_res + 
                        ifelse(sum * 100 / ncol(cov.subset) >= usePC, 1, 0)
                }
                n_TRUE = ifelse(!is.na(suppressWarnings(as.numeric(filterVars$minCond))), 
                    as.numeric(filterVars$minCond), -1)
                if(n_TRUE == -1) n_TRUE = length(cond_vars)
                res = (sum_res >= n_TRUE)
            } else {
                sum = rowSums(cov >= minimum / 100)
                res = ifelse(sum * 100 / ncol(cov) >= usePC, TRUE, FALSE)
            }
            if("EventTypes" %in% names(filterVars)) {
                res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
            }
            return(res)
        } else if(filterType == "Consistency") {    # requires: 
            message("Running Consistency filter")
            maximum = ifelse("maximum" %in% names(filterVars),
                filterVars$maximum, 1)

            if(use_cond == TRUE) {
                cond_vec = unlist(colData[, which(colnames(colData) == filterVars$condition)])
                cond_vars = unique(cond_vec)
            }
            Up_Inc = as.matrix(S4Vectors::metadata(filterObject)$Up_Inc)[
                rowData(filterObject)$EventName[
                    rowData(filterObject)$EventType %in% c("IR", "MXE", "SE")
                ],
            ]
            Down_Inc = as.matrix(S4Vectors::metadata(filterObject)$Down_Inc)[
                rowData(filterObject)$EventName[
                    rowData(filterObject)$EventType %in% c("IR", "MXE", "SE")
                ],
            ]
            IntronDepth = as.matrix(assay(filterObject, "Included"))
            IntronDepth = IntronDepth[rowData$EventType %in% c("IR", "MXE", "SE"),]
            minDepth.Inc = Up_Inc + Down_Inc
            # do not test if depth below threshold
            Up_Inc[minDepth.Inc < minDepth] = IntronDepth[minDepth.Inc < minDepth]
            Down_Inc[minDepth.Inc < minDepth] = IntronDepth[minDepth.Inc < minDepth]

            Excluded = as.matrix(assay(filterObject, "Excluded"))
            Excluded = Excluded[rowData$EventType %in% c("MXE"),]
            Up_Exc = as.matrix(S4Vectors::metadata(filterObject)$Up_Exc)
            Down_Exc = as.matrix(S4Vectors::metadata(filterObject)$Down_Exc)
            minDepth.Exc = Up_Exc + Down_Exc
            # do not test if depth below threshold
            Up_Exc[minDepth.Exc < minDepth] = Excluded[minDepth.Exc < minDepth]    
            Down_Exc[minDepth.Exc < minDepth] = Excluded[minDepth.Exc < minDepth]

            sum_res = rep(0, nrow(filterObject))
            if(use_cond == TRUE) {
                for(cond in cond_vars) {
                    Up_Inc.subset = Up_Inc[, which(cond_vec == cond)]
                    Down_Inc.subset = Down_Inc[, which(cond_vec == cond)]
                    IntronDepth.subset = IntronDepth[, which(cond_vec == cond)]
                    Up_Exc.subset = Up_Exc[, which(cond_vec == cond)]
                    Down_Exc.subset = Down_Exc[, which(cond_vec == cond)]
                    Excluded.subset = Excluded[, which(cond_vec == cond)]

                    sum_inc = rowSums(
                        abs(log2(Up_Inc.subset + 1) - log2(IntronDepth.subset +1)) < maximum &
                        abs(log2(Down_Inc.subset + 1) - log2(IntronDepth.subset +1)) < maximum
                    )
                    sum_exc = rowSums(
                        abs(log2(Up_Exc.subset + 1) - log2(Excluded.subset +1)) < maximum &
                        abs(log2(Down_Exc.subset + 1) - log2(Excluded.subset +1)) < maximum
                    )
                    sum_inc = c(sum_inc, rep(ncol(Up_Inc.subset), 
                        sum(!(rowData$EventType %in% c("IR", "MXE", "SE")))))
                    sum_exc = c(rep(ncol(Up_Inc.subset), 
                        sum(rowData$EventType == "IR")),
                    sum_exc, rep(ncol(Up_Inc.subset), 
                        sum(!(rowData$EventType %in% c("IR", "MXE")))))
                    sum = 0.5 * (sum_inc + sum_exc)
                    sum_res = sum_res + 
                        ifelse(sum * 100 / ncol(Up_Inc.subset) >= usePC, 1, 0)
                }
                n_TRUE = ifelse(!is.na(suppressWarnings(as.numeric(filterVars$minCond))), 
                    as.numeric(filterVars$minCond), -1)
                if(n_TRUE == -1) n_TRUE = length(cond_vars)
                res = (sum_res >= n_TRUE)
            } else {
                sum_inc = rowSums(
                    abs(log2(Up_Inc + 1) - log2(IntronDepth +1)) < filterVars$maximum &
                    abs(log2(Down_Inc + 1) - log2(IntronDepth +1)) < filterVars$maximum
                )
                sum_exc = rowSums(
                    abs(log2(Up_Exc + 1) - log2(Excluded +1)) < filterVars$maximum &
                    abs(log2(Down_Exc + 1) - log2(Excluded +1)) < filterVars$maximum
                )
                sum_inc = c(sum_inc, rep(ncol(Up_Inc), 
                    sum(!(rowData$EventType %in% c("IR", "MXE", "SE")))))
                sum_exc = c(
                    rep(ncol(Up_Inc), sum(rowData$EventType == "IR")),
                    sum_exc, 
                    rep(
                        ncol(Up_Inc), 
                        sum(!(rowData$EventType %in% c("IR", "MXE")))
                    )
                )
                sum = 0.5 * (sum_inc + sum_exc)
                res = ifelse(sum * 100 / ncol(Up_Inc) >= usePC, TRUE, FALSE)
            }
            if("EventTypes" %in% names(filterVars)) {
                res[!(rowData(filterObject)$EventType %in% filterVars$EventTypes)] = TRUE
            }
            return(res)
        }
    } else if(filterClass == "Annotation") {
        if(filterType == "Protein_Coding") {
            # returns if any of included or excluded is protein_coding
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[
                get("Inc_Is_Protein_Coding") == TRUE | 
                get("Exc_Is_Protein_Coding") == TRUE]
            rowSelected = rowSelected[get("EventType") != "IR" | 
                get("Inc_Is_Protein_Coding") == TRUE] # filter for CDS introns
            res = rowData$EventName %in% rowSelected$EventName
        } else if(filterType == "NMD_Switching") {
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[!is.na(get("Inc_Is_NMD")) & !is.na(get("Exc_Is_NMD"))]
            rowSelected = rowSelected[get("Inc_Is_NMD") != get("Exc_Is_NMD")]
            res = rowData$EventName %in% rowSelected$EventName
        } else if(filterType == "Transcript_Support_Level") {
            if(!("minimum" %in% names(filterVars))) {
                minimum = 1
            } else {
                minimum = filterVars$minimum
            }
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[get("Inc_TSL") != "NA" & get("Exc_TSL") != "NA"]
            rowSelected[, c("Inc_TSL") := as.numeric(get("Inc_TSL"))]
            rowSelected[, c("Exc_TSL") := as.numeric(get("Exc_TSL"))]
            rowSelected = rowSelected[get("Inc_TSL") <= minimum & get("Exc_TSL") <= minimum]
            res = rowData$EventName %in% rowSelected$EventName
        }
        if("EventTypes" %in% names(filterVars)) {
            res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
        }
        return(res)
    } else {
        return(filterResult)
    }
}

#' Convenience function to apply a list of filters to a SummarizedExperiment object
#'
#' See `?runFilter` for details regarding filters
#' 
#' @param se A SummarizedExperiment object created by `MakeSE()`
#' @param filters A list of filters to apply. Each filter must contain the elements
#'   `filterClass`, `filterType` and `filterVars`. See `?runFilter` for details
#' @return A vector of logicals, with `TRUE` indicating events to be retained, and
#'   `FALSE` for events to be filtered out
#' @md
#' @export
apply_filters <- function(se, filters) {
    # filters are a list of filters to apply on se
    # returns a vector of TRUE / FALSE
    # a filtered se can be made using:
    #       se.filtered = se[apply_filters(se, filters),]
   
    assert_that(is(filters, "list"), msg = "filters must be a list")
    for(i in seq_len(length(filters))) {
        assert_that("filterVars" %in% names(filters[[i]]),
            msg = paste("filterVars is missing from filters @ index #", i))
        assert_that("filterClass" %in% names(filters[[i]]),
            msg = paste("filterClass is missing from filters @ index #", i))
        assert_that("filterType" %in% names(filters[[i]]),
            msg = paste("filterType is missing from filters @ index #", i))
    }
    assert_that(is(se, "NxtSE"), 
        msg = "se must be a NxtSE object")

    filterSummary = rep(TRUE, nrow(se))
    for(i in seq_len(length(filters))) {
        filterSummary = filterSummary & runFilter(
            filters[[i]]$filterClass,
            filters[[i]]$filterType,
            filters[[i]]$filterVars,
            se
        )
    }
    
    return(filterSummary)
}