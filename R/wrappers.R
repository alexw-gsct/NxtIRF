# wrappers to R/C++

run_IRFinder_multithreaded = function(
    reference_path = "./Reference", 
    bamfiles = "Unsorted.bam", 
    output_files = "./Sample",
    max_threads = max(parallel::detectCores() - 2, 1),
    Use_OpenMP = TRUE,
    run_featureCounts = FALSE,
    localHub = FALSE,
    ah = AnnotationHub(localHub = localHub)
        ) {

    s_bam = normalizePath(bamfiles)
    s_ref = normalizePath(reference_path)
    
    assert_that(max_threads == 1 || max_threads <= parallel::detectCores() - 1,
        msg = paste(max_threads, " threads is not allowed for this system"))
    
    assert_that(all(file.exists(s_bam)),
        msg = paste(
            paste(
                unique(s_bam[!file.exists(s_bam)]),
                collapse = ""
            ),
        " - files not found"))

    assert_that(all(dir.exists(dirname(output_files))),
        msg = paste(
            paste(
                unique(dirname(output_files[!dir.exists(dirname(output_files))])),
                collapse = ""
            ),
        " - directories not found"))   
       
    assert_that(file.exists(file.path(s_ref, "settings.Rds")),
        msg = paste(file.path(s_ref, "settings.Rds"), "not found"))   

    assert_that(file.exists(file.path(s_ref, "IRFinder.ref.gz")),
        msg = paste(file.path(s_ref, "IRFinder.ref.gz"), "not found"))

    ref_file = normalizePath(file.path(s_ref, "IRFinder.ref.gz"))
    
    assert_that(length(s_bam) == length(output_files),
        msg = "Number of output files and bam files must be the same")
    
    if(Has_OpenMP() > 0 & Use_OpenMP) {
        n_threads = floor(max_threads)
        n_threads = min(n_threads, length(s_bam))
        IRF_main_multithreaded(ref_file, s_bam, output_files, n_threads)
    } else {
        # Use BiocParallel
        n_rounds = ceiling(length(s_bam) / floor(max_threads))
        n_threads = ceiling(length(s_bam) / n_rounds)

        BPPARAM = BiocParallel::bpparam()
        if(Sys.info()["sysname"] == "Windows") {
          BPPARAM_mod = BiocParallel::SnowParam(n_threads)
          message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
        } else {
          BPPARAM_mod = BiocParallel::MulticoreParam(n_threads)
          message(paste("Using MulticoreParam", BPPARAM_mod$workers, "threads"))
        }

        row_starts = seq(1, by = n_threads,
            length.out = n_rounds)
            
        for(i in seq_len(n_rounds)) {
            selected_rows_subset = seq(row_starts[i], 
                min(length(s_bam), row_starts[i] + n_threads - 1)
            )
            BiocParallel::bplapply(selected_rows_subset,
                function(i, IRF_main, s_bam, reference_file, output_files) {
                    IRF_main(s_bam[i], reference_file, output_files[i])
                }, 

                s_bam = s_bam,
                output_files = output_files,
                IRF_main = IRF_main, 
                reference_file = ref_file,
                BPPARAM = BPPARAM_mod
            )
        }       
    }
    
    if(run_featureCounts == TRUE) {

        NxtIRF.CheckPackageInstalled("Rsubread", "2.4.0")
        settings = readRDS(file.path(s_ref, "settings.Rds"))
        ah_transcriptome = settings$ah_transcriptome
        gtf_file = settings$gtf_file

        if(ah_transcriptome != "") {
            assert_that(substr(ah_transcriptome,1,2) == "AH",
                msg = "Given transcriptome AnnotationHub reference is incorrect")
            assert_that(ah_transcriptome %in% names(ah),
                msg = paste(ah_transcriptome, "is not a record in given AnnotationHub object")
            )
            gtf_file = AnnotationHub::cache(ah[names(ah) == ah_transcriptome])
        } else {
            assert_that(file.exists(gtf_file),
                msg = paste("Given transcriptome gtf file", gtf_file,
                    "not found"))
        }
        
        assert_that(file.exists(gtf_file),
            msg = paste("Given gtf reference file", gtf_file, "does not exist. Could not run featureCounts")
        )
        
        # determine paired_ness, strandedness 
        data.list = get_multi_DT_from_gz(
            normalizePath(paste0(output_files[1], ".txt.gz")), 
            c("BAM", "Directionality")
        )
        stats = data.list$BAM
        direct = data.list$Directionality

        if(stats$Value[3] == 0 & stats$Value[4] > 0) {
            paired = TRUE
        } else if(stats$Value[3] > 0 && 
                stats$Value[4] / stats$Value[3] / 1000) {
            paired = TRUE
        } else {
            paired = FALSE
        }
        strand = direct$Value[9]
        if(strand == -1) strand = 2
        
        res = Rsubread::featureCounts(
            s_bam,
            annot.ext = gtf_file,
            isGTFAnnotationFile = TRUE,
            strandSpecific = strand,
            isPairedEnd = paired,
            requireBothEndsMapped = paired,
            nthreads = n_threads
        )

        # Append to existing main.FC.Rds if exists:
        
        if(file.exists(file.path(dirname(output_files[1]), "main.FC.Rds"))) {
            res.old = readRDS(file.path(dirname(output_files[1]), "main.FC.Rds"))

            # Check md5 of annotation to show same reference was used
            md5.old = openssl::md5(paste(
                res.old$annotation$GeneID, res.old$annotation$Chr,
                res.old$annotation$Start, res.old$annotation$End, 
                res.old$annotation$Strand, collapse=" "
                ))
            md5 = openssl::md5(paste(
                res$annotation$GeneID, res$annotation$Chr,
                res$annotation$Start, res$annotation$End, 
                res$annotation$Strand, collapse=" "
                ))
            md5.old.stat = openssl::md5(paste(
                res.old$stat$Status, collapse=" "
                ))
            md5.stat = openssl::md5(paste(
                res$stat$Status, collapse=" "
                ))
            if(md5 == md5.old & md5.stat == md5.old.stat) {
                # cbind stats
                new_samples = res$targets[!(res$targets %in% res.old$targets)]
                res$targets = c(res.old$targets, new_samples)

                res$stat = cbind(res.old$stat, res$stat[,new_samples])
                # cbind counts            
                res$counts = cbind(res.old$counts, res$counts[,new_samples])
            }
        }
        
        if(all(c("counts", "annotation", "targets", "stat") %in% names(res))) {
            saveRDS(res, file.path(dirname(output_files[1]), "main.FC.Rds"))
        }
        
    }   

}    

run_IRFinder = function(
            bamfile = "Unsorted.bam", 
            ref_file = "./IRFinder.ref.gz", 
            output_file = "./Sample",
            run_featureCounts = FALSE,
            gtf = "transcripts.gtf",
            ah_transcriptome = "",
            localHub = FALSE,
            ah = AnnotationHub(localHub = localHub)
        ) {
    s_bam = normalizePath(bamfile)
    s_ref = normalizePath(ref_file)
    system.time({
        IRF_main(s_bam, s_ref, output_file)
    })

    assert_that(file.exists(paste0(output_file, ".txt.gz")),
        msg = paste("IRFinder appears to have failed.",
            paste0(output_file, ".txt.gz"), "does not exist")
    )
    
    if(run_featureCounts == TRUE) {

        NxtIRF.CheckPackageInstalled("Rsubread", "2.4.0")
        if(ah_transcriptome != "") {
            assert_that(substr(ah_transcriptome,1,2) == "AH",
                msg = "Given transcriptome AnnotationHub reference is incorrect")
            assert_that(ah_transcriptome %in% names(ah),
                msg = paste(ah_transcriptome, "is not a record in given AnnotationHub object")
            )
            gtf_file = AnnotationHub::cache(ah[names(ah) == ah_transcriptome])
        } else {
            assert_that(file.exists(gtf),
                msg = paste("Given transcriptome gtf file", gtf,
                    "not found"))
            gtf_file = gtf
        }
        assert_that(file.exists(gtf_file),
            msg = paste("Given gtf reference file", gtf_file, "does not exist. Could not run featureCounts")
        )
        
        # determine paired_ness, strandedness 
        data.list = get_multi_DT_from_gz(
            normalizePath(paste0(output_file, ".txt.gz")), 
            c("BAM", "Directionality")
        )
        stats = data.list$BAM
        direct = data.list$Directionality

        if(stats$Value[3] == 0 & stats$Value[4] > 0) {
            paired = TRUE
        } else if(stats$Value[3] > 0 && 
                stats$Value[4] / stats$Value[3] / 1000) {
            paired = TRUE
        } else {
            paired = FALSE
        }
        strand = direct$Value[9]
        if(strand == -1) strand = 2
        
        res = Rsubread::featureCounts(
            s_bam,
            annot.ext = gtf_file,
            isGTFAnnotationFile = TRUE,
            strandSpecific = strand,
            isPairedEnd = paired,
            requireBothEndsMapped = paired
        )
        
        # Append to existing main.FC.Rds if exists:
        
        if(file.exists(file.path(dirname(output_files[1]), "main.FC.Rds"))) {
            res.old = readRDS(file.path(dirname(output_files[1]), "main.FC.Rds"))

            # Check md5 of annotation to show same reference was used
            md5.old = openssl::md5(paste(
                res.old$annotation$GeneID, res.old$annotation$Chr,
                res.old$annotation$Start, res.old$annotation$End, 
                res.old$annotation$Strand, collapse=" "
                ))
            md5 = openssl::md5(paste(
                res$annotation$GeneID, res$annotation$Chr,
                res$annotation$Start, res$annotation$End, 
                res$annotation$Strand, collapse=" "
                ))
            md5.old.stat = openssl::md5(paste(
                res.old$stat$Status, collapse=" "
                ))
            md5.stat = openssl::md5(paste(
                res$stat$Status, collapse=" "
                ))
            if(md5 == md5.old & md5.stat == md5.old.stat) {
                # cbind stats
                new_samples = res$targets[!(res$targets %in% res.old$targets)]
                res$target = c(res.old$targets, new_samples)

                res$stat = cbind(res.old$stat, res$stat[,new_samples])
                # cbind counts            
                res$counts = cbind(res.old$counts, res$counts[,new_samples])
            }
        }

        if(all(c("counts", "annotation", "targets", "stat") %in% names(res))) {
            saveRDS(res, file.path(dirname(output_files[1]), "main.FC.Rds"))
        }
        
    }
    
}

run_FeatureCounts = function(s_bam, sample_names, s_output, reference_path, strand, paired,
        localHub = FALSE,
        ah = AnnotationHub(localHub = localHub),       
        max_threads = max(parallel::detectCores() - 2, 1)) {
    assert_that(length(s_bam) > 0 && length(s_bam) == length(sample_names),
        msg = "Mismatch between s_bam and sample_names. Make sure same length")
    assert_that(all(file.exists(s_bam)),
        msg = "Some bam files do not exist")
    assert_that(all(dir.exists(dirname(s_output))),
        msg = "Root directory of s_output does not exist")
    assert_that(max_threads == 1 || max_threads <= parallel::detectCores() - 1,
        msg = paste(max_threads, " threads is not allowed for this system"))

    s_ref = normalizePath(reference_path)

    assert_that(file.exists(file.path(s_ref, "settings.Rds")),
        msg = paste(file.path(s_ref, "settings.Rds"), "not found"))   

    settings = readRDS(file.path(s_ref, "settings.Rds"))
    ah_transcriptome = settings$ah_transcriptome
    gtf_file = settings$gtf_file

    if(ah_transcriptome != "") {
        assert_that(substr(ah_transcriptome,1,2) == "AH",
            msg = "Given transcriptome AnnotationHub reference is incorrect")
        assert_that(ah_transcriptome %in% names(ah),
            msg = paste(ah_transcriptome, "is not a record in given AnnotationHub object")
        )
        gtf_file = AnnotationHub::cache(ah[names(ah) == ah_transcriptome])
    } else {
        assert_that(file.exists(gtf_file),
            msg = paste("Given transcriptome gtf file", gtf_file,
                "not found"))
    }

    res = Rsubread::featureCounts(
        s_bam,
        annot.ext = gtf_file,
        isGTFAnnotationFile = TRUE,
        strandSpecific = strand,
        isPairedEnd = paired,
        requireBothEndsMapped = paired
    )
    
    res$targets = sample_names
    colnames(res$counts) = sample_names
    colnames(res$stat)[-1] = sample_names

    saveRDS(res, paste0(s_output, ".FC.Rds"))
}

#' @export
merge_FeatureCounts = function(record_1 = "main.FC.Rds", 
        record_2 = "addit.FC.Rds", overwrite = FALSE) {
    assert_that(all(file.exists(c(record_1, record_2))),
        msg = "Some records do not exist")
        
    res1 = readRDS(record_1)
    res2 = readRDS(record_2)
    
    # Check md5 of annotation to show same reference was used
    md5.1 = openssl::md5(paste(
        res1$annotation$GeneID, res1$annotation$Chr,
        res1$annotation$Start, res1$annotation$End, 
        res1$annotation$Strand, collapse=" "
    ))
    md5.2 = openssl::md5(paste(
        res2$annotation$GeneID, res2$annotation$Chr,
        res2$annotation$Start, res2$annotation$End, 
        res2$annotation$Strand, collapse=" "
    ))
    md5.1.stat = openssl::md5(paste(
        res1$stat$Status, collapse=" "
    ))
    md5.2.stat = openssl::md5(paste(
        res2$stat$Status, collapse=" "
    ))
    res = list()
    assert_that(md5.1 == md5.2 & md5.1.stat == md5.2.stat, 
        msg = paste("Record annotations do not match.",
            "Perhaps they were generated under different gene annotations"))
    # cbind stats
    if(overwrite) {
        old_samples = res1$targets[!(res1$targets %in% res2$targets)]
        res$targets = c(old_samples, res2$targets)
        res$stat = cbind(res1$stat[,1], res1$stat[,old_samples],
            res2$stat[-1])
        res$counts = cbind(res1$counts[,old_samples], res2$counts)
    } else {
        new_samples = res2$targets[!(res2$targets %in% res1$targets)]
        res$targets = c(res1$targets, new_samples)
        res$stat = cbind(res1$stat, res2$stat[new_samples])
        res$counts = cbind(res1$counts, res2$counts[, new_samples])
    }
    res$annotation = res1$annotation

    return(res)
}

#' @export
run_IRFinder_GenerateMapReads = function(genome.fa = "", out.fa, read_len = 70, read_stride = 10, error_pos = 35) {
  return(
    IRF_GenerateMappabilityReads(normalizePath(genome.fa), 
        file.path(normalizePath(dirname(out.fa)), out.fa),
        read_len = read_len, 
        read_stride = read_stride, 
        error_pos = error_pos)
  )
}

#' @export
run_IRFinder_MapExclusionRegions = function(bamfile = "", output_file, threshold = 4, includeCov = FALSE) {
  s_bam = normalizePath(bamfile)
  assert_that(file.exists(s_bam),
    msg = paste(s_bam, "does not exist"))
  return(
    IRF_GenerateMappabilityRegions(s_bam, 
        normalizePath(output_file),
        threshold = threshold,
        includeCov = includeCov)
  )
}


#' @export
run_Gunzip = function(infile = "", outfile) {
  file_to_read = normalizePath(infile)
  assert_that(file.exists(file_to_read),
    msg = paste(file_to_read, "does not exist"))
    
  assert_that(dir.exists(dirname(outfile)),
    msg = paste(dirname(outfile), "does not exist"))
    
  IRF_gunzip(file_to_read, outfile)
}

#' @export
get_multi_DT_from_gz = function(infile = "", block_headers = c("Header1", "Header2")) {
  file_to_read = normalizePath(infile)
  assert_that(file.exists(file_to_read),
    msg = paste(file_to_read, "does not exist"))


  df.list = IRF_gunzip_DF(file_to_read, block_headers)
  for(i in seq_len(length(df.list))) {
    for(j in seq_len(length(df.list[[i]]))) {
      suppressWarnings({
        if(all(df.list[[i]][[j]] == "NA" | !is.na(as.numeric(df.list[[i]][[j]])))) {
          df.list[[i]][[j]] = as.numeric(df.list[[i]][[j]])
        }
      })
    }
    df.list[[i]] = as.data.table(df.list[[i]])
  }
  return(df.list)
}

#' @export
GetCoverage = function(file, seqname = "", start = 0, end = 0, strand = 2) {
  assert_that(as.numeric(strand) %in% c(0,1,2),
                          msg = "Invalid strand. Must be either 0 (+), 1 (-) or 2(*)")
  assert_that(as.numeric(start) <= as.numeric(end) | end == 0,
                          msg = "Null or negative regions not allowed")
  message(paste("Fetching file=",file,"coords", seqname, start, end, strand))                       
  if(seqname == "") {
    raw_list = IRF_RLEList_From_Cov(normalizePath(file), strand)
    final_list = list()
    if(length(raw_list) > 0) {
      for(i in 1:length(raw_list)) {
        final_list[[i]] = S4Vectors::Rle(raw_list[[i]]$values, raw_list[[i]]$length)
      }
    } else {
      return(NULL)
    }
    final_RLE = as(final_list, "RleList")
    names(final_RLE) = names(raw_list)
    return(final_RLE)
  } else if(end == 0) {
    raw_RLE = IRF_RLE_From_Cov(normalizePath(file), as.character(seqname), 0,0, as.numeric(strand))
		final_RLE = S4Vectors::Rle(raw_RLE$values, raw_RLE$length)
  } else {
    raw_RLE = IRF_RLE_From_Cov(normalizePath(file), as.character(seqname), 
			round(as.numeric(start)), round(as.numeric(end)), as.numeric(strand))
		final_RLE = S4Vectors::Rle(raw_RLE$values, raw_RLE$length)													 
  }
  
}

#' @export
IsCOV = function(coverage_files) {
    for(i in coverage_files) {
        if(file.exists(i) && IRF_Check_Cov(normalizePath(i))) {
            # do nothing
        } else {
            return(FALSE)
        }
    }    
    return(TRUE)
}