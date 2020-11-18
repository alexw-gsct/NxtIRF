# wrappers to R/C++

#' @export
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
        msg = paste("IRFinder appears to have failed")
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

        if(all(c("counts", "annotation", "targets", "stat")) %in% names(res)) {
            saveRDS(res, paste(output_file, "fc.Rds", sep="."))
        }
        
    }
    
}

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
  if(seqname == "") {
    raw_list = IRF_RLEList_From_Cov(normalizePath(file), strand)
    final_list = list()
    if(length(raw_list) > 0) {
      for(i in 1:length(raw_list)) {
        final_list[[i]] = S4Vectors::Rle(raw_list[[i]]$values, raw_list[[i]]$length)
      }
    }
    # final_RLE = S4Vectors:::new_SimpleList_from_list("SimpleRleList", final_list)
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
