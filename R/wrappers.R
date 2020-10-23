# wrappers to R/C++

#' @export
run_IRFinder = function(bamfile = "Unsorted.bam", ref_file = "./IRFinder.ref", output_file = "./Output.txt") {
  s_bam = normalizePath(bamfile)
  s_ref = normalizePath(ref_file)
#  s_output = normalizePath(output_file)
  system.time({
    IRF_main(s_bam, s_ref, output_file)
  })
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
  assertthat::assert_that(file.exists(file_to_read),
    msg = paste(file_to_read, "does not exist"))
    
  assertthat::assert_that(dir.exists(dirname(outfile)),
    msg = paste(dirname(outfile), "does not exist"))
    
  IRF_gunzip(file_to_read, outfile)
}

#' @export
get_multi_DT_from_gz = function(infile = "", block_headers = c("Header1", "Header2")) {
  file_to_read = normalizePath(infile)
  assertthat::assert_that(file.exists(file_to_read),
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
  assertthat::assert_that(as.numeric(strand) %in% c(0,1,2),
                          msg = "Invalid strand. Must be either 0 (+), 1 (-) or 2(*)")
  assertthat::assert_that(as.numeric(start) <= as.numeric(end) | end == 0,
                          msg = "Null or negative regions not allowed")
  if(seqname == "") {
    raw_list = IRF_RLEList_From_Cov(normalizePath(file), strand)
    final_list = list()
    if(length(raw_list) > 0) {
      for(i in 1:length(raw_list)) {
        final_list[[i]] = S4Vectors::Rle(raw_list[[i]]$values, raw_list[[i]]$length)
      }
    }
    final_RLE = S4Vectors:::new_SimpleList_from_list("SimpleRleList", final_list)
    names(final_RLE) = names(raw_list)
    return(final_RLE)
  } else if(end == 0) {
    raw_RLE = IRF_RLE_From_Cov(normalizePath(file), as.character(seqname), 0,0, as.numeric(strand))
  } else {
    raw_RLE = IRF_RLE_From_Cov(normalizePath(file), as.character(seqname), 
                               round(as.numeric(start)), round(as.numeric(end)), as.numeric(strand))
  }
  final_RLE = S4Vectors::Rle(raw_RLE$values, raw_RLE$length)
}
