# wrappers to R/C++

#' @export
run_IRFinder = function(bamfile = "Unsorted.bam", ref_path = "./REF", output_path = "./Output") {
  s_bam = normalizePath(bamfile)
  s_ref = normalizePath(ref_path)
  if(substr(output_path, nchar(output_path), nchar(output_path)) == "/") {
    s_output = normalizePath(paste0(output_path, ".."))
  } else {
    s_output = normalizePath(paste(output_path, "..", sep="/"))
  }
  system.time({
    IRF_main(s_bam, s_ref, s_output)
  })
}

#' @export
run_IRFinder_mappability = function(bamfile = "Unsorted.bam", ref_path = "./REF") {
  s_bam = normalizePath(bamfile)
  s_ref = normalizePath(ref_path)
  system.time({
    IRF_genmap(s_bam, s_ref)
  })
}