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

#' @export
run_IRFinder_mappability = function(bamfile = "Unsorted.bam", ref_path = "./REF") {
  s_bam = normalizePath(bamfile)
  s_ref = normalizePath(ref_path)
  system.time({
    IRF_genmap(s_bam, s_ref)
  })
}

#' @export
run_IRFinder_GenerateMapReads = function(genome.fa = "", out.fa, read_len = 70, read_stride = 10, error_pos = 35) {
  return(
    IRF_SupplyMappaReads(genome.fa, out.fa, read_len = read_len, read_stride = read_stride, error_pos = error_pos)
  )
}

#' @export
run_IRFinder_GenerateMapRegionReads = function(genome.fa = "", out.fa, region.txt, read_len = 70, read_stride = 10, error_pos = 35) {
  return(
    IRF_SupplyMappaRegionReads(genome.fa, region.txt, out.fa, read_len = read_len, read_stride = read_stride, error_pos = error_pos)
  )
}


#' @export
run_IRFinder_gunzip = function(infile, outfile) {
  IRF_gunzip(normalizePath(infile), outfile)
}