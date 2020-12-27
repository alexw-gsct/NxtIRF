globalVariables(c(":=","."))

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

NxtIRF.gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

#' @export
NxtIRF.CoordToGR = function(coordinates) {
	temp = tstrsplit(coordinates,split="/")
	strand = as.character(temp[[2]])
	temp2 = tstrsplit(temp[[1]],split=":")
	seqnames = temp2[[1]]
	temp3 = tstrsplit(temp2[[2]],split="-")
	start = temp3[[1]]
	end = temp3[[2]]
	return(GRanges(seqnames = seqnames, ranges = IRanges(
		start = as.numeric(start), end = as.numeric(end)),
		strand = strand))
}

#' @export
NxtIRF.CoordToInt = function(coordinates, type = "") {
	temp = tstrsplit(coordinates,split="/")
	strand = as.character(temp[[2]])
	temp2 = tstrsplit(temp[[1]],split=":")
	seqnames = temp2[[1]]
	temp3 = tstrsplit(temp2[[2]],split="-")
	start = temp3[[1]]
	end = temp3[[2]]
	if(type == "start") {
		return(as.numeric(start))
	} else if(type == "end") {
		return(as.numeric(end))
	} else if(type == "strand") {
		return(as.character(strand))
	} else if(type == "seqnames") {
		return(as.character(seqnames))
	} else {
		final = list(as.character(seqnames),
			as.numeric(start),as.numeric(end),
			as.character(strand))
		names(final) = c("seqnames","start","end","strand")
		return(final)
	}
}

NxtIRF.StartBenchmarkTimer <- function(irf = NULL, handle = "1") {

	assert_that(is(irf, "NxtProject"),
		msg = "Valid NxtProject object required")

	irf@benchmarks[[handle]] <- Sys.time()
	return(irf)

}

NxtIRF.Benchmark <- function(irf = NULL, handle = "1", msg = "Time elapsed:", units = "auto") {

	assert_that(is(irf, "NxtProject"),
		msg = "Valid NxtProject object required")

	end_time <- Sys.time()
	time_diff = round(difftime(end_time, irf@benchmarks[[handle]], units = units),2)
	message(handle, ":", msg,
		time_diff,
		units(time_diff), "\n")
	irf@benchmarks[[handle]] = NULL
	return(irf)	
	
}

NxtIRF.CheckPackageInstalled <- function(package = "DESeq2", version = "1.0.0") {
	assert_that(
		tryCatch(ifelse(packageVersion(package)>=version, TRUE, FALSE),
		error = function(e) FALSE),
		msg = paste(package, "package is not installed; and is required for this function")
	)
}

NxtIRF.SplitVector <- function(vector = "", n_workers = 1) {
	assert_that(n_workers >= 1,
		msg = "n_workers must be at least  1")
	n_workers_use = as.integer(n_workers)
	
	assert_that(length(vector) >= 1,
		msg = "vector to split must be of length at least 1")
  
  if(n_workers_use > length(vector)) n_workers_use = length(vector)
	vector_starts = round(seq(1, length(vector) + 1, length.out = n_workers_use + 1))
	vector_starts = unique(vector_starts)
	
	return_val = list()
	for(i in seq_len(length(vector_starts) - 1)) {
		return_val[[i]] = vector[seq(vector_starts[i], vector_starts[i+1] - 1)]
	}
	return(return_val)
}

semi_join.DT = function(A, B, by, nomatch = 0) {
	A[A[B, on = by, which = TRUE, nomatch = nomatch]]
}

is_valid <- function(x) {
    !is.null(x) && length(x) > 0 && 
        (isS4(x) || !is.na(x)) && 
        (!is.character(x) || (x != "" && x != "(none)"))
}

make.path.relative = function(base, target) {
    common = sub('^([^|]*)[^|]*(?:\\|\\1[^|]*)$', '^\\1/?', paste0(base, '|', target))
    
    paste0(gsub('[^/]+/?', '../', sub(common, '', base)),
           sub(common, '', target))
}