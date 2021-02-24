globalVariables(c(":=","."))

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

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
    if(Sys.info()["sysname"] == "Windows") {
        base = normalizePath(base, winslash = "/")
    }
    common = sub('^([^|]*)[^|]*(?:\\|\\1[^|]*)$', '^\\1/?', paste0(base, '|', target))
    
    paste0(gsub('[^/]+/?', '../', sub(common, '', base)),
           sub(common, '', target))
}

# GGPLOT themes

#' @export
theme_white = theme(axis.line.x = element_line(colour = "black"),
			panel.grid.major = element_line(size = rel(0.5), colour="grey"),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			panel.background = element_blank(),
			legend.position = "none",
			axis.title.x.top = element_blank(),
			# axis.text.x.bottom = element_blank(),
			# axis.text.y = element_blank(),
			axis.line.x.bottom = element_blank(),
			axis.text=element_text(size=rel(1.0)),
            plot.title = element_text(hjust = 0.5),
			# axis.title.x=element_blank(),
			# axis.title.y=element_blank()
            )

#' @export
theme_white_legend = theme(axis.line.x = element_line(colour = "black"),
			panel.grid.major = element_line(size = rel(0.5), colour="grey"),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			panel.background = element_blank(),
			# legend.position = "none",
			axis.title.x.top = element_blank(),
			# axis.text.x.bottom = element_blank(),
			# axis.text.y = element_blank(),
			axis.line.x.bottom = element_blank(),
			axis.text=element_text(size=rel(1.0)),
            plot.title = element_text(hjust = 0.5),
			# axis.title.x=element_blank(),
			# axis.title.y=element_blank()
            )

#' @export
download_NxtIRF_example <- function(destination_dir = tempdir()) {
    assert_that(dir.exists(dirname(destination_dir)),
        msg = paste(dirname(destination_dir), "must exist"))
        
        resource_path = "https://raw.github.com/alexw-gsct/NxtIRF_resources/main/example"
        
    files = c("UT_1.bam", "UT_2.bam", "UT_3.bam", "D2_1.bam", "D2_2.bam", "D2_3.bam")
    
    for(sample_file in files) {
        cache_file = parse_valid_file(file.path(resource_path, sample_file))
        if(file.exists(cache_file)) {
            file.copy(cache_file, file.path(destination_dir, sample_file))
            message(paste(sample_file, "downloaded to", destination_dir))
        } else {
            message(paste(sample_file, "could not be downloaded from", resource_path))
        }
    }
    
    if(all(file.exists(file.path(destination_dir, files)))) {
        message(paste("All files successfully downloaded to", destination_dir))
    } else {
        message(paste("One or more example files could not be downloaded to", destination_dir))        
    }
}

dash_progress <- function(message = "", total_items = 1) {
    assert_that(total_items = round(total_items) & total_items > 0,
        msg = "dash_progress needs at least 1 item")
    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/total_items, message = message)
    }
}