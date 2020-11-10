#' @export
FetchAH <- function(ah_record_name, localHub = FALSE, 
    ah = AnnotationHub(localHub = localHub), verbose = FALSE) {
################################################################################
    assert_that(substr(ah_record_name, 1, 2) == "AH", msg = 
        paste(ah_record_name,
        "does not appear to be a valid AnnotationHub record name"))
    assert_that(ah_record_name %in% names(ah), msg = 
        paste(ah_record_name, 
            "is not found in AnnotationHub index.",
            "Perhaps check online connection or record name"))

    ah.record = ah[names(ah) == ah_record_name]

	if(verbose) message("Downloading asset from AnnotationHub, if required...", 
        appendLF = FALSE)
    cache_loc = AnnotationHub::cache(ah.record)
    if(verbose) message("done")
    
    assert_that(file.exists(cache_loc), 
        msg = "AnnotationHub cache error - asset not found")

    if(ah.record$rdataclass == "GRanges") {
        if(verbose) message("Importing to memory as GRanges object...", 
            appendLF = FALSE)
        gtf = rtracklayer::import(cache_loc, "gtf", genome = ah.record$genome)
        if(verbose) message("done")
        return(gtf)
    } else if(ah.record$rdataclass == "TwoBitFile") {
        if(verbose) message("Importing to memory as GRanges object...", 
            appendLF = FALSE)
        twobit = rtracklayer::TwoBitFile(cache_loc)
        if(verbose) message("done")
        return(twobit)
    } else {
		message(paste(ah_record_name, 
            "does not appear to be a GRanges object or TwoBitFile,",
            "returning NULL"))
		return(NULL)
	}
}

################################################################################
parse_valid_file <- function(file, msg) {
    if(is_valid(file) && !file.exists(file)) {
        message(paste(file, "not found.",
            "Reference generated without", msg))
        return("")
    } else if(!is_valid(file) == 0) {
         return("")  
    } else if (file.exists(file)) {
        return(file)
    } else {
        return("")
    }
}

fetch_genome_defaults <- function(genome_type, nonPolyARef, MappabilityRef,
        BlacklistRef) {
    if(genome_type == "hg38") {
        nonPolyAFile = system.file("extra-input-files/Human_hg38_nonPolyA_ROI.bed",
            package = "NxtIRF")
        MappabilityFile = 
            system.file("extra-input-files/Mappability_Regions_hg38_v94.txt.gz", 
            package = "NxtIRF")
    } else if(genome_type == "hg19")  {
        nonPolyAFile = system.file("extra-input-files/Human_hg19_nonPolyA_ROI.bed", 
            package = "NxtIRF")
        MappabilityFile = 
            system.file("extra-input-files/Mappability_Regions_hg19_v75.txt.gz", 
            package = "NxtIRF")
    } else if(genome_type == "mm10")  {
        nonPolyAFile = system.file("extra-input-files/Mouse_mm10_nonPolyA_ROI.bed", 
            package = "NxtIRF")
        MappabilityFile =
            system.file("extra-input-files/Mappability_Regions_mm10_v94.txt.gz",
            package = "NxtIRF")
    } else if(genome_type == "mm9")  {
        nonPolyAFile = system.file("extra-input-files/Mouse_mm9_nonPolyA_ROI.bed",
            package = "NxtIRF")
        MappabilityFile = 
            system.file("extra-input-files/Mappability_Regions_mm9_v67.txt.gz",
            package = "NxtIRF")
    } else {
        nonPolyAFile = 
            parse_valid_file(nonPolyARef, "non-polyA reference")
        MappabilityFile = 
            parse_valid_file(MappabilityRef, "Mappability reference")
    }

    BlacklistFile =  
        parse_valid_file(BlacklistRef, "Blacklist exclusion")
 
    final = list(nonPolyAFile = nonPolyAFile, MappabilityFile = MappabilityFile,
        BlacklistFile = BlacklistFile)
    return(final)
}

################################################################################
prep_ref_path <- function(reference_path) {
    base_output_path = normalizePath(dirname(reference_path))
    if(!dir.exists(file.path(base_output_path, basename(reference_path)))) {
        dir.create(file.path(base_output_path, basename(reference_path)))
    }
    if(!dir.exists(file.path(base_output_path, basename(reference_path), "fst"))) {
        dir.create(file.path(base_output_path, basename(reference_path), "fst"))
    }
}

################################################################################
process_gtf <- function(gtf.gr, reference_path) {

    message("Processing gtf file...", appendLF = FALSE)

    # fix gene / transcript names with '/' (which breaks IRFinder code)
    gtf.gr$gene_name = gsub("/","_",gtf.gr$gene_name)
    gtf.gr$transcript_name = gsub("/","_",gtf.gr$transcript_name)

    # Extracting and saving Genes, Transcripts, Exons, Proteins and saving as .fst files for faster random access    
    Genes = gtf.gr[gtf.gr$type == "gene"]
    Genes <- GenomeInfoDb::sortSeqlevels(Genes)
    Genes <- sort(Genes)
    Genes$gene_display_name = paste0(Genes$gene_name, " (", Genes$gene_id, ")")

    # Annotate gene_groups_stranded / unstranded
    Genes.Group.stranded = as.data.table(reduce(Genes))
    # Order by seqnames then by start coordinate then by strand
    setorder(Genes.Group.stranded, seqnames, start, strand)
		gene_group_stranded <- gene_group_unstranded <- NULL
    Genes.Group.stranded[, gene_group_stranded := .I]
    OL = findOverlaps(
        Genes, makeGRangesFromDataFrame(as.data.frame(Genes.Group.stranded)))
    Genes$gene_group_stranded[OL@from] =
        Genes.Group.stranded$gene_group_stranded[OL@to]

    Genes.Group.unstranded = as.data.table(reduce(Genes, ignore.strand = TRUE))
    # Order by seqnames then by start coordinate then by strand
    setorder(Genes.Group.unstranded, seqnames, start)
    Genes.Group.unstranded[, gene_group_unstranded := .I]
    OL = findOverlaps(
        Genes, makeGRangesFromDataFrame(as.data.frame(Genes.Group.unstranded),
        ignore.strand = TRUE))
    Genes$gene_group_unstranded[OL@from] =
        Genes.Group.unstranded$gene_group_unstranded[OL@to]
    
    write.fst(as.data.frame(Genes), file.path(reference_path,"fst","Genes.fst"))
        
    Transcripts = gtf.gr[gtf.gr$type == "transcript"]
    Transcripts <- GenomeInfoDb::sortSeqlevels(Transcripts)
    Transcripts <- sort(Transcripts)
    
    if("transcript_biotype" %in% names(mcols(Transcripts))) {
        # do nothing
    } else if("transcript_type" %in% names(mcols(Exons))) {
        colnames(mcols(Transcripts))[which(colnames(mcols(Transcripts)) == 
                "transcript_type")] = "transcript_biotype"
    } else {
        mcols(Transcripts)$transcript_biotype = "protein_coding"
    }

    if("transcript_support_level" %in% names(mcols(Transcripts))) {
        Transcripts$transcript_support_level =
            tstrsplit(Transcripts$transcript_support_level, split=" ")[[1]]
        Transcripts$transcript_support_level[
            is.na(Transcripts$transcript_support_level)] = "NA"
    }
   
    write.fst(as.data.frame(Transcripts),
        file.path(reference_path,"fst","Transcripts.fst"))

    Exons = gtf.gr[gtf.gr$type == "exon"]
    Exons <- GenomeInfoDb::sortSeqlevels(Exons)
    Exons <- sort(Exons)

    # transcript_biotype is very important field. If Gencode, this is transcript_type.
    #    In rare case we do not have this field
    # This next bit ensures transcript_biotype exists.

    if("transcript_biotype" %in% names(mcols(Exons))) {

    } else if("transcript_type" %in% names(mcols(Exons))) {
        colnames(mcols(Exons))[
            which(colnames(mcols(Exons)) == "transcript_type")] = 
            "transcript_biotype"
    } else {
        mcols(Exons)$transcript_biotype = "protein_coding"
    }
    tmp.exons.exclude =  Exons[!grepl("intron", Exons$transcript_biotype)]

    # Assign gene groups then bake exon-groups into Exons
    tmp.Exons.Group.stranded = as.data.table(reduce(tmp.exons.exclude))
    
    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(tmp.Exons.Group.stranded)), 
        makeGRangesFromDataFrame(as.data.frame(Genes.Group.stranded)))
    tmp.Exons.Group.stranded$gene_group[OL@from] =
        Genes.Group.stranded$gene_group_stranded[OL@to]

    # Some retained_intron transcripts have terminal exons lying outside to that of main transcripts. Include these also
    tmp.exons.exclude.span = split(
      makeGRangesFromDataFrame(as.data.frame(tmp.Exons.Group.stranded)),
      tmp.Exons.Group.stranded$gene_group)
    tmp.exons.exclude.span = unlist(range(tmp.exons.exclude.span),use.names=TRUE)
 
    tmp.exons.RI =  Exons[grepl("intron", Exons$transcript_biotype)]
    if(length(tmp.exons.RI) > 0) {
      OL = findOverlaps(
          tmp.exons.RI, 
          tmp.exons.exclude.span
      )
      tmp.exons.RI = tmp.exons.RI[-OL@from]
      tmp.exons.exclude = c(tmp.exons.exclude, tmp.exons.RI)
    # Reassign everything
      tmp.Exons.Group.stranded = as.data.table(reduce(tmp.exons.exclude))
      
      OL = findOverlaps(
          makeGRangesFromDataFrame(as.data.frame(tmp.Exons.Group.stranded)), 
          makeGRangesFromDataFrame(as.data.frame(Genes.Group.stranded)))
      tmp.Exons.Group.stranded$gene_group[OL@from] =
        Genes.Group.stranded$gene_group_stranded[OL@to]
    }

    setorder(tmp.Exons.Group.stranded, seqnames, start, strand)
		gene_group <- exon_group <- NULL
    tmp.Exons.Group.stranded[, exon_group := data.table::rowid(gene_group)]
    tmp.Exons.Group.stranded[strand == "-", 
        exon_group := max(exon_group) + 1 - exon_group, by = "gene_group"]

    tmp.Exons.Group.unstranded = as.data.table(reduce(tmp.exons.exclude,
        ignore.strand = TRUE))
    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(tmp.Exons.Group.unstranded)), 
        makeGRangesFromDataFrame(as.data.frame(Genes.Group.unstranded)),
        , ignore.strand = TRUE)
    tmp.Exons.Group.unstranded$gene_group[OL@from] =
        Genes.Group.unstranded$gene_group_unstranded[OL@to]
    setorder(tmp.Exons.Group.unstranded, seqnames, start, strand)
    tmp.Exons.Group.unstranded[, exon_group := data.table::rowid(gene_group)]
    tmp.Exons.Group.unstranded[strand == "-", 
        exon_group := max(exon_group) + 1 - exon_group, by = "gene_group"]

    # Now annotate all exons in Exons with the gene and exon groups
    OL = findOverlaps(
        Exons, 
        makeGRangesFromDataFrame(as.data.frame(tmp.Exons.Group.stranded)))
    Exons$gene_group_stranded[OL@from] =
        tmp.Exons.Group.stranded$gene_group[OL@to]
    Exons$exon_group_stranded[OL@from] =
        tmp.Exons.Group.stranded$exon_group[OL@to]

    OL = findOverlaps(
        Exons, 
        makeGRangesFromDataFrame(as.data.frame(tmp.Exons.Group.unstranded)),
        ignore.strand = TRUE)
    Exons$gene_group_unstranded[OL@from] =
        tmp.Exons.Group.unstranded$gene_group[OL@to]
    Exons$exon_group_unstranded[OL@from] =
        tmp.Exons.Group.unstranded$exon_group[OL@to]

#   Why filter by protein_coding or processed_transcript? Make this an option

    
    # Finally write to disk
    write.fst(as.data.frame(Exons), file.path(reference_path,"fst","Exons.fst"))
    # Also write tmp.Exon groups
    write.fst(rbind(tmp.Exons.Group.stranded, tmp.Exons.Group.unstranded), 
        file.path(reference_path,"fst","Exons.Group.fst"))
        
    Proteins = gtf.gr[gtf.gr$type == "CDS"]
    Proteins <- GenomeInfoDb::sortSeqlevels(Proteins)
    Proteins <- sort(Proteins)
    write.fst(as.data.frame(Proteins),
        file.path(reference_path,"fst","Proteins.fst"))
    
    gtf.misc = gtf.gr[!gtf.gr$type %in% c("gene", "transcript", "exon", "CDS")]
    gtf.misc <- GenomeInfoDb::sortSeqlevels(gtf.misc)
    gtf.misc <- sort(gtf.misc)
    write.fst(as.data.frame(gtf.misc),
        file.path(reference_path,"fst","Misc.fst"))
    
    message("done\n")
}

################################################################################
process_introns <- function(reference_path, genome, UseExtendedTranscripts) {
    message("Processing introns...", appendLF = FALSE)    

    Exons = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path,"fst","Exons.fst")),
        keep.extra.columns = TRUE
    )
    Transcripts = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path,"fst","Transcripts.fst")),
        keep.extra.columns = TRUE
    )    
    Proteins = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path,"fst","Proteins.fst")),
        keep.extra.columns = TRUE
    )
    Exons.Group = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path,"fst","Exons.Group.fst")),
        keep.extra.columns = TRUE
    )
    Exons.Group.stranded = Exons.Group[strand(Exons.Group) != "*"]
    Exons.Group.unstranded = Exons.Group[strand(Exons.Group) == "*"]
    
    if(UseExtendedTranscripts == FALSE) {
        candidate.transcripts = Exons[Exons$transcript_biotype %in% 
            c("processed_transcript", "protein_coding")]
    } else {
        candidate.transcripts = Exons   
    }
    candidate.introns = grlGaps(
        split(candidate.transcripts, candidate.transcripts$transcript_id)
    )
    candidate.introns = as.data.table(candidate.introns)    
    candidate.introns$group = NULL
    colnames(candidate.introns)[1] = "transcript_id"
    setorder(candidate.introns, seqnames, start, end, strand)

# Annotating Introns:
    intron_number <- transcript_id <- intron_id <- i.gene_name <-i.gene_id <- 
    i.transcript_name <- i.transcript_biotype <- i.transcript_support_level <-
    transcript_support_level <- NULL
        
    candidate.introns[,intron_number := data.table::rowid(transcript_id)]
    candidate.introns[strand == "-", intron_number := 
        max(intron_number) + 1 - intron_number, by = "transcript_id"]
    candidate.introns[,intron_id := 
        paste0(transcript_id, "_Intron", intron_number)]
    candidate.introns[as.data.table(Transcripts), on = "transcript_id", 
        c("gene_name", "gene_id", "transcript_name", "transcript_biotype") := 
        list(i.gene_name, i.gene_id, i.transcript_name, i.transcript_biotype)]
    
    # Grab splice motifs at this point; filter by valid splice motifs
    donor.introns = data.frame(seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+", 
            candidate.introns$start, candidate.introns$end - 1),
        stop = ifelse(candidate.introns$strand == "+",
            candidate.introns$start + 1, candidate.introns$end),
        strand = candidate.introns$strand)
    donor.seq = getSeq(genome, makeGRangesFromDataFrame(donor.introns))
    acceptor.introns = data.frame(seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+", 
            candidate.introns$end - 1, candidate.introns$start),
        stop = ifelse(candidate.introns$strand == "+", 
            candidate.introns$end, candidate.introns$start + 1),
        strand = candidate.introns$strand)
    acceptor.seq = getSeq(genome, makeGRangesFromDataFrame(acceptor.introns))
    candidate.introns$splice_motif = paste0(donor.seq, acceptor.seq)

# Do other annotations here:
    i.protein_id <- i.ccds_id <- intron_start <- intron_end <-
    exon_group_stranded_upstream <- exon_group_stranded_downstream <- 
    exon_group_unstranded_upstream <- exon_group_unstranded_downstream <- NULL
    
    candidate.introns[as.data.table(Transcripts), on = "transcript_id", 
        c("gene_name", "gene_id", "transcript_name") :=
        list(i.gene_name, i.gene_id, i.transcript_name)]
    if("transcript_support_level" %in% names(mcols(Transcripts))) {
        candidate.introns[as.data.table(Transcripts), on = "transcript_id", 
        c("transcript_support_level") := list(i.transcript_support_level)]
        candidate.introns[, transcript_support_level := 
            tstrsplit(transcript_support_level, split=" ")[[1]]]
        candidate.introns[is.na(transcript_support_level),
            transcript_support_level := "NA"]
    }
    if("protein_id" %in% names(mcols(Proteins))) {
        candidate.introns[as.data.table(Proteins), on = "transcript_id", 
        c("protein_id") := list(i.protein_id)]
    }
    if("ccds_id" %in% names(mcols(Exons))) {
        candidate.introns[as.data.table(Exons), on = "transcript_id", 
        c("ccds_id") := list(i.ccds_id)]
    }

    # Cannot annotate candidate introns by min and max exon_groups
		# because retained introns will overlap one or more exon groups
		# need to walk start -=1, end += 1, then do the overlap thing
    candidate.introns[,intron_start := start]
    candidate.introns[,intron_end := end]

    candidate.introns$gene_group_stranded = NA
    candidate.introns$exon_group_stranded_upstream = NA
    candidate.introns$exon_group_stranded_downstream = NA        
				
    candidate.introns[strand == "+", start := intron_start - 1]
    candidate.introns[strand == "+", end := intron_start]
    candidate.introns[strand == "-", start := intron_end]
    candidate.introns[strand == "-", end := intron_end + 1]		
    
    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns)), 
        makeGRangesFromDataFrame(as.data.frame(Exons.Group.stranded)))				
    candidate.introns$gene_group_stranded[OL@from] =
        Exons.Group.stranded$gene_group[OL@to]
    candidate.introns$exon_group_stranded_upstream[OL@from] =
        Exons.Group.stranded$exon_group[OL@to]		

    candidate.introns[strand == "-", start := intron_start - 1]
    candidate.introns[strand == "-", end := intron_start]
    candidate.introns[strand == "+", start := intron_end]
    candidate.introns[strand == "+", end := intron_end + 1]		
    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns)), 
        makeGRangesFromDataFrame(as.data.frame(Exons.Group.stranded)))			
    candidate.introns$exon_group_stranded_downstream[OL@from] =
        Exons.Group.stranded$exon_group[OL@to]		

################################################################################

    # Need fix for retained_introns or sense_intronic where junction extends into the obligate introns
    tmp = makeGRangesFromDataFrame(as.data.frame(Exons.Group.stranded),
        keep.extra.columns = TRUE)
    tmp.Introns.Group.stranded = grlGaps(
        split(tmp, tmp$gene_group)
    )
    tmp.Introns.Group.stranded = as.data.table(tmp.Introns.Group.stranded)
    data.table::setnames(tmp.Introns.Group.stranded, "group_name", "gene_group")
    tmp.Introns.Group.stranded[,intron_number := data.table::rowid(gene_group)]
    tmp.Introns.Group.stranded[strand == "-", 
        intron_number := max(intron_number) + 1 - intron_number,
        by = "gene_group"]

    candidate.introns.subset =
        candidate.introns[is.na(exon_group_stranded_upstream)]
    candidate.introns = candidate.introns[!is.na(exon_group_stranded_upstream)]
		candidate.introns.subset[strand == "+", start := intron_start - 1]
		candidate.introns.subset[strand == "+", end := intron_start]
		candidate.introns.subset[strand == "-", start := intron_end]
		candidate.introns.subset[strand == "-", end := intron_end + 1]		

    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns.subset)), 
        makeGRangesFromDataFrame(as.data.frame(tmp.Introns.Group.stranded)))			
    candidate.introns.subset$gene_group_stranded[OL@from] =
        tmp.Introns.Group.stranded$gene_group[OL@to]
    candidate.introns.subset$exon_group_stranded_upstream[OL@from] =
        tmp.Introns.Group.stranded$intron_number[OL@to]		

    candidate.introns = rbind(candidate.introns, candidate.introns.subset)

    candidate.introns.subset =
        candidate.introns[is.na(exon_group_stranded_downstream)]
    candidate.introns =
        candidate.introns[!is.na(exon_group_stranded_downstream)]
    
    candidate.introns.subset[strand == "-", start := intron_start - 1]
    candidate.introns.subset[strand == "-", end := intron_start]
    candidate.introns.subset[strand == "+", start := intron_end]
    candidate.introns.subset[strand == "+", end := intron_end + 1]		

    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns.subset)), 
        makeGRangesFromDataFrame(as.data.frame(tmp.Introns.Group.stranded)))				
    candidate.introns.subset$exon_group_stranded_downstream[OL@from] =
        tmp.Introns.Group.stranded$intron_number[OL@to] + 1		

    candidate.introns = rbind(candidate.introns, candidate.introns.subset)
    
################################################################################

# Now repeat the same for unstranded condition
    candidate.introns$gene_group_unstranded = NA
    candidate.introns$exon_group_unstranded_upstream = NA
    candidate.introns$exon_group_unstranded_downstream = NA

    candidate.introns[strand == "+", start := intron_start - 1]
    candidate.introns[strand == "+", end := intron_start]
    candidate.introns[strand == "-", start := intron_end]
    candidate.introns[strand == "-", end := intron_end + 1]		
    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns)), 
        makeGRangesFromDataFrame(as.data.frame(Exons.Group.unstranded)))			
    candidate.introns$gene_group_unstranded[OL@from] =
        Exons.Group.unstranded$gene_group[OL@to]
    candidate.introns$exon_group_unstranded_upstream[OL@from] =
        Exons.Group.unstranded$exon_group[OL@to]		

    candidate.introns[strand == "-", start := intron_start - 1]
    candidate.introns[strand == "-", end := intron_start]
    candidate.introns[strand == "+", start := intron_end]
    candidate.introns[strand == "+", end := intron_end + 1]		
    
    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns)), 
        makeGRangesFromDataFrame(as.data.frame(Exons.Group.unstranded)))			
    candidate.introns$exon_group_unstranded_downstream[OL@from] =
        Exons.Group.unstranded$exon_group[OL@to]		

    # Need fix for retained_introns or sense_intronic where junction extends into the obligate introns
    tmp = makeGRangesFromDataFrame(
        as.data.frame(Exons.Group.unstranded), keep.extra.columns = TRUE)
    tmp.Introns.Group.unstranded = grlGaps(
        split(tmp, tmp$gene_group)
    )
    tmp.Introns.Group.unstranded = as.data.table(tmp.Introns.Group.unstranded)
    data.table::setnames(tmp.Introns.Group.unstranded,
        "group_name", "gene_group")
    tmp.Introns.Group.unstranded[,
        intron_number := data.table::rowid(gene_group)]

    candidate.introns.subset =
        candidate.introns[is.na(exon_group_unstranded_upstream)]
    candidate.introns =
        candidate.introns[!is.na(exon_group_unstranded_upstream)]

    candidate.introns.subset[strand == "+", start := intron_start - 1]
    candidate.introns.subset[strand == "+", end := intron_start]
    candidate.introns.subset[strand == "-", start := intron_end]
    candidate.introns.subset[strand == "-", end := intron_end + 1]		

    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns.subset)), 
        makeGRangesFromDataFrame(as.data.frame(tmp.Introns.Group.unstranded)))			
		
    candidate.introns.subset$gene_group_unstranded[OL@from] =
        tmp.Introns.Group.unstranded$gene_group[OL@to]
    candidate.introns.subset$exon_group_unstranded_upstream[OL@from] =
        tmp.Introns.Group.unstranded$intron_number[OL@to]		

    candidate.introns = rbind(candidate.introns, candidate.introns.subset)

    candidate.introns.subset =
        candidate.introns[is.na(exon_group_unstranded_downstream)]
    candidate.introns =
        candidate.introns[!is.na(exon_group_unstranded_downstream)]
    
    candidate.introns.subset[strand == "-", start := intron_start - 1]
    candidate.introns.subset[strand == "-", end := intron_start]
    candidate.introns.subset[strand == "+", start := intron_end]
    candidate.introns.subset[strand == "+", end := intron_end + 1]		

    OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(candidate.introns.subset)), 
        makeGRangesFromDataFrame(as.data.frame(tmp.Introns.Group.unstranded)))			
		
    candidate.introns.subset$exon_group_unstranded_downstream[OL@from] =
        tmp.Introns.Group.unstranded$intron_number[OL@to] + 1		

    candidate.introns = rbind(candidate.introns, candidate.introns.subset)

	# reset
    candidate.introns[,start := intron_start]
    candidate.introns[,end := intron_end]
    candidate.introns[,Event :=  paste0(seqnames, ":", intron_start, "-", intron_end, "/", strand)]
    
    write.fst(candidate.introns, file.path(reference_path,"fst","junctions.fst"))
    
    message("done\n")    
}

################################################################################
gen_irf <- function(reference_path, extra_files, genome) {
    message("Generating ref-cover.bed ...", appendLF = FALSE)    

    # Generating IRFinder-base references
    Genes = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path,"fst","Genes.fst"))
    )
    Genes.rev = Genes
    strand(Genes.rev) = ifelse(strand(Genes.rev) == "+", "-", 
        ifelse(strand(Genes.rev) == "-", "+", "*")) # Invert strand
    Genes.Extended = reduce(c(flank(Genes.rev, 5000), 
        flank(Genes.rev, 1000, start = FALSE)))

    candidate.introns = as.data.table(
        read.fst(file.path(reference_path,"fst","junctions.fst"))
    )
    Exons = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path,"fst","Exons.fst"))
    )
    Transcripts = makeGRangesFromDataFrame(
        read.fst(file.path(reference_path,"fst","Transcripts.fst"))
    )
    
    tmp.exons.exclude =  Exons[!grepl("intron", Exons$transcript_biotype)]

    candidate.introns = candidate.introns[transcript_biotype %in% 
        c("protein_coding", "processed_transcript",
        "lincRNA", "antisense", "nonsense_mediated_decay")]

    candidate.introns[, transcript_biotype := factor(transcript_biotype,
        c("protein_coding", "processed_transcript",
        "lincRNA", "antisense", "nonsense_mediated_decay"), ordered = TRUE)]
        
    if("transcript_support_level" %in% colnames(candidate.introns)) {
        # Sort by tsl first, then reverse later
        setorder(candidate.introns, transcript_biotype, transcript_support_level)    
    } else {
        setorder(candidate.introns, transcript_biotype)
    }

    introns.unique = unique(candidate.introns,
        by = c("seqnames", "start", "end", "width", "strand"))
    setorder(introns.unique, seqnames, start, end, strand)
    introns.unique = makeGRangesFromDataFrame(
        as.data.frame(introns.unique), keep.extra.columns=TRUE)
        
    setorder(candidate.introns, seqnames, start, end, strand)
    
    exclude.directional = as.data.table(tmp.exons.exclude)
    exclude.directional = unique(exclude.directional,
        by = c("seqnames", "start", "end", "width", "strand"))
    exclude.directional[, start := start - 5]
    exclude.directional[, end := end + 5]
    
    exclude.directional.reverse = copy(exclude.directional)
    exclude.directional.reverse[strand == "-", strand:= "P"]
    exclude.directional.reverse[strand == "+", strand:= "-"]
    exclude.directional.reverse[strand == "P", strand:= "+"]

    exclude.omnidirectional = GenomicRanges::GRanges(NULL)
    
    if(extra_files$MappabilityFile != "") {
        exclude.omnidirectional = c(exclude.omnidirectional,
            rtracklayer::import(extra_files$MappabilityFile, "bed")
        )
    }
    if(extra_files$BlacklistFile != "") {
        exclude.omnidirectional = c(exclude.omnidirectional,
            rtracklayer::import(extra_files$BlacklistFile, "bed")
        )
    }
    exclude.omnidirectional =
        reduce(exclude.omnidirectional, min.gapwidth = 9) # merge with any gaps <= 9
    if(length(exclude.omnidirectional) > 0) {
        introns.unique.blacklisted = findOverlaps(introns.unique,
            exclude.omnidirectional, type = "within")
        introns.unique = introns.unique[-introns.unique.blacklisted@from] 
        # clean introns by those lying completely within blacklist regions
    }
    
    introns.unique.exon.dir = findOverlaps(introns.unique, 
        makeGRangesFromDataFrame(exclude.directional), type = "within")
    introns.unique.exon.nd = findOverlaps(introns.unique, 
        makeGRangesFromDataFrame(exclude.directional),
        type = "within", ignore.strand=TRUE)

    introns.unique$known_exon_dir =
        ( seq_len(length(introns.unique)) %in% introns.unique.exon.dir@from )
    introns.unique$known_exon_nd =
        ( seq_len(length(introns.unique)) %in% introns.unique.exon.nd@from )

    introns.unique.antiover = findOverlaps(introns.unique, Genes.rev)
    introns.unique.antinear = findOverlaps(introns.unique, Genes.Extended)

    introns.unique$antiover =
        ( seq_len(length(introns.unique)) %in% introns.unique.antiover@from )
    introns.unique$antinear =
        ( seq_len(length(introns.unique)) %in% introns.unique.antinear@from )

# Now subset introns by punching holes using blacklist regions

    introns.unique$intron_width = BiocGenerics::width(introns.unique)

# Remove introns less than 50 bp:
    introns.unique = introns.unique[BiocGenerics::width(introns.unique) > 50]

# remove 5 bases from start & end
    BiocGenerics::start(introns.unique) =
        BiocGenerics::start(introns.unique) + 5
    BiocGenerics::end(introns.unique) =
        BiocGenerics::end(introns.unique) - 5

    introns.unique.dir = introns.unique
    introns.unique.nd = introns.unique
################################################################################

# Dir
    if(length(exclude.omnidirectional) > 0) {
        introns.intersect.dir = GenomicRanges::intersect(introns.unique.dir, 
            c(exclude.omnidirectional,
            makeGRangesFromDataFrame(exclude.directional)))
    } else {
        introns.intersect.dir = GenomicRanges::intersect(introns.unique.dir, 
            c(makeGRangesFromDataFrame(exclude.directional)))    
    }
    introns.intersect.ol =
        findOverlaps(introns.unique.dir, introns.intersect.dir)
    # make a GRanges same size as the number of intersections
    introns.intersect.dir.final =
        introns.intersect.dir[introns.intersect.ol@to]
    introns.intersect.dir.final$intron_id =
        introns.unique$intron_id[introns.intersect.ol@from]

    introns.unique.dir.ID = 
        split(introns.unique.dir, introns.unique.dir$intron_id)
    introns.intersect.dir.ID = 
        split(introns.intersect.dir.final,
            introns.intersect.dir.final$intron_id)
    introns.unique.dir.ID.compare = introns.unique.dir.ID[
        names(introns.unique.dir.ID) %in% names(introns.intersect.dir.ID)
    ]

# nd
    if(exists("exclude.omnidirectional")) {
        introns.intersect.nd = GenomicRanges::intersect(
            introns.unique.nd, c(
                exclude.omnidirectional, 
                makeGRangesFromDataFrame(exclude.directional), 
                makeGRangesFromDataFrame(exclude.directional.reverse)
            )
        )
    } else {
        introns.intersect.nd = GenomicRanges::intersect(
            introns.unique.nd, c(
                makeGRangesFromDataFrame(exclude.directional), 
                makeGRangesFromDataFrame(exclude.directional.reverse)
            )
        )    
    }
    introns.intersect.ol =
        findOverlaps(introns.unique.nd, introns.intersect.nd)
    # make a GRanges same size as the number of intersections
    introns.intersect.nd.final =
        introns.intersect.nd[introns.intersect.ol@to]
    introns.intersect.nd.final$intron_id =
        introns.unique$intron_id[introns.intersect.ol@from]

    introns.unique.nd.ID =
        split(introns.unique.nd, introns.unique.nd$intron_id)
    introns.intersect.nd.ID =
        split(introns.intersect.nd.final, introns.intersect.nd.final$intron_id)
    introns.unique.nd.ID.compare =
        introns.unique.nd.ID[
            names(introns.unique.nd.ID) %in% names(introns.intersect.nd.ID)
        ]

# Dir setdiff
		i.seqnames <- i.intron_start <- i.intron_end <- i.intron_width <- i.width <- i.strand <- NULL
	i.transcript_id <- i.known_exon_dir <- i.gene_group_stranded <- i.exon_group_stranded_upstream <- NULL
	i.exon_group_stranded_downstream <- exclbases <- intron_width <- inclbases <- IRFname <- NULL
  gene_name <- num_blocks <- known_exon_dir <- i.known_exon_nd <- i.antiover <- i.antinear <- NULL
  i.gene_group_unstranded <- i.exon_group_unstranded_upstream <- i.exon_group_unstranded_downstream <- NULL
  known_exon_nd <- antiover <- antinear <- V1 <- V2 <- V3 <- V6 <- V9 <- NULL
  
################################################################################  
    tmpdir.IntronCover =
        setdiff(introns.unique.dir.ID.compare, introns.intersect.dir.ID)
    # now add back introns that did not require intersection (or would have been excluded as known-exons
    tmpdir.IntronCover = c(tmpdir.IntronCover, 
        introns.unique.dir.ID[
            !(names(introns.unique.dir.ID) %in% names(introns.intersect.dir.ID))
        ],
        introns.unique.dir.ID[
            names(introns.unique.dir.ID) %in%
                introns.unique.dir$intron_id[
                    introns.unique.dir$known_exon_dir == TRUE
                ]
        ]
    )

    tmpdir.IntronCover = as.data.table(tmpdir.IntronCover)
    tmpdir.IntronCover = tmpdir.IntronCover[,
        c("seqnames", "start", "end", "strand", "width", "group_name")]
    colnames(tmpdir.IntronCover)[6] = "intron_id"

    tmpdir.IntronCover.summa = tmpdir.IntronCover
    tmpdir.IntronCover.summa[, 
        c("num_blocks", "inclbases") := list(.N, sum(width)), by = "intron_id"]
    tmpdir.IntronCover.summa = unique(
        tmpdir.IntronCover.summa[,c("intron_id", "num_blocks", "inclbases")],
        by = "intron_id")
    tmpdir.IntronCover.summa[as.data.table(introns.unique), 
        on = "intron_id", c("seqnames", "intron_start", "intron_end",
            "intron_width", "width", "strand", "gene_name", "transcript_id",
            "known_exon_dir", "GG", "EG_up", "EG_down")
        := list(i.seqnames, i.intron_start, i.intron_end, i.intron_width,
            i.width, i.strand, i.gene_name, i.transcript_id, i.known_exon_dir,
            i.gene_group_stranded, i.exon_group_stranded_upstream,
            i.exon_group_stranded_downstream)]
    tmpdir.IntronCover.summa[, exclbases := intron_width - inclbases]
        # Exclude exclbases / width > 0.3
    tmpdir.IntronCover.summa = 
        tmpdir.IntronCover.summa[exclbases / intron_width < 0.3]

    tmpdir.IntronCover = semi_join.DT(
        tmpdir.IntronCover, tmpdir.IntronCover.summa,
        by = "intron_id")
 
    tmpdir.IntronCover.summa[, IRFname := paste("dir", gene_name, intron_id,
        strand, num_blocks, 
        sprintf("%.f", intron_start - 1), sprintf("%.f", intron_end),
        inclbases, exclbases,
        ifelse(known_exon_dir, "known-exon","clean"), sep="/")]

    tmpdir.IntronCover = makeGRangesFromDataFrame(
        tmpdir.IntronCover, keep.extra.columns=TRUE)
    tmpdir.IntronCover = split(tmpdir.IntronCover, tmpdir.IntronCover$intron_id)

    names(tmpdir.IntronCover) = tmpdir.IntronCover.summa$IRFname[match(
        names(tmpdir.IntronCover), tmpdir.IntronCover.summa$intron_id)]

################################################################################
# Nondir setdiff
    tmpnd.IntronCover = 
        setdiff(introns.unique.nd.ID.compare, introns.intersect.nd.ID)
    # now add back introns that did not require intersection
    tmpnd.IntronCover = c(tmpnd.IntronCover, 
        introns.unique.nd.ID[
            !(names(introns.unique.nd.ID) %in% names(introns.intersect.nd.ID))
        ],
        introns.unique.nd.ID[
            names(introns.unique.nd.ID) %in% 
                introns.unique.nd$intron_id[
                    introns.unique.dir$known_exon_nd == TRUE
                ]
        ]
    )

    tmpnd.IntronCover = as.data.table(tmpnd.IntronCover)
    tmpnd.IntronCover = tmpnd.IntronCover[,
        c("seqnames", "start", "end", "strand", "width", "group_name")]
    colnames(tmpnd.IntronCover)[6] = "intron_id"

    tmpnd.IntronCover.summa = tmpnd.IntronCover
    tmpnd.IntronCover.summa[, 
        c("num_blocks", "inclbases") := list(.N, sum(width)), by = "intron_id"]
    tmpnd.IntronCover.summa = unique(
        tmpnd.IntronCover.summa[
            ,c("intron_id", "num_blocks", "inclbases")
        ], by = "intron_id")
        
    tmpnd.IntronCover.summa[as.data.table(introns.unique), 
        on = "intron_id", c("seqnames","intron_start", "intron_end",
            "intron_width", "width", "strand", "gene_name", "transcript_id",
            "known_exon_nd", "antiover", "antinear", "GG", "EG_up", "EG_down")
          := list(i.seqnames, i.intron_start, i.intron_end, i.intron_width,
          i.width, i.strand, i.gene_name, i.transcript_id, 
          i.known_exon_nd, i.antiover, i.antinear, i.gene_group_unstranded,
          i.exon_group_unstranded_upstream,
          i.exon_group_unstranded_downstream)]
    tmpnd.IntronCover.summa[, exclbases := intron_width - inclbases]
        # Exclude exclbases / width > 0.3
    tmpnd.IntronCover.summa = 
        tmpnd.IntronCover.summa[exclbases / intron_width < 0.3]

    tmpnd.IntronCover = semi_join.DT(
        tmpnd.IntronCover, tmpnd.IntronCover.summa,
        by = "intron_id")
 
    tmpnd.IntronCover.summa[, IRFname := paste("nd", gene_name, intron_id,
        strand, num_blocks,
        sprintf("%.f", intron_start - 1), sprintf("%.f", intron_end),
        inclbases, exclbases, sep="/")]
    # casewise naming of last condition
    tmpnd.IntronCover.summa[known_exon_nd & antiover & antinear, IRFname := 
        paste(IRFname, "known-exon+anti-over+anti-near",sep="/")]
    tmpnd.IntronCover.summa[known_exon_nd & antiover & !antinear, IRFname := 
        paste(IRFname, "known-exon+anti-over",sep="/")]
    tmpnd.IntronCover.summa[known_exon_nd & !antiover & antinear, IRFname := 
        paste(IRFname, "known-exon+anti-near",sep="/")]
    tmpnd.IntronCover.summa[!known_exon_nd & antiover & antinear, IRFname := 
        paste(IRFname, "anti-over+anti-near",sep="/")]
    tmpnd.IntronCover.summa[!known_exon_nd & !antiover & antinear, IRFname := 
        paste(IRFname, "anti-near",sep="/")]
    tmpnd.IntronCover.summa[!known_exon_nd & antiover & !antinear, IRFname := 
        paste(IRFname, "anti-over",sep="/")]
    tmpnd.IntronCover.summa[known_exon_nd & !antiover & !antinear, IRFname := 
        paste(IRFname, "known-exon",sep="/")]
    tmpnd.IntronCover.summa[!known_exon_nd & !antiover & !antinear, IRFname := 
        paste(IRFname, "clean",sep="/")]

    tmpnd.IntronCover = makeGRangesFromDataFrame(
        tmpnd.IntronCover, keep.extra.columns=TRUE)
    tmpnd.IntronCover = split(tmpnd.IntronCover, tmpnd.IntronCover$intron_id)

    names(tmpnd.IntronCover) = tmpnd.IntronCover.summa$IRFname[match(
        names(tmpnd.IntronCover), tmpnd.IntronCover.summa$intron_id)]
################################################################################
# Sort and Export out
	setorder(tmpnd.IntronCover.summa, seqnames,
        intron_start, intron_end, strand)
    tmpnd.IntronCover = tmpnd.IntronCover[tmpnd.IntronCover.summa$IRFname]
    
	setorder(tmpdir.IntronCover.summa, seqnames,
        intron_start, intron_end, strand)
    tmpdir.IntronCover = tmpdir.IntronCover[tmpdir.IntronCover.summa$IRFname]

    rtracklayer::export(tmpdir.IntronCover, file.path(reference_path,
        "tmpdir.IntronCover.bed"))
    rtracklayer::export(tmpnd.IntronCover, file.path(reference_path,
        "tmpnd.IntronCover.bed"))

# Generate final ref-cover.bed

    tmpdir.IntronCover = fread(file.path(reference_path,
        "tmpdir.IntronCover.bed"), sep="\t")
    tmpdir.IntronCover[,cat := "dir"]
    tmpnd.IntronCover = fread(file.path(reference_path,
        "tmpnd.IntronCover.bed"), sep="\t")
    tmpnd.IntronCover[,cat := "nd"]

    ref.cover = rbind(tmpdir.IntronCover, tmpnd.IntronCover)
    setorder(ref.cover, V1, V2, V3, V6, cat)
    ref.cover$cat = NULL
    ref.cover[, V9 := as.character(V9)]
    ref.cover[, V9 := "255,0,0"]

    # fwrite(ref.cover, file.path(reference_path, "ref-cover.bed"), sep="\t", col.names = F, scipen = 50)

    # cleanup
    if(file.exists(file.path(reference_path, "tmpdir.IntronCover.bed"))) {
        file.remove(file.path(reference_path, "tmpdir.IntronCover.bed"))
    }
    if(file.exists(file.path(reference_path, "tmpnd.IntronCover.bed"))){
        file.remove(file.path(reference_path, "tmpnd.IntronCover.bed"))
    }
    
# Now compile list of IRFinder introns here
    write.fst(tmpnd.IntronCover.summa, file.path(reference_path, "fst",
        "Introns.ND.fst"))
    write.fst(tmpdir.IntronCover.summa, file.path(reference_path, "fst",
        "Introns.Dir.fst"))

    message("done\n")

    message("Generating ref-ROI.bed ...", appendLF = FALSE)    

################################################################################
# ROI
    if("gene_biotype" %in% names(mcols(Transcripts))) {
        rRNA = as.data.frame(Transcripts[grepl("rRNA",
            Transcripts$gene_biotype)])
        rRNA$start = rRNA$start - 1
        rRNA$name = with(rRNA, paste("rRNA", seqnames, start, end, strand,
            transcript_id, gene_biotype, gene_id, gene_name, sep="/"))
        rRNA = rRNA[, c("seqnames", "start", "end", "name")]
    } else if("gene_type" %in% names(mcols(Transcripts))) {
        rRNA = as.data.frame(Transcripts[grepl("rRNA", Transcripts$gene_type)])
        rRNA$start = rRNA$start - 1
        rRNA$name = with(rRNA, paste("rRNA", seqnames, start, end, strand,
            transcript_id, gene_type, gene_id, gene_name, sep="/"))
        rRNA = rRNA[, c("seqnames", "start", "end", "name")]    
    } else {
        rRNA = c()
    }
    
    nonPolyA = GenomicRanges::GRanges(NULL)
    if(extra_files$nonPolyAFile != "") {
        nonPolyA = c(nonPolyA,
            rtracklayer::import(extra_files$nonPolyAFile, "bed")
        )
        nonPolyA = as.data.frame(nonPolyA)
        nonPolyA = nonPolyA[, c("seqnames", "start", "end")]
        nonPolyA$name = "NonPolyA"
    } else {
        nonPolyA = c()
    }
    
    AllChr = data.frame(seqnames = names(seqinfo(genome)),
        start = 1, end = seqinfo(genome)@seqlengths,
        names = names(seqinfo(genome)))
    AllChr = makeGRangesListFromDataFrame(AllChr, split.field = "names")
    Genes.chr = c(Genes, flank(Genes, 10000),
        flank(Genes, 10000, start = FALSE))
    Genes.chr = reduce(Genes.chr, min.gapwidth = 1000)
    Genes.chr$chr = seqnames(Genes.chr)
    Genes.chr = split(Genes.chr, Genes.chr$chr)
    
    AllChr = AllChr[names(Genes.chr)]
    AllChr.split = setdiff(AllChr, Genes.chr, ignore.strand = TRUE)
    Intergenic = unlist(AllChr.split)
    if(length(Intergenic) > 0) {
        names(Intergenic) = seq_len(length(Intergenic))
        Intergenic = as.data.frame(Intergenic)
        Intergenic = Intergenic[,c("seqnames", "start", "end")]
        Intergenic$name = paste("Intergenic", Intergenic$seqnames, sep="/")
        Intergenic$start = Intergenic$start - 1
    } else {
        Intergenic = c()
    }
    ref.ROI = rbind(rRNA, nonPolyA, Intergenic) %>%
        dplyr::arrange(seqnames, start)
    
    ref.ROI$start = ref.ROI$start - 1   # convert back to 0-based
    message("done\n")

    message("Generating ref-read-continues.ref ...", appendLF = FALSE)    
    
################################################################################
# ref-read-continues.ref
    introns.unique.readcons = rbind(
        tmpdir.IntronCover.summa[,
            c("seqnames", "intron_start", "intron_end", "strand")],
        tmpnd.IntronCover.summa[,
            c("seqnames", "intron_start", "intron_end", "strand")])
    introns.unique.readcons[, intron_start := intron_start - 1]     # 0-based
    readcons.left = introns.unique.readcons[,
        c("seqnames", "intron_start", "strand")]
    readcons.right = introns.unique.readcons[,
        c("seqnames", "intron_end", "strand")]
    colnames(readcons.left) = c("V1", "V2", "V3")
    colnames(readcons.right) = c("V1", "V2", "V3")
    readcons = rbind(readcons.left, readcons.right)
    setorder(readcons, V1, V2, V3)
    readcons = unique(readcons)
    
    message("done\n")

    message("Generating ref-sj.ref ...", appendLF = FALSE)    
    
# ref-sj.ref
    # Reload candidate introns here, as we've filtered this before
    candidate.introns = as.data.table(
        read.fst(file.path(reference_path,"fst","junctions.fst")))

    ref.sj = candidate.introns[,c("seqnames", "start", "end", "strand")]
    ref.sj = unique(ref.sj)
    ref.sj[,start := start - 1]

    message("done\n")
    
    IRF_file = file.path(reference_path, "IRFinder.ref.gz")
# Concatenate all 4 reference files into one file
    fwrite(list(">ref-cover.bed"), IRF_file, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)
    fwrite(ref.cover, IRF_file, append = TRUE, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)
    fwrite(list(">ref-read-continues.ref"), IRF_file, append = TRUE, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)
    fwrite(readcons, IRF_file, append = TRUE, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)
    fwrite(list(">ref-ROI.bed"), IRF_file, append = TRUE, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)
    fwrite(ref.ROI, IRF_file, append = TRUE, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)
    fwrite(list(">ref-sj.ref"), IRF_file, append = TRUE, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)
    fwrite(ref.sj, IRF_file, append = TRUE, 
        sep="\t", eol = "\n", col.names = F, scipen = 50)    

}

################################################################################
gen_nmd = function(reference_path, genome) {
    Exons.tr = as.data.table(
        read.fst(file.path(reference_path,"fst","Exons.fst"))
    )
    candidate.introns = as.data.table(
        read.fst(file.path(reference_path,"fst","junctions.fst"))
    )
    Misc = as.data.table(
        read.fst(file.path(reference_path,"fst","Misc.fst"))
    )
    start.DT = Misc[type == "start_codon"]

    Exons.tr = Exons.tr[transcript_id %in% start.DT$transcript_id]
    
    Exons.tr[start.DT, on = c("transcript_id"), 
      c("sc_start", "sc_end") := list(i.start, i.end)]    
    Exons.tr[start < sc_start & strand == "+", start := sc_start]
    Exons.tr[end < sc_start & strand == "+", end := sc_start]
    Exons.tr[start > sc_end & strand == "-", start := sc_end]
    Exons.tr[end > sc_end & strand == "-", end := sc_end]
    Exons.tr = Exons.tr[start < end]
    
    protein.introns = candidate.introns[
        transcript_id %in% Exons.tr$transcript_id]
    # determine here whether protein introns are CDS, 5' or 3' UTR introns
    UTR5 = Misc[type == "five_prime_utr"]
    UTR5.introns = grlGaps(
        split(makeGRangesFromDataFrame(as.data.frame(UTR5)),
            UTR5$transcript_id)
    )
    UTR5.introns = as.data.table(UTR5.introns)
    UTR3 = Misc[type == "three_prime_utr"]
    UTR3.introns = grlGaps(
        split(makeGRangesFromDataFrame(as.data.frame(UTR3)),
            UTR3$transcript_id)
    )
    UTR3.introns = as.data.table(UTR3.introns)
    
    CDS.introns = grlGaps(
        split(makeGRangesFromDataFrame(as.data.frame(Exons.tr)),
            Exons.tr$transcript_id)
    )
    CDS.introns = as.data.table(CDS.introns)
    
    protein.introns[UTR5.introns, 
        on = c("seqnames", "start", "end", "strand"),
            intron_type := "UTR5"]
    protein.introns[UTR3.introns,
        on = c("seqnames", "start", "end", "strand"),
            intron_type := "UTR3"]    
    protein.introns[CDS.introns, 
        on = c("seqnames", "start", "end", "strand"),
            intron_type := "CDS"]
    
    NMD.Table = DetermineNMD(Exons.tr, protein.introns, genome, 50)
    NMD.Table[protein.introns, on = "intron_id", 
        intron_type := i.intron_type]
    
    write.fst(NMD.Table, file.path(reference_path, "fst", "IR.NMD.fst"))
}

################################################################################
gen_splice <- function(reference_path, genome) {

    message("Annotating Splice Events\n")

    GeneOrder = as.data.table(
        read.fst(file.path(reference_path,"fst","Genes.fst"))
    )
    setorder(GeneOrder, seqnames, start, end, strand)
    candidate.introns = as.data.table(
        read.fst(file.path(reference_path,"fst","junctions.fst"))
    )
    introns.skipcoord = copy(candidate.introns)
    setorder(introns.skipcoord, gene_id, transcript_name, intron_number)
    
	introns.skipcoord[strand == "+", 
        skip_coord := ifelse(intron_number == max(intron_number), NA, 
		paste0(seqnames,":",intron_start,
            "-", data.table::shift(intron_end, 1, NA, "lead"),"/",strand)),
        by = transcript_id]
	introns.skipcoord[strand == "-", 
        skip_coord := ifelse(intron_number == max(intron_number), NA, 
		paste0(seqnames,":",data.table::shift(intron_start, 1, NA, "lead"),
		"-", intron_end,"/",strand)), by = transcript_id]
	introns.skipcoord[, 
        skip_coord_2 := data.table::shift(skip_coord, 1, NA, "lag")]

	introns.skippedJn = introns.skipcoord[skip_coord %in% Event,
        c("gene_id","gene_name","skip_coord")]
	introns.skippedJn = unique(introns.skippedJn)
 
    chrOrder = names(seqinfo(genome))
 
    message("Annotating Mutually-Exclusive-Exon Splice Events...",
        appendLF = FALSE)

	introns.search.MXE = introns.skipcoord[introns.skipcoord[,
		.I[intron_number < max(intron_number)], by = transcript_id]$V1]
	introns.search.MXE = introns.search.MXE[introns.search.MXE[
        ,.N,by = c("gene_id", "skip_coord")],
		on = c("gene_id", "skip_coord"), N := i.N]
	introns.search.MXE = introns.search.MXE[N > 1]
	introns.search.MXE.pos = introns.search.MXE[strand == "+"]
	setorder(introns.search.MXE.pos, seqnames, intron_start, intron_end)
	introns.search.MXE.neg = introns.search.MXE[strand == "-"]
	setorder(introns.search.MXE.neg, seqnames, intron_end, -intron_start)
	introns.search.MXE = rbindlist(
        list(introns.search.MXE.pos, introns.search.MXE.neg)
    )
	
	introns.search.MXE = introns.search.MXE[, 
        c("skip_coord", "gene_id", "Event", "transcript_id",
		"transcript_name", "intron_number")]
	data.table::setnames(introns.search.MXE, old = "Event", new = "Event1")
	
	introns.search.MXE2 = introns.skipcoord[,
        c("skip_coord_2","gene_id","Event","transcript_id","transcript_name")]
	data.table::setnames(introns.search.MXE2, 
        old = c("skip_coord_2","Event"), 
            new = c("skip_coord", "Event2"))
	
	introns.search.MXE[introns.search.MXE2, 
        on = c("gene_id","transcript_id","transcript_name","skip_coord"),
		Event2 := i.Event2]
	
	introns.search.MXE = unique(introns.search.MXE, 
        by= c("gene_id","skip_coord","Event1"))
	introns.search.MXE = unique(introns.search.MXE, 
        by= c("gene_id","skip_coord","Event2"))
	introns.search.MXE = introns.search.MXE[, if(.N>1) .SD, 
        by = c("gene_id","skip_coord")]
	
    if(nrow(introns.search.MXE) > 0) {
        introns.found.MXE = introns.search.MXE[ , {
            edge1 = rep(seq_len(.N), (.N:1) - 1L)
            i = 2L:(.N * (.N - 1L) / 2L + 1L)
                o = cumsum(c(0, (.N-2L):1))
            edge2 = i - o[edge1]
            .(
                gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
                Event1a = Event1[edge1], Event1b = Event1[edge2],
                Event2a = Event2[edge1], Event2b = Event2[edge2],
                transcript_id_a = transcript_id[edge1],
                transcript_id_b = transcript_id[edge2],
                transcript_name_a = transcript_name[edge1],
                transcript_name_b = transcript_name[edge2],
                intron_number_a = intron_number[edge1],
                intron_number_b = intron_number[edge2]
            )
            }, by = skip_coord]

        introns.found.MXE[, 
            gene_id := factor(gene_id,GeneOrder$gene_id,ordered=TRUE)]
        setorder(introns.found.MXE, gene_id, transcript_name_a)
        introns.found.MXE[, 
            EventName := paste0("MXE:", transcript_name_a, "-exon",
                (1 + intron_number_a),";",
                transcript_name_b,"-exon",(1 + intron_number_b))
        ]
        introns.found.MXE[, EventID := paste0("MXE#", seq_len(.N))]
        data.table::setnames(introns.found.MXE, 
            old = "skip_coord", new = "EventRegion")
        introns.found.MXE[, EventType := "MXE"]
        introns.found.MXE = introns.found.MXE[, 
            c(
                "EventType","EventID","EventName","Event1a","Event1b",
                "Event2a","Event2b", "gene_id","gene_id_b","EventRegion",
                "transcript_id_a","transcript_name_a","intron_number_a",
                "transcript_id_b","transcript_name_b","intron_number_b"
            )
        ]

        introns.found.MXE = unique(introns.found.MXE, 
            by = c("Event1a", "Event1b", "Event2a", "Event2b"))
    } else {
        introns.found.MXE = c()
    }
	message("done\n")
    gc() 
    
################################################################################
# annotate skipped junctions with two included junctions
    message("Annotating Skipped-Exon Splice Events...", appendLF = FALSE)

	introns.found.SE = introns.skippedJn[,"skip_coord"]
	introns.search.SE = candidate.introns[,
        c("gene_id","Event","transcript_id", 
            "transcript_name","intron_number")
    ]
	setnames(introns.search.SE, 
        old = c("Event", "transcript_id", "transcript_name", "intron_number"),
		new = c("skip_coord", "skip_transcript_id", "skip_transcript_name", 
            "skip_intron_number"))
	introns.found.SE[introns.search.SE, on = "skip_coord", 
		c("gene_id_b", "skip_transcript_id",
            "skip_transcript_name", "skip_intron_number") :=
		list(i.gene_id, i.skip_transcript_id,
            i.skip_transcript_name, i.skip_intron_number)]
	introns.found.SE = unique(introns.found.SE,
        by = c("gene_id_b","skip_coord"))	

	introns.search.SE2 = introns.skipcoord[,
        c("skip_coord","gene_id","Event","transcript_id","transcript_name",
			"intron_number")]
	introns.found.SE[introns.search.SE2, on = "skip_coord",
        c("gene_id", "inc_coord_upst", "inc_transcript_id",
            "inc_transcript_name", "inc_intron_number") :=
        list(i.gene_id, i.Event, i.transcript_id,
            i.transcript_name, i.intron_number)]

	introns.search.SE3 = introns.skipcoord[,
        c("skip_coord_2","gene_id","Event",
            "transcript_id","transcript_name")]
	setnames(introns.search.SE3, 
        old = c("skip_coord_2", "transcript_id"),
        new = c("skip_coord", "inc_transcript_id"))
	introns.found.SE[introns.search.SE3, 
        on = c("skip_coord", "inc_transcript_id"),
        inc_coord_downst := i.Event]

	introns.found.SE = unique(introns.found.SE, 
        by = c("gene_id","skip_coord",
            "inc_transcript_id", "inc_transcript_name"))
	introns.found.SE[, 
        gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
	setorder(introns.found.SE, gene_id, inc_transcript_name)
	introns.found.SE[, EventName := paste0(
			"SE:",inc_transcript_name,"-exon",(1 + inc_intron_number),";",
			skip_transcript_name,"-int",skip_intron_number
			)]
	introns.found.SE[, EventID := paste0("SE#",seq_len(.N))]
	introns.found.SE[, EventType := "SE"]
	introns.found.SE[, Event2b := NA]
	introns.found.SE[, EventRegion := skip_coord]
	introns.found.SE = introns.found.SE[, 
        c(
            "EventType","EventID","EventName","inc_coord_upst",
			"skip_coord","inc_coord_downst","Event2b",
			"gene_id","gene_id_b", "EventRegion",
			"inc_transcript_id","inc_transcript_name","inc_intron_number",
			"skip_transcript_id","skip_transcript_name","skip_intron_number"
        )
    ]
	setnames(introns.found.SE, 
        old = c("inc_coord_upst", "inc_coord_downst", "skip_coord"),
		new = c("Event1a", "Event2a", "Event1b"))
	setnames(introns.found.SE, 
        old = c("inc_transcript_id", "inc_transcript_name",
            "skip_transcript_id","skip_transcript_name"),
		new = c("transcript_id_a", "transcript_name_a",
            "transcript_id_b","transcript_name_b"))
	setnames(introns.found.SE, 
        new = c("intron_number_a","intron_number_b"),
        old = c("inc_intron_number", "skip_intron_number")
    )

	message("done\n")
    gc()

################################################################################
    message("Annotating Alternate First / Last Exon Splice Events...",
        appendLF = FALSE)

	# AFE/ALE

	introns.search.AFE = candidate.introns[intron_number == 1]
	introns.search.AFE[,
        seqnames := factor(seqnames, chrOrder, ordered = TRUE)]

	introns.search.AFE.pos = introns.search.AFE[strand == "+"]
	setorder(introns.search.AFE.pos, seqnames, intron_end, -intron_start)
	introns.search.AFE.pos = introns.search.AFE.pos[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number")]
	setnames(introns.search.AFE.pos, old = "intron_end", new = "intron_coord")

	introns.search.AFE.neg = introns.search.AFE[strand == "-"]
	setorder(introns.search.AFE.neg, seqnames, intron_start, intron_end)
	introns.search.AFE.neg = introns.search.AFE.neg[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number")]
	setnames(introns.search.AFE.neg, old = "intron_start", new = "intron_coord")

	introns.search.AFE = rbindlist(
        list(introns.search.AFE.pos, introns.search.AFE.neg))
	introns.search.AFE = unique(introns.search.AFE, by = "Event")

	introns.search.ALE = candidate.introns[candidate.introns[,
			.I[intron_number == max(intron_number)], by = transcript_id]$V1]
	introns.search.ALE[, 
        seqnames := factor(seqnames, chrOrder, ordered = TRUE)]

	introns.search.ALE.pos = introns.search.ALE[strand == "+"]
	setorder(introns.search.ALE.pos, seqnames, intron_start, intron_end)
	introns.search.ALE.pos = introns.search.ALE.pos[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number")]
	setnames(introns.search.ALE.pos, 
        old = "intron_start", new = "intron_coord")

	introns.search.ALE.neg = introns.search.ALE[strand == "-"]
	setorder(introns.search.ALE.neg, seqnames, intron_end, -intron_start)
	introns.search.ALE.neg = introns.search.ALE.neg[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number")]
	setnames(introns.search.ALE.neg, old = "intron_end", new = "intron_coord")

	introns.search.ALE = rbindlist(
        list(introns.search.ALE.pos, introns.search.ALE.neg))
	introns.search.ALE = unique(introns.search.ALE, by = "Event")

	introns.found.AFE = introns.search.AFE[ , {
		edge1 = rep(seq_len(.N), (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1],
            transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1],
            transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1],
            intron_number_b = intron_number[edge2]
		)
	}, by = c("seqnames", "intron_coord")]
	introns.found.AFE = introns.found.AFE[!is.na(gene_id)]
	
	introns.found.ALE = introns.search.ALE[ , {
		edge1 = rep(seq_len(.N), (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1],
            transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1],
            transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1],
            intron_number_b = intron_number[edge2]
		)
	}, by = c("seqnames", "intron_coord")]
	introns.found.ALE = introns.found.ALE[!is.na(gene_id)]

################################################################################
	introns.found.AFE = unique(introns.found.AFE, by = c("Event1a", "Event1b"))
	introns.found.AFE = introns.found.AFE[, 
        gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
	setorder(introns.found.AFE, gene_id, transcript_name_a)
	introns.found.AFE = introns.found.AFE[, 
        EventName := paste0("AFE:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.AFE = introns.found.AFE[, EventType := "AFE"]
	introns.found.AFE = introns.found.AFE[, EventRegion := Event1b]
	introns.found.AFE[, EventID := paste0("AFE#", seq_len(.N))]
	introns.found.AFE = introns.found.AFE[, 
        c("EventType", "EventID", "EventName",
		"Event1a", "Event1b", "Event2a", "Event2b",
		"gene_id", "gene_id_b", "EventRegion",
		"transcript_id_a", "transcript_name_a", "intron_number_a",
		"transcript_id_b", "transcript_name_b", "intron_number_b")]
	
	introns.found.ALE = unique(introns.found.ALE, by = c("Event1a", "Event1b"))
	introns.found.ALE = introns.found.ALE[, 
        gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
	setorder(introns.found.ALE, gene_id, transcript_name_a)
	introns.found.ALE = introns.found.ALE[, 
        EventName := paste0("ALE:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.ALE = introns.found.ALE[, EventType := "ALE"]
	introns.found.ALE = introns.found.ALE[, EventRegion := Event1b]
	introns.found.ALE[, EventID := paste0("ALE#", seq_len(.N))]
	introns.found.ALE = introns.found.ALE[,
        c("EventType", "EventID", "EventName",
		"Event1a", "Event1b", "Event2a", "Event2b",
		"gene_id", "gene_id_b", "EventRegion",
		"transcript_id_a", "transcript_name_a", "intron_number_a",
		"transcript_id_b", "transcript_name_b", "intron_number_b")]

	message("done\n")
    gc()
################################################################################
	message("Annotating Alternate 5' / 3' Splice Site Splice Events...",
        appendLF = FALSE)

    candidate.introns.ASS = candidate.introns[
        !is.na(exon_group_stranded_upstream) &
        !is.na(exon_group_stranded_downstream)]
    
    setnames(candidate.introns.ASS, 
        old = c("exon_group_stranded_upstream",
            "exon_group_stranded_downstream"),
        new = c("exon_groups_start","exon_groups_end"))
    
	introns.search.A5SS = copy(candidate.introns.ASS)
	introns.search.A5SS[, 
        seqnames := factor(seqnames, chrOrder, ordered = TRUE)]

	introns.search.A5SS.pos = introns.search.A5SS[strand == "+"]
	setorder(introns.search.A5SS.pos, seqnames, intron_end, -intron_start)
	introns.search.A5SS.pos = introns.search.A5SS.pos[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start","exon_groups_end")]
	setnames(introns.search.A5SS.pos, 
        old = "intron_end", new = "intron_coord")

	introns.search.A5SS.neg = introns.search.A5SS[strand == "-"]
	setorder(introns.search.A5SS.neg, seqnames, intron_start, intron_end)
	introns.search.A5SS.neg = introns.search.A5SS.neg[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start","exon_groups_end")]
	setnames(introns.search.A5SS.neg, 
        old = "intron_start", new = "intron_coord")

	introns.search.A5SS = rbindlist(
        list(introns.search.A5SS.pos, introns.search.A5SS.neg))
	introns.search.A5SS = unique(introns.search.A5SS, by = "Event")

	introns.search.A3SS = copy(candidate.introns.ASS)
	introns.search.A3SS[, 
        seqnames := factor(seqnames, chrOrder, ordered = TRUE)]

	introns.search.A3SS.pos = introns.search.A3SS[strand == "+"]
	setorder(introns.search.A3SS.pos, seqnames, intron_start, intron_end)
	introns.search.A3SS.pos = introns.search.A3SS.pos[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start","exon_groups_end")]
    setnames(introns.search.A3SS.pos, 
        old = "intron_start", new = "intron_coord")

	introns.search.A3SS.neg = introns.search.A3SS[strand == "-"]
	setorder(introns.search.A3SS.neg, seqnames, intron_end, -intron_start)
	introns.search.A3SS.neg = introns.search.A3SS.neg[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id", 
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start","exon_groups_end")]
	data.table::setnames(introns.search.A3SS.neg, 
        old = "intron_end", new = "intron_coord")

	introns.search.A3SS = rbindlist(
        list(introns.search.A3SS.pos, introns.search.A3SS.neg))
	introns.search.A3SS = unique(introns.search.A3SS, by = "Event")

################################################################################
	introns.found.A5SS = introns.search.A5SS[ , {
		edge1 = rep(seq_len(.N), (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1],
            transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1],
            transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1],
            intron_number_b = intron_number[edge2],
            exon_groups_start_a = exon_groups_start[edge1],
            exon_groups_start_b = exon_groups_start[edge2], 
            exon_groups_end_a = exon_groups_end[edge1],
            exon_groups_end_b = exon_groups_end[edge2]
		)
		}, by = intron_coord]
	introns.found.A5SS = introns.found.A5SS[!is.na(gene_id)]
	
	introns.found.A3SS = introns.search.A3SS[ , {
		edge1 = rep(seq_len(.N), (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1],
            transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1],
            transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1],
            intron_number_b = intron_number[edge2],
            exon_groups_start_a = exon_groups_start[edge1],
            exon_groups_start_b = exon_groups_start[edge2], 
            exon_groups_end_a = exon_groups_end[edge1],
            exon_groups_end_b = exon_groups_end[edge2]
		)
		}, by = intron_coord]
	introns.found.A3SS = introns.found.A3SS[!is.na(gene_id)]

	introns.found.A5SS = unique(introns.found.A5SS,
        by = c("Event1a", "Event1b"))
    # filter by same exon group starts and ends:
   	introns.found.A5SS = introns.found.A5SS[
        exon_groups_start_a == exon_groups_start_b]
   	introns.found.A5SS = introns.found.A5SS[
        exon_groups_end_a == exon_groups_end_b]

	introns.found.A5SS = introns.found.A5SS[, 
        gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
	setorder(introns.found.A5SS, gene_id, transcript_name_a)
	introns.found.A5SS = introns.found.A5SS[, 
        EventName := paste0("A5SS:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.A5SS = introns.found.A5SS[, EventType := "A5SS"]
	introns.found.A5SS = introns.found.A5SS[, EventRegion := Event1b]
	introns.found.A5SS = introns.found.A5SS[!introns.found.AFE,
        on = c("Event1a", "Event1b")]
	introns.found.A5SS[, EventID := paste0("A5SS#", seq_len(.N))]
	introns.found.A5SS = introns.found.A5SS[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]

	introns.found.A3SS = unique(introns.found.A3SS,
        by = c("Event1a", "Event1b"))
    # filter by same exon group starts and ends:
   	introns.found.A3SS = introns.found.A3SS[
        exon_groups_start_a == exon_groups_start_b]
   	introns.found.A3SS = introns.found.A3SS[
        exon_groups_end_a == exon_groups_end_b]
	introns.found.A3SS = introns.found.A3SS[,
        gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
	setorder(introns.found.A3SS, gene_id, transcript_name_a)
	introns.found.A3SS = introns.found.A3SS[,
        EventName := paste0("A3SS:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.A3SS = introns.found.A3SS[, EventType := "A3SS"]
	introns.found.A3SS = introns.found.A3SS[, EventRegion := Event1b]
	introns.found.A3SS = introns.found.A3SS[!introns.found.ALE, 
        on = c("Event1a", "Event1b")]
	introns.found.A3SS[, EventID := paste0("A3SS#", seq_len(.N))]
	introns.found.A3SS = introns.found.A3SS[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    
    is_valid_splice_type <- function(x) !is.null(x) && nrow(x) > 0
    
    tmp_AS = list(introns.found.MXE, introns.found.SE,
            introns.found.AFE, introns.found.ALE,
            introns.found.A5SS, introns.found.A3SS)
    tmp_AS <- base::Filter(is_valid_splice_type, tmp_AS)
    AS_Table = rbindlist(tmp_AS)    
    
################################################################################
    if(nrow(AS_Table) > 0) {
        candidate.introns.order = copy(candidate.introns)
        if(!("transcript_support_level" %in% colnames(candidate.introns))) {
          candidate.introns.order[, transcript_support_level := "NA"]
        }
        candidate.introns.order[, is_protein_coding := !is.na(protein_id)]
        candidate.introns.order[, by = "transcript_id",
            is_last_intron := (intron_number == max(intron_number))]
        
        AS_Table.search.a = AS_Table[,
            c("EventType", "EventID", "Event1a", "Event2a")]
        AS_Table.search.a[,Event := Event1a]
        AS_Table.search.a = candidate.introns.order[AS_Table.search.a,
            on = "Event", 
            c("EventType","EventID", "Event1a", "Event2a", "transcript_id",
                "transcript_support_level", "is_protein_coding",
                "is_last_intron", "intron_number")]
        setnames(AS_Table.search.a, "intron_number", "in_1a")
        AS_Table.search.a = 
            AS_Table.search.a[EventType !=  "AFE" | in_1a == 1]
        AS_Table.search.a = 
            AS_Table.search.a[EventType !=  "ALE" | is_last_intron]
        AS_Table.search.a[,Event := Event2a]
        AS_Table.search.a[is.na(Event),Event := Event1a]
        AS_Table.search.a = candidate.introns.order[AS_Table.search.a, 
            on = c("Event", "transcript_id", "transcript_support_level"),
            c("EventType","EventID", "Event1a", "Event2a", "transcript_id",
                "transcript_support_level", "is_protein_coding",
                "is_last_intron","in_1a", "intron_number")]
        AS_Table.search.a = AS_Table.search.a[!is.na(intron_number)]
        setnames(AS_Table.search.a, "intron_number", "in_2a")
        
        AS_Table.search.b = AS_Table[,
            c("EventType", "EventID", "Event1b", "Event2b")]
        AS_Table.search.b[,Event := Event1b]
        AS_Table.search.b = candidate.introns.order[AS_Table.search.b,
            on = "Event", 
            c("EventType","EventID", "Event1b", "Event2b", "transcript_id",
                "transcript_support_level", "is_protein_coding",
                "is_last_intron","intron_number")]
        setnames(AS_Table.search.b, "intron_number", "in_1b")
        AS_Table.search.b = 
            AS_Table.search.b[EventType !=  "AFE" | in_1b == 1]
        AS_Table.search.b = 
            AS_Table.search.b[EventType !=  "ALE" | is_last_intron]
        AS_Table.search.b[,Event := Event2b]
        AS_Table.search.b[is.na(Event),Event := Event1b]
        AS_Table.search.b = candidate.introns.order[AS_Table.search.b, 
            on = c("Event", "transcript_id", "transcript_support_level"),
            c("EventType","EventID", "Event1b", "Event2b", "transcript_id",
                "transcript_support_level", "is_protein_coding",
                "is_last_intron", "in_1b", "intron_number")]
        AS_Table.search.b = AS_Table.search.b[!is.na(intron_number)]
        setnames(AS_Table.search.b, "intron_number", "in_2b")

        AS_Table.search.a[candidate.introns.order, 
            on = "transcript_id", transcript_name := i.transcript_name]
        AS_Table.search.b[candidate.introns.order, 
            on = "transcript_id", transcript_name := i.transcript_name]
        setorder(AS_Table.search.a, 
            transcript_support_level, -is_protein_coding, transcript_name)
        setorder(AS_Table.search.b, 
            transcript_support_level, -is_protein_coding, transcript_name)
        AS_Table.find.a = unique(AS_Table.search.a, by = "EventID")
        AS_Table.find.a = AS_Table.find.a[AS_Table[, "EventID"], on = "EventID"]
        AS_Table.find.b = unique(AS_Table.search.b, by = "EventID")
        AS_Table.find.b = AS_Table.find.b[AS_Table[, "EventID"], on = "EventID"]
################################################################################        
        AS_Table$transcript_id_a = AS_Table.find.a$transcript_id
        AS_Table$transcript_name_a = AS_Table.find.a$transcript_name
        AS_Table$intron_number_a = AS_Table.find.a$in_1a
        AS_Table$transcript_id_b = AS_Table.find.b$transcript_id
        AS_Table$transcript_name_b = AS_Table.find.b$transcript_name
        AS_Table$intron_number_b = AS_Table.find.b$in_1b
        
        AS_Table[EventType == "MXE", 
            EventName := paste0("MXE:", transcript_name_a,"-exon",
                as.character(as.numeric(intron_number_a) + 1), ";", 
                transcript_name_b,"-exon",
                as.character(as.numeric(intron_number_b) + 1)
            )
        ]
        AS_Table[EventType == "SE", 
            EventName := paste0("SE:", transcript_name_a,"-exon",
                as.character(as.numeric(intron_number_a) + 1), ";",
                transcript_name_b,"-int",
                as.character(as.numeric(intron_number_b))
            )
        ]
        AS_Table[EventType == "AFE", 
            EventName := paste0("AFE:", transcript_name_a,"-exon1;", 
            transcript_name_b,"-exon1")]
        AS_Table[EventType == "ALE", 
            EventName := paste0("ALE:", transcript_name_a, "-exon", 
                as.character(as.numeric(intron_number_a) + 1), ";", 
                transcript_name_b, "-exon",
                as.character(as.numeric(intron_number_b) + 1)
            )
        ]
        AS_Table[EventType == "A5SS", 
            EventName := paste0("A5SS:", transcript_name_a,"-exon", 
                as.character(as.numeric(intron_number_a)), ";", 
                transcript_name_b,"-exon",
                as.character(as.numeric(intron_number_b))
            )
        ]
        AS_Table[EventType == "A3SS", 
            EventName := paste0("A3SS:", transcript_name_a,"-exon", 
                as.character(as.numeric(intron_number_a + 1)), ";",
                transcript_name_b,"-exon",
                as.character(as.numeric(intron_number_b + 1))
            )
        ]
        write.fst(as.data.frame(AS_Table), file.path(reference_path,"fst","Splice.fst"))
        
        setnames(AS_Table.search.a, 
            old = c("Event1a", "Event2a", "in_1a", "in_2a"),
            new = c("Event1", "Event2", "in_1", "in_2"))
        AS_Table.search.a[, isoform := "A"]
        setnames(AS_Table.search.b,
            old = c("Event1b", "Event2b", "in_1b", "in_2b"),
            new = c("Event1", "Event2", "in_1", "in_2"))
        AS_Table.search.b[, isoform := "B"]
        
        write.fst(as.data.frame(rbind(AS_Table.search.a, AS_Table.search.b)), 
          file.path(reference_path,"fst","Splice.options.fst"))
        message("done\n")
    } else {
        message("no splice events found\n")
    }
    gc()

}

################################################################################
gen_splice_proteins <- function(reference_path, genome) {
    message("Translating Alternate Splice Peptides...", appendLF = FALSE)

    AS_Table <- as.data.table(
        read.fst(file.path(reference_path,"fst","Splice.fst"))
    )
    Proteins.Splice = as.data.table(
        read.fst(file.path(reference_path,"fst","Proteins.fst"))
    )

    AS_Table.Extended = copy(AS_Table)
    
    Proteins.Splice$exon_number = as.numeric(Proteins.Splice$exon_number)
    # make phase easier for me to understand
    Proteins.Splice[, phase := -phase %% 3]
    # Upstream applicable for MXE, SE, ALE, A3SS
    Upstream = AS_Table[EventType %in% c("MXE", "SE", "ALE", "A3SS")]

    # Do A
    Upstream.A = Upstream[,
        c("EventID", "transcript_id_a", "intron_number_a")]
    Upstream.A[, c("transcript_id", "exon_number") := 
        list(transcript_id_a, intron_number_a)]
    # left_join with Exons
    Upstream.A = Proteins.Splice[Upstream.A, 
        on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Upstream.A.gr = makeGRangesFromDataFrame(as.data.frame(
        na.omit(Upstream.A)), keep.extra.columns = TRUE)
    Upstream.A.seq = getSeq(genome, Upstream.A.gr)
    Upstream.A[!is.na(seqnames),seq := as.character(Upstream.A.seq)]
    # Trim sequence by phase
    seq = substr(Upstream.A$seq[!is.na(Upstream.A$seqnames)],
        1 + (3 - Upstream.A$phase[!is.na(Upstream.A$seqnames)]) %% 3,
        nchar(Upstream.A$seq[!is.na(Upstream.A$seqnames)]))        
    # trim last n bases
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    # translate
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Upstream.A[!is.na(seqnames), AA_seq := as.character(prot)]
    AS_Table.Extended[EventType %in% c("MXE", "SE", "ALE", "A3SS"),
        AA_upstr.A := Upstream.A$AA_seq]

    # repeat for B:
    Upstream.B = Upstream[,
        c("EventID", "transcript_id_b", "intron_number_b")]
    Upstream.B[, c("transcript_id", "exon_number") :=
        list(transcript_id_b, intron_number_b)]
    # left_join with Exons
    Upstream.B = Proteins.Splice[Upstream.B, 
        on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Upstream.B.gr = makeGRangesFromDataFrame(as.data.frame(
        na.omit(Upstream.B)), keep.extra.columns = TRUE)
    Upstream.B.seq = getSeq(genome, Upstream.B.gr)
    Upstream.B[!is.na(seqnames),seq := as.character(Upstream.B.seq)]
    # Trim sequence by phase
    seq = substr(Upstream.B$seq[!is.na(Upstream.B$seqnames)],
        1 + (3 - Upstream.B$phase[!is.na(Upstream.B$seqnames)]) %% 3,
        nchar(Upstream.B$seq[!is.na(Upstream.B$seqnames)]))        
    # trim last n bases
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    # translate
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Upstream.B[!is.na(seqnames), AA_seq := as.character(prot)]
    AS_Table.Extended[EventType %in% c("MXE", "SE", "ALE", "A3SS"),
        AA_upstr.B := Upstream.B$AA_seq]
        
################################################################################    
    # Do downstream seq before casette:
    Downstream = AS_Table[EventType %in% c("MXE", "SE", "AFE", "A5SS")]
    # Add EventType as exon_number is conditional on this
    Downstream.A = Downstream[, 
        c("EventType", "EventID", "transcript_id_a", "intron_number_a")]
    Downstream.A[, c("transcript_id", "exon_number") := 
        list(transcript_id_a, intron_number_a)]
    # Modify downstream exon number
    Downstream.A[EventType %in% c("MXE", "SE"), 
        exon_number := exon_number + 2]
    Downstream.A[EventType %in% c("AFE", "A5SS"), 
        exon_number := exon_number + 1]
    # left_join with Exons
    Downstream.A = Proteins.Splice[Downstream.A, 
        on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Downstream.A.gr = makeGRangesFromDataFrame(as.data.frame(
        na.omit(Downstream.A)), keep.extra.columns = TRUE)
    Downstream.A.seq = getSeq(genome, Downstream.A.gr)
    Downstream.A[!is.na(seqnames),seq := as.character(Downstream.A.seq)]
    # Trim sequence by phase
    seq = substr(Downstream.A$seq[!is.na(Downstream.A$seqnames)],
        1 + (3 - Downstream.A$phase[!is.na(Downstream.A$seqnames)]) %% 3,
        nchar(Downstream.A$seq[!is.na(Downstream.A$seqnames)]))        
    # trim last n bases
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    # translate
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Downstream.A[!is.na(seqnames), AA_seq := as.character(prot)]
    AS_Table.Extended[EventType %in% c("MXE", "SE", "AFE", "A5SS"),
        AA_downstr.A := Downstream.A$AA_seq]
    # B:
    Downstream.B = Downstream[, 
        c("EventType", "EventID", "transcript_id_b", "intron_number_b")]
    Downstream.B[, c("transcript_id", "exon_number") := 
        list(transcript_id_b, intron_number_b)]
    # Modify downstream exon number: Note SE is different for B
    Downstream.B[EventType %in% c("MXE"), exon_number := exon_number + 2]
    Downstream.B[EventType %in% c("SE", "AFE", "A5SS"), 
        exon_number := exon_number + 1]
    # left_join with Exons
    Downstream.B = Proteins.Splice[Downstream.B, 
        on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Downstream.B.gr = makeGRangesFromDataFrame(as.data.frame(
        na.omit(Downstream.B)), keep.extra.columns = TRUE)
    Downstream.B.seq = getSeq(genome, Downstream.B.gr)
    Downstream.B[!is.na(seqnames),seq := as.character(Downstream.B.seq)]
    # Trim sequence by phase
    seq = substr(Downstream.B$seq[!is.na(Downstream.B$seqnames)],
        1 + (3 - Downstream.B$phase[!is.na(Downstream.B$seqnames)]) %% 3,
        nchar(Downstream.B$seq[!is.na(Downstream.B$seqnames)]))        
    # trim last n bases
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    # translate
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Downstream.B[!is.na(seqnames), AA_seq := as.character(prot)]
    AS_Table.Extended[EventType %in% c("MXE", "SE", "AFE", "A5SS"),
        AA_downstr.B := Downstream.B$AA_seq]
        
################################################################################
    # Casette A
    Casette.A = AS_Table[, 
        c("EventType", "EventID", "transcript_id_a", "intron_number_a")]
    Casette.A[, c("transcript_id", "exon_number") := 
        list(transcript_id_a, intron_number_a)]
    Casette.A[EventType %in% c("MXE", "SE", "ALE", "A3SS"), 
        exon_number := exon_number + 1]

    Casette.A = Proteins.Splice[Casette.A, 
        on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Casette.A.gr = makeGRangesFromDataFrame(as.data.frame(
        na.omit(Casette.A)), keep.extra.columns = TRUE)
    Casette.A.seq = getSeq(genome, Casette.A.gr)
    Casette.A[!is.na(seqnames),casette_seq := as.character(Casette.A.seq)]

    data.table::setnames(Casette.A, "phase", "phase_casette")
# Add nucleotides from upstream and downstream
    Casette.A = Upstream.A[Casette.A, on = "EventID", 
        c("EventID", "phase_casette", "casette_seq", "seq")]
    data.table::setnames(Casette.A, "seq", "upstr_seq")
    Casette.A = Downstream.A[Casette.A, on = "EventID", 
        c("EventID", "phase_casette", "casette_seq", "upstr_seq", "seq")]
    data.table::setnames(Casette.A, "seq", "Downstr_seq")
    
# Construct extended casette sequence:
    Casette.A[, casette_seq_extended := casette_seq]
    # Trim casette_seq_extended if upstream sequence does not exists
    Casette.A[!is.na(phase_casette) & is.na(upstr_seq),
        casette_seq_extended := substr(
            casette_seq_extended, 
            phase_casette + 1, 
            nchar(casette_seq_extended)
        )
    ]    
    Casette.A[!is.na(phase_casette) & phase_casette > 0 & !is.na(upstr_seq), 
        casette_seq_extended := paste0(
            substr(upstr_seq, 
                nchar(upstr_seq) + 1 - phase_casette, 
                nchar(upstr_seq)),
            casette_seq_extended
        )
    ]
    Casette.A[nchar(casette_seq_extended) %% 3 > 0 & !is.na(Downstr_seq), 
        casette_seq_extended := paste0(
            casette_seq_extended,
            substr(Downstr_seq, 1, 3 - (nchar(casette_seq_extended) %% 3))
        )
    ]
# Translate:
    seq = Casette.A$casette_seq_extended[
        !is.na(Casette.A$casette_seq_extended)]
    # trim out-of-phase to be tidy:
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Casette.A[!is.na(casette_seq_extended), 
        AA_seq := as.character(prot)]
    AS_Table.Extended[, AA_casette.A := Casette.A$AA_seq]
    
################################################################################
    # Casette B
    Casette.B = AS_Table[EventType != "SE", 
        c("EventType", "EventID", "transcript_id_b", "intron_number_b")]
    Casette.B[, c("transcript_id", "exon_number") := 
        list(transcript_id_b, intron_number_b)]
    Casette.B[EventType %in% c("MXE", "ALE", "A3SS"), 
        exon_number := exon_number + 1]

    Casette.B = Proteins.Splice[Casette.B, 
        on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Casette.B.gr = makeGRangesFromDataFrame(as.data.frame(
        na.omit(Casette.B)), keep.extra.columns = TRUE)
    Casette.B.seq = getSeq(genome, Casette.B.gr)
    Casette.B[!is.na(seqnames),casette_seq := as.character(Casette.B.seq)]

    data.table::setnames(Casette.B, "phase", "phase_casette")
# Add nucleotides from upstream and downstream
    Casette.B = Upstream.B[Casette.B, on = "EventID", 
        c("EventID", "phase_casette", "casette_seq", "seq")]
    data.table::setnames(Casette.B, "seq", "upstr_seq")
    Casette.B = Downstream.B[Casette.B, on = "EventID", 
        c("EventID", "phase_casette", "casette_seq", "upstr_seq", "seq")]
    data.table::setnames(Casette.B, "seq", "Downstr_seq")
    
# Construct extended casette sequence:
    Casette.B[, casette_seq_extended := casette_seq]
    # Trim casette_seq_extended if upstream sequence does not exists
    Casette.B[!is.na(phase_casette) & is.na(upstr_seq),
        casette_seq_extended := substr(
            casette_seq_extended, phase_casette + 1, 
            nchar(casette_seq_extended)
        )
    ]    
    Casette.B[!is.na(phase_casette) & phase_casette > 0 & !is.na(upstr_seq), 
        casette_seq_extended := paste0(
            substr(upstr_seq, nchar(upstr_seq) + 1 - phase_casette, 
                nchar(upstr_seq)
            ),
            casette_seq_extended
        )
    ]
    Casette.B[nchar(casette_seq_extended) %% 3 > 0 & !is.na(Downstr_seq), 
        casette_seq_extended := paste0(casette_seq_extended,
            substr(Downstr_seq, 1, 3 - (nchar(casette_seq_extended) %% 3))
        )
    ]
################################################################################    
# Translate:
    seq = Casette.B$casette_seq_extended[
        !is.na(Casette.B$casette_seq_extended)]
    # trim out-of-phase to be tidy:
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Casette.B[!is.na(casette_seq_extended), AA_seq := as.character(prot)]
    AS_Table.Extended[EventType != "SE", AA_casette.B := Casette.B$AA_seq]

    AS_Table.Extended[, AA_full.A := ""]
    AS_Table.Extended[!is.na(AA_upstr.A), 
        AA_full.A := paste0(AA_full.A, AA_upstr.A)]
    AS_Table.Extended[!is.na(AA_casette.A), 
        AA_full.A := paste0(AA_full.A, AA_casette.A)]
    AS_Table.Extended[!is.na(AA_downstr.A), 
        AA_full.A := paste0(AA_full.A, AA_downstr.A)]
    AS_Table.Extended[, AA_full.B := ""]
    AS_Table.Extended[!is.na(AA_upstr.B), 
        AA_full.B := paste0(AA_full.B, AA_upstr.B)]
    AS_Table.Extended[!is.na(AA_casette.B), 
        AA_full.B := paste0(AA_full.B, AA_casette.B)]
    AS_Table.Extended[!is.na(AA_downstr.B), 
        AA_full.B := paste0(AA_full.B, AA_downstr.B)]
    
    write.fst(as.data.frame(AS_Table.Extended), 
        file.path(reference_path,"fst","Splice.Extended.fst"))

    message("done\n")
}

# FetchAH_FTP <- function(ah_record, reference_path, localHub = FALSE, ah = AnnotationHub(localHub = localHub)) {
  # ah.record = ah[names(ah) == ah_record]
  
  # Best to fetch from URL
  # if(!dir.exists(file.path(reference_path, "resource"))) dir.create(file.path(reference_path, "resource"))

  # transcripts.gtf = file.path(normalizePath(reference_path), "resource", "transcripts.gtf")       
  # if(!file.exists(transcripts.gtf)) {
      # url = ah.record$sourceurl
      # assert_that(substr(url,1,3) == "ftp",
        # msg = paste("ftp site not found for", ah_record))
      # urlfile = basename(url)
      # if(substr(urlfile, nchar(urlfile) -6, nchar(urlfile)) == ".gtf.gz") {
        # download.file(url, destfile = paste(transcripts.gtf, "gz", sep="."))
        # GEOquery::gunzip(paste(transcripts.gtf, "gz", sep="."))
        # transcripts.gtf = paste0(transcripts.gtf, ".gz")
      # } else if(substr(urlfile, nchar(urlfile) - 3, nchar(urlfile)) == ".gtf") {
        # download.file(url, destfile = transcripts.gtf)
      # } else {
        # warning("sourceurl entry for AnnotationHub resource is not a valid gtf.gz or gtf file")
      # }
  # }
  # rtracklayer::import(transcripts.gtf, "gtf")
# }

#' Fetch genome / transcriptome reference from AnnotationHub and writes to reference_path
#'
#' @export
GenerateMappabilityReads <- function(fasta = "genome.fa", ah_genome = "", 
    reference_path, read_len = 70, read_stride = 10, error_pos = 35,
    verbose = FALSE, localHub = FALSE, 
    ah = AnnotationHub(localHub = localHub)) {

  if(ah_genome != "") {
    assert_that(substr(ah_genome,1,2) == "AH",
        msg = paste(ah_genome, "- Given AnnotationHub reference is incorrect"))

    genome = FetchAH(ah_genome, ah = ah, verbose = verbose)

    if(!dir.exists(file.path(reference_path, "resource"))) {
        dir.create(file.path(reference_path, "resource"))
    }
    fasta_file = file.path(reference_path, "resource", 
        paste(ah_genome, "fa", sep="."))
        
    if(verbose) message(paste("Preparing to export as", fasta_file))
    # Write fasta to file
    genome.DNA = rtracklayer::import(genome)
    gc()
    genome = Biostrings::replaceAmbiguities(genome)
    gc()
    
    rtracklayer::export(genome, fasta_file, "fasta")
    if(verbose) message(paste("Successful export of", fasta_file))
  } else {
    assert_that(file.exists(normalizePath(fasta)),
        msg = paste("Given genome file", normalizePath(fasta), "not found"))
    fasta_file = fasta
  }

# Run map read generator:
    run_IRFinder_GenerateMapReads(normalizePath(fasta_file), 
      file.path(normalizePath(reference_path), paste("MappabilityReads", 
        ifelse(ah_genome == "", "genome", ah_genome),"fa", sep=".")), 
        read_len, read_stride, error_pos)
}

#' @export
GenerateMappabilityBED = function(BAM = "", out.bed, threshold = 4) {
  assert_that(file.exists(BAM),
    msg = paste(BAM, "BAM file does not exist"))
  assert_that(dir.exists(dirname(out.bed)),
    msg = paste(dirname(out.bed), "directory does not exist"))
    
  return(
    IRF_GenerateMappabilityRegions(normalizePath(BAM), 
        file.path(normalizePath(dirname(out.bed)), out.bed),
        threshold = threshold)
  )
}


#' Builds reference files used by IRFinder / NxtIRF.
#'
#' @description
#' This function builds the reference required by the IRFinder engine, as well
#' as access-ready refined splice annotation data for NxtIRF. This reference
#' can be created using either:
#' 1. AnnotationHub genome and gene annotation (Ensembl): supply the names of
#'    the genome sequence `ah_genome` and gene annotations `ah_transcriptome`
#'    (see example below), or
#' 2. User-supplied FASTA and GTF file.
#' 
#' @param fasta: The path to the user-supplied genome fasta file
#' @param gtf: The path to the user-supplied transcript gtf file
#' @param ah_genome: The name of the AnnotationHub record containing the genome 2bit file.
#'   Leave blank to use user-supplied `fasta` file.
#' @param ah_transcriptome: The name of the AnnotationHub record containing the transcript gtf file
#'   Leave blank to use user-supplied `gtf` file.
#' @param reference_path: The directory to store the reference files
#' @param genome_type: Allows `BuildReference()` to select default `nonPolyARef` and `MappabilityRef`
#'   for selected genomes. Allowed options include: 'hg38', 'hg19', 'mm9', 'mm10'. Leave blank to
#'   supply custom `nonPolyARef` and `MappabilityRef` files
#' @param nonPolyARef: A BED file (3 unnamed columns containing chromosome, start and end 
#'   coordinates) of regions defining known non-polyadenylated transcripts. This file is used for
#'   QC analysis of IRFinder-processed files to measure Poly-A enrichment quality of samples.
#'   Leave blank to not use a `nonPolyARef` file (or to use default - see `genome_type`)
#' @param MappabilityRef: A BED file (3 unnamed columns containing chromosome, start and end 
#'   coordinates) of poorly-mapped regions due to repeat elements in the genome. We recommend
#'   using the default Mappability files supplied (see `genome_type`). Alternately, this
#'   reference can be generated by running `GenerateMappabilityReads()` on the genome sequence,
#'   followed by alignment of the produced fasta file to an aligner of choice (e.g. STAR, HISAT2).
#'   The aligned sequences (as BAM file) should then be analysed using `GenerateMappabilityBED()`,
#'   which will provide the Mappability file to be used here.
#' @param BlackListRef: A BED file (3 unnamed columns containing chromosome, start and end 
#'   coordinates) of regions to be otherwise excluded from IR analysis. Leave blank to not use a 
#'   `BlackListRef` file.
#' @param UseExtendedTranscripts: Should IRFinder include non-protein-coding transcripts such as
#'   anti-sense and lincRNAs? Setting `FALSE` (default IRFinder) will exclude transcripts other than 
#'   `protein_coding` and `processed_transcript` transcripts from IR analysis.
#' @param localHub: See `?AnnotationHub::AnnotationHub()`. Setting `TRUE` will run `AnnotationHub()`
#'   in offline mode
#' @param ah: An AnnotationHub object containing the records `ah_genome` and/or `ah_transcriptome`
#'   records to be used.
#' @param BPPARAM: A BPPARAM object for BiocParallel multi-threading.
#' @return Nothing. The created reference will be written to the given directory. This includes:
#' * `reference_path`/IRFinder.ref.gz: A gzipped text file containing collated IRFinder references 
#'   to be used as input for the IRFinder analysis
#' * `reference_path`/fst/: Contains fst files for subsequent easy access to NxtIRF generated
#'   references
#' * `reference_path`/resource/: Contains a TwoBitFile generated by this function for easy
#'   subsequent access to the genome. Only applies if user-supplied `fasta` file is used.
#' @examples
#' ## create an AnnotationHub object
#' ah = AnnotationHub::AnnotationHub()
#' 
#' ## filter for entries containing c("Homo Sapiens", "Ensembl", "release-94")
#' records = AnnotationHub::query(ah, c("Homo Sapiens", "Ensembl", "release-94"))
#' 
#' ## display summary of filtered records
#' print(records)
#' @seealso [GenerateMappabilityReads()], [GenerateMappabilityBED()], [AnnotationHub::AnnotationHub())]
#' @md
#' @export
BuildReference <- function(fasta = "genome.fa", gtf = "transcripts.gtf", 
    ah_genome = "", ah_transcriptome = "", reference_path = "./Reference",
    genome_type = "", nonPolyARef = "", MappabilityRef = "", BlacklistRef = "",
    UseExtendedTranscripts = TRUE,
    localHub = FALSE, ah = AnnotationHub(localHub = localHub), 
    BPPARAM = BiocParallel::bpparam()) {

################################################################################
    assert_that(genome_type != "",
        msg = paste("genome_type not specified.",
        "This should be either one of 'hg38', 'hg19', 'mm10', 'mm9', or",
        "'other'. If 'other', please provide a nonPolyARef file or leave",
        "blank to omit polyA profiling.")
    )

    extra_files = fetch_genome_defaults(genome_type, nonPolyARef, 
        MappabilityRef, BlacklistRef)
        
    assert_that(
        tryCatch(ifelse(normalizePath(dirname(reference_path)) != "",TRUE, TRUE),
            error = function(e) FALSE),
            msg = paste("Base path of ", reference_path, " does not exist"))

    prep_ref_path(reference_path)
    
    N = 9
    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(1/N, message = "Reading Reference Files")
    }
    
    if(ah_genome != "") {
        assert_that(substr(ah_genome,1,2) == "AH",
            msg = "Given genome AnnotationHub reference is incorrect")
        genome = FetchAH(ah_genome, ah = ah)
        message("done\n")
        fasta_file = ""
    } else {
        assert_that(file.exists(normalizePath(fasta)),
            msg = paste("Given genome fasta file", normalizePath(fasta),
                "not found"))

        # Convert genome to TwoBitFile for easy access:
        genome = Biostrings::readDNAStringSet(fasta)
        # convert to local 2bit for better memory management
        if(!dir.exists(file.path(reference_path, "resource"))) {
            dir.create(file.path(reference_path, "resource"))
        }
        rtracklayer::export(genome, file.path(reference_path, "resource", 
            "genome.2bit"), "2bit")
        message("Genome converted to Twobit file\n")
    
        message("Connecting to genome TwoBitFile...", appendLF = FALSE)
        genome = rtracklayer::TwoBitFile(file.path(reference_path, "resource",
            "genome.2bit"))
        gc()
        message("done\n")
        fasta_file = fasta
    }
    if(ah_transcriptome != "") {
        assert_that(substr(ah_transcriptome,1,2) == "AH",
            msg = "Given transcriptome AnnotationHub reference is incorrect")
        gtf.gr = FetchAH(ah_transcriptome, ah = ah)
        message("done\n")
        gtf_file = ""
    } else {
        assert_that(file.exists(gtf),
            msg = paste("Given transcriptome gtf file", gtf,
                "not found"))
        gtf_file = gtf

        message("Reading source GTF file...", appendLF = FALSE)
        gtf.gr = rtracklayer::import(gtf, "gtf")
        message("done\n")
    }
    gc()

    chrOrder = names(seqinfo(genome))

    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Processing gtf file")
    }
    
    process_gtf(gtf.gr, reference_path)
   
    # To save memory, remove original gtf
    rm(gtf.gr)
    gc()

    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Processing introns")
    }

    process_introns(reference_path, genome, UseExtendedTranscripts)
    gc()


    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Generating IRFinder Reference")
    }

# Finished annotating introns, now use it to build reference:
    gen_irf(reference_path, extra_files, genome)

# Annotate IR-NMD


    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Annotating IR-NMD")
    }
    gen_nmd(reference_path, genome)
    gc()
    
# Annotating Alternative Splicing Events
        

    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Annotating Splice Events")
    }
    
    gen_splice(reference_path, genome)
    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Finalising Splice Event Annotations")
    }
    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Translating AS Peptides")
    }	
    gen_splice_proteins(reference_path, genome)

    message("Splice Annotations finished\n")

	message("Reference build finished")
    if(!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::incProgress(1/N, message = "Reference build finished")
    }
  
  # create settings.csv only after everything is finalised
	settings.list = list(
        fasta_file = fasta_file, gtf_file = gtf_file, 
        ah_genome = ah_genome, ah_transcriptome = ah_transcriptome,
		reference_path = reference_path, genome_type = genome_type, 
        nonPolyARef = nonPolyARef, MappabilityRef = MappabilityRef,
        BlacklistRef = BlacklistRef,
		UseExtendedTranscripts = UseExtendedTranscripts
    )
	saveRDS(settings.list, file.path(reference_path, "settings.Rds"))
	
}


DetermineNMD <- function(exon_list, intron_list, genome, threshold = 50) {
  # transcript_list can be a GRanges, data.frame, or data.table coerce-able to a GRanges object
  # All members of transcript must have the same transcript-id and must not completely overlap
  
  # Also it is presumed the first nucleotide is the beginning of the stop codon
  
  assert_that(is(exon_list, "GRanges") | is(exon_list, "data.frame") | 
    is(exon_list, "data.table"), 
      msg = "exon_list must be a GRanges, data.frame, or data.table coerce-able to a GRanges object")

  assert_that(is(exon_list, "GRanges") || 
    all(c("seqnames", "start", "end", "strand") %in% colnames(exon_list)),
    msg = "exon_list must be a GRanges, data.frame, or data.table coerce-able to a GRanges object")
  
  assert_that(is(intron_list, "GRanges") | is(intron_list, "data.frame") | 
    is(intron_list, "data.table"), 
      msg = "intron_list must be a GRanges, data.frame, or data.table coerce-able to a GRanges object")

  assert_that(is(intron_list, "GRanges") || 
    all(c("seqnames", "start", "end", "strand") %in% colnames(intron_list)),
    msg = "intron_list must be a GRanges, data.frame, or data.table coerce-able to a GRanges object")
  
  message("Calculating IR-NMD")
  elem_number <- transcript_id <- is_last_elem <- AA <- stop_pos <- splice_len <- stop_to_EJ <- 
	i.stop_pos <- i.splice_len <- i.stop_to_EJ <- splice_is_NMD <- splice_start_to_last_EJ <- 
	splice_stop_to_last_EJ <- intron_id <- type <- first_exon_id <- i.first_exon_id <- intron_pos <- 
	i.elem_number <- i.seq <- IRT_len <- use_short <- IRT_is_NMD <- i.IRT_len <- i.use_short <- 
	i.IRT_is_NMD <- NULL
	
  exon.DT = as.data.table(exon_list)
  intron.DT = as.data.table(intron_list)

  exon.DT = exon.DT[, c("seqnames", "start", "end", "strand", "transcript_id")]

  exon.gr = makeGRangesFromDataFrame(as.data.frame(exon.DT))
  exon.DT[, seq := as.character(getSeq(genome, exon.gr))]

  gc()
  
  # Easy bit: test whether spliced transcript is NMD-inducing

  exon.MLE.DT = copy(exon.DT)

  setorder(exon.MLE.DT, start)
  exon.MLE.DT[, elem_number := data.table::rowid(transcript_id)]
  exon.MLE.DT[strand == "-", elem_number := max(elem_number) + 1 - elem_number, by = "transcript_id"]
  exon.MLE.DT[, by = "transcript_id", is_last_elem := (elem_number == max(elem_number))]

  exon.MLE.DT = exon.MLE.DT[is_last_elem == FALSE]
  
  # sort by order
  setorder(exon.MLE.DT, transcript_id, elem_number)

  exon.MLE.DT = exon.MLE.DT[, c("transcript_id", "seq")]
  # construct sequences:
  splice = exon.MLE.DT[, lapply(.SD, paste0, collapse = ""), by = "transcript_id"]
  splice[is.na(rowSums(stringr::str_locate(seq, "N"))),
    AA := as.character(suppressWarnings(Biostrings::translate(as(seq, "DNAStringSet"))))]
  
  # Find nucleotide position of first stop codon
  splice[, stop_pos := stringr::str_locate(AA, "\\*")[,1] * 3 - 2]
  splice[, splice_len := nchar(seq)]
  splice[!is.na(AA), stop_to_EJ := splice_len - stop_pos]

  # Hard bit: test whether IR transcript is NMD-inducing

  intron.DT = intron.DT[, c("seqnames", "start", "end", "strand", "transcript_id", "intron_id")]  
  final = intron.DT[, c("intron_id", "transcript_id")]

  final[splice, on = "transcript_id", c("splice_stop_pos", "splice_start_to_last_EJ", "splice_stop_to_last_EJ") := 
    list(i.stop_pos, i.splice_len, i.stop_to_EJ)]
  
  final[, splice_is_NMD := ifelse(splice_start_to_last_EJ - splice_stop_to_last_EJ >= threshold, TRUE, FALSE)]
  final[is.na(splice_stop_to_last_EJ), splice_is_NMD := FALSE]
  final[is.na(splice_start_to_last_EJ), splice_is_NMD := NA]
  
  i_partition = seq(1, nrow(intron.DT), by = 10000)
  i_partition = append(i_partition, nrow(intron.DT) + 1)
    
  pb = txtProgressBar(max = length(i_partition) - 1, style = 3)	
  for(i in seq_len(length(i_partition) - 1)) {
    setTxtProgressBar(pb, i)
    intron.part = intron.DT[seq(i_partition[i], i_partition[i+1] - 1)]
    intron.part$type = "intron"
# join exons with introns to determine phase of intron    
    exon.DT.skinny = exon.DT[, -("seq")]
    intron.part.upstream = as.data.table(rbind(
      as.data.frame(intron.part) %>% dplyr::select("transcript_id", "intron_id",
          "seqnames", "start", "end", "strand", "type"),
      dplyr::left_join(as.data.frame(intron.part) %>% dplyr::select("transcript_id", "intron_id"), 
          as.data.frame(exon.DT.skinny), by = "transcript_id") %>% dplyr::mutate(type = "exon")
    ))

    setorder(intron.part.upstream, start)
    intron.part.upstream[, elem_number := data.table::rowid(intron_id)]
    intron.part.upstream[strand == "-", elem_number := max(elem_number) + 1 - elem_number, by = "intron_id"]
    
    # remove introns upstream to first available exon
    intron.part.upstream[intron.part.upstream[type == "exon", first_exon_id := min(elem_number), by = "intron_id"], 
      on = "intron_id", first_exon_id := i.first_exon_id]    
    intron.part.upstream = intron.part.upstream[type == "exon" | elem_number > first_exon_id]
    # trim exons downstream of intron
    intron.part.upstream[intron.part.upstream[type == "intron"], 
      on = "intron_id", intron_pos := i.elem_number]    

    intron.part.upstream = intron.part.upstream[!is.na(intron_pos)]
    intron.part.upstream = intron.part.upstream[elem_number < intron_pos | type == "intron" ]
    
    # Retrieve first 1000 bases of intron sequence. Most introns may be screened using this method. Save overhead
    
    intron.part.short = intron.part.upstream[type == "intron"]
    intron.part.short[ strand == "+" & end - start > 1000, end := start + 1000 ]
    intron.part.short[ strand == "-" & end - start > 1000, start := end - 1000 ]
    
    intron.short.gr = makeGRangesFromDataFrame(as.data.frame(intron.part.short))
    intron.part.short[, seq := as.character(getSeq(genome, intron.short.gr))]    
    intron.part.upstream[exon.DT, on = c("transcript_id", "seqnames", "start", "end", "strand"), seq := i.seq]
    intron.part.upstream[intron.part.short, on = c("intron_id", "type"), seq := i.seq]

    # Test introns by translation
    setorder(intron.part.upstream, transcript_id, elem_number)
    intron.part.upstream = intron.part.upstream[, c("intron_id", "seq")]

    IRT = intron.part.upstream[, lapply(.SD, paste0, collapse = ""), by = "intron_id"]
    IRT[is.na(rowSums(stringr::str_locate(seq, "N"))),
      AA := as.character(suppressWarnings(Biostrings::translate(as(seq, "DNAStringSet"))))]
    
    # Find nucleotide position of first stop codon
    IRT[, stop_pos := stringr::str_locate(AA, "\\*")[,1] * 3 - 2]
    IRT[, IRT_len := nchar(seq)]
    IRT[!is.na(AA), stop_to_EJ := IRT_len - stop_pos]
    IRT[, use_short := TRUE]
    
    IRT[, IRT_is_NMD := ifelse(stop_to_EJ >= threshold, TRUE, FALSE)]
    IRT[is.na(stop_pos), IRT_is_NMD := FALSE]
    IRT[is.na(IRT_len), IRT_is_NMD := NA]
    
    # Annotate into final
    final[IRT, on = "intron_id", c("IRT_stop_pos", "IRT_start_to_last_EJ", "IRT_stop_to_last_EJ", 
      "IRT_use_short", "IRT_is_NMD") := 
      list(i.stop_pos, i.IRT_len, i.stop_to_EJ, i.use_short, i.IRT_is_NMD)]
    
    # Now exclude NMD-TRUE introns from full analysis
    
    intron.part = intron.part[!(intron_id %in% IRT$intron_id[IRT$IRT_is_NMD == TRUE])]
    
    # Exclude introns preceding any ORF exons:
    intron.part = intron.part[intron_id %in% IRT$intron_id]
    
    intron.gr = makeGRangesFromDataFrame(as.data.frame(intron.part))
    intron.part[, seq := as.character(getSeq(genome, intron.gr))]

    intron.MLE.DT = as.data.table(
      rbind(
        as.data.frame(intron.part) %>% dplyr::select("transcript_id", "intron_id",
          "seqnames", "start", "end", "strand", "seq"),
        dplyr::left_join(as.data.frame(intron.part) %>% dplyr::select("transcript_id", "intron_id"), 
          as.data.frame(exon.DT), by = "transcript_id")
      )
    )

    # sort
    setorder(intron.MLE.DT, start)
    intron.MLE.DT[, elem_number := data.table::rowid(transcript_id)]
    intron.MLE.DT[strand == "-", elem_number := max(elem_number) + 1 - elem_number, by = "transcript_id"]
    intron.MLE.DT[, by = "transcript_id", is_last_elem := (elem_number == max(elem_number))]

    intron.MLE.DT = intron.MLE.DT[is_last_elem == FALSE]
    
    # sort by order
    setorder(intron.MLE.DT, transcript_id, elem_number)

    intron.MLE.DT = intron.MLE.DT[, c("intron_id", "seq")]

    IRT = intron.MLE.DT[, lapply(.SD, paste0, collapse = ""), by = "intron_id"]
    IRT[is.na(rowSums(stringr::str_locate(seq, "N"))),
      AA := as.character(suppressWarnings(Biostrings::translate(as(seq, "DNAStringSet"))))]
    
    # Find nucleotide position of first stop codon
    IRT[, stop_pos := stringr::str_locate(AA, "\\*")[,1] * 3 - 2]
    IRT[, IRT_len := nchar(seq)]
    IRT[!is.na(AA), stop_to_EJ := IRT_len - stop_pos]
    IRT[, use_short := FALSE]
    
    IRT[, IRT_is_NMD := ifelse(stop_to_EJ >= threshold, TRUE, FALSE)]
    IRT[is.na(stop_pos), IRT_is_NMD := FALSE]
    IRT[is.na(IRT_len), IRT_is_NMD := NA]
    
    final[IRT, on = "intron_id", c("IRT_stop_pos", "IRT_start_to_last_EJ", "IRT_stop_to_last_EJ", 
      "IRT_use_short", "IRT_is_NMD") := 
      list(i.stop_pos, i.IRT_len, i.stop_to_EJ, i.use_short, i.IRT_is_NMD)]
    gc()
  }
	setTxtProgressBar(pb,i)
	close(pb)
  message("done\n")
  
  return(final)
}

annotateIntronGroups <- function(candidate.introns, Exons.Group, stranded = TRUE) {

}

grlGaps<-function(grl) {
	psetdiff(unlist(range(grl),use.names=TRUE),grl)
}
