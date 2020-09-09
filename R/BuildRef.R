#' Builds NxtIRF Reference
#'
#' Builds reference files used by NxtIRF. Requires Ensembl GTF file 'transcripts.gtf' and genome file 'genome.fa'
#'
#' @param fasta: Either the path to the genome fasta file, or the name of the AnnotationHub object containing the TwoBitFile of the genome (beginning with "AH...")
#' @param gtf: Either the path to the transcript gtf file, or the name of the AnnotationHub object containing the GRanges object of the transcript reference (beginning with "AH...")
#' @param RefPath: The directory in which genome.fa and transcripts.gtf is found. If not provided, a temp directory will be generated.
#' @param stopcodons: A vector containing stop codons for your species of interest. Default is c("TAA","TAG","TGA")
#' @param PolyA_length: (Default = 15) Minimum length of repeat adenosine sequence detection.
#' @param PolyA_mismatch: (Default = 1) Maximum number of allowed mismatch of repeat adenosine sequence detection.
#' @return Nothing. The created NxtAnnotation object is saved in 'RefPath/Ref_Filename'
#' @export
NxtIRF.BuildReference <- function(fasta = "genome.fa", gtf = "transcripts.gtf",
    RefPath = "", stopcodons = c("TAA","TAG","TGA"), PolyA_length = 15, PolyA_mismatch = 1) {

    # Use temp directory if RefPath not provided
    if(RefPath == "") {
        RefPath = tempdir()
    }

	# Benchmarking
	irf = new("NxtProject")

	irf <- NxtIRF.StartBenchmarkTimer(irf, "BuildReference")

	irf <- NxtIRF.StartBenchmarkTimer(irf)

    if(substr(fasta,1,2) != "AH" | substr(gtf,1,2) != "AH") {
        UseAH = FALSE
        assertthat::assert_that(file.exists(paste(RefPath,fasta,sep="/")),
            msg = paste(paste(RefPath,fasta,sep="/"), " not found. Aborting reference build\n"))
        
        assertthat::assert_that(file.exists(paste(RefPath,gtf,sep="/")),
            msg = paste(paste(RefPath,gtf,sep="/"), " not found. Aborting reference build\n"))

        # Require Samtools
        # if(!file.exists(paste(RefPath,"genome.fa.fai",sep="/"))) {
            # message("Building genome index", appendLF = F)
            # indexFa(paste(RefPath,"genome.fa",sep="/"))
            # message("\u2713\n")
        # } else {
            # message("genome.fa.fai detected\n")
        # }
        
        # GenomeIndex = fread(paste(RefPath,"genome.fa.fai",sep="/"))
        # chrOrder <- GenomeIndex$V1

        genome = FaFile(paste(RefPath,fasta,sep="/"))

        # Import gtf
        message("Reading source GTF file...", appendLF = F)
        gtf.gr = rtracklayer::import(paste(RefPath,gtf,sep="/"))
        
    } else {
        UseAH = TRUE
        ah = AnnotationHub::AnnotationHub()

        genome = ah[[fasta]]        
        gtf.gr = ah[[gtf]]
    }
    
    chrOrder = levels(seqnames(gtf.gr))
    GenomeIndex = as.data.frame(gtf.gr) %>% dplyr::group_by(seqnames) %>% dplyr::summarize(V2 = max(end)) %>%
        dplyr::ungroup() %>% dplyr::rename(V1 = seqnames)   
    
	# Convert Gencode to Ensembl format (if detected)
	if(all(c("gene_type", "transcript_type", "ccdsid") %in% names(values(gtf.gr)))) {
		names(values(gtf.gr))[names(values(gtf.gr)) == "gene_type"] = "gene_biotype"
		names(values(gtf.gr))[names(values(gtf.gr)) == "transcript_type"] = "transcript_biotype"
		names(values(gtf.gr))[names(values(gtf.gr)) == "ccdsid"] = "ccds_id"
	}

	Supposed.Columns = c("source", "type", "score", "phase", 
		"gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype",
		"transcript_support_level", "exon_number", "exon_id", "ccds_id", "tag")
	
	if(all(Supposed.Columns %in% names(values(gtf.gr)))) {			
		gtf.gr = gtf.gr[,Supposed.Columns]
	} else {
		Missing.Columns = Supposed.Columns[!(Supposed.Columns %in% names(values(gtf.gr)))]
		message("transcripts.gtf is missing the following required columns:", Missing.Columns)
		return(irf)
	}
	message("\u2713\n")

	message("Generating gene annotations...", appendLF = F)

	annogenes = gtf.gr[gtf.gr$type == "gene",c("gene_id","gene_name","gene_biotype", "source")]
		annogenes.unstranded = annogenes
		strand(annogenes.unstranded) = "*"
		annogenes$unstranded_group = CreateGroups.GR(annogenes.unstranded)
		rm(annogenes.unstranded)
		annogenes$stranded_group = CreateGroups.GR(annogenes)

	saveRDS(annogenes, paste(RefPath,"annogenes.Rds",sep="/"))
	message("\u2713\n")

	annogenes = readRDS(paste(RefPath,"annogenes.Rds",sep="/"))

	message("Generating transcript annotations...", appendLF = F)

		annotrans = gtf.gr[gtf.gr$type == "transcript", c("source",
		"gene_id","transcript_id", 
			"gene_name","gene_biotype", "transcript_name", "transcript_biotype",
			"ccds_id", "tag", "transcript_support_level")
		]
		annotrans <- sortSeqlevels(annotrans)
		annotrans <- sort(annotrans)

		annotrans = as.data.table(annotrans)

		annotrans[, tsl_priority := suppressWarnings(as.numeric(substr(annotrans$transcript_support_level,1,1)))]
		annotrans[is.na(tsl_priority), tsl_priority := Inf]
		annotrans[tsl_priority < 1, tsl_priority := Inf]

		annotrans[, cds_priority := ifelse(is.na(ccds_id),Inf,
				ifelse(nchar(ccds_id) == 0,0,as.numeric(substr(ccds_id,5,nchar(ccds_id)))
				))]

		annotrans[, orf_priority := ifelse(transcript_biotype == "protein_coding",1,
				ifelse(transcript_biotype == "processed_transcript",2,3))]
		
		StartCodon = gtf.gr[gtf.gr$type == "start_codon", "transcript_id"]
		StartCodon = as.data.table(StartCodon)
		StartCodon[strand == "+", ATGcoord := min(start), by = transcript_id]
		StartCodon[strand == "-", ATGcoord := max(end), by = transcript_id]

		StopCodon = gtf.gr[gtf.gr$type == "stop_codon", "transcript_id"]
		StopCodon = as.data.table(StopCodon)
		StopCodon[strand == "+", StopCodon := min(start), by = transcript_id]
		StopCodon[strand == "-", StopCodon := max(end), by = transcript_id]
		
		annotrans = annotrans[StartCodon, on = "transcript_id", ATGcoord := i.ATGcoord]
		annotrans = annotrans[StopCodon, on = "transcript_id", StopCodon := i.StopCodon]

		setorder(annotrans,  tsl_priority, cds_priority, seqnames, start, end, strand)

		annotrans = makeGRangesFromDataFrame(setDF(annotrans), keep.extra.columns=T)

	saveRDS(annotrans, paste(RefPath,"annotrans.Rds",sep="/"))

	annotrans = readRDS(paste(RefPath,"annotrans.Rds",sep="/"))

	message("\u2713\n")

	message("Generating exon annotations...", appendLF = F)
	
	annoexons = gtf.gr[gtf.gr$type == "exon", c("source",
	"gene_id","transcript_id", "exon_number",
		"gene_name","gene_biotype", "transcript_name", "transcript_biotype",
		"exon_id", "ccds_id")
	]
		annoexons <- sortSeqlevels(annoexons)
		annoexons <- sort(annoexons)
		annoexons$exon_groups = CreateGroups.GRL(annoexons, annoexons$gene_id)
		
		annoexons_protein_coding = which(annoexons$transcript_biotype == "protein_coding")	
		annoexons$exon_groups_prot = 0
		annoexons$exon_groups_prot[annoexons_protein_coding] = 
			CreateGroups.GRL(annoexons[annoexons_protein_coding,], annoexons$gene_id[annoexons_protein_coding])
		
		annoexons_CCDS = which(!is.na(annoexons$ccds_id))
		annoexons$exon_groups_ccds = 0
		annoexons$exon_groups_ccds[annoexons_CCDS] = 
			CreateGroups.GRL(annoexons[annoexons_CCDS,], annoexons$gene_id[annoexons_CCDS])
	
# ORF data
		
	saveRDS(annoexons, paste(RefPath,"annoexons.Rds",sep="/"))
		
	message("\u2713\n")
	
	message("Exon Annotations finished. ", appendLF = F)

	irf <- NxtIRF.Benchmark(irf)
	irf <- NxtIRF.StartBenchmarkTimer(irf)

	annoexons = readRDS(paste(RefPath,"annoexons.Rds",sep="/"))

	# message("Generating intron annotations...", appendLF = F)
	message("Generating intron annotations: Initializing...", appendLF = F)

	annointrons = grlGaps(
		split(annoexons, annoexons$transcript_id) # annoexons.split
	)
		annointrons = unlist(annointrons)
		annointrons$transcript_id = names(annointrons)
		annointrons = sortSeqlevels(annointrons)
		annointrons = sort(annointrons)
	
	annotrans.summary = as.data.table(annotrans)

	annointrons = as.data.table(annointrons)
	annointrons[annotrans.summary, on = "transcript_id",
		c("gene_id", "gene_name", "transcript_name") := list(i.gene_id, i.gene_name, i.transcript_name)]

	annointrons[,intron_number := rowid(transcript_id)]
	annointrons[strand == "-", intron_number := max(intron_number) + 1 - intron_number, by = "transcript_id"]
	annointrons[,c("intron_start", "intron_end") := list(start - 1, end + 1)]

	annoexons.summary = as.data.table(annoexons)
	annoexons.summary = annoexons.summary[!duplicated(annoexons.summary[,c("seqnames", "start", "end", "strand", "gene_id")])]	
	annoexons.summary[, c("intron_start", "intron_end") := list(end, start)]
	
	annointrons[annoexons.summary, on = c("seqnames", "intron_start", "strand", "gene_id"), 
		c("exon_groups_start","exon_groups_prot_start","exon_groups_ccds_start") := 
		list(i.exon_groups, i.exon_groups_prot, i.exon_groups_ccds)]

	annointrons[annoexons.summary, on = c("seqnames", "intron_end", "strand", "gene_id"), 
		c("exon_groups_end","exon_groups_prot_end","exon_groups_ccds_end") := 
		list(i.exon_groups, i.exon_groups, i.exon_groups_ccds)]		

	annointrons[, last_intron := (intron_number == max(intron_number)), by = "transcript_id"]
	annointrons[, first_intron := (intron_number == 1)]
	annointrons[, only_intron := (last_intron == TRUE & intron_number == 1)]


	# Revised UTR annotation
	
	gtf.DT = as.data.table(gtf.gr)
	gtf.transcript = gtf.DT[type == "exon"]
	gtf.start = gtf.DT[type == "start_codon"]
	gtf.stop = gtf.DT[type == "stop_codon"]
	gtf.transcript[, c("min_left", "max_right") := 
		list(min(start), max(end)), by = c("transcript_id", "strand")]
	gtf.start[, c("ATG_start", "ATG_end") := 
		list(min(start), max(end)), by = c("transcript_id", "strand")]
	gtf.stop[, c("stopcodon_start", "stopcodon_end") := 
		list(min(start), max(end)), by = c("transcript_id", "strand")]

	gtf.transcript = gtf.transcript[, c("transcript_id", "strand", "min_left", "max_right")]
	gtf.transcript = gtf.transcript[!duplicated(gtf.transcript)]

	gtf.transcript[gtf.start, on = c("transcript_id","strand"),
		c("ATG_start", "ATG_end") := list(i.ATG_start, i.ATG_end)]
	gtf.transcript[gtf.stop, on = c("transcript_id","strand"),
		c("stopcodon_start", "stopcodon_end") := list(i.stopcodon_start, i.stopcodon_end)]

	gtf.transcript[, utr5_start := ifelse(!is.finite(ATG_start) | !is.finite(ATG_end),
		as.double(NA), ifelse(strand == "+", as.double(min_left), as.double(ATG_start + 1)))]
	gtf.transcript[, utr5_end := ifelse(!is.finite(ATG_start) | !is.finite(ATG_end),
		as.double(NA), ifelse(strand == "+", as.double(ATG_start - 1), as.double(max_right)))]
		
	gtf.transcript[, utr3_start := ifelse(!is.finite(stopcodon_start) | !is.finite(stopcodon_end),
		as.double(NA), ifelse(strand == "+", as.double(stopcodon_end + 1), as.double(min_left)))]
	gtf.transcript[, utr3_end := ifelse(!is.finite(stopcodon_start) | !is.finite(stopcodon_end),
		as.double(NA), ifelse(strand == "+", as.double(max_right), as.double(stopcodon_start - 1)))]
		
	annointrons[gtf.transcript, on = c("transcript_id"),
		c("ATG_start", "ATG_end", "stopcodon_start", "stopcodon_end",
			"utr5_start", "utr5_end", "utr3_start", "utr3_end") := list(
			ATG_start, ATG_end, stopcodon_start, stopcodon_end,
			utr5_start, utr5_end, utr3_start, utr3_end)]

	annointrons[, intron_protein_coding := !is.na(ATG_start) & !is.na(stopcodon_start)]

	annointrons[, intron_type := ifelse(intron_protein_coding == TRUE, 
		ifelse(strand == "+", 
		ifelse(intron_start > ATG_end & intron_end < stopcodon_start,"intron_CDS",
			ifelse(intron_start < ATG_end,"intron_utr5","intron_utr3")),
		# condition negative strand
		ifelse(intron_start < ATG_start & intron_end > stopcodon_end,"intron_CDS",
			ifelse(intron_start > ATG_start,"intron_utr5","intron_utr3"))),
		"intron_non_protein_coding_transcript")]

	message("\u2713")

# Intron Frames
	message("Generating intron annotations: Determining intron reading frames...", appendLF = F)
	
	annointrons.frame = annointrons[intron_type == "intron_CDS"]
	setorder(annointrons.frame, transcript_id, intron_number)

	annointrons.frame[strand == "+", 
		intron_frame := (intron_start - ATG_end - cumsum(data.table::shift(width, 1, fill = 0, type = "lag"))
			) %% 3 + 1, by = "transcript_id"]

	annointrons.frame[strand == "-", 
		intron_frame := (ATG_start - intron_end - cumsum(data.table::shift(width, 1, fill = 0, type = "lag"))
			) %% 3 + 1, by = "transcript_id"]			
	
	# annointrons.frame$intron_5pframe_start = 0
	# annointrons.frame$intron_5pframe_end = 0
	
	annointrons.frame[strand == "+", intron_5pframe_start := 
		intron_start + 2 - intron_frame]
	annointrons.frame[strand == "+", intron_5pframe_end := 
		intron_end + 1 - (intron_end - intron_5pframe_start + 2) %% 3]

	annointrons.frame[strand == "-", intron_5pframe_end := 
		intron_end - 2 + intron_frame]
	annointrons.frame[strand == "-", intron_5pframe_start := 
		intron_start - 1 + (intron_5pframe_end - intron_start + 2) %% 3]		

	message("\u2713")
			
	message("Retrieving intron sequences to determine stop codons...")

	annointrons.frame.condensed = annointrons.frame[, 
		c("seqnames", "intron_start", "intron_end", "intron_5pframe_start", "intron_5pframe_end", "strand")]
	annointrons.frame.condensed = unique(annointrons.frame.condensed)
    
	gr = with(annointrons.frame.condensed,
		makeGRangesFromDataFrame(data.frame(
		seqnames = seqnames, start = intron_5pframe_start,
			end = intron_5pframe_end,strand = strand)))

	# get ~1000 introns at a time to save memory
	iter_num = 1 + round(nrow(annointrons.frame.condensed) / 10000)
		# ~ 45 Mb per segment
			
	object = NxtIRF.SplitVector(seq_len(nrow(annointrons.frame.condensed)), iter_num)

	# stopcodons = c("TAA","TAG","TGA")
	pb = txtProgressBar(max = length(object), style = 3)
	for(q in seq_len(length(object))) {
		setTxtProgressBar(pb,q)
		indices = object[[q]]
		intron.seq = DNAStringSet(getSeq(genome,gr[indices]))
		
		matches = list()
		for(i in seq_len(length(stopcodons))) {
			matches[[i]] = startIndex(vmatchPattern(stopcodons[i],intron.seq))
			matches[[i]] = sapply(matches[[i]], FindMinInFrame)
		}
		
		annointrons.frame.condensed[indices,
			stopcodon_inframe := matrixStats::rowMins(do.call(cbind,matches), na.rm = T)]
		annointrons.frame.condensed[indices,
			frameshift := ifelse((intron_end - intron_start - 1) %% 3 == 0,"in_frame","frameshift") ]
		annointrons.frame.condensed[indices,
			readthrough := ifelse(is.infinite(stopcodon_inframe) & frameshift == "in_frame","in_frame_readthrough","-") ]
			
	}
	
	annointrons.frame[annointrons.frame.condensed, 
		on =  c("seqnames", "intron_start", "intron_end", "intron_5pframe_start", "intron_5pframe_end", "strand"),
		c("stopcodon_inframe", "frameshift", "readthrough") := list(i.stopcodon_inframe, i.frameshift, i.readthrough)]

	annointrons[annointrons.frame, on = c("seqnames", "intron_start", "intron_end", "strand", "transcript_id"),
		c("intron_frame","intron_5pframe_start","intron_5pframe_end",
		"stopcodon_inframe","frameshift","readthrough") :=
		list(i.intron_frame, i.intron_5pframe_start, i.intron_5pframe_end, 
		i.stopcodon_inframe, i.frameshift, i.readthrough)]

	annointrons$intron_type = factor(annointrons$intron_type, c("intron_CDS", "intron_utr5", "intron_utr3",
		"intron_non_protein_coding_transcript", ordered = TRUE))

	# Set transcript order
	
	annointrons$transcript_id = factor(annointrons$transcript_id, annotrans$transcript_id, ordered = TRUE)
	setorder(annointrons, intron_type, transcript_id, intron_number)

	annointrons[, EventName := paste0("IR:",
		transcript_name,"-","Int",intron_number)]
	annointrons[, Event := paste0(seqnames,":",
		intron_start,"-",intron_end, "/", strand)]
	annointrons = annointrons[, c("EventName",              "Event",         "seqnames",
		"start",                  "end",                   
		"width",                  "strand",
		"intron_start",          "intron_end",                "transcript_id",        
		"gene_id",               "gene_name",             "transcript_name",      
		"intron_number",                   
		"exon_groups_start",     "exon_groups_prot_start","exon_groups_ccds_start",
		"exon_groups_end",       "exon_groups_prot_end",  "exon_groups_ccds_end", 
		"last_intron",           "first_intron",          "only_intron",          
		"ATG_start",             "ATG_end",               "stopcodon_start",      
		"stopcodon_end",         "utr5_start",            "utr5_end",             
		"utr3_start",            "utr3_end",              "intron_protein_coding",
		"intron_type",
		"intron_frame", "stopcodon_inframe","frameshift","readthrough"
		)]
	message("\u2713\n")
	
	# Annotate intron retention overhangs for single cell
	annointrons[strand == "+", Overhang_5prime := paste0(seqnames,":", intron_start,"/",strand)]
	annointrons[strand == "-", Overhang_5prime := paste0(seqnames,":", intron_end - 1,"/",strand)]
	
	annointrons[strand == "-", Overhang_3prime := paste0(seqnames,":", intron_start,"/",strand)]
	annointrons[strand == "+", Overhang_3prime := paste0(seqnames,":", intron_end - 1,"/",strand)]
	
	# PolyA detection and exclusion
    
	PolyN = PolyA_length
	PolyMM = PolyA_mismatch
	message("\nDetecting repeat PolyA sequences (at least ",PolyN," nt stretches, allow ",PolyMM," mismatch)...\n")
	Genome.gr = makeGRangesFromDataFrame(data.frame(
		seqnames = GenomeIndex$V1, start = 1, end = GenomeIndex$V2))

	polyA_blocks = list()

	message("Processing positive strands\n")
	pb = txtProgressBar(max = length(Genome.gr), style = 3)
	for(i in 1:length(Genome.gr)) {
		setTxtProgressBar(pb,i)

		seq = getSeq(genome, Genome.gr[i])

		matches = vmatchPattern(paste(rep("A",PolyN), collapse = ""), seq, max.mismatch=PolyMM)
		if(elementNROWS(matches) > 0) {
			matches.gr = makeGRangesFromDataFrame(data.frame(
				seqnames = as.character(Genome.gr[i]@seqnames),
				start = unlist(startIndex(matches)), end = unlist(endIndex(matches)), 
				strand = "+"))
			matches.gr = reduce(matches.gr)
			
			polyA_blocks[[length(polyA_blocks) + 1]] = as.data.frame(matches.gr)
		}
	}
	setTxtProgressBar(pb,i)
	close(pb)

	message("Processing negative strands\n")
	pb = txtProgressBar(max = length(Genome.gr), style = 3)
	for(i in 1:length(Genome.gr)) {
		setTxtProgressBar(pb,i)

		seq = getSeq(genome, Genome.gr[i])

		matches = vmatchPattern(paste(rep("T",PolyN), collapse = ""), seq, max.mismatch=PolyMM)
		if(elementNROWS(matches) > 0) {
			matches.gr = makeGRangesFromDataFrame(data.frame(
				seqnames = as.character(Genome.gr[i]@seqnames),
				start = unlist(startIndex(matches)), end = unlist(endIndex(matches)), 
				strand = "-"))
			matches.gr = reduce(matches.gr)			
			polyA_blocks[[length(polyA_blocks) + 1]] = as.data.frame(matches.gr)
		}
	}
	setTxtProgressBar(pb,i)
	close(pb)
	rm(seq)
	gc()
	polyA_blocks = rbindlist(polyA_blocks)

	Genes.Gr = annogenes
	PolyA.Gr = makeGRangesFromDataFrame(polyA_blocks)
	
	Gene_PolyA = suppressWarnings(findOverlaps(PolyA.Gr, Genes.Gr))

	Gene_PolyA.DT = data.table(block_id = Gene_PolyA@from, gene_index = Gene_PolyA@to)
	Gene_PolyA.DT[, gene_id := annogenes$gene_id[gene_index]]
	Gene_PolyA.DT[, strand := as.vector(strand(annogenes))[gene_index]]
	Gene_PolyA.DT[, start := polyA_blocks$start[block_id]]
	Gene_PolyA.DT[, end := polyA_blocks$end[block_id]]
	
	Gene_PolyA.summary = Gene_PolyA.DT[, .(min(start), max(end)), 
		by = c("gene_id", "strand")]

	annogenes$Last_PolyA = 0
	annogenes$Last_PolyA[match(Gene_PolyA.summary$gene_id, annogenes$gene_id)] = ifelse(
		Gene_PolyA.summary$strand == "+", Gene_PolyA.summary$V2, Gene_PolyA.summary$V1)

	annointrons[strand == "+", Post_PolyA := (start > annogenes$Last_PolyA[match(gene_id, annogenes$gene_id)])]
	annointrons[strand == "-", Post_PolyA := (end < annogenes$Last_PolyA[match(gene_id, annogenes$gene_id)] |
		annogenes$Last_PolyA[match(gene_id, annogenes$gene_id)] == 0)]
	
	annointrons_dedup = annointrons[!duplicated(annointrons[,
				c("seqnames","intron_start","intron_end","strand")])]	


	# Interrogate splice site motifs
	message("Collating splice site motifs...", appendLF = F)

	introns_ss_left.gr = with(annointrons_dedup, makeGRangesFromDataFrame(data.frame(
		seqnames = seqnames, start = intron_start + 1,
			end = intron_start + 2,strand = strand)))
	introns_ss_right.gr = with(annointrons_dedup, makeGRangesFromDataFrame(data.frame(
		seqnames = seqnames, start = intron_end - 2,
			end = intron_end - 1, strand = strand)))
	intron.seq.left = DNAStringSet(getSeq(genome,introns_ss_left.gr))
	intron.seq.right = DNAStringSet(getSeq(genome,introns_ss_right.gr))
	
	annointrons_dedup$SS_motif = with(annointrons_dedup, ifelse(strand == "+",
		paste0(intron.seq.left, intron.seq.right),
		paste0(intron.seq.right, intron.seq.left)))

	annointrons[annointrons_dedup, on = c("seqnames", "intron_start", "intron_end", "strand"),
		SS_motif := i.SS_motif]

	saveRDS(annointrons, paste(RefPath, "annointrons.Rds", sep="/"))
	saveRDS(annointrons_dedup, paste(RefPath, "annointrons_dedup.Rds", sep="/"))

	message("\u2713\n")

	message("Intron Annotations finished\n")
	
	irf <- NxtIRF.Benchmark(irf)
	irf <- NxtIRF.StartBenchmarkTimer(irf)

# Annotating Alternative Splicing Events


	message("Annotating Splice Events\n")

	annointrons = readRDS(paste(RefPath, "annointrons.Rds", sep="/"))

	introns.skipcoord = annointrons
	setorder(introns.skipcoord, gene_id, transcript_name, intron_number)
	introns.skipcoord[strand == "+", skip_coord := ifelse(intron_number == max(intron_number), NA, 
		paste0(seqnames,":",intron_start,
		"-", data.table::shift(intron_end, 1, NA, "lead"),"/",strand)), by = transcript_id]
	introns.skipcoord[strand == "-", skip_coord := ifelse(intron_number == max(intron_number), NA, 
		paste0(seqnames,":",data.table::shift(intron_start, 1, NA, "lead"),
		"-", intron_end,"/",strand)), by = transcript_id]
	introns.skipcoord[, skip_coord_2 := data.table::shift(skip_coord, 1, NA, "lag")]

	introns.skippedJn = introns.skipcoord[skip_coord %in% Event, c("gene_id","gene_name","skip_coord")]
	introns.skippedJn = unique(introns.skippedJn)
	
message("Annotating Mutually-Exclusive-Exon AS Events...", appendLF = F)

	introns.search.MXE = introns.skipcoord[introns.skipcoord[,
		.I[intron_number < max(intron_number)], by = transcript_id]$V1]
	introns.search.MXE = introns.search.MXE[introns.search.MXE[,.N,by = c("gene_id", "skip_coord")],
		on = c("gene_id", "skip_coord"), N := i.N]
	introns.search.MXE = introns.search.MXE[N > 1]
	introns.search.MXE.pos = introns.search.MXE[strand == "+"]
	setorder(introns.search.MXE.pos, seqnames, intron_start, intron_end)
	introns.search.MXE.neg = introns.search.MXE[strand == "-"]
	setorder(introns.search.MXE.neg, seqnames, intron_end, -intron_start)
	introns.search.MXE = rbindlist(list(introns.search.MXE.pos, introns.search.MXE.neg))
	
	introns.search.MXE = introns.search.MXE[, c("skip_coord", "gene_id", "Event", "transcript_id",
		"transcript_name", "intron_number")]
	setnames(introns.search.MXE, old = "Event", new = "Event1")
	
	introns.search.MXE2 = introns.skipcoord[,c("skip_coord_2","gene_id","Event","transcript_id","transcript_name")]
	setnames(introns.search.MXE2, old = c("skip_coord_2","Event"), new = c("skip_coord", "Event2"))
	
	introns.search.MXE[introns.search.MXE2, on =  c("gene_id","transcript_id","transcript_name","skip_coord"),
		Event2 := i.Event2]
	
	introns.search.MXE = unique(introns.search.MXE, by= c("gene_id","skip_coord","Event1"))
	introns.search.MXE = unique(introns.search.MXE, by= c("gene_id","skip_coord","Event2"))
	introns.search.MXE = introns.search.MXE[, if(.N>1) .SD, by = c("gene_id","skip_coord")]
	
	introns.found.MXE = introns.search.MXE[ , {
		edge1 = rep(1:.N, (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event1[edge1], Event1b = Event1[edge2],
			Event2a = Event2[edge1], Event2b = Event2[edge2],
			transcript_id_a = transcript_id[edge1], transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1], transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1], intron_number_b = intron_number[edge2]
		)
		}, by = skip_coord]

	introns.found.MXE[, gene_id := factor(gene_id,annogenes$gene_id,ordered=TRUE)]
	setorder(introns.found.MXE, gene_id, transcript_name_a)
	introns.found.MXE[, EventName := paste0("MXE:", transcript_name_a, "-exon",(1 + intron_number_a),";",
			transcript_name_b,"-exon",(1 + intron_number_b))]
	introns.found.MXE[, EventID := paste0("MXE#", 1:.N)]
	setnames(introns.found.MXE, old = "skip_coord", new = "EventRegion")
	introns.found.MXE[, EventType := "MXE"]
	introns.found.MXE = introns.found.MXE[, c("EventType","EventID","EventName","Event1a","Event1b","Event2a","Event2b",
		"gene_id","gene_id_b","EventRegion",
		"transcript_id_a","transcript_name_a","intron_number_a",
		"transcript_id_b","transcript_name_b","intron_number_b")]

	message("\u2713\n")

message("Annotating Skipped-Exon AS Events...", appendLF = F)

# annotate skipped junctions with two included junctions

	introns.found.SE = introns.skippedJn[,"skip_coord"]
	introns.search.SE = annointrons[, c("gene_id","Event","transcript_id", 
				"transcript_name","intron_number")]
	setnames(introns.search.SE, old = c("Event", "transcript_id", "transcript_name", "intron_number"),
		new = c("skip_coord", "skip_transcript_id", "skip_transcript_name", "skip_intron_number"))
	introns.found.SE[introns.search.SE, on = "skip_coord", 
		c("gene_id_b", "skip_transcript_id", "skip_transcript_name", "skip_intron_number") :=
		list(i.gene_id, i.skip_transcript_id, i.skip_transcript_name, i.skip_intron_number)]
	introns.found.SE = unique(introns.found.SE, by = c("gene_id_b","skip_coord"))	

	introns.search.SE2 = introns.skipcoord[,c("skip_coord","gene_id","Event","transcript_id","transcript_name",
				"intron_number")]
	introns.found.SE[introns.search.SE2, on = "skip_coord",
			c("gene_id", "inc_coord_upst", "inc_transcript_id", "inc_transcript_name", "inc_intron_number") :=
			list(i.gene_id, i.Event, i.transcript_id, i.transcript_name, i.intron_number)]

	introns.search.SE3 = introns.skipcoord[,c("skip_coord_2","gene_id","Event","transcript_id","transcript_name")]
	setnames(introns.search.SE3, old = c("skip_coord_2", "transcript_id"), new = c("skip_coord", "inc_transcript_id"))
	introns.found.SE[introns.search.SE3, on = c("skip_coord", "inc_transcript_id"),
			inc_coord_downst := i.Event]

	introns.found.SE = unique(introns.found.SE, by = c("gene_id","skip_coord","inc_transcript_id", "inc_transcript_name"))
	introns.found.SE[, gene_id := factor(gene_id, annogenes$gene_id, ordered = TRUE)]
	setorder(introns.found.SE, gene_id, inc_transcript_name)
	introns.found.SE[, EventName := paste0(
			"SE:",inc_transcript_name,"-exon",(1 + inc_intron_number),";",
			skip_transcript_name,"-int",skip_intron_number
			)]
	introns.found.SE[, EventID := paste0("SE#",seq_len(.N))]
	introns.found.SE[, EventType := "SE"]
	introns.found.SE[, Event2b := NA]
	introns.found.SE[, EventRegion := skip_coord]
	introns.found.SE = introns.found.SE[, c("EventType","EventID","EventName","inc_coord_upst",
			"skip_coord","inc_coord_downst","Event2b",
			"gene_id","gene_id_b", "EventRegion",
			"inc_transcript_id","inc_transcript_name","inc_intron_number",
			"skip_transcript_id","skip_transcript_name","skip_intron_number")]
	setnames(introns.found.SE, old = c("inc_coord_upst", "inc_coord_downst", "skip_coord"),
		new = c("Event1a", "Event2a", "Event1b"))
	setnames(introns.found.SE, old = c("inc_transcript_id", "inc_transcript_name", "skip_transcript_id","skip_transcript_name"),
		new = c("transcript_id_a", "transcript_name_a", "transcript_id_b","transcript_name_b"))
	setnames(introns.found.SE, new = c("intron_number_a","intron_number_b"), old = c("inc_intron_number", "skip_intron_number"))

	message("\u2713\n")

message("Annotating Alternate First / Last Exon AS Events...", appendLF = F)

	# AFE/ALE

	introns.search.AFE = annointrons[intron_number == 1]
	introns.search.AFE[, seqnames := factor(seqnames, chrOrder, ordered = TRUE)]

	introns.search.AFE.pos = introns.search.AFE[strand == "+"]
	setorder(introns.search.AFE.pos, seqnames, intron_end, -intron_start)
	introns.search.AFE.pos = introns.search.AFE.pos[,c("seqnames", "intron_end", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number")]
	# introns.search.AFE.pos = introns.search.AFE.pos[introns.search.AFE.pos[, .I[.N > 1], by = c("seqnames", "intron_end")]$V1]
	setnames(introns.search.AFE.pos, old = "intron_end", new = "intron_coord")

	introns.search.AFE.neg = introns.search.AFE[strand == "-"]
	setorder(introns.search.AFE.neg, seqnames, intron_start, intron_end)
	introns.search.AFE.neg = introns.search.AFE.neg[,c("seqnames", "intron_start", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number")]
	# introns.search.AFE.neg = introns.search.AFE.neg[introns.search.AFE.neg[, .I[.N > 1], by = c("seqnames", "intron_start")]$V1]
	setnames(introns.search.AFE.neg, old = "intron_start", new = "intron_coord")

	introns.search.AFE = rbindlist(list(introns.search.AFE.pos, introns.search.AFE.neg))
	introns.search.AFE = unique(introns.search.AFE, by = "Event")

	introns.search.ALE = annointrons[annointrons[,
			.I[intron_number == max(intron_number)], by = transcript_id]$V1]
	introns.search.ALE[, seqnames := factor(seqnames, chrOrder, ordered = TRUE)]

	introns.search.ALE.pos = introns.search.ALE[strand == "+"]
	setorder(introns.search.ALE.pos, seqnames, intron_start, intron_end)
	introns.search.ALE.pos = introns.search.ALE.pos[,c("seqnames", "intron_start", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number")]
	# introns.search.ALE.pos = introns.search.ALE.pos[introns.search.ALE.pos[, .I[.N > 1], by = c("seqnames", "intron_start")]$V1]
	setnames(introns.search.ALE.pos, old = "intron_start", new = "intron_coord")

	introns.search.ALE.neg = introns.search.ALE[strand == "-"]
	setorder(introns.search.ALE.neg, seqnames, intron_end, -intron_start)
	introns.search.ALE.neg = introns.search.ALE.neg[,c("seqnames", "intron_end", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number")]
	# introns.search.ALE.neg = introns.search.ALE.neg[introns.search.ALE.neg[, .I[.N > 1], by = c("seqnames", "intron_end")]$V1]
	setnames(introns.search.ALE.neg, old = "intron_end", new = "intron_coord")

	introns.search.ALE = rbindlist(list(introns.search.ALE.pos, introns.search.ALE.neg))
	introns.search.ALE = unique(introns.search.ALE, by = "Event")

	introns.found.AFE = introns.search.AFE[ , {
		edge1 = rep(1:.N, (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1], transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1], transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1], intron_number_b = intron_number[edge2]
		)
	}, by = c("seqnames", "intron_coord")]
	introns.found.AFE = introns.found.AFE[!is.na(gene_id)]
	
	introns.found.ALE = introns.search.ALE[ , {
		edge1 = rep(1:.N, (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1], transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1], transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1], intron_number_b = intron_number[edge2]
		)
	}, by = c("seqnames", "intron_coord")]
	introns.found.ALE = introns.found.ALE[!is.na(gene_id)]

	introns.found.AFE = unique(introns.found.AFE, by = c("Event1a", "Event1b"))
	introns.found.AFE = introns.found.AFE[, gene_id := factor(gene_id, annogenes$gene_id, ordered = TRUE)]
	setorder(introns.found.AFE, gene_id, transcript_name_a)
	introns.found.AFE = introns.found.AFE[, EventName := paste0("AFE:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.AFE = introns.found.AFE[, EventType := "AFE"]
	introns.found.AFE = introns.found.AFE[, EventRegion := Event1b]
	introns.found.AFE[, EventID := paste0("AFE#", 1:.N)]
	introns.found.AFE = introns.found.AFE[, c("EventType", "EventID", "EventName",
		"Event1a", "Event1b", "Event2a", "Event2b",
		"gene_id", "gene_id_b", "EventRegion",
		"transcript_id_a", "transcript_name_a", "intron_number_a",
		"transcript_id_b", "transcript_name_b", "intron_number_b")]
	
	introns.found.ALE = unique(introns.found.ALE, by = c("Event1a", "Event1b"))
	introns.found.ALE = introns.found.ALE[, gene_id := factor(gene_id, annogenes$gene_id, ordered = TRUE)]
	setorder(introns.found.ALE, gene_id, transcript_name_a)
	introns.found.ALE = introns.found.ALE[, EventName := paste0("ALE:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.ALE = introns.found.ALE[, EventType := "ALE"]
	introns.found.ALE = introns.found.ALE[, EventRegion := Event1b]
	introns.found.ALE[, EventID := paste0("ALE#", 1:.N)]
	introns.found.ALE = introns.found.ALE[, c("EventType", "EventID", "EventName",
		"Event1a", "Event1b", "Event2a", "Event2b",
		"gene_id", "gene_id_b", "EventRegion",
		"transcript_id_a", "transcript_name_a", "intron_number_a",
		"transcript_id_b", "transcript_name_b", "intron_number_b")]

	message("\u2713\n")

	message("Annotating Alternate 5' / 3' Splice Site AS Events...", appendLF = F)

	introns.search.A5SS = annointrons
	introns.search.A5SS[, seqnames := factor(seqnames, chrOrder, ordered = TRUE)]
	# introns.search.A5SS = unique(introns.search.A5SS, by = c("gene_id","Event"))

	introns.search.A5SS.pos = introns.search.A5SS[strand == "+"]
	setorder(introns.search.A5SS.pos, seqnames, intron_end, -intron_start)
	introns.search.A5SS.pos = introns.search.A5SS.pos[,c("seqnames", "intron_end", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number", "exon_groups_start","exon_groups_end")]
	# introns.search.A5SS.pos = introns.search.A5SS.pos[introns.search.A5SS.pos[, .I[.N > 1], 
		# by = c("intron_end", "exon_groups_start","exon_groups_end")]$V1]
	setnames(introns.search.A5SS.pos, old = "intron_end", new = "intron_coord")

	introns.search.A5SS.neg = introns.search.A5SS[strand == "-"]
	setorder(introns.search.A5SS.neg, seqnames, intron_start, intron_end)
	introns.search.A5SS.neg = introns.search.A5SS.neg[,c("seqnames", "intron_start", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number","exon_groups_start","exon_groups_end")]
	# introns.search.A5SS.neg = introns.search.A5SS.neg[introns.search.A5SS.neg[, .I[.N > 1], 
		# by = c("intron_start","exon_groups_start","exon_groups_end")]$V1]
	setnames(introns.search.A5SS.neg, old = "intron_start", new = "intron_coord")

	introns.search.A5SS = rbindlist(list(introns.search.A5SS.pos, introns.search.A5SS.neg))
	introns.search.A5SS = unique(introns.search.A5SS, by = "Event")

	introns.search.A3SS = annointrons
	introns.search.A3SS[, seqnames := factor(seqnames, chrOrder, ordered = TRUE)]
	# introns.search.A3SS = unique(introns.search.A3SS, by = c("gene_id","Event"))

	introns.search.A3SS.pos = introns.search.A3SS[strand == "+"]
	setorder(introns.search.A3SS.pos, seqnames, intron_start, intron_end)
	introns.search.A3SS.pos = introns.search.A3SS.pos[,c("seqnames", "intron_start", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number","exon_groups_start","exon_groups_end")]
	# introns.search.A3SS.pos = introns.search.A3SS.pos[introns.search.A3SS.pos[, .I[.N > 1], 
		# by = c("intron_start","exon_groups_start","exon_groups_end")]$V1]
	setnames(introns.search.A3SS.pos, old = "intron_start", new = "intron_coord")

	introns.search.A3SS.neg = introns.search.A3SS[strand == "-"]
	setorder(introns.search.A3SS.neg, seqnames, intron_end, -intron_start)
	introns.search.A3SS.neg = introns.search.A3SS.neg[,c("seqnames", "intron_end", "strand", "Event", "gene_id", 
		"transcript_id", "transcript_name", "intron_number","exon_groups_start","exon_groups_end")]
	# introns.search.A3SS.neg = introns.search.A3SS.neg[introns.search.A3SS.neg[, .I[.N > 1], 
		# by = c("intron_end", "exon_groups_start","exon_groups_end")]$V1]
	setnames(introns.search.A3SS.neg, old = "intron_end", new = "intron_coord")

	introns.search.A3SS = rbindlist(list(introns.search.A3SS.pos, introns.search.A3SS.neg))
	introns.search.A3SS = unique(introns.search.A3SS, by = "Event")

	introns.found.A5SS = introns.search.A5SS[ , {
		edge1 = rep(1:.N, (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1], transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1], transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1], intron_number_b = intron_number[edge2]
		)
		}, by = intron_coord]
	introns.found.A5SS = introns.found.A5SS[!is.na(gene_id)]
	
	introns.found.A3SS = introns.search.A3SS[ , {
		edge1 = rep(1:.N, (.N:1) - 1L)
		i = 2L:(.N * (.N - 1L) / 2L + 1L)
			o = cumsum(c(0, (.N-2L):1))
		edge2 = i - o[edge1]
		.(
			gene_id = gene_id[edge1], gene_id_b = gene_id[edge2],
			Event1a = Event[edge1], Event1b = Event[edge2],
			Event2a = NA, Event2b = NA, EventRegion = Event[edge2],
			transcript_id_a = transcript_id[edge1], transcript_id_b = transcript_id[edge2],
			transcript_name_a = transcript_name[edge1], transcript_name_b = transcript_name[edge2],
			intron_number_a = intron_number[edge1], intron_number_b = intron_number[edge2]
		)
		}, by = intron_coord]
	introns.found.A3SS = introns.found.A3SS[!is.na(gene_id)]

	introns.found.A5SS = unique(introns.found.A5SS, by = c("Event1a", "Event1b"))
	introns.found.A5SS = introns.found.A5SS[, gene_id := factor(gene_id, annogenes$gene_id, ordered = TRUE)]
	setorder(introns.found.A5SS, gene_id, transcript_name_a)
	introns.found.A5SS = introns.found.A5SS[, EventName := paste0("A5SS:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.A5SS = introns.found.A5SS[, EventType := "A5SS"]
	introns.found.A5SS = introns.found.A5SS[, EventRegion := Event1b]
	introns.found.A5SS = introns.found.A5SS[!introns.found.AFE, on = c("Event1a", "Event1b")]
	introns.found.A5SS[, EventID := paste0("A5SS#", 1:.N)]
	introns.found.A5SS = introns.found.A5SS[, c("EventType", "EventID", "EventName",
		"Event1a", "Event1b", "Event2a", "Event2b",
		"gene_id", "gene_id_b", "EventRegion",
		"transcript_id_a", "transcript_name_a", "intron_number_a",
		"transcript_id_b", "transcript_name_b", "intron_number_b")]

	introns.found.A3SS = unique(introns.found.A3SS, by = c("Event1a", "Event1b"))
	introns.found.A3SS = introns.found.A3SS[, gene_id := factor(gene_id, annogenes$gene_id, ordered = TRUE)]
	setorder(introns.found.A3SS, gene_id, transcript_name_a)
	introns.found.A3SS = introns.found.A3SS[, EventName := paste0("A3SS:",transcript_name_a,"-exon",intron_number_a,
			";",transcript_name_b,"-exon",intron_number_b)]
	introns.found.A3SS = introns.found.A3SS[, EventType := "A3SS"]
	introns.found.A3SS = introns.found.A3SS[, EventRegion := Event1b]
	introns.found.A3SS = introns.found.A3SS[!introns.found.ALE, on = c("Event1a", "Event1b")]
	introns.found.A3SS[, EventID := paste0("A3SS#", 1:.N)]
	introns.found.A3SS = introns.found.A3SS[, c("EventType", "EventID", "EventName",
		"Event1a", "Event1b", "Event2a", "Event2b",
		"gene_id", "gene_id_b", "EventRegion",
		"transcript_id_a", "transcript_name_a", "intron_number_a",
		"transcript_id_b", "transcript_name_b", "intron_number_b")]

	AS_Table = rbindlist(list(introns.found.MXE, introns.found.SE,
		introns.found.AFE, introns.found.ALE,
		introns.found.A5SS, introns.found.A3SS))
	
	# Annotate 5' and 3' splice events for single cell
	
	AS_Table[, Splice_5prime_A := Event1a]
	AS_Table[, Splice_3prime_A := Event1a]
	AS_Table[EventType %in% c("SE", "MXE"), Splice_3prime_A := Event2a]
	AS_Table[, Splice_5prime_B := Event1b]
	AS_Table[, Splice_3prime_B := Event1b]
	AS_Table[EventType == "MXE", Splice_3prime_B := Event2b]
	
	saveRDS(AS_Table, paste(RefPath, "annosplice.Rds", sep="/"))
	
	# fwrite(AS_Table,paste(RefPath,"Splice.Annotation.csv",sep="/"))

	message("\u2713\n")

	message("Splice Annotations finished. ", appendLF = F)
	
	irf <- NxtIRF.Benchmark(irf)
	
	# Merge all partial references into single reference object
	
	Annotation = new("NxtAnnotation")
	
	Annotation@chrOrder <- chrOrder
	
	Annotation@Genes = readRDS(paste(RefPath,"annogenes.Rds",sep="/"))
	Annotation@Transcripts = readRDS(paste(RefPath,"annotrans.Rds",sep="/"))
	Annotation@Exons = readRDS(paste(RefPath,"annoexons.Rds",sep="/"))
	
	Annotation@Introns = readRDS(paste(RefPath,"annointrons.Rds",sep="/"))
	Annotation@Splice = readRDS(paste(RefPath,"annosplice.Rds",sep="/"))
	
	# deduplicate introns table
	Annotation@Introns_dedup = readRDS(paste(RefPath,"annointrons_dedup.Rds",sep="/"))

	file.remove(paste(RefPath, c("annogenes.Rds", "annotrans.Rds","annoexons.Rds",
		"annointrons.Rds", "annosplice.Rds", "annointrons_dedup.Rds"), sep="/"))

	if(UseAH) {
        Annotation@RefPath = "AnnotationHub"
    } else {
        Annotation@RefPath = normalizePath(RefPath)
    }
	Annotation@RefFasta = fasta
    Annotation@RefGTF = gtf
    
	message("Reference build finished", appendLF = F)
	irf <- NxtIRF.Benchmark(irf, "BuildReference")
	return(Annotation)
}

NxtIRF.FinaliseReference <- function(RefPath = "Reference") {

# Orientation
	assertthat::assert_that(dir.exists(RefPath),
		msg = "Reference Directory not found")

	assertthat::assert_that(all(c("genome.fa.fai", "annogenes.Rds", "annotrans.Rds","annoexons.Rds",
		"annointrons.Rds", "annosplice.Rds", "annointrons_dedup.Rds") %in% list.files(RefPath)),
		msg = "Reference in current directory not built. Use NxtIRF.BuildReference(RefPath) prior to using NxtIRF")
		
	Annotation = new("NxtAnnotation")
	
	GenomeIndex = fread(paste(Annotation@RefPath,"genome.fa.fai",sep="/"))
	Annotation@chrOrder <- GenomeIndex$V1
	
	Annotation@Genes = readRDS(paste(Annotation@RefPath,"annogenes.Rds",sep="/"))
	Annotation@Transcripts = readRDS(paste(Annotation@RefPath,"annotrans.Rds",sep="/"))
	Annotation@Exons = readRDS(paste(Annotation@RefPath,"annoexons.Rds",sep="/"))
	
	Annotation@Introns = readRDS(paste(Annotation@RefPath,"annointrons.Rds",sep="/"))
	Annotation@Splice = readRDS(paste(Annotation@RefPath,"annosplice.Rds",sep="/"))
	
	# deduplicate introns table
	Annotation@Introns_dedup = readRDS(paste(Annotation@RefPath,"annointrons_dedup.Rds",sep="/"))

	# saveRDS(Annotation, paste(Annotation@RefPath, Ref_Filename, sep="/"))
	# message(paste("Annotation saved as", paste(Annotation@RefPath, Ref_Filename, sep="/")))
	
	# Remove intermediate files
	file.remove(paste(Annotation@RefPath, c("annogenes.Rds", "annotrans.Rds","annoexons.Rds",
		"annointrons.Rds", "annosplice.Rds", "annointrons_dedup.Rds"), sep="/"))

	return(Annotation)
}

FindMinInFrame <- function(matchList) {
	if(any(matchList %% 3 == 1)) {
		return(min(matchList[matchList %% 3 == 1]))
	}  else {return(NA)}
}

CreateGroups.GR <- function(gr) {
	gr.group = reduce(gr, with.revmap = TRUE)
	gr$vector = 0
	gr$vector[unlist(gr.group$revmap)] = rep(seq_len(length(gr.group$revmap)),
		lengths(gr.group$revmap))
	return(gr$vector)
}

CreateGroups.GRL <- function(gr, split.vector) {
	gr$order = seq_len(length(gr))
	grl = split(gr, split.vector)
	grl.group = reduce(grl, with.revmap = TRUE)
	gr.group = unlist(grl.group)
	add_to_revmap = unname(rep(c(0, cumsum(lengths(grl))), c(lengths(grl), 0)))
	gr.DT = as.data.table(grl)
	gr.DT[, vector := 0]
	gr.DT$vector[unlist(gr.group$revmap) + add_to_revmap] = rep(
		seq_len(length(gr.group$revmap)),
		lengths(gr.group$revmap))
	gr.DT[, vector := vector - min(vector) + 1, by = "gene_id"]
	gr.DT[, vector_final := ifelse((strand == "-"), max(vector) - vector + 1, vector), by = "gene_id"]
	return_val = gr.DT$vector_final[match(seq_len(length(gr)), gr.DT$order)]
}

grlGaps<-function(grl) {
	psetdiff(unlist(range(grl),use.names=TRUE),grl)
}


CreateGroups <- function(FeatureList) {
	FeatureList = FeatureList %>% group_by(seqnames) %>%
		mutate(Split.Vector = ifelse(
			start > cummax(ifelse(is.na(lag(end)),0,lag(end))),
			1,0)) %>% 
		ungroup() %>%
		mutate(Group = cumsum(Split.Vector)) 
	return(FeatureList$Group)
}

CreateExonGroups <- function(ExonTable, chrOrder) {
	
	returnval = c()
	
	for(i in seq_len(length(chrOrder))) {
		temp = ExonTable %>% filter(seqnames == chrOrder[i]) %>%
			group_by(gene_id) %>%
			mutate(Split.Vector = ifelse(
				start > cummax(ifelse(is.na(lag(end)),0,lag(end))),1,0)) %>%
			mutate(ExonGroup = cumsum(Split.Vector)) %>%
			mutate(ExonGroup = ifelse(strand == "-",max(ExonGroup) - ExonGroup + 1,ExonGroup))
		returnval = c(returnval,temp$ExonGroup)
	}
	return(returnval)
		
}
