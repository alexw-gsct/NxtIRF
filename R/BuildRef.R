#' Fetch genome / transcriptome reference from AnnotationHub and writes to reference
#'
#' @export
FetchAnnotation <- function(ah_genome = "AH65745", ah_transcriptome = "AH64631", reference_path = "./Reference", 
    verbose = FALSE) {
    NxtIRF.CheckPackageInstalled("AnnotationHub", "2.20.0")
    if (verbose) message("Initializing AnnotationHub")
    
    ah = AnnotationHub::AnnotationHub()

    genome.fasta.tmp = paste0(normalizePath(reference_path), "/", ah_genome, ".genome.fa.tmp")
    genome.fasta = paste0(normalizePath(reference_path), "/", ah_genome, ".genome.fa")

    if (file.exists(genome.fasta)) {
        message("Genome file already exists, skipping...")
    } else {
        if (verbose) message("Fetching genome from AnnotationHub")
        
        genome = ah[[ah_genome]]
        chrOrder = names(seqinfo(genome))
        gr.df = data.frame(seqnames = chrOrder,
            start = 1, end = seqinfo(genome)[chrOrder]@seqlengths)
        genome.seq = getSeq(genome, GenomicRanges::makeGRangesFromDataFrame(gr.df))
        names(genome.seq) = chrOrder
        rtracklayer::export(genome.seq, genome.fasta, "fasta")
        # Rsamtools::indexFa(genome.fasta)
        # for(chr_i in seq_len(length(chrOrder))) {
            # chr = chrOrder[chr_i]
            # maxlen = seqinfo(genome)[chrOrder[chr_i]]@seqlengths
            # gr = GRanges(seqnames = chr, ranges = IRanges(1, maxlen))
            # seq.df = data.frame(
                # name = paste0(">",chr), seq = Biostrings::DNAStringSet(getSeq(genome, gr)))
            # data.table::fwrite(seq.df, file = genome.fasta.tmp, append = file.exists(genome.fasta.tmp),
                # sep= "\n", eol = "\n", col.names = FALSE)
        # }
        # message(paste("Cleaning and processing ",genome.fasta))
        # run_IRFinder_CleanFasta(genome.fasta.tmp, genome.fasta)
        # file.remove(genome.fasta.tmp)
        message(paste("Genome fasta file written to ",genome.fasta))
    }
    
    transcripts.gtf = paste0(normalizePath(reference_path), "/", ah_transcriptome, ".transcripts.gtf")
    if(file.exists(transcripts.gtf)) {
        message("Transcriptome file already exists, skipping...")    
    } else {

        if (verbose) message("Fetching transcript gtf file from AnnotationHub")
        gtf = ah[[ah_transcriptome]]
        genome(gtf) = genome(gtf)[1]

        rtracklayer::export(gtf, transcripts.gtf, "gtf")

        message(paste("Transcriptome gtf file written to ", transcripts.gtf))
    }
}

#' Builds IRFinder Reference
#'
#' Builds reference files used by IRFinder / rIRFinder. Requires genome fasta and transcriptome gtf file
#'
#' @param fasta: The path to the genome fasta file
#' @param gtf: The path to the transcript gtf file
#' @param reference_path: The directory to write the (r)/IRFinder reference files
#' @return Nothing. The created reference will be written to the given directory
#' @export
BuildReference <- function(fasta = "genome.fa", gtf = "transcripts.gtf", fasta_ah = "", gtf_ah = "",
    reference_path = "./Reference",
    genome_type = "", nonPolyARef = "", MappabilityRef = "", BlacklistRef = "",
    FilterIRByProcessedTranscript = FALSE) {

    # genome_type = match.arg(genome_type)
    # if(genome_type != "") message(paste(genome_type, "specified as genome type. Using corresponding nonPolyA reference"))
    
    assertthat::assert_that(genome_type != "",
            msg = "genome_type not specified. This should be either one of 'hg38', 'hg19', 'mm10', 'mm9', or 'other'. If 'other', please provide a nonPolyARef file or leave blank to omit polyA profiling.")
    
    nonPolyAFile = ""
    if(genome_type == "hg38") {
        nonPolyAFile = system.file("extra-input-files/Human_hg38_nonPolyA_ROI.bed", package = "rIRFinder")
    } else if(genome_type == "hg19")  {
        nonPolyAFile = system.file("extra-input-files/Human_hg19_nonPolyA_ROI.bed", package = "rIRFinder")
    } else if(genome_type == "mm10")  {
        nonPolyAFile = system.file("extra-input-files/Mouse_mm10_nonPolyA_ROI.bed", package = "rIRFinder")
    } else if(genome_type == "mm9")  {
        nonPolyAFile = system.file("extra-input-files/Mouse_mm9_nonPolyA_ROI.bed", package = "rIRFinder")
    } else {
        nonPolyAFile = nonPolyARef
    }
    if(nonPolyAFile != "" & !file.exists(nonPolyAFile)) {
        message(paste(nonPolyARef, "not found. Reference generated without non-polyA reference"))
        nonPolyAFile = ""
    }
    
    MappabilityFile = ""
    if (MappabilityRef == "") {
        message("Mappability table not provided. IRFinder reference will be generated without mappability exclusion")        
    } else if(MappabilityRef != "" & !file.exists(MappabilityRef)) {
        message(paste(MappabilityRef, "not found. Reference generated without mappability exclusion"))
    } else {
        MappabilityFile = MappabilityRef
    }

    BlacklistFile = ""
    if (BlacklistRef == "") {
        message("Blacklist table not provided. IRFinder reference will be generated without Blacklist exclusion")        
    } else if(BlacklistRef != "" & !file.exists(BlacklistRef)) {
        message(paste(BlacklistRef, "not found. Reference generated without Blacklist exclusion"))
    } else {
        BlacklistFile = BlacklistRef
    }

    if(fasta_ah != "") {
        assertthat::assert_that(substr(fasta_ah,1,2) == "AH",
            msg = "Given genome AnnotationHub reference is incorrect")
        message("Loading AnnotationHub")            
        ah = AnnotationHub::AnnotationHub()
        message("Connecting to AnnotationHub genome...", appendLF = F)
        genome = ah[[fasta_ah]]
        genome_ah = TRUE
        message("done\n")
    } else {
        assertthat::assert_that(file.exists(normalizePath(fasta)),
            msg = paste("Given genome fasta file", normalizePath(fasta), "not found"))
        message("Connecting to genome file...", appendLF = F)
        genome = Biostrings::readDNAStringSet(fasta)
        genome_ah = FALSE
        message("done\n")
    }
    if(gtf_ah != "") {
        assertthat::assert_that(substr(gtf_ah,1,2) == "AH",
            msg = "Given transcriptome AnnotationHub reference is incorrect")
            if(!exists("ah")) {
                message("Loading AnnotationHub")            
                ah = AnnotationHub::AnnotationHub()
            }
            message("Reading AnnotationHub GTF file...", appendLF = F)
            gtf.gr = ah[[gtf_ah]]
            message("done\n")
    } else {
        assertthat::assert_that(file.exists(normalizePath(gtf)),
            msg = paste("Given transcriptome gtf file", normalizePath(gtf), "not found"))
        message("Reading source GTF file...", appendLF = F)
        gtf.gr = rtracklayer::import(gtf)
        message("done\n")
    }
 
    if(genome_ah == TRUE) {
        chrOrder = names(rtracklayer::seqinfo(genome))    
    } else {
        chrOrder = names(BSgenome::seqinfo(genome))
    }

    message("Processing gtf file...", appendLF = F)
   
    # fix gene / transcript names with '/' (which breaks IRFinder code)
    gtf.gr$gene_name = gsub("/","_",gtf.gr$gene_name)
    gtf.gr$transcript_name = gsub("/","_",gtf.gr$transcript_name)

    # Extracting and saving Genes, Transcripts, Exons, Proteins and saving as .fst files for faster random access
    Genes = gtf.gr[gtf.gr$type == "gene"]
    Genes <- GenomeInfoDb::sortSeqlevels(Genes)
    Genes <- sort(Genes)
        fst::write.fst(as.data.frame(Genes), paste(reference_path,"Genes.fst", sep="/"))

    Transcripts = gtf.gr[gtf.gr$type == "transcript"]
    Transcripts <- GenomeInfoDb::sortSeqlevels(Transcripts)
    Transcripts <- sort(Transcripts)
        fst::write.fst(as.data.frame(Transcripts), paste(reference_path,"Transcripts.fst", sep="/"))

    Exons = gtf.gr[gtf.gr$type == "exon"]
    Exons <- GenomeInfoDb::sortSeqlevels(Exons)
    Exons <- sort(Exons)
        fst::write.fst(as.data.frame(Exons), paste(reference_path,"Exons.fst", sep="/"))

    Proteins = gtf.gr[gtf.gr$type == "CDS"]
    Proteins <- GenomeInfoDb::sortSeqlevels(Proteins)
    Proteins <- sort(Proteins)
        fst::write.fst(as.data.frame(Proteins), paste(reference_path,"Proteins.fst", sep="/"))
    
    gtf.misc = gtf.gr[!gtf.gr$type %in% c("gene", "transcript", "exon", "CDS")]
    gtf.misc <- GenomeInfoDb::sortSeqlevels(gtf.misc)
    gtf.misc <- sort(gtf.misc)
        fst::write.fst(as.data.frame(gtf.misc), paste(reference_path,"Misc.fst", sep="/"))
    
    # Generating IRFinder-base references
    Genes.rev = Genes
    GenomicRanges::strand(Genes.rev) = ifelse(GenomicRanges::strand(Genes.rev) == "+", "-", 
        ifelse(GenomicRanges::strand(Genes.rev) == "-", "+", "*")) # Invert strand
    Genes.Extended = GenomicRanges::reduce(c(GenomicRanges::flank(Genes.rev, 5000), 
        GenomicRanges::flank(Genes.rev, 1000, start = F)))

#   Why filter by protein_coding or processed_transcript? Make this an option
    if(FilterIRByProcessedTranscript == TRUE) {
        if("transcript_biotype" %in% names(GenomicRanges::mcols(Exons))) {
            candidate.transcripts = Exons[Exons$transcript_biotype %in% c("processed_transcript", "protein_coding")]
        } else if("transcript_type" %in% names(GenomicRanges::mcols(Exons))) {
            candidate.transcripts = Exons[Exons$transcript_type %in% c("processed_transcript", "protein_coding")]    
        } else {
            candidate.transcripts = Exons
        }
    } else {
        candidate.transcripts = Exons    
    }

    if("transcript_biotype" %in% names(GenomicRanges::mcols(Exons))) {
        tmp.exons.exclude =  Exons[!grepl("intron", Exons$transcript_biotype)]
    } else if("transcript_type" %in% names(GenomicRanges::mcols(Exons))) {
        tmp.exons.exclude =  Exons[!grepl("intron", Exons$transcript_type)]
    } else {
        tmp.exons.exclude =  Exons
    }
    message("done\n")

    message("Processing introns...", appendLF = F)    

    candidate.introns = grlGaps(
        GenomicRanges::split(candidate.transcripts, candidate.transcripts$transcript_id)
    )
    candidate.introns = data.table::as.data.table(candidate.introns)    
    candidate.introns$group = NULL
    colnames(candidate.introns)[1] = "transcript_id"
    data.table::setorder(candidate.introns, seqnames, start, end, strand)
    
    # Grab splice motifs at this point; filter by valid splice motifs
    donor.introns = data.frame(seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+", candidate.introns$start, candidate.introns$end - 1),
        stop = ifelse(candidate.introns$strand == "+", candidate.introns$start + 1, candidate.introns$end),
        strand = candidate.introns$strand)
    donor.seq = BSgenome::getSeq(genome, GenomicRanges::makeGRangesFromDataFrame(donor.introns))
    acceptor.introns = data.frame(seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+", candidate.introns$end - 1, candidate.introns$start),
        stop = ifelse(candidate.introns$strand == "+", candidate.introns$end, candidate.introns$start + 1),
        strand = candidate.introns$strand)
    acceptor.seq = BSgenome::getSeq(genome, GenomicRanges::makeGRangesFromDataFrame(acceptor.introns))
    candidate.introns$splice_motif = paste0(donor.seq, acceptor.seq)

    # Acceptable splice motifs: GT-AG, GC-AG, AT-AC, AT-AG, GT-AC
    # novel.introns = candidate.introns[!(splice_motif %in% c("GTAG", "GCAG", "ATAC", "ATAG", "GTAC"))]
    candidate.introns = candidate.introns[splice_motif %in% c("GTAG", "GCAG", "ATAC", "ATAG", "GTAC")]

# Annotating Introns:
    candidate.introns[,intron_number := data.table::rowid(transcript_id)]
    candidate.introns[strand == "-", intron_number := max(intron_number) + 1 - intron_number, by = "transcript_id"]
    candidate.introns[,intron_id := paste0(transcript_id, "/Intron", intron_number)]
    candidate.introns[data.table::as.data.table(Transcripts), on = "transcript_id", 
        c("gene_name", "gene_id", "transcript_name") := list(i.gene_name, i.gene_id, i.transcript_name)]
    
# Do other annotations here:
    candidate.introns[data.table::as.data.table(Transcripts), on = "transcript_id", 
        c("gene_name", "gene_id", "transcript_name") := list(i.gene_name, i.gene_id, i.transcript_name)]
    if("transcript_support_level" %in% names(GenomicRanges::mcols(Transcripts))) {
        candidate.introns[data.table::as.data.table(Transcripts), on = "transcript_id", 
        c("transcript_support_level") := list(i.transcript_support_level)]
        candidate.introns[, transcript_support_level := 
            data.table::tstrsplit(transcript_support_level, split=" ")[[1]]]
        candidate.introns[is.na(transcript_support_level), transcript_support_level := "NA"]
    }
    if("protein_id" %in% names(GenomicRanges::mcols(Proteins))) {
        candidate.introns[data.table::as.data.table(Proteins), on = "transcript_id", 
        c("protein_id") := list(i.protein_id)]
    }
    if("ccds_id" %in% names(GenomicRanges::mcols(Exons))) {
        candidate.introns[data.table::as.data.table(Exons), on = "transcript_id", 
        c("ccds_id") := list(i.ccds_id)]
    }
        fst::write.fst(candidate.introns, paste(reference_path,"candidate.introns.fst", sep="/"))

    message("done\n")

    message("Generating ref-cover.bed ...", appendLF = F)    

# Finished annotating introns, now use it to build reference:
    # Sort by tsl first, then reverse later
    
    if("transcript_support_level" %in% colnames(candidate.introns)) setorder(candidate.introns, transcript_support_level)    
        introns.unique = unique(candidate.introns, by = c("seqnames", "start", "end", "width", "strand"))
        setorder(introns.unique, seqnames, start, end, strand)
        introns.unique = GenomicRanges::makeGRangesFromDataFrame(as.data.frame(introns.unique), keep.extra.columns=TRUE)
    if("transcript_support_level" %in% colnames(candidate.introns)) data.table::setorder(candidate.introns, seqnames, start, end, strand)
    
    exclude.directional = data.table::as.data.table(tmp.exons.exclude)
    exclude.directional = unique(exclude.directional, by = c("seqnames", "start", "end", "width", "strand"))
    exclude.directional[, start := start - 5]
    exclude.directional[, end := end + 5]
    
    exclude.directional.reverse = exclude.directional[strand == "-", strand:= "P"]
    exclude.directional.reverse[strand == "+", strand:= "-"]
    exclude.directional.reverse[strand == "P", strand:= "+"]

    # mappability = data.table::fread("c:/alex/IRFinder_data/REF/mappa.csv")
    # TODO: check mappability and blacklist are valid Bed3 formats (chr, start, end)
    if(MappabilityFile != "") {
        mappability = data.table::fread(MappabilityFile)
        mappability = mappability[,1:3]
        colnames(mappability) = c("seqnames", "start", "end")
        exclude.omnidirectional = GenomicRanges::makeGRangesFromDataFrame(mappability) # + merge with any blacklists
        exclude.omnidirectional = GenomicRanges::reduce(exclude.omnidirectional, min.gapwidth = 9) # merge with any gaps <= 9
        if(BlacklistFile != "") {
            blacklist = data.table::fread(BlacklistFile)
            blacklist = blacklist[,1:3]
            colnames(blacklist) = c("seqnames", "start", "end")
            exclude.omnidirectional = c(exclude.omnidirectional,
                GenomicRanges::makeGRangesFromDataFrame(blacklist))
        }
    } else if(BlacklistFile != "") {
        blacklist = data.table::fread(BlacklistFile)
        blacklist = blacklist[,1:3]
        colnames(blacklist) = c("seqnames", "start", "end")
        exclude.omnidirectional = GenomicRanges::makeGRangesFromDataFrame(blacklist)
    }

    introns.unique.blacklisted = GenomicRanges::findOverlaps(introns.unique, exclude.omnidirectional, type = "within")
    # length(introns.unique.blacklisted@from)
    # [1] 3071
    introns.unique = introns.unique[-introns.unique.blacklisted@from] # clean introns by those lying completely within blacklist regions

    # Now label introns as "known-exon", "anti-over", or "anti-near"
    introns.unique.exon.dir = GenomicRanges::findOverlaps(introns.unique, 
        GenomicRanges::makeGRangesFromDataFrame(exclude.directional), type = "within")
    introns.unique.exon.nd = GenomicRanges::findOverlaps(introns.unique, 
        GenomicRanges::makeGRangesFromDataFrame(exclude.directional), type = "within", ignore.strand=TRUE)

    introns.unique$known_exon_dir = ( seq(length(introns.unique)) %in% introns.unique.exon.dir@from )
    introns.unique$known_exon_nd = ( seq(length(introns.unique)) %in% introns.unique.exon.nd@from )

    introns.unique.antiover = GenomicRanges::findOverlaps(introns.unique, Genes.rev)
    introns.unique.antinear = GenomicRanges::findOverlaps(introns.unique, Genes.rev)

    introns.unique$antiover = ( seq(length(introns.unique)) %in% introns.unique.antiover@from )
    introns.unique$antinear = ( seq(length(introns.unique)) %in% introns.unique.antinear@from )

# Now subset introns by punching holes using blacklist regions

    introns.unique$intron_start = BiocGenerics::start( introns.unique ) - 1
    introns.unique$intron_end = BiocGenerics::end( introns.unique )
    introns.unique$intron_width = BiocGenerics::width(introns.unique)

# Conditions for final intron inclusion
  # if ($newlen > 40 && ($newlen/$len) >= 0.7) {
# Thus, Remove introns less than 50 bp:
    introns.unique = introns.unique[BiocGenerics::width(introns.unique) > 50]

# remove 5 bases from start & end
    BiocGenerics::start(introns.unique) = BiocGenerics::start(introns.unique) + 5
    BiocGenerics::end(introns.unique) = BiocGenerics::end(introns.unique) - 5

    introns.unique.dir = introns.unique
    introns.unique.nd = introns.unique

# Dir
    introns.intersect.dir = GenomicRanges::intersect(introns.unique.dir, c(exclude.omnidirectional, 
        GenomicRanges::makeGRangesFromDataFrame(exclude.directional)))

    introns.intersect.ol = GenomicRanges::findOverlaps(introns.unique, introns.intersect.dir)
    # make a GRanges same size as the number of intersections
    introns.intersect.dir.final = introns.intersect.dir[introns.intersect.ol@to]
    introns.intersect.dir.final$intron_id = introns.unique$intron_id[introns.intersect.ol@from]

    introns.unique.dir.ID = GenomicRanges::split(introns.unique.dir, introns.unique.dir$intron_id)
    introns.intersect.dir.ID = GenomicRanges::split(introns.intersect.dir.final, introns.intersect.dir.final$intron_id)
    introns.unique.dir.ID.compare = introns.unique.dir.ID[names(introns.unique.dir.ID) %in% names(introns.intersect.dir.ID)]

# nd
    introns.intersect.nd = GenomicRanges::intersect(introns.unique.nd, c(exclude.omnidirectional, 
        GenomicRanges::makeGRangesFromDataFrame(exclude.directional), 
        GenomicRanges::makeGRangesFromDataFrame(exclude.directional.reverse)))

    introns.intersect.ol = GenomicRanges::findOverlaps(introns.unique, introns.intersect.nd)
    # make a GRanges same size as the number of intersections
    introns.intersect.nd.final = introns.intersect.nd[introns.intersect.ol@to]
    introns.intersect.nd.final$intron_id = introns.unique$intron_id[introns.intersect.ol@from]

    introns.unique.nd.ID = GenomicRanges::split(introns.unique.nd, introns.unique.nd$intron_id)
    introns.intersect.nd.ID = GenomicRanges::split(introns.intersect.nd.final, introns.intersect.nd.final$intron_id)
    introns.unique.nd.ID.compare = introns.unique.nd.ID[names(introns.unique.nd.ID) %in% names(introns.intersect.nd.ID)]

# Dir setdiff
    tmpdir.IntronCover = GenomicRanges::setdiff(introns.unique.dir.ID.compare, introns.intersect.dir.ID)
    # now add back introns that did not require intersection
    tmpdir.IntronCover = c(tmpdir.IntronCover, 
        introns.unique.dir.ID[!(names(introns.unique.dir.ID) %in% names(introns.intersect.dir.ID))])

    tmpdir.IntronCover = data.table::as.data.table(tmpdir.IntronCover)
    tmpdir.IntronCover = tmpdir.IntronCover[,c("seqnames", "start", "end", "strand", "width", "group_name")]
    colnames(tmpdir.IntronCover)[6] = "intron_id"

    tmpdir.IntronCover.summa = tmpdir.IntronCover
    tmpdir.IntronCover.summa[, c("num_blocks", "inclbases") := list(.N, sum(width)), by = "intron_id"]
    tmpdir.IntronCover.summa = unique(tmpdir.IntronCover.summa[,c("intron_id", "num_blocks", "inclbases")]
        , by = "intron_id")
    tmpdir.IntronCover.summa[data.table::as.data.table(introns.unique), 
        on = "intron_id", c("seqnames", "intron_start", "intron_end", "intron_width", "width", "strand", "gene_name", "transcript_id", "known_exon_dir")
          := list(i.seqnames, i.intron_start, i.intron_end, i.intron_width, i.width, i.strand, i.gene_name, i.transcript_id, i.known_exon_dir)]
    tmpdir.IntronCover.summa[, exclbases := width - inclbases]
        # Exclude exclbases / width > 0.3
    tmpdir.IntronCover.summa = tmpdir.IntronCover.summa[exclbases / intron_width < 0.3]

    tmpdir.IntronCover = semi_join.DT(tmpdir.IntronCover, tmpdir.IntronCover.summa, by = "intron_id")
 
    tmpdir.IntronCover.summa[, IRFname := paste("dir", gene_name, transcript_id, strand, num_blocks, 
        intron_start, intron_end, inclbases, exclbases,
        ifelse(known_exon_dir, "known-exon","clean"), sep="/")]

    tmpdir.IntronCover = GenomicRanges::makeGRangesFromDataFrame(tmpdir.IntronCover, keep.extra.columns=TRUE)
    tmpdir.IntronCover = GenomicRanges::split(tmpdir.IntronCover, tmpdir.IntronCover$intron_id)

    names(tmpdir.IntronCover) = tmpdir.IntronCover.summa$IRFname[match(
        names(tmpdir.IntronCover), tmpdir.IntronCover.summa$intron_id)]

# Nondir setdiff
    tmpnd.IntronCover = GenomicRanges::setdiff(introns.unique.nd.ID.compare, introns.intersect.nd.ID)
    # now add back introns that did not require intersection
    tmpnd.IntronCover = c(tmpnd.IntronCover, 
        introns.unique.nd.ID[!(names(introns.unique.nd.ID) %in% names(introns.intersect.nd.ID))])

    tmpnd.IntronCover = data.table::as.data.table(tmpnd.IntronCover)
    tmpnd.IntronCover = tmpnd.IntronCover[,c("seqnames", "start", "end", "strand", "width", "group_name")]
    colnames(tmpnd.IntronCover)[6] = "intron_id"

    tmpnd.IntronCover.summa = tmpnd.IntronCover
    tmpnd.IntronCover.summa[, c("num_blocks", "inclbases") := list(.N, sum(width)), by = "intron_id"]
    tmpnd.IntronCover.summa = unique(tmpnd.IntronCover.summa[,c("intron_id", "num_blocks", "inclbases")]
        , by = "intron_id")
    tmpnd.IntronCover.summa[data.table::as.data.table(introns.unique), 
        on = "intron_id", c("seqnames","intron_start", "intron_end", "intron_width", "width", "strand", "gene_name", "transcript_id", 
        "known_exon_nd", "antiover", "antinear")
          := list(i.seqnames, i.intron_start, i.intron_end, i.intron_width, i.width, i.strand, i.gene_name, i.transcript_id, 
          i.known_exon_nd, i.antiover, i.antinear)]
    tmpnd.IntronCover.summa[, exclbases := width - inclbases]
        # Exclude exclbases / width > 0.3
    tmpnd.IntronCover.summa = tmpnd.IntronCover.summa[exclbases / intron_width < 0.3]

    tmpnd.IntronCover = semi_join.DT(tmpnd.IntronCover, tmpnd.IntronCover.summa, by = "intron_id")
 
    tmpnd.IntronCover.summa[, IRFname := paste("nd", gene_name, transcript_id, strand, num_blocks, 
        intron_start, intron_end, inclbases, exclbases, sep="/")]
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

    tmpnd.IntronCover = GenomicRanges::makeGRangesFromDataFrame(tmpnd.IntronCover, keep.extra.columns=TRUE)
    tmpnd.IntronCover = GenomicRanges::split(tmpnd.IntronCover, tmpnd.IntronCover$intron_id)

    names(tmpnd.IntronCover) = tmpnd.IntronCover.summa$IRFname[match(
        names(tmpnd.IntronCover), tmpnd.IntronCover.summa$intron_id)]

# Sort and Export out
	data.table::setorder(tmpnd.IntronCover.summa, seqnames, intron_start, intron_end, strand)
    tmpnd.IntronCover = tmpnd.IntronCover[tmpnd.IntronCover.summa$IRFname]
    
	data.table::setorder(tmpdir.IntronCover.summa, seqnames, intron_start, intron_end, strand)
    tmpdir.IntronCover = tmpdir.IntronCover[tmpdir.IntronCover.summa$IRFname]

    rtracklayer::export(tmpdir.IntronCover, paste(reference_path, "tmpdir.IntronCover.bed", sep="/"))
    rtracklayer::export(tmpnd.IntronCover, paste(reference_path, "tmpnd.IntronCover.bed", sep="/"))

# Generate final ref-cover.bed

    tmpdir.IntronCover = data.table::fread(paste(reference_path, "tmpdir.IntronCover.bed", sep="/"), sep="\t")
    tmpdir.IntronCover[,cat := "dir"]
    tmpnd.IntronCover = data.table::fread(paste(reference_path, "tmpnd.IntronCover.bed", sep="/"), sep="\t")
    tmpnd.IntronCover[,cat := "nd"]

    ref.cover = rbind(tmpdir.IntronCover, tmpnd.IntronCover)
    data.table::setorder(ref.cover, V1, V2, V3, V6, cat)
    ref.cover$cat = NULL
    ref.cover[, V9 := as.character(V9)]
    ref.cover[, V9 := "255,0,0"]

    data.table::fwrite(ref.cover, paste(reference_path, "ref-cover.bed", sep="/"), sep="\t", col.names = F)

    message("done\n")

    message("Generating ref-ROI.bed ...", appendLF = F)    
 
# ROI
    if("gene_biotype" %in% names(GenomicRanges::mcols(Transcripts))) {
        rRNA = as.data.frame(Transcripts[grepl("rRNA", Transcripts$gene_biotype)])
        rRNA$start = rRNA$start - 1
        rRNA$name = with(rRNA, paste("rRNA", seqnames, start, end, strand,
            transcript_id, gene_biotype, gene_id, gene_name, sep="/"))
        rRNA = rRNA[, c("seqnames", "start", "end", "name")]
    } else if("gene_type" %in% names(GenomicRanges::mcols(Transcripts))) {
        rRNA = as.data.frame(Transcripts[grepl("rRNA", Transcripts$gene_type)])
        rRNA$start = rRNA$start - 1
        rRNA$name = with(rRNA, paste("rRNA", seqnames, start, end, strand,
            transcript_id, gene_type, gene_id, gene_name, sep="/"))
        rRNA = rRNA[, c("seqnames", "start", "end", "name")]    
    } else {
        rRNA = c()
    }
    
    if(nonPolyAFile != "") {
        nonPolyA = data.table::fread(nonPolyAFile, sep="\t")
        nonPolyA = nonPolyA[,1:3]
        colnames(nonPolyA) = c("seqnames", "start", "end")
        nonPolyA$name = with(nonPolyA,paste("NonPolyA", seqnames, start, end, sep="/"))
    } else {
        nonPolyA = c()
    }
    
    AllChr = data.frame(seqnames = names(BSgenome::seqinfo(genome)),
        start = 1, end = BSgenome::seqinfo(genome)@seqlengths, names = names(BSgenome::seqinfo(genome)))
    AllChr = GenomicRanges::makeGRangesListFromDataFrame(AllChr, split.field = "names")
    Genes.chr = c(Genes, GenomicRanges::flank(Genes, 10000), GenomicRanges::flank(Genes, 10000, start = F))
    Genes.chr = GenomicRanges::reduce(Genes.chr, min.gapwidth = 1000)
    Genes.chr$chr = GenomicRanges::seqnames(Genes.chr)
    Genes.chr = GenomicRanges::split(Genes.chr, Genes.chr$chr)
    
    AllChr = AllChr[names(Genes.chr)]
    AllChr.split = GenomicRanges::setdiff(AllChr, Genes.chr, ignore.strand = TRUE)
    Intergenic = unlist(AllChr.split)
    names(Intergenic) = 1:length(Intergenic)
    Intergenic = as.data.frame(Intergenic)
    Intergenic = Intergenic[,c("seqnames", "start", "end")]
    Intergenic$name = paste("Intergenic", Intergenic$seqnames, sep="/")
    Intergenic$start = Intergenic$start - 1

    ref.ROI = rbind(rRNA, nonPolyA, Intergenic) %>% dplyr::arrange(seqnames, start)
    
    data.table::fwrite(ref.ROI, paste(reference_path, "ref-ROI.bed", sep="/"), sep="\t", col.names = F)

    message("done\n")

    message("Generating ref-read-continues.ref ...", appendLF = F)    

# ref-read-continues.ref
    introns.unique.readcons = rbind(tmpdir.IntronCover.summa[, c("seqnames", "intron_start", "intron_end", "strand")],
        tmpnd.IntronCover.summa[, c("seqnames", "intron_start", "intron_end", "strand")])
    readcons.left = introns.unique.readcons[,c("seqnames", "intron_start", "strand")]
    readcons.right = introns.unique.readcons[,c("seqnames", "intron_end", "strand")]
    colnames(readcons.left) = c("V1", "V2", "V3")
    colnames(readcons.right) = c("V1", "V2", "V3")
    readcons = rbind(readcons.left, readcons.right) %>% dplyr::arrange(V1, V2, V3) %>% 
        dplyr::filter(!duplicated(.))
    
    data.table::fwrite(readcons, paste(reference_path, "ref-read-continues.ref", sep="/"), sep="\t", col.names = F)

    message("done\n")

    message("Generating ref-sj.ref ...", appendLF = F)    
    
# ref-sj.ref
    ref.sj = candidate.introns[,c("seqnames", "start", "end", "strand")]
    ref.sj = unique(ref.sj)
    ref.sj[,start := start - 1]
    data.table::fwrite(ref.sj, paste(reference_path, "ref-sj.ref", sep="/"), sep="\t", col.names = F)

    message("done\n")
    
# Concatenate all 4 reference files into one file
    data.table::fwrite(list(">ref-cover.bed"), paste(reference_path, "IRFinder.ref", sep="/"), 
        sep="\t", eol = "\n", col.names = F)
    data.table::fwrite(ref.cover, paste(reference_path, "IRFinder.ref", sep="/"), append = TRUE, 
        sep="\t", eol = "\n", col.names = F)
    data.table::fwrite(list(">ref-read-continues.ref"), paste(reference_path, "IRFinder.ref", sep="/"), append = TRUE, 
        sep="\t", eol = "\n", col.names = F)
    data.table::fwrite(readcons, paste(reference_path, "IRFinder.ref", sep="/"), append = TRUE, 
        sep="\t", eol = "\n", col.names = F)
    data.table::fwrite(list(">ref-ROI.bed"), paste(reference_path, "IRFinder.ref", sep="/"), append = TRUE, 
        sep="\t", eol = "\n", col.names = F)
    data.table::fwrite(ref.ROI, paste(reference_path, "IRFinder.ref", sep="/"), append = TRUE, 
        sep="\t", eol = "\n", col.names = F)
    data.table::fwrite(list(">ref-sj.ref"), paste(reference_path, "IRFinder.ref", sep="/"), append = TRUE, 
        sep="\t", eol = "\n", col.names = F)
    data.table::fwrite(ref.sj, paste(reference_path, "IRFinder.ref", sep="/"), append = TRUE, 
        sep="\t", eol = "\n", col.names = F)
    
# Annotating Alternative Splicing Events

    message("Annotating Splice Events\n")

    candidate.introns[,intron_start := start]
    candidate.introns[,intron_end := end]
    candidate.introns[,Event :=  paste0(seqnames, ":", intron_start, "-", intron_end, "/", strand)]
    
    GeneOrder = data.table::as.data.table(Genes)
    setorder(GeneOrder, seqnames, start, end, strand)
    
	introns.skipcoord = copy(candidate.introns)
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
	
message("Annotating Mutually-Exclusive-Exon Splice Events...", appendLF = F)

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

	introns.found.MXE[, gene_id := factor(gene_id,GeneOrder$gene_id,ordered=TRUE)]
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

	message("done\n")

message("Annotating Skipped-Exon Splice Events...", appendLF = F)

# annotate skipped junctions with two included junctions

	introns.found.SE = introns.skippedJn[,"skip_coord"]
	introns.search.SE = candidate.introns[, c("gene_id","Event","transcript_id", 
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
	introns.found.SE[, gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
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

	message("done\n")

message("Annotating Alternate First / Last Exon Splice Events...", appendLF = F)

	# AFE/ALE

	introns.search.AFE = candidate.introns[intron_number == 1]
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

	introns.search.ALE = candidate.introns[candidate.introns[,
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
	introns.found.AFE = introns.found.AFE[, gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
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
	introns.found.ALE = introns.found.ALE[, gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
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

	message("done\n")

	message("Annotating Alternate 5' / 3' Splice Site Splice Events...", appendLF = F)

    # This section requires "exon groups"
    #   First we group genes into strand-specific contiguous groups of genes
    #   Then all transcripts are combined into list separated by gene groups
    #   Then each GRangesList is collapsed to derive exon groups
    #   All exons are labelled by exon groups
    #   A5SS and A3SS require both splice alternatives occur in exons belonging to same exon group
    
    Genes.Group = GenomicRanges::reduce(Genes)
    Genes.Group$gene_group = 1:length(Genes.Group)
    Exons.Group = GenomicRanges::reduce(candidate.transcripts)
    Exons.Group.GG.ol = GenomicRanges::findOverlaps(
        GenomicRanges::makeGRangesFromDataFrame(Exons.Group), Genes.Group)
    Exons.Group$gene_group[Exons.Group.GG.ol@from] = Genes.Group$gene_group[Exons.Group.GG.ol@to]
    Exons.Group = data.table::as.data.table(Exons.Group)
    setorder(Exons.Group, seqnames, start, end)
    Exons.Group[, exon_group := data.table::rowid(gene_group)]
    Exons.Group[strand == "-", exon_group := max(exon_group) + 1 - exon_group, by = "gene_group"]

    candidate.introns.ASS = copy(candidate.introns)
    candidate.introns.ASS[, c("start", "end") := list(start - 1, end + 1)]
    candidate.introns.ASS.ol = GenomicRanges::findOverlaps(
        GenomicRanges::makeGRangesFromDataFrame(candidate.introns.ASS),
        GenomicRanges::makeGRangesFromDataFrame(Exons.Group))
    candidate.introns.ASS.summa = data.table::data.table(
        intron_id = candidate.introns.ASS$intron_id[candidate.introns.ASS.ol@from],
        exon_group = Exons.Group$exon_group[candidate.introns.ASS.ol@to])
    candidate.introns.ASS.summa[, c("exon_groups_start", "exon_groups_end") := 
        list(min(exon_group), max(exon_group)), by = "intron_id"]
    candidate.introns.ASS.summa = unique(candidate.introns.ASS.summa, by = "intron_id")
    candidate.introns.ASS[candidate.introns.ASS.summa, on = "intron_id",
        c("exon_groups_start", "exon_groups_end") := list(i.exon_groups_start, i.exon_groups_end)]
    
    
	introns.search.A5SS = copy(candidate.introns.ASS)
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

	introns.search.A3SS = copy(candidate.introns.ASS)
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
			intron_number_a = intron_number[edge1], intron_number_b = intron_number[edge2],
            exon_groups_start_a = exon_groups_start[edge1], exon_groups_start_b = exon_groups_start[edge2], 
            exon_groups_end_a = exon_groups_end[edge1], exon_groups_end_b = exon_groups_end[edge2]
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
			intron_number_a = intron_number[edge1], intron_number_b = intron_number[edge2],
            exon_groups_start_a = exon_groups_start[edge1], exon_groups_start_b = exon_groups_start[edge2], 
            exon_groups_end_a = exon_groups_end[edge1], exon_groups_end_b = exon_groups_end[edge2]
		)
		}, by = intron_coord]
	introns.found.A3SS = introns.found.A3SS[!is.na(gene_id)]

	introns.found.A5SS = unique(introns.found.A5SS, by = c("Event1a", "Event1b"))
    # filter by same exon group starts and ends:
   	introns.found.A5SS = introns.found.A5SS[exon_groups_start_a == exon_groups_start_b]
   	introns.found.A5SS = introns.found.A5SS[exon_groups_end_a == exon_groups_end_b]

	introns.found.A5SS = introns.found.A5SS[, gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
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
    # filter by same exon group starts and ends:
   	introns.found.A3SS = introns.found.A3SS[exon_groups_start_a == exon_groups_start_b]
   	introns.found.A3SS = introns.found.A3SS[exon_groups_end_a == exon_groups_end_b]
	introns.found.A3SS = introns.found.A3SS[, gene_id := factor(gene_id, GeneOrder$gene_id, ordered = TRUE)]
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
	
    # Rename based on tsl and protein-coding ability, if applicable
    if("transcript_support_level" %in% colnames(candidate.introns)) {
        candidate.introns.order = copy(candidate.introns)
        candidate.introns.order[, is_protein_coding := !is.na(protein_id)]
        candidate.introns.order[, by = "transcript_id", is_last_intron := (intron_number == max(intron_number))]
        
        AS_Table.search.a = AS_Table[, c("EventType", "EventID", "Event1a", "Event2a")]
        AS_Table.search.a[,Event := Event1a]
        AS_Table.search.a = candidate.introns.order[AS_Table.search.a, on = "Event", 
            c("EventType","EventID", "Event1a", "Event2a", "transcript_id", "transcript_support_level", "is_protein_coding", "is_last_intron", "intron_number")]
        setnames(AS_Table.search.a, "intron_number", "in_1a")
        AS_Table.search.a = AS_Table.search.a[EventType !=  "AFE" | in_1a == 1]
        AS_Table.search.a = AS_Table.search.a[EventType !=  "ALE" | is_last_intron]
        AS_Table.search.a[,Event := Event2a]
        AS_Table.search.a[is.na(Event),Event := Event1a]
        AS_Table.search.a = candidate.introns.order[AS_Table.search.a, 
            on = c("Event",  "transcript_id", "transcript_support_level"),
            c("EventType","EventID", "Event1a", "Event2a", "transcript_id", "transcript_support_level", "is_protein_coding", "is_last_intron","in_1a", "intron_number")]
        AS_Table.search.a = AS_Table.search.a[!is.na(intron_number)]
        setnames(AS_Table.search.a, "intron_number", "in_2a")
        
        AS_Table.search.b = AS_Table[, c("EventType", "EventID", "Event1b", "Event2b")]
        AS_Table.search.b[,Event := Event1b]
        AS_Table.search.b = candidate.introns.order[AS_Table.search.b, on = "Event", 
            c("EventType","EventID", "Event1b", "Event2b", "transcript_id", "transcript_support_level", "is_protein_coding", "is_last_intron","intron_number")]
        setnames(AS_Table.search.b, "intron_number", "in_1b")
        AS_Table.search.b = AS_Table.search.b[EventType !=  "AFE" | in_1b == 1]
        AS_Table.search.b = AS_Table.search.b[EventType !=  "ALE" | is_last_intron]
        AS_Table.search.b[,Event := Event2b]
        AS_Table.search.b[is.na(Event),Event := Event1b]
        AS_Table.search.b = candidate.introns.order[AS_Table.search.b, 
            on = c("Event",  "transcript_id", "transcript_support_level"),
            c("EventType","EventID", "Event1b", "Event2b", "transcript_id", "transcript_support_level", "is_protein_coding", "is_last_intron", "in_1b", "intron_number")]
        AS_Table.search.b = AS_Table.search.b[!is.na(intron_number)]
        setnames(AS_Table.search.b, "intron_number", "in_2b")

        AS_Table.search.a[candidate.introns.order, on = "transcript_id", transcript_name := i.transcript_name]
        AS_Table.search.b[candidate.introns.order, on = "transcript_id", transcript_name := i.transcript_name]
        setorder(AS_Table.search.a, transcript_support_level, -is_protein_coding, transcript_name)
        setorder(AS_Table.search.b, transcript_support_level, -is_protein_coding, transcript_name)
        AS_Table.search.a = unique(AS_Table.search.a, by = "EventID")
        AS_Table.search.a = AS_Table.search.a[AS_Table[, "EventID"], on = "EventID"]
        AS_Table.search.b = unique(AS_Table.search.b, by = "EventID")
        AS_Table.search.b = AS_Table.search.b[AS_Table[, "EventID"], on = "EventID"]
        
        AS_Table$transcript_id_a = AS_Table.search.a$transcript_id
        AS_Table$transcript_name_a = AS_Table.search.a$transcript_name
        AS_Table$intron_number_a = AS_Table.search.a$in_1a
        AS_Table$transcript_id_b = AS_Table.search.b$transcript_id
        AS_Table$transcript_name_b = AS_Table.search.b$transcript_name
        AS_Table$intron_number_b = AS_Table.search.b$in_1b
        
        AS_Table[EventType == "MXE", EventName := paste0("MXE:", transcript_name_a,"-exon",
            as.character(as.numeric(intron_number_a) + 1), ";", transcript_name_b,"-exon",
            as.character(as.numeric(intron_number_b) + 1))]
        AS_Table[EventType == "SE", EventName := paste0("SE:", transcript_name_a,"-exon",
            as.character(as.numeric(intron_number_a) + 1), ";", transcript_name_b,"-int",
            as.character(as.numeric(intron_number_b)))]
        AS_Table[EventType == "AFE", EventName := paste0("AFE:", transcript_name_a,"-exon1;", 
            transcript_name_b,"-exon1")]
        AS_Table[EventType == "ALE", EventName := paste0("ALE:", transcript_name_a, "-exon", 
            as.character(as.numeric(intron_number_a) + 1), ";", transcript_name_b, "-exon",
            as.character(as.numeric(intron_number_b) + 1))]
        AS_Table[EventType == "A5SS", EventName := paste0("A5SS:", transcript_name_a,"-exon", 
            as.character(as.numeric(intron_number_a)), ";", transcript_name_b,"-exon",
            as.character(as.numeric(intron_number_b)))]
        AS_Table[EventType == "A3SS", EventName := paste0("A3SS:", transcript_name_a,"-exon", 
            as.character(as.numeric(intron_number_a + 1)), ";", transcript_name_b,"-exon",
            as.character(as.numeric(intron_number_b + 1)))]
    }
	
        fst::write.fst(as.data.frame(AS_Table), paste(reference_path,"Splice.fst", sep="/"))

	message("done\n")

	message("Translating Alternate Splice Peptides", appendLF = F)

    # Proteomic Consequences of Splicing
    AS_Table.Extended = copy(AS_Table)
    Proteins.Splice = as.data.table(Proteins)
    Proteins.Splice$exon_number = as.numeric(Proteins.Splice$exon_number)
    Proteins.Splice[, phase := -phase %% 3]    # make phase easier for me to understand
    # Upstream applicable for MXE, SE, ALE, A3SS
    Upstream = AS_Table[EventType %in% c("MXE", "SE", "ALE", "A3SS")]

    # Do A
    Upstream.A = Upstream[, c("EventID", "transcript_id_a", "intron_number_a")]
    Upstream.A[, c("transcript_id", "exon_number") := list(transcript_id_a, intron_number_a)]
    # left_join with Exons
    Upstream.A = Proteins.Splice[Upstream.A, on = c("transcript_id", "exon_number"), c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Upstream.A.gr = makeGRangesFromDataFrame(as.data.frame(na.omit(Upstream.A)), keep.extra.columns = T)
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
    AS_Table.Extended[EventType %in% c("MXE", "SE", "ALE", "A3SS"), AA_upstr.A := Upstream.A$AA_seq]

    # repeat for B:
    Upstream.B = Upstream[, c("EventID", "transcript_id_b", "intron_number_b")]
    Upstream.B[, c("transcript_id", "exon_number") := list(transcript_id_b, intron_number_b)]
    # left_join with Exons
    Upstream.B = Proteins.Splice[Upstream.B, on = c("transcript_id", "exon_number"), c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Upstream.B.gr = makeGRangesFromDataFrame(as.data.frame(na.omit(Upstream.B)), keep.extra.columns = T)
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
    AS_Table.Extended[EventType %in% c("MXE", "SE", "ALE", "A3SS"), AA_upstr.B := Upstream.B$AA_seq]
    
    # Do downstream seq before casette:
    Downstream = AS_Table[EventType %in% c("MXE", "SE", "AFE", "A5SS")]
    # Add EventType as exon_number is conditional on this
    Downstream.A = Downstream[, c("EventType", "EventID", "transcript_id_a", "intron_number_a")]
    Downstream.A[, c("transcript_id", "exon_number") := list(transcript_id_a, intron_number_a)]
    # Modify downstream exon number
    Downstream.A[EventType %in% c("MXE", "SE"), exon_number := exon_number + 2]
    Downstream.A[EventType %in% c("AFE", "A5SS"), exon_number := exon_number + 1]
    # left_join with Exons
    Downstream.A = Proteins.Splice[Downstream.A, on = c("transcript_id", "exon_number"), c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Downstream.A.gr = makeGRangesFromDataFrame(as.data.frame(na.omit(Downstream.A)), keep.extra.columns = T)
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
    AS_Table.Extended[EventType %in% c("MXE", "SE", "AFE", "A5SS"), AA_downstr.A := Downstream.A$AA_seq]
    # B:
    Downstream.B = Downstream[, c("EventType", "EventID", "transcript_id_b", "intron_number_b")]
    Downstream.B[, c("transcript_id", "exon_number") := list(transcript_id_b, intron_number_b)]
    # Modify downstream exon number: Note SE is different for B
    Downstream.B[EventType %in% c("MXE"), exon_number := exon_number + 2]
    Downstream.B[EventType %in% c("SE", "AFE", "A5SS"), exon_number := exon_number + 1]
    # left_join with Exons
    Downstream.B = Proteins.Splice[Downstream.B, on = c("transcript_id", "exon_number"), c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Downstream.B.gr = makeGRangesFromDataFrame(as.data.frame(na.omit(Downstream.B)), keep.extra.columns = T)
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
    AS_Table.Extended[EventType %in% c("MXE", "SE", "AFE", "A5SS"), AA_downstr.B := Downstream.B$AA_seq]

    # Casette A
    Casette.A = AS_Table[, c("EventType", "EventID", "transcript_id_a", "intron_number_a")]
    Casette.A[, c("transcript_id", "exon_number") := list(transcript_id_a, intron_number_a)]
    Casette.A[EventType %in% c("MXE", "SE", "ALE", "A3SS"), exon_number := exon_number + 1]

    Casette.A = Proteins.Splice[Casette.A, on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Casette.A.gr = makeGRangesFromDataFrame(as.data.frame(na.omit(Casette.A)), keep.extra.columns = T)
    Casette.A.seq = getSeq(genome, Casette.A.gr)
    Casette.A[!is.na(seqnames),casette_seq := as.character(Casette.A.seq)]

    setnames(Casette.A, "phase", "phase_casette")
# Add nucleotides from upstream and downstream
    Casette.A = Upstream.A[Casette.A, on = "EventID", c("EventID", "phase_casette", "casette_seq", "seq")]
    setnames(Casette.A, "seq", "upstr_seq")
    Casette.A = Downstream.A[Casette.A, on = "EventID", c("EventID", "phase_casette", "casette_seq", "upstr_seq", "seq")]
    setnames(Casette.A, "seq", "Downstr_seq")
    
# Construct extended casette sequence:
    Casette.A[, casette_seq_extended := casette_seq]
    # Trim casette_seq_extended if upstream sequence does not exists
    Casette.A[!is.na(phase_casette) & is.na(upstr_seq), casette_seq_extended := 
        substr(casette_seq_extended, phase_casette + 1, nchar(casette_seq_extended)) ]    
    Casette.A[!is.na(phase_casette) & phase_casette > 0 & !is.na(upstr_seq), casette_seq_extended := paste0(
        substr(upstr_seq, nchar(upstr_seq) + 1 - phase_casette, nchar(upstr_seq)), casette_seq_extended)]
    Casette.A[nchar(casette_seq_extended) %% 3 > 0 & !is.na(Downstr_seq), casette_seq_extended := paste0(casette_seq_extended,
        substr(Downstr_seq, 1, 3 - (nchar(casette_seq_extended) %% 3)))]
# Translate:
    seq = Casette.A$casette_seq_extended[!is.na(Casette.A$casette_seq_extended)]
    # trim out-of-phase to be tidy:
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Casette.A[!is.na(casette_seq_extended), AA_seq := as.character(prot)]
    AS_Table.Extended[, AA_casette.A := Casette.A$AA_seq]
    
    # Casette B
    Casette.B = AS_Table[EventType != "SE", c("EventType", "EventID", "transcript_id_b", "intron_number_b")]
    Casette.B[, c("transcript_id", "exon_number") := list(transcript_id_b, intron_number_b)]
    Casette.B[EventType %in% c("MXE", "ALE", "A3SS"), exon_number := exon_number + 1]

    Casette.B = Proteins.Splice[Casette.B, on = c("transcript_id", "exon_number"), 
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Casette.B.gr = makeGRangesFromDataFrame(as.data.frame(na.omit(Casette.B)), keep.extra.columns = T)
    Casette.B.seq = getSeq(genome, Casette.B.gr)
    Casette.B[!is.na(seqnames),casette_seq := as.character(Casette.B.seq)]

    setnames(Casette.B, "phase", "phase_casette")
# Add nucleotides from upstream and downstream
    Casette.B = Upstream.B[Casette.B, on = "EventID", c("EventID", "phase_casette", "casette_seq", "seq")]
    setnames(Casette.B, "seq", "upstr_seq")
    Casette.B = Downstream.B[Casette.B, on = "EventID", c("EventID", "phase_casette", "casette_seq", "upstr_seq", "seq")]
    setnames(Casette.B, "seq", "Downstr_seq")
    
# Construct extended casette sequence:
    Casette.B[, casette_seq_extended := casette_seq]
    # Trim casette_seq_extended if upstream sequence does not exists
    Casette.B[!is.na(phase_casette) & is.na(upstr_seq), casette_seq_extended := 
        substr(casette_seq_extended, phase_casette + 1, nchar(casette_seq_extended)) ]    
    Casette.B[!is.na(phase_casette) & phase_casette > 0 & !is.na(upstr_seq), casette_seq_extended := paste0(
        substr(upstr_seq, nchar(upstr_seq) + 1 - phase_casette, nchar(upstr_seq)), casette_seq_extended)]
    Casette.B[nchar(casette_seq_extended) %% 3 > 0 & !is.na(Downstr_seq), casette_seq_extended := paste0(casette_seq_extended,
        substr(Downstr_seq, 1, 3 - (nchar(casette_seq_extended) %% 3)))]
# Translate:
    seq = Casette.B$casette_seq_extended[!is.na(Casette.B$casette_seq_extended)]
    # trim out-of-phase to be tidy:
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    prot = Biostrings::translate(as(seq, "DNAStringSet"))
    Casette.B[!is.na(casette_seq_extended), AA_seq := as.character(prot)]
    AS_Table.Extended[EventType != "SE", AA_casette.B := Casette.B$AA_seq]

    AS_Table.Extended[, AA_full.A := ""]
    AS_Table.Extended[!is.na(AA_upstr.A), AA_full.A := paste0(AA_full.A, AA_upstr.A)]
    AS_Table.Extended[!is.na(AA_casette.A), AA_full.A := paste0(AA_full.A, AA_casette.A)]
    AS_Table.Extended[!is.na(AA_downstr.A), AA_full.A := paste0(AA_full.A, AA_downstr.A)]
    AS_Table.Extended[, AA_full.B := ""]
    AS_Table.Extended[!is.na(AA_upstr.B), AA_full.B := paste0(AA_full.B, AA_upstr.B)]
    AS_Table.Extended[!is.na(AA_casette.B), AA_full.B := paste0(AA_full.B, AA_casette.B)]
    AS_Table.Extended[!is.na(AA_downstr.B), AA_full.B := paste0(AA_full.B, AA_downstr.B)]
    
    fst::write.fst(as.data.frame(AS_Table.Extended), paste(reference_path,"Splice.Extended.fst", sep="/"))

	message("Splice Annotations finished\n")
    
	message("Reference build finished")
}

grlGaps<-function(grl) {
	GenomicRanges::psetdiff(unlist(range(grl),use.names=TRUE),grl)
}

