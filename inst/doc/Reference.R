## ---- eval = FALSE------------------------------------------------------------
#  library(NxtIRF)

## ---- eval = FALSE------------------------------------------------------------
#  if(!dir.exists("Reference")) dir.create("Reference")
#  reference_path = "./Reference"

## ---- eval = FALSE------------------------------------------------------------
#  library(AnnotationHub)
#  ah = AnnotationHub()
#  ah.options = query(ah, c("Homo Sapiens", "Ensembl", "release-94"))

## ---- eval = FALSE------------------------------------------------------------
#  ah.options

## ---- eval = FALSE------------------------------------------------------------
#  ah_genome = "AH65745"
#  ah_transcriptome = "AH64631"

## ---- eval = FALSE------------------------------------------------------------
#  BuildReference(
#      ah_genome = ah_genome,
#      ah_transcriptome = ah_transcriptome,
#      reference_path = reference_path,
#      genome_type = "hg38")

## ---- eval = FALSE------------------------------------------------------------
#  fasta = "genome.fa"
#  gtf = "transcripts.gtf"
#  
#  BuildReference(
#      fasta = fasta,
#      gtf = gtf,
#      genome_type = "mm10"
#  )

## ---- eval = FALSE------------------------------------------------------------
#  run_IRFinder_GenerateMapReads(
#      genome.fa = "genome.fa",
#      out.fa = "mappability_reads.fa",
#      read_len = 70,
#      read_stride = 10,
#      error_pos = 35
#  )

## ---- eval = FALSE------------------------------------------------------------
#  run_IRFinder_MapExclusionRegions(
#      bamfile = "Aligned.bam",
#      output_file = "MappabilityExclusion.bed",
#      threshold = 4,
#      includeCov = FALSE
#  )

## ---- eval = FALSE------------------------------------------------------------
#  fasta = "genome.fa"
#  gtf = "transcripts.gtf"
#  MappabilityRef = "MappabilityExclusion.bed"
#  
#  BuildReference(
#      fasta = fasta,
#      gtf = gtf,
#      genome_type = "",
#      MappabilityRef = MappabilityRef
#  )

## ---- eval = FALSE------------------------------------------------------------
#  library(NxtIRF)
#  nxtIRF()

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Title page"----
knitr::include_graphics("img/Ref_1_title.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Reference"----
knitr::include_graphics("img/Ref_2_ref_empty.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Create a new folder for the Reference"----
knitr::include_graphics("img/Ref_3_path_create.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select Reference path dialog box"----
knitr::include_graphics("img/Ref_3_path.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select species"----
knitr::include_graphics("img/Ref_4_species.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select Ensembl release"----
knitr::include_graphics("img/Ref_5_version.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select GTF file"----
knitr::include_graphics("img/Ref_6_GTF.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select FASTA file"----
knitr::include_graphics("img/Ref_7_FASTA.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Supply custom FASTA/GTF files"----
knitr::include_graphics("img/Ref_8_USER.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select defaults for hg38"----
knitr::include_graphics("img/Ref_9_genome_type.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Ready to run Reference"----
knitr::include_graphics("img/Ref_10_final.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Running BuildReference"----
knitr::include_graphics("img/Ref_11_progress.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Reference build complete"----
knitr::include_graphics("img/Ref_12_complete.png")

