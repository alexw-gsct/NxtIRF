## ---- eval = FALSE------------------------------------------------------------
#  library(NxtIRF)

## ---- eval = FALSE------------------------------------------------------------
#  dir.create("bams")      # Create the directory to contain the BAM files
#  download_NxtIRF_example(destination_dir = "./bams")

## ---- eval = FALSE------------------------------------------------------------
#  bamfiles.df = Find_Bams("./bams")

## ---- eval = FALSE------------------------------------------------------------
#  bamfiles.df = Find_Bams("./bams", use_subdir = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  bamfiles.df = Find_Bams("./bams", use_subdir = FALSE).
#  View(bamfiles.df)   # To view the output of this function in RStudio

## ---- eval = FALSE------------------------------------------------------------
#  IRFinder(
#      bamfiles = bamfiles.df$BAM,
#      sample_names = bamfiles.df$sample,
#      reference_path = "./Reference",
#      output_path = "./IRFinder_Output",
#      n_threads = 1
#  )

## ---- eval = FALSE------------------------------------------------------------
#  expr.df = Find_IRFinder_Output("./IRFinder_Output")
#  View(expr.df)

## ---- eval = FALSE------------------------------------------------------------
#  CollateData(
#      Experiment = expr.df,
#      reference_path = "./Reference",
#      output_path = "./NxtIRF_Output",
#      IRMode = "SpliceOverMax",       # Use IRMode = "SpliceMax" for legacy IRFinder calculation of IR-Ratio
#      low_memory_mode = TRUE,         # Only set to FALSE in large servers with lots of memory
#      n_threads = 1
#  )

## ---- eval = FALSE------------------------------------------------------------
#  colData = data.frame(sample = expr.df$sample)
#  colData$Treatment = data.table::tstrsplit(colData$sample, split="_")[[1]]
#  View(colData)

## ---- eval = FALSE------------------------------------------------------------
#  se = MakeSE(
#      fst_path = "./NxtIRF_Output",
#      colData = colData
#  )

## ---- eval = FALSE------------------------------------------------------------
#  filters = get_default_filters()
#  
#  se.filtered = se[apply_filters(se, filters),]

## ---- eval = FALSE------------------------------------------------------------
#  res.limma = limma_ASE(
#      se = se.filtered,
#      test_factor = "Treatment",
#      test_nom = "D2",
#      test_denom = "UT",
#  )

## ---- eval = FALSE------------------------------------------------------------
#  library(NxtIRF)
#  nxtIRF()

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Title page"----
knitr::include_graphics("img/Ref_1_title.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="System"----
knitr::include_graphics("img/System_single_thread.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Experiment"----
knitr::include_graphics("img/Expr_1_empty.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select Reference Path"----
knitr::include_graphics("img/Expr_2_ref_select.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Details of the Loaded Reference"----
knitr::include_graphics("img/Expr_2a_ref_display.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Select BAM Path"----
knitr::include_graphics("img/Expr_3_bam_select.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap = "List of BAM files"----
knitr::include_graphics("img/Expr_3a_bam_display.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Selecting BAM files to run IRFinder"----
knitr::include_graphics("img/Expr_4_IRF_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Confirm to run IRFinder"----
knitr::include_graphics("img/Expr_4_IRF_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="IRFinder is complete"----
knitr::include_graphics("img/Expr_4_IRF_3.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Add an annotation column"----
knitr::include_graphics("img/Expr_5_anno_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Annotate your samples"----
knitr::include_graphics("img/Expr_5_anno_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Collating the IRFinder output"----
knitr::include_graphics("img/Expr_6_Collate.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="NxtIRF collation is done!"----
knitr::include_graphics("img/Expr_6_Collate_Done.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Load the SummarizedExperiment"----
knitr::include_graphics("img/Expr_7_SE.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="QC tab"----
knitr::include_graphics("img/Expr_8_QC_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Read depth per sample"----
knitr::include_graphics("img/Expr_8_QC_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Comparisons between two QC parameters"----
knitr::include_graphics("img/Expr_8_QC_3.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Filters tab"----
knitr::include_graphics("img/Expr_9_Filter_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Loading and running default filters"----
knitr::include_graphics("img/Expr_9_Filter_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Selecing log10 scale"----
knitr::include_graphics("img/Expr_9_Filter_3.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Differential Expression Analysis tab"----
knitr::include_graphics("img/Expr_10_DE_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Settings to perform differential analysis on D2 vs UT samples"----
knitr::include_graphics("img/Expr_10_DE_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Differential analysis results"----
knitr::include_graphics("img/Expr_10_DE_3.png")

