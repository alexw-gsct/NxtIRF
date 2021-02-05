## ---- eval = FALSE------------------------------------------------------------
#  library(NxtIRF)

## ---- eval = FALSE------------------------------------------------------------
#  nxtIRF()

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Loading the Reference and Experiment"----
knitr::include_graphics("img/Vis_1_intro.png")

## ---- eval = FALSE------------------------------------------------------------
#  setwd("/path/to/project")
#  expr.df = Find_IRFinder_Output("./IRFinder_output")
#  
#  colData = data.frame(sample = expr.df$sample)
#  colData$Treatment = data.table::tstrsplit(colData$sample, split="_")[[1]]
#  View(colData)

## ---- eval = FALSE------------------------------------------------------------
#  se = MakeSE(
#      fst_path = "./NxtIRF_Output",
#      colData = colData
#  )

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Loading the Differential Expression Analysis"----
knitr::include_graphics("img/Vis_2_DE.png")

## ---- eval = FALSE------------------------------------------------------------
#  res.limma = data.table::fread("res.limma.csv", data.table = FALSE)

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Volcano Plot"----
knitr::include_graphics("img/Vis_3_Volc_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Volcano Plot of Top 100 events"----
knitr::include_graphics("img/Vis_3_Volc_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Volcano Plot by Type of Alternative Splicing"----
knitr::include_graphics("img/Vis_3_Volc_3.png")

## ---- eval = FALSE------------------------------------------------------------
#  library(ggplot2)
#  
#  ggplot(res.limma, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point()

## ---- eval = FALSE------------------------------------------------------------
#  ggplot(res.limma, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point() +
#      facet_wrap(vars(EventType))

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Scatter Plot settings"----
knitr::include_graphics("img/Vis_4_Diag_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Scatter Plot"----
knitr::include_graphics("img/Vis_4_Diag_2.png")

## ---- eval = FALSE------------------------------------------------------------
#  ggplot(res.limma, aes(x = AvgPSI_D2, y = AvgPSI_UT)) + geom_point()

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Selecting Events"----
knitr::include_graphics("img/Vis_5_Gate_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Selected Top Events"----
knitr::include_graphics("img/Vis_5_Gate_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Scatter Plot highlighting Selected Events"----
knitr::include_graphics("img/Vis_5_Gate_3.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Heatmap of Selected Events"----
knitr::include_graphics("img/Vis_6_Heat_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Unselecting Selected Events"----
knitr::include_graphics("img/Vis_7_Reset_1.png")

## ---- eval = FALSE------------------------------------------------------------
#  mat = make_matrix(
#      se = se.filtered,
#      event_list = res.limma$EventName[1:20],
#      method = "PSI"      # use "logit" for logit-transformed values, or "Z-score" for Z-score transformed values
#  )

## ---- eval = FALSE------------------------------------------------------------
#  library("genefilter")
#  library("pheatmap")
#  library("RColorBrewer")
#  
#  color = colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(100)
#  
#  pheatmap(mat, color = color, breaks = seq(0, 1, length.out = 101),
#      main = "Top 20 ASE Events in sample files")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Coverage Plot Tab - Selecting Events"----
knitr::include_graphics("img/Vis_8_Cov_1.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Setting Tracks and Generating Plots"----
knitr::include_graphics("img/Vis_8_Cov_2.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Coverage Plot of Individual Tracks"----
knitr::include_graphics("img/Vis_8_Cov_3.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Coverage Plot by Condition"----
knitr::include_graphics("img/Vis_8_Cov_4.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Stacked Coverage Plot"----
knitr::include_graphics("img/Vis_8_Cov_5.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Stacked Coverage Plot with T-test"----
knitr::include_graphics("img/Vis_8_Cov_6.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Plot by Gene"----
knitr::include_graphics("img/Vis_8_Cov_7.png")

## ---- echo=FALSE, out.width='640pt', fig.align = 'center', fig.cap="Normalise by 3-prime most intron"----
knitr::include_graphics("img/Vis_8_Cov_8.png")

## ---- eval = FALSE------------------------------------------------------------
#  cov_data = prepare_covplot_data("./Reference")

## ---- eval = FALSE------------------------------------------------------------
#  i = 1
#  pl = Plot_Coverage(
#      se = se,
#      Event = res.limma$EventName[i],
#      cov_data = cov_data,
#      tracks = c("UT_1", "D2_1"),
#      )

## ---- eval = FALSE------------------------------------------------------------
#  pl$final_plot

## ---- eval = FALSE------------------------------------------------------------
#  egg::ggarrange(pl$ggplot[[1]], pl$ggplot[[2]], pl$ggplot[[6]], ncol = 1)

## ---- eval = FALSE------------------------------------------------------------
#  i = 1
#  pl = Plot_Coverage(
#      se = se,
#      Event = res.limma$EventName[i],
#      cov_data = cov_data,
#      tracks = c("UT", "D2"),
#      condition = "Treatment",
#      stack_tracks = TRUE,
#      t_test = TRUE
#      )

## ---- eval = FALSE------------------------------------------------------------
#  egg::ggarrange(pl$ggplot[[1]], pl$ggplot[[5]], pl$ggplot[[6]], ncol = 1)

## ---- eval = FALSE------------------------------------------------------------
#  i = 1
#  pl = Plot_Coverage(
#      se = se,
#      Event = res.limma$EventName[i],
#      cov_data = cov_data,
#      tracks = c("D2", "UT"),
#      condition = "Treatment",
#      stack_tracks = TRUE, t_test = TRUE,
#      Gene = "ANAPC2",
#      zoom_factor = 0,
#      strand = "-")
#  pl$final_plot

## ---- eval = FALSE------------------------------------------------------------
#  res.limma$EventName[which(grepl("ANAPC2", res.limma$EventName))]

## ---- eval = FALSE------------------------------------------------------------
#  pl = Plot_Coverage(
#      se = se,
#      Event = "ANAPC2/ENST00000323927_Intron12/clean",
#      cov_data = cov_data,
#      tracks = c("D2", "UT"),
#      condition = "Treatment",
#      stack_tracks = TRUE, t_test = TRUE,
#      Gene = "ANAPC2",
#      zoom_factor = 0,
#      strand = "-")
#  pl$final_plot

