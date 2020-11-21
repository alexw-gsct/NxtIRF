
#' @export
nxtIRF <- function(offline = FALSE, BPPARAM = BiocParallel::bpparam()) {

	assert_that(interactive(),
		msg = "NxtIRF App can only be run in interactive mode (i.e. RStudio).")

    # ah = AnnotationHub(localHub = offline)
	
	ui_dash <- dashboardPage(
		dashboardHeader(title = "NxtIRF"),
		ui_sidebar(),
		dashboardBody(
			tabItems(
                ui_tab_title(),
                ui_tab_system(),
                ui_tab_ref_new(),			

                ui_tab_expr(),
				
				tabItem(tabName = "navQC",
                    div(style = 'overflow-x: scroll',  DT::dataTableOutput('DT_QC'))
				),

                ui_tab_filter(),

                ui_tab_analyse(),

                ui_tab_diag(),

                ui_tab_volcano(),

                ui_tab_heatmap(),
				
                ui_tab_coverage()
			)
		)
	)
	
	runApp(shinyApp(ui_dash, dash_server))
		
}