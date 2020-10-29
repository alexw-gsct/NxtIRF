is_valid <- function(.) !is.null(.) && . != "" && . != "(none)"

#' @export
startNxtIRF <- function(offline = FALSE, BPPARAM = BiocParallel::bpparam()) {

	assertthat::assert_that(interactive(),
		msg = "NxtIRF App can only be run in interactive mode (i.e. RStudio).")

  ah = AnnotationHub::AnnotationHub(localHub = offline)

  # data.table::setDTthreads()

	filterModule_UI <- function(id, label = "Counter") {
		ns <- NS(id)
		wellPanel(
			h5(label),	# e.g. "Filter #1"
			selectInput(ns("filterClass"), "Filter Class", width = '100%', choices = c("", "Annotation", "Data")),
			selectInput(ns("filterType"), "Filter Type", width = '100%', choices = c("")),
      conditionalPanel(ns = ns,
        condition = "['Consistency'].indexOf(input.filterType) >= 0",
        shinyWidgets::sliderTextInput(ns("d1_1c"), "log-fold maximum", choices = seq(0.2, 5, by = 0.2), selected = 1)
      ),
      conditionalPanel(ns = ns,
        condition = "['Coverage'].indexOf(input.filterType) >= 0",
        sliderInput(ns("d1_1d"), "Percent Coverage", min = 0, max = 100, value = 80)
      ),
      conditionalPanel(ns = ns,
        condition = "['Depth'].indexOf(input.filterType) >= 0",
        shinyWidgets::sliderTextInput(ns("d1_1"), "Minimum", choices = c(1,2,3,5,10,20,30,50,100,200,300,500), selected = 20),
      ),
      conditionalPanel(ns = ns,
        condition = "['Depth', 'Coverage'].indexOf(input.filterType) >= 0",
        tagList(
          shinyWidgets::sliderTextInput(ns("d1_2"), "Minimum Conditions Satisfy Criteria (-1 = ALL)", choices = seq(-1,8), selected = -1),
          selectInput(ns("cond1"), "Condition", width = '100%',
            choices = c("")),
          sliderInput(ns("d1_3"), "Percent samples per condition satisfying criteria", min = 0, max = 100, value = 80)
        )
      ),
      conditionalPanel(ns = ns,
        condition = "['Coverage', 'Consistency'].indexOf(input.filterType) >= 0",
        shinyWidgets::sliderTextInput(ns("d1_1b"), "Signal Threshold to apply criteria", 
          choices = c(1,2,3,5,10,20,30,50,100,200,300,500), selected = 20),
      ),
      conditionalPanel(ns = ns,
        condition = "['Data'].indexOf(input.filterClass) >= 0",
        selectInput(ns("EventType"), "Splice Type", width = '100%', multiple = TRUE,
          choices = c("IR", "MXE", "SE", "AFE", "ALE", "A5SS", "A3SS"))
      )
		)
	}

	filterModule_server <- function(id, conditionList) {
		moduleServer(
			id,
			function(input, output, session) {
			# returns a list of filter options
				observeEvent(input$filterClass, {
					if(input$filterClass == "Annotation") {
						updateSelectInput(session = session, inputId = "filterType", 
							choices = c("Filter A", "Filter B"))
					} else if(input$filterClass == "Data") {
						updateSelectInput(session = session, inputId = "filterType", 
							choices = c("Depth", "Coverage", "Consistency"))
					} else {
						updateSelectInput(session = session, inputId = "filterType", 
							choices = c(""))						
					}
				})


				final <- reactiveValues(
          filterClass = "",
          filterType = "",
          filterVars = list()
        )

        observeEvent(input$filterType, {
          message("Update being triggered")
          if(input$filterType %in% c("Depth", "Coverage")) {
          message("Update triggered")
            updateSelectInput(session = session, inputId = "cond1", choices = conditionList())
          }
        })
        
        toListen <- reactive({
          list(input$filterClass, input$filterType,
             input$d1_1, input$d1_1b, input$d1_1c, input$d1_1d, 
             input$d1_2, input$d1_3,
             input$cond1, input$EventType
          )
        })
        
        observeEvent(toListen(), {
          final$filterClass = input$filterClass
          final$filterType = input$filterType
					if(final$filterType == "Depth") {
						final$filterVars$minimum = input$d1_1
					} else {
						final$filterVars$minimum = input$d1_1d
					}
					final$filterVars$maximum = input$d1_1c
          final$filterVars$minDepth = input$d1_1b
          final$filterVars$minCond = input$d1_2
          final$filterVars$condition = input$cond1
          final$filterVars$pcTRUE = input$d1_3
					final$filterVars$EventTypes = input$EventType
          if(length(final$filterVars) > 0 && all(sapply(final$filterVars, is_valid))) {
            final$trigger = runif(1)
          } else {
            final$trigger = NULL
          }
        })
        
				# Returns filter list from module
				return(final)
			}
		)
	}

	ui <- navbarPage("NxtIRF", id = "navSelection",
# Title Page
		tabPanel("About", value = "navTitle",
			img(src="https://pbs.twimg.com/profile_images/1310789966293655553/7HawCItY_400x400.jpg")
		),
		tabPanel("Multithreading", value = "navThreads",
			wellPanel(
				selectInput('expr_Cores', 'Number of Processors to Use', width = '100%',
					choices = 1, selected = 1),					
			),
		),		
# Reference
		navbarMenu("Reference",
		# New reference
			tabPanel("New", value = "navRef_New",
				fluidRow(
					column(5,
						h4("Select Reference Directory"),
						shinyDirButton("dir_reference_path", label = "Choose reference path", title = "Choose reference path"),
							verbatimTextOutput("txt_reference_path"),
						br(),
						h4("Select Ensembl Reference from AnnotationHub"),
						selectInput('newrefAH_Species', 'Select Species', width = '100%',
							choices = c("")),
						selectInput('newrefAH_Version_Trans', 'Select Transcriptome Version', width = '100%',
							choices = c("")),
						selectInput('newrefAH_Trans', 'Select Transcriptome Reference', width = '100%',
							choices = c("")),
						selectInput('newrefAH_Assembly', 'Select Genome Assembly', width = '100%',
							choices = c("")),
						selectInput('newrefAH_Version_Genome', 'Select Genome Version', width = '100%',
							choices = c("")),
						selectInput('newrefAH_Genome', 'Select Genome Reference', width = '100%',
							choices = c("")),
						br(),
						h4("or select Reference from File"),
						br(),
						shinyFilesButton("file_genome", label = "Choose genome FASTA File", title = "Choose genome FASTA File", multiple = FALSE),
						verbatimTextOutput("txt_genome"),
						shinyFilesButton("file_gtf", label = "Choose transcriptome GTF File", title = "Choose transcriptome GTF File", multiple = FALSE),
						verbatimTextOutput("txt_gtf"),
					),
					column(5,
						selectInput('newref_genome_type', 'Select Genome Type to set Mappability and non-PolyA files (leave empty to reset)', 
							c("hg38", "mm10", "hg19", "mm9", "other", "custom")),
						shinyFilesButton("file_mappa", label = "Choose Mappability Exclusion file", 
							title = "Choose Mappability Exclusion file", multiple = FALSE),
						verbatimTextOutput("txt_mappa"),
						shinyFilesButton("file_NPA", label = "Choose non-PolyA BED file", title = "Choose non-PolyA BED file", multiple = FALSE),
						verbatimTextOutput("txt_NPA"),
						shinyFilesButton("file_bl", label = "Choose blacklist BED file", title = "Choose blacklist BED file", multiple = FALSE),
						verbatimTextOutput("txt_bl"),
						br(),
						actionButton("buildRef", "Build Reference"),
						actionButton("clearNewRef", "Clear settings"),
						uiOutput("refStatus")
					)
				)
			),
		# Load reference
			tabPanel("Load", value = "navRef_Load",
				fluidRow(
					column(9,
						h4("Select Reference Directory"),
						shinyDirButton("dir_reference_path_load", label = "Choose reference path", title = "Choose reference path"),
							textOutput("txt_reference_path_load"),
						br(),
						actionButton("clearLoadRef", "Clear settings"), # TODO
						br(),
						textOutput("loadRef_field1"), br(),
						textOutput("loadRef_field2"), br(),
						textOutput("loadRef_field3"), br(),
						textOutput("loadRef_field4"), br(),
						textOutput("loadRef_field5"), br(),
						textOutput("loadRef_field6"), br(),
						textOutput("loadRef_field7"),						
					)
				)
			),
			tabPanel("View", value = "navRef_View",
				fluidRow(style='height:20vh',
					column(4, 
            verbatimTextOutput("warning_view_ref"),
						selectizeInput('genes_view', 'Genes', choices = "(none)"),
						selectInput('events_view', 'Events', multiple = TRUE,
							c("(none)"))
          ),
					column(8,
						div(style="display: inline-block;vertical-align:top; width: 80px;",
							selectInput("chr_view_ref", label = "Chr", c("(none)"), selected = "(none)")),
						div(style="display: inline-block;vertical-align:top; width: 120px;",
							textInput("start_view_ref", label = "Left", c(""))),
						div(style="display: inline-block;vertical-align:top; width: 120px;",
							textInput("end_view_ref", label = "Right", c(""))),						
						div(style="display: inline-block;vertical-align:top; width: 250px;",
							sliderInput("zoom_view_ref", label = "Zoom", value = 0, min = 0, max = 12)),
						div(style="display: inline-block;vertical-align:top;",
              shinyWidgets::materialSwitch("move_labels_view_ref", label = "Move transcript labels"))
					)
				),
				fluidRow(
						plotlyOutput("plot_view_ref", height = "800px"),
            textOutput("plotly_event")
				)
      )
		),
		# Experiment
		tabPanel("Experiment", value = "navExpr",
			fluidRow(
				column(4,
					textOutput("txt_reference_path_expr"), # done
					br(),
					
					shinyDirButton("dir_bam_path_load", 
						label = "Choose BAM path", title = "Choose BAM path"), # done
					textOutput("txt_bam_path_expr"),
					br(),
					
					shinyDirButton("dir_irf_path_load", 
						label = "Choose IRFinder output path", title = "Choose IRFinder output path"), # done					
					textOutput("txt_irf_path_expr"),
					br(),
          actionButton("run_irf_expr", "Run IRFinder on selected bam files"), # TODO
					textOutput("txt_run_irf_expr"),
					br(),
					conditionalPanel(
						condition = "['Annotations'].indexOf(input.hot_switch_expr) >= 0",
						wellPanel(
							uiOutput("newcol_expr"), # done
							div(class='row',
								div(class= "col-sm-6",
									radioButtons("type_newcol_expr", "Type", c("character", "integer", "double"))
								),
								div(class = "col-sm-6", 
									actionButton("addcolumn_expr", "Add"), br(),  # done
									actionButton("removecolumn_expr", "Remove") # done
								)
							)
						)
					),

					shinyDirButton("dir_collate_path_load", 
						label = "Choose NxtIRF FST output path", title = "Choose NxtIRF FST output path"), # done
					textOutput("txt_collate_path_expr"), # done
					br(),
					
					actionButton("run_collate_expr", "Compile NxtIRF FST files"),
					textOutput("txt_run_col_expr")
				),
				column(8,
					shinyFilesButton("loadexpr_expr", label = "Load Experiment", 
						title = "Load Experiment Data Frame", multiple = FALSE), # done
					shinySaveButton("saveexpr_expr", "Save Experiment", "Save Experiment as...", 
						filetype = list(dataframe = "csv")),
					actionButton("clear_expr", "Clear Experiment"),
					textOutput("txt_run_save_expr"),
          br(),
					shinyFilesButton("file_expr_path_load", label = "Choose Sample Annotation Table", 
						title = "Choose Sample Annotation Table", multiple = FALSE), # done
					textOutput("txt_sample_anno_expr"), # done
          br(),
					actionButton("build_expr", "Build SummarizedExperiment"),

					shinyWidgets::radioGroupButtons(
						 inputId = "hot_switch_expr",
						 label = "Experiment Display",
						 choices = c("Files", "Annotations"),
						 selected = "Files"
					),
					conditionalPanel(
						condition = "['Files'].indexOf(input.hot_switch_expr) >= 0",
						rHandsontableOutput("hot_files_expr")
					),
					conditionalPanel(
						condition = "['Annotations'].indexOf(input.hot_switch_expr) >= 0",
						rHandsontableOutput("hot_anno_expr")
					)
				)	# last column
			),
		),
		# Analyze
		navbarMenu("Analyse",
			tabPanel("View QC", value = "navQC"),
			tabPanel("Filters", value = "navFilter",
		# Takes experimental data frame, sets filters, then constructs SummarizedExperiment object
				# Current Experiment
        fluidRow(
					column(4,	
						wellPanel(style = "overflow-y:scroll; max-height: 800px",
							filterModule_UI("filter1", "Filter #1"),
							filterModule_UI("filter2", "Filter #2"),
							filterModule_UI("filter3", "Filter #3"),
							filterModule_UI("filter4", "Filter #4")
						)
					),
					column(4,	
						wellPanel(style = "overflow-y:scroll; max-height: 800px",
							filterModule_UI("filter5", "Filter #5"),
							filterModule_UI("filter6", "Filter #6"),
							filterModule_UI("filter7", "Filter #7"),
							filterModule_UI("filter8", "Filter #8")
						)
					),
					column(4,	
						textOutput("current_expr_Filters"), br(),
						textOutput("current_ref_Filters"), br(),
						# actionButton("load_filterdata_Filters", "Load Data"),
						actionButton("refresh_filters_Filters", "Refresh Filters"),
						plotlyOutput("plot_filtered_Events"),
						selectInput('graphscale_Filters', 'Y-axis Scale', width = '100%',
							choices = c("linear", "log10")), 
						shinySaveButton("saveAnalysis_Filters", "Save SummarizedExperiment", "Save SummarizedExperiment as...", 
							filetype = list(RDS = "Rds")),
					)
        )
      ),
			tabPanel("Differential Analysis", value = "navAnalyse",
        fluidRow(
					column(4,	
						textOutput("warning_DE"),
            selectInput('method_DE', 'Method', 
							c("DESeq2", "limma", "DSS")),
            selectInput('variable_DE', 'Variable', 
							c("(none)")),
            selectInput('nom_DE', 'Nominator', 
							c("(none)")),
            selectInput('denom_DE', 'Denominator', 
							c("(none)")),
            selectInput('batch1_DE', 'Batch Factor 1', 
							c("(none)")),
            selectInput('batch2_DE', 'Batch Factor 2', 
							c("(none)")),
						actionButton("perform_DE", "Perform DE"),
					shinyFilesButton("load_DE", label = "Load DE", 
							title = "Load DE from RDS", multiple = FALSE),
						shinySaveButton("save_DE", "Save DE", "Save DE as...", 
							filetype = list(RDS = "Rds")),
					),
					column(8,	
            DT::dataTableOutput("DT_DE")
          )
        )
      )	# DESeq2 or DSS

		),

		navbarMenu("Display"

		) # last navbar

	)



	server = function(input, output, session) {

		settings_newref <- shiny::reactiveValues(
			newref_path = "",
			newref_fasta = "",
			newref_gtf = "",
			newref_AH_fasta = "",
			newref_AH_gtf = "",
			newref_mappa = "",
			newref_NPA = "",
			newref_bl = ""
		)

    default_volumes <- c("Working Directory" = getwd(), "Home" = fs::path_home(), getVolumes()())
    addit_volume = c()

	# tabEvent Observer
		observeEvent(input$navSelection, {
			if(input$navSelection == "navRef_New") {
				# autopopulate if previous settings detectedrend
				if(settings_newref$newref_path != "") output$txt_reference_path <- renderText(settings_newref$newref_path)
				if(input$newrefAH_Species != "") {
				} else {
					ah.filtered = ah[ah$dataprovider == "Ensembl"]
					ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
					ah.filtered = ah.filtered[ah.filtered$sourcetype == "GTF"]
					updateSelectInput(session = session, inputId = "newrefAH_Species", 
						choices = c("", sort(unique(ah.filtered$species))))				
				}
			} else if(input$navSelection == "navRef_Load") {
        if(settings_loadref$loadref_path != "") {
          load_ref()
        }
        output$txt_reference_path_load <- renderText({
            validate(need(settings_loadref$loadref_path, "Please select reference path"))
            settings_loadref$loadref_path
        })
			} else if(input$navSelection == "navRef_View") {
        if(settings_loadref$loadref_path != "") {
          load_ref()
        }
        output$warning_view_ref <- renderText({
            validate(need(settings_loadref$loadref_path, "Please select reference path"))            
            settings_loadref$loadref_path
        })
        req(settings_loadref$loadref_path)
        settings_ViewRef$gene_list <- getGeneList()
				
        settings_ViewRef$elem.DT <- loadViewRef()
        settings_ViewRef$transcripts.DT <- loadTranscripts()
        
				if(!is.null(settings_ViewRef$gene_list)) {
          message(paste("Populating drop-down box with", 
            length(unique(settings_ViewRef$gene_list$gene_display_name)),"genes"))
					updateSelectInput(session = session, inputId = "chr_view_ref", 
						choices = c("(none)", sort(unique(settings_ViewRef$gene_list$seqnames))), selected = "(none)")    								          
					updateSelectizeInput(session = session, inputId = "genes_view", server = TRUE,
						choices = c("(none)", settings_ViewRef$gene_list$gene_display_name), selected = "(none)")    								
				} else {
					updateSelectInput(session = session, inputId = "chr_view_ref", 
						choices = c("(none)", settings_ViewRef$gene_list$gene_display_name), selected = "(none)")    								
					updateSelectizeInput(session = session, inputId = "genes_view", server = TRUE,
						choices = c("(none)"), selected = "(none)") 
				}
        
			} else if(input$navSelection == "navThreads") {	
				max_cores = parallel::detectCores() - 2
				updateSelectInput(session = session, inputId = "expr_Cores", 
					choices = seq(max_cores, 1), selected = max_cores)        								
			} else if(input$navSelection == "navExpr") {
				# Determine IRFinder cores

				output$txt_reference_path <- renderText({
					validate(
						need(settings_loadref$loadref_path, "Please Reference->Load and select reference path")
					)
					settings_loadref$loadref_path
				})
			} else if(input$navSelection == "navFilter") {
				if(!is.null(settings_SE$se)) {
					output$current_expr_Filters = renderText("SummarizedExperiment loaded")
				} else {
					output$current_expr_Filters = renderText("Please load SummarizedExperiment first")
				}
				if(settings_loadref$loadref_path != "") {
					output$current_ref_Filters = renderText("Reference loaded")
				} else {
					output$current_ref_Filters = renderText("Please load reference first")
				}
        if(!is.null(settings_SE$se) & settings_loadref$loadref_path != "") {
          processFilters()
        }
			} else if(input$navSelection == "navAnalyse") {
        output$warning_DE = renderText({
          validate(need(settings_SE$se, "Please load experiment via 'Experiment' tab"))
          "Experiment Loaded"
        })
        req(settings_SE$se)
        colData = SummarizedExperiment::colData(settings_SE$se)
        updateSelectInput(session = session, inputId = "variable_DE", 
          choices = c("(none)", colnames(colData)), selected = "(none)")
        updateSelectInput(session = session, inputId = "batch1_DE", 
          choices = c("(none)", colnames(colData)), selected = "(none)")						
        updateSelectInput(session = session, inputId = "batch2_DE", 
          choices = c("(none)", colnames(colData)), selected = "(none)")
			}
		})
    
		observe({  
			shinyDirChoose(input, "dir_reference_path", roots = c(default_volumes, addit_volume), session = session)
			output$txt_reference_path <- renderText({
					validate(need(input$dir_reference_path, "Please select reference path"))
					settings_newref$newref_path = parseDirPath(c(default_volumes, addit_volume), input$dir_reference_path)
			})
		})
		observe({
			shinyFileChoose(input, "file_genome", roots = c(default_volumes, addit_volume), 
				session = session, filetypes = c("fa", "fasta", "gz"))
			if(!is.null(input$file_genome)){
			 file_selected<-parseFilePaths(c(default_volumes, addit_volume), input$file_genome)
			 settings_newref$newref_fasta = as.character(file_selected$datapath)
			 output$txt_genome <- renderText(as.character(file_selected$datapath))
			}
		})
		observe({  
			shinyFileChoose(input, "file_gtf", roots = c(default_volumes, addit_volume), 
				session = session, filetypes = c("gtf", "gz"))
			if(!is.null(input$file_gtf)){
			 file_selected<-parseFilePaths(c(default_volumes, addit_volume), input$file_gtf)
			 settings_newref$newref_gtf = as.character(file_selected$datapath)
			 output$txt_gtf <- renderText(as.character(file_selected$datapath))
			}
		})
		observe({  
		shinyFileChoose(input, "file_mappa", roots = c(default_volumes, addit_volume), 
			session = session, filetypes = c("txt", "gz"))
			if(!is.null(input$file_mappa)){
			 file_selected<-parseFilePaths(c(default_volumes, addit_volume), input$file_mappa)
			 settings_newref$newref_mappa = as.character(file_selected$datapath)
			}
		})
    observeEvent(settings_newref$newref_mappa, {
      req(settings_newref$newref_mappa)
			output$txt_mappa <- renderText(settings_newref$newref_mappa)
    })
		observe({  
			shinyFileChoose(input, "file_NPA", roots = c(default_volumes, addit_volume), 
				session = session, filetypes = c("bed", "txt", "gz"))
			if(!is.null(input$file_NPA)){
				file_selected<-parseFilePaths(c(default_volumes, addit_volume), input$file_NPA)
				settings_newref$newref_NPA = as.character(file_selected$datapath)
			}
		})
    observeEvent(settings_newref$newref_NPA, {
      req(settings_newref$newref_NPA)
			output$txt_NPA <- renderText(settings_newref$newref_NPA)    
    })
		observe({  
			shinyFileChoose(input, "file_bl", roots = c(default_volumes, addit_volume), 
				session = session, filetypes = c("bed", "txt", "gz"))
			if(!is.null(input$file_bl)){
			 file_selected<-parseFilePaths(c(default_volumes, addit_volume), input$file_bl)
			 settings_newref$newref_bl = as.character(file_selected$datapath)
			}
		})
    observeEvent(settings_newref$newref_bl, {
      req(settings_newref$newref_bl)
			output$txt_bl <- renderText(settings_newref$newref_bl)
    })
		observeEvent(input$newref_genome_type, {
			if(input$newref_genome_type == "hg38") {
				settings_newref$newref_NPA = system.file("extra-input-files/Human_hg38_nonPolyA_ROI.bed", package = "NxtIRF")
        settings_newref$newref_mappa = system.file("extra-input-files/Mappability_Regions_hg38_v94.txt.gz", package = "NxtIRF")
			} else if(input$newref_genome_type == "hg19")  {
				settings_newref$newref_NPA = system.file("extra-input-files/Human_hg19_nonPolyA_ROI.bed", package = "NxtIRF")
        settings_newref$newref_mappa = system.file("extra-input-files/Mappability_Regions_hg19_v75.txt.gz", package = "NxtIRF")
			} else if(input$newref_genome_type == "mm10")  {
				settings_newref$newref_NPA = system.file("extra-input-files/Mouse_mm10_nonPolyA_ROI.bed", package = "NxtIRF")
        settings_newref$newref_mappa = system.file("extra-input-files/Mappability_Regions_mm10_v94.txt.gz", package = "NxtIRF")
			} else if(input$newref_genome_type == "mm9")  {
				settings_newref$newref_NPA = system.file("extra-input-files/Mouse_mm9_nonPolyA_ROI.bed", package = "NxtIRF")
        settings_newref$newref_mappa = system.file("extra-input-files/Mappability_Regions_mm9_v67.txt.gz", package = "NxtIRF")
			} else if(input$newref_genome_type == "other") {
        # do nothing. This allows user to first select the default and then change to user-defined files
			} else {
				settings_newref$newref_NPA = ""
        settings_newref$newref_mappa = ""
			}
		})

		observeEvent(input$newrefAH_Species, {
			if(input$newrefAH_Species != "") {
				ah.filtered = ah[ah$dataprovider == "Ensembl"]
				ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
				ah.filtered = ah.filtered[ah.filtered$sourcetype == "GTF"]
				ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]

				queryfirst = data.table::tstrsplit(ah.filtered[1]$sourceurl, split="/", fixed=TRUE)
				query_index = which(sapply(queryfirst, function(x) grepl("release", x)))
				choices = unlist(unique(data.table::tstrsplit(ah.filtered$sourceurl, split="/", fixed=TRUE)[query_index]))
				updateSelectInput(session = session, inputId = "newrefAH_Version_Trans", 
					choices = c("", sort(choices)))
			} else {
				updateSelectInput(session = session, inputId = "newrefAH_Version_Trans", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Trans", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Assembly", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Version_Genome", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Genome", 
					choices = c(""))        
			}
		})
		observeEvent(input$newrefAH_Version_Trans, {
			if(input$newrefAH_Version_Trans != "") {
				ah.filtered = ah[ah$dataprovider == "Ensembl"]
				ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
				ah.filtered = ah.filtered[ah.filtered$sourcetype == "GTF"]
				ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]
				ah.filtered = ah.filtered[grepl(input$newrefAH_Version_Trans, ah.filtered$sourceurl)]

				updateSelectInput(session = session, inputId = "newrefAH_Trans", 
					choices = c("", paste(names(ah.filtered), basename(ah.filtered$sourceurl), sep=": ")))
									
				# Also search for compatible genome
				genomes_avail = ah.filtered$genome
				ah.filtered = ah[ah$dataprovider == "Ensembl"]
				ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
				ah.filtered = ah.filtered[ah.filtered$sourcetype == "FASTA"]

				if(any(ah.filtered$genome %in% genomes_avail)) {
						ah.filtered = ah.filtered[ah.filtered$genome %in% genomes_avail]                
						updateSelectInput(session = session, inputId = "newrefAH_Assembly", 
								choices = c("", unique(ah.filtered$genome)))
				} else {
						updateSelectInput(session = session, inputId = "newrefAH_Assembly", 
								choices = c(""))            
				}
			} else {
				updateSelectInput(session = session, inputId = "newrefAH_Trans", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Assembly", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Version_Genome", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Genome", 
					choices = c(""))        
			}
		})
		observeEvent(input$newrefAH_Assembly, {
			if(input$newrefAH_Assembly != "") {
				ah.filtered = ah[ah$dataprovider == "Ensembl"]
				ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
				ah.filtered = ah.filtered[ah.filtered$sourcetype == "FASTA"]
				ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]
				ah.filtered = ah.filtered[ah.filtered$genome == input$newrefAH_Assembly]

				queryfirst = data.table::tstrsplit(ah.filtered[1]$sourceurl, split="/", fixed=TRUE)
				query_index = which(sapply(queryfirst, function(x) grepl("release", x)))
				choices = unlist(unique(data.table::tstrsplit(ah.filtered$sourceurl, split="/", fixed=TRUE)[query_index]))
				updateSelectInput(session = session, inputId = "newrefAH_Version_Genome", 
					choices = c("", sort(choices)))
			} else {
				updateSelectInput(session = session, inputId = "newrefAH_Version_Genome", 
					choices = c(""))
				updateSelectInput(session = session, inputId = "newrefAH_Genome", 
					choices = c(""))        
			}
		})
		observeEvent(input$newrefAH_Version_Genome, {
			if(input$newrefAH_Version_Genome != "") {
				ah.filtered = ah[ah$dataprovider == "Ensembl"]
				ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
				ah.filtered = ah.filtered[ah.filtered$sourcetype == "FASTA"]
				ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]
				ah.filtered = ah.filtered[ah.filtered$genome == input$newrefAH_Assembly]
				ah.filtered = ah.filtered[grepl(input$newrefAH_Version_Genome, ah.filtered$sourceurl)]

				updateSelectInput(session = session, inputId = "newrefAH_Genome", 
					choices = c("", paste(names(ah.filtered), basename(ah.filtered$sourceurl), sep=": ")))
			} else {
				updateSelectInput(session = session, inputId = "newrefAH_Genome", 
					choices = c(""))        
			}
		})
		observeEvent(input$newrefAH_Genome, {
			if(input$newrefAH_Genome != "") {
				settings_newref$newref_AH_fasta = input$newrefAH_Genome
			} else {
				settings_newref$newref_AH_fasta = ""
			}
		})
		observeEvent(input$newrefAH_Trans, {
			if(input$newrefAH_Trans != "") {
				settings_newref$newref_AH_gtf = input$newrefAH_Trans
			} else {
				settings_newref$newref_AH_gtf = ""
			}
		})
		
		# buildRef Button
		observeEvent(input$buildRef, {
			args <- list(reference_path = settings_newref$newref_path,
				ah_genome_tmp = settings_newref$newref_AH_fasta, ah_gtf_tmp = settings_newref$newref_AH_gtf, 
				fasta = settings_newref$newref_fasta, gtf = settings_newref$newref_gtf,
				genome_type = "Interactive", nonPolyARef = settings_newref$newref_NPA, MappabilityRef = settings_newref$newref_mappa,
				BlacklistRef = settings_newref$newref_bl)

			args <- Filter(is_valid, args)
		
			if(!("reference_path" %in% names(args))) {
				output$refStatus = renderText({ "Reference path not set" })
			} else if(!any(c("fasta", "ah_genome_tmp") %in% names(args))) {
				output$refStatus = renderText({ "Genome not provided" })        
			} else if(!any(c("gtf", "ah_gtf_tmp") %in% names(args))) {
				output$refStatus = renderText({ "Gene annotations not provided" })
			} else {        
				args.df = as.data.frame(t(as.data.frame(args)))
				colnames(args.df) = "value"
				# data.table::fwrite(args.df, file.path(args$reference_path, "settings_newref.csv"), row.names = TRUE)
				if("ah_genome_tmp" %in% names(args)) {
					args$ah_genome = data.table::tstrsplit(args$ah_genome_tmp, split=":", fixed=TRUE)[[1]]
					args$ah_genome_tmp = NULL
				}
				if("ah_gtf_tmp" %in% names(args)) {
					args$ah_transcriptome = data.table::tstrsplit(args$ah_gtf_tmp, split=":", fixed=TRUE)[[1]]
					args$ah_gtf_tmp = NULL
				}
				do.call(BuildReference, args)
				# If successfully created, load this reference automatically
				if(file.exists(file.path(settings_newref$newref_path, "settings.Rds"))) {
					settings_loadref$loadref_path = settings_newref$newref_path
				}
			}
		})
		
		# clearNewRef Button
		observeEvent(input$clearNewRef, {
			settings_newref$newref_path = getwd()
			settings_newref$newref_fasta = ""
			settings_newref$newref_gtf = ""
			settings_newref$newref_AH_fasta = ""
			settings_newref$newref_AH_gtf = ""
			settings_newref$newref_mappa = ""
			settings_newref$newref_NPA = ""
			settings_newref$newref_bl = ""
			output$txt_reference_path <- renderText("Please select reference path")
			output$txt_genome <- renderText("")
			output$txt_gtf <- renderText("")
			output$txt_mappa <- renderText("")
			output$txt_NPA <- renderText("")
			output$txt_bl <- renderText("")
      updateSelectInput(session = session, inputId = "newrefAH_Species", 
          choices = c(""))
      updateSelectInput(session = session, inputId = "newrefAH_Version_Trans", 
          choices = c(""))
      updateSelectInput(session = session, inputId = "newrefAH_Trans", 
          choices = c(""))
      updateSelectInput(session = session, inputId = "newrefAH_Assembly", 
          choices = c(""))
      updateSelectInput(session = session, inputId = "newrefAH_Version_Genome", 
          choices = c(""))
      updateSelectInput(session = session, inputId = "newrefAH_Genome", 
          choices = c(""))        
		})

# Load Reference Page
		settings_loadref <- shiny::reactiveValues(
			loadref_path = "",
			settings = c()
		)
		shinyDirChoose(input, "dir_reference_path_load", roots = c(default_volumes, addit_volume), session = session)
		observeEvent(input$dir_reference_path_load,{  
      req(input$dir_reference_path_load)
      settings_loadref$loadref_path = parseDirPath(c(default_volumes, addit_volume), 
				input$dir_reference_path_load)
    })
		load_ref = function() {
      req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds")))
			settings_loadref$settings = readRDS(file.path(settings_loadref$loadref_path, "settings.Rds"))
      if("reference_path" %in% names(settings_loadref$settings)) {
        if("ah_genome" %in% names(settings_loadref$settings)) {
          output$loadRef_field1 <- renderText({
            paste("AnnotationHub genome:",
              settings_loadref$settings$ah_genome, "\n",
              ah$description[which(names(ah) == settings_loadref$settings$ah_genome)], "\n",
              ah$sourceurl[which(names(ah) == settings_loadref$settings$ah_genome)]
            )
          })
        }
        if("ah_transcriptome" %in% names(settings_loadref$settings)) {
          output$loadRef_field2 <- renderText({
            paste("AnnotationHub gene annotations:",
              settings_loadref$settings$ah_transcriptome, "\n",
              ah$description[which(names(ah) == settings_loadref$settings$ah_transcriptome)], "\n",
              ah$sourceurl[which(names(ah) == settings_loadref$settings$ah_transcriptome)]
            )
          })
        }
        if("fasta_file" %in% names(settings_loadref$settings)) {
          output$loadRef_field3 <- renderText({
            paste("Genome FASTA file (user):",
              settings_loadref$settings$fasta_file
            )
          })					
        }
        if("gtf_file" %in% names(settings_loadref$settings)) {
          output$loadRef_field4 <- renderText({
            paste("Annotation GTF file (user):",
              settings_loadref$settings$gtf_file
            )
          })					
        }
        if("MappabilityRef" %in% names(settings_loadref$settings)) {
          output$loadRef_field5 <- renderText({
            paste("Mappability Exclusion file:",
              settings_loadref$settings$MappabilityRef
            )
          })					
        }
        if("nonPolyARef" %in% names(settings_loadref$settings)) {
          output$loadRef_field6 <- renderText({
            paste("Non-PolyA file:",
              settings_loadref$settings$nonPolyARef
            )
          })					
        }
        if("BlacklistRef" %in% names(settings_loadref$settings)) {
          output$loadRef_field7 <- renderText({
            paste("Blacklist Exclusion file:",
              settings_loadref$settings$BlacklistRef
            )
          })					
        }
      } else {
        settings_loadref$loadref_path = ""
      }
		}

		observeEvent(settings_loadref$loadref_path,{
      req(settings_loadref$loadref_path)
			if(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) {
				load_ref()
			}
			output$txt_reference_path_load <- renderText({
					validate(need(settings_loadref$loadref_path, "Please select reference path"))
					settings_loadref$loadref_path
			})
		})
# View Ref page

    settings_ViewRef <- shiny::reactiveValues(
			gene_list = NULL,
      elem.DT = NULL,
      transcripts.DT = NULL,
      view_chr = "",
      view_start = "",
      view_end = "",
      data_start = 0,
      data_end = 0,
      blockView = FALSE,
      repan = NULL,
      plot_ini = FALSE
		)
		
    loadTranscripts <- function() {
      req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 
			file_path = file.path(settings_loadref$loadref_path, "fst", "Transcripts.fst")
      
      Transcripts.DT = as.data.table(fst::read.fst(file_path))

      if("transcript_support_level" %in% colnames(Transcripts.DT)) {
          Transcripts.DT$transcript_support_level = tstrsplit(Transcripts.DT$transcript_support_level, split=" ")[[1]]
          Transcripts.DT$transcript_support_level[is.na(Transcripts.DT$transcript_support_level)] = "NA"
      } else {
        Transcripts.DT$transcript_support_level = 1
      }
      
      return(Transcripts.DT)
    }



    loadViewRef <- function() {
      req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 
			dir_path = file.path(settings_loadref$loadref_path, "fst")

      exons.DT = as.data.table(fst::read.fst(file.path(dir_path, "Exons.fst")))
      exons.DT = exons.DT[transcript_id != "protein_coding"]

      protein.DT = as.data.table(fst::read.fst(file.path(dir_path, "Proteins.fst")))
      misc.DT = as.data.table(fst::read.fst(file.path(dir_path, "Misc.fst")))
      introns.DT = as.data.table(fst::read.fst(file.path(dir_path, "junctions.fst")))
      introns.DT[, type := "intron"]
      
      total.DT = rbindlist(list(
        exons.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")],
        protein.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")],
        misc.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")],
        introns.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")]
      ))
      return(total.DT)
    }
    
		getGeneList <- function() {
      req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 
			file_path = file.path(settings_loadref$loadref_path, "fst", "Genes.fst")
			if(!file.exists(file_path)) return(NULL)
			
			df = as.data.table(fst::read.fst(file_path))
			return(df)
		}

    plot_view_ref_fn = function(chr, start, end) {
      view_chr = input$chr_view_ref
      view_start = as.numeric(input$start_view_ref)
      view_end = as.numeric(input$end_view_ref)

      req(view_chr)
      req(view_start)
      req(view_end)
      req(view_chr != settings_ViewRef$view_chr | view_start != settings_ViewRef$view_start | 
        view_end != settings_ViewRef$view_end)

      settings_ViewRef$view_chr = view_chr
      settings_ViewRef$view_start = view_start
      settings_ViewRef$view_end = view_end

      # set data_start and data_end at 3X zoom outside the view start / end
      view_center = round(0.5 * (view_start + view_end))
      
      message("plot_view_ref_fn")

      data_start = view_start - (view_end - view_start)
      data_end = view_end + (view_end - view_start)

      settings_ViewRef$data_start = data_start
      settings_ViewRef$data_end = data_end


      if(is.null(settings_ViewRef$elem.DT)) settings_ViewRef$elem.DT <- loadViewRef()
      if(is.null(settings_ViewRef$transcripts.DT)) settings_ViewRef$transcripts.DT <- loadTranscripts()

      req(settings_ViewRef$elem.DT)
      req(settings_ViewRef$transcripts.DT)
      message("stuff loaded")
      
      screen.DT = settings_ViewRef$elem.DT[seqnames == view_chr]
      screen.DT = screen.DT[start <= data_end & end >= data_start]
      
      tr_list = unique(screen.DT$transcript_id)
      message(paste(length(tr_list), " transcripts"))
      
      # limit transcript list
      req(length(tr_list) < 70)
      
      transcripts.DT = settings_ViewRef$transcripts.DT[transcript_id %in% tr_list]
      transcripts.DT = transcripts.DT
      setorder(transcripts.DT, transcript_support_level, width)

      transcripts.DT[, FOV_start := start]
      transcripts.DT[, FOV_end := end]   
      
      # apply plot_order on transcripts.DT
      transcripts.DT$plot_level = 0
      i = 0
      while(any(transcripts.DT$plot_level == 0)) {
        i = i + 1
        remaining = which(transcripts.DT$plot_level == 0)
        while(length(remaining) > 0) {
          transcripts.DT$plot_level[remaining[1]] = i
          ol.gr = GenomicRanges::reduce(
            GenomicRanges::makeGRangesFromDataFrame(as.data.frame(transcripts.DT[plot_level == i])),
            ignore.strand = TRUE
          )
          OL = GenomicRanges::findOverlaps(
            ol.gr,
            GenomicRanges::makeGRangesFromDataFrame(as.data.frame(transcripts.DT)),
            ignore.strand = TRUE
          )
          remaining = which(transcripts.DT$plot_level == 0)
          remaining = remaining[which(!(remaining %in% OL@to))]
        }
      }
      transcripts.DT[strand == "+", display_name := paste(transcript_name, "-", transcript_biotype, " ->>")]
      transcripts.DT[strand == "-", display_name := paste("<-- ", transcript_name, "-", transcript_biotype)]
      transcripts.DT[, disp_x := 0.5 * (start + end)]
      transcripts.DT[start < view_start, disp_x := 0.5 * (view_start + end)]
      transcripts.DT[end > view_end, disp_x := 0.5 * (start + view_end)]
      transcripts.DT[start < view_start & end > view_end, disp_x := 0.5 * (view_start + view_end)]
      transcripts.DT = transcripts.DT[end > view_start & start < view_end]

      plot.DT = settings_ViewRef$elem.DT[transcript_id %in% transcripts.DT$transcript_id]
      plot.DT$transcript_id = factor(plot.DT$transcript_id, unique(transcripts.DT$transcript_id), ordered = TRUE)
      plot.DT[transcripts.DT, on = "transcript_id", 
        c("transcript_name", "transcript_biotype", "transcript_support_level", "display_name", "plot_level") := 
        list(i.transcript_name, i.transcript_biotype, i.transcript_support_level, i.display_name, i.plot_level)]
      
      p = ggplot(plot.DT, aes(text = display_name))
      # text
      # p = p + geom_text(data = transcripts.DT, aes(
          # x = 0.5 * (max(start, view_start) + min(end, view_end)), y = plot_level - 0.4, label = display_name
        # )
      # )
      
      if(nrow(subset(as.data.frame(plot.DT), type = "intron")) > 0) {
        p = p + geom_segment(data = subset(as.data.frame(plot.DT), type = "intron"), 
          aes(x = start, xend = end, y = plot_level, yend = plot_level))
      }
      if(nrow(subset(as.data.frame(plot.DT), type != "intron")) > 0) {
        p = p + 
          geom_rect(data = subset(as.data.frame(plot.DT), type != "intron"), 
            aes(xmin = start, xmax = end, 
              ymin = plot_level - 0.1 - ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0), 
              ymax = plot_level + 0.1 + ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0)
            )
          )
      }
      # p = p + xlim(view_start, view_end)
      
      anno = list(
        x = transcripts.DT$disp_x,
        y = transcripts.DT$plot_level - 0.4,
        text = transcripts.DT$display_name,
        xref = "x", yref = "y", showarrow = FALSE)
        
      pl = ggplotly(p, tooltip = "text") %>% layout(
        annotations = anno,
        dragmode = "pan",
        xaxis = list(range = c(view_start, view_end)),
        yaxis = list(fixedrange = TRUE)
      )
      if(input$move_labels_view_ref == TRUE) {
        pl = pl %>% plotly::config(editable = TRUE)
      }
      output$plot_view_ref <- renderPlotly({
        print(pl)
      })
      settings_ViewRef$plot_ini = TRUE      
    }
    
    observeEvent(input$move_labels_view_ref, {
      settings_ViewRef$blockView = TRUE
      settings_ViewRef$blockView = FALSE      
    })

    settings_ViewRef$repan = reactive({
      req(settings_ViewRef$plot_ini == TRUE)
      event_data("plotly_relayout")
    })

    observeEvent(settings_ViewRef$repan(), {
      repan = settings_ViewRef$repan()
      message(names(repan))
      req(length(repan) == 2)
      req(all(c("xaxis.range[0]", "xaxis.range[1]") %in% names(repan)))
      settings_ViewRef$blockView = TRUE

      updateTextInput(session = session, inputId = "start_view_ref", 
        value = max(1, round(repan[["xaxis.range[0]"]])))
      updateTextInput(session = session, inputId = "end_view_ref", 
        value = round(repan[["xaxis.range[1]"]]))
        
      settings_ViewRef$blockView = FALSE
      
    })


		observeEvent(input$genes_view, {
      req(input$genes_view)
      req(input$genes_view != "(none)")

      gene_id_view = settings_ViewRef$gene_list[gene_display_name == input$genes_view]

      settings_ViewRef$blockView = TRUE

      # change settings on input based on this
      updateSelectInput(session = session, inputId = "chr_view_ref", 
        selected = gene_id_view$seqnames[1])
      updateTextInput(session = session, inputId = "start_view_ref", 
        value = gene_id_view$start[1])
      updateTextInput(session = session, inputId = "end_view_ref", 
        value = gene_id_view$end[1])
        
      settings_ViewRef$blockView = FALSE
		})

    observeEvent(settings_ViewRef$blockView, {
      req(settings_ViewRef$blockView == FALSE)
      plot_view_ref_fn()      
    })
    observeEvent(input$chr_view_ref, {
      req(settings_ViewRef$blockView == FALSE)
      req(input$chr_view_ref != "(none)")
      plot_view_ref_fn()
    })
    observeEvent(input$start_view_ref, {
      req(settings_ViewRef$blockView == FALSE)
      plot_view_ref_fn()
    })
    observeEvent(input$end_view_ref, {
      req(settings_ViewRef$blockView == FALSE)
      plot_view_ref_fn()
    })
    
# Design Experiment page
		settings_expr <- shiny::reactiveValues(
			expr_path = "",
			bam_path = "",
			irf_path = "",
			anno_file = "",
			collate_path = "",
			df = c(),
			df.files = c(),
			df.anno = c()
		)
		files_header = c("bam_file", "irf_file", "cov_file", "junc_file")

    # Make sure df.anno exists when df.files exist:
    observeEvent(settings_expr$df.files, {
      req(settings_expr$df.files)
      if(!is_valid(settings_expr$df.anno)) {
        DT = as.data.table(settings_expr$df.files)
        settings_expr$df.anno = DT[, "sample"]
      }
    })
    
		## Handsontable auto-updates settings_expr$df on user edit
    observeEvent(input$hot_files_expr,{
      req(input$hot_files_expr)
      settings_expr$df.files = hot_to_r(input$hot_files_expr) 
    })
    observeEvent(input$hot_anno_expr,{
      req(input$hot_anno_expr)
      settings_expr$df.anno = hot_to_r(input$hot_anno_expr)
    })
		output$hot_files_expr <- renderRHandsontable({
			if (!is.null(settings_expr$df.files)) {     
				rhandsontable(settings_expr$df.files, useTypes = TRUE, stretchH = "all")
			} else {
        NULL
      }
		})
		output$hot_anno_expr <- renderRHandsontable({
			if (!is.null(settings_expr$df.anno)) {     
				rhandsontable(settings_expr$df.anno, useTypes = TRUE, stretchH = "all")
			} else {
        NULL
      }
		})
		
		observeEvent(settings_expr$expr_path, {
      output$txt_reference_path_expr <- renderText({
				paste("Root folder:", settings_expr$expr_path)
			})
		})
    observe({  
      shinyDirChoose(input, "dir_bam_path_load", roots = c(default_volumes, addit_volume), session = session)
			output$txt_bam_path_expr <- renderText({
					validate(need(input$dir_bam_path_load, "Please select path where BAMs are kept"))
          settings_expr$expr_path = dirname(parseDirPath(c(default_volumes, addit_volume), input$dir_bam_path_load))
					settings_expr$bam_path = parseDirPath(c(default_volumes, addit_volume), input$dir_bam_path_load)
			})
    })
    Expr_Load_BAMs = function() {
		# First assume bams are named by subdirectory names
			temp.DT = FindSamples(settings_expr$bam_path, suffix = ".bam", use_subdir = TRUE)
			if(!is.null(temp.DT)) {
        temp.DT = as.data.table(temp.DT)
				if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
					# Assume subdirectory names designate sample names
				} else {
					temp.DT = as.data.table(FindSamples(
						settings_expr$bam_path, suffix = ".bam", use_subdir = FALSE))
					if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
				# Else assume bam names designate sample names					
					} else {
						output$txt_bam_path_expr <- renderText("BAM file names (or its path names) must be unique")							
						settings_expr$bam_path = ""
						temp.DT = NULL
					}
				}
			} else {
				output$txt_bam_path_expr <- renderText("No bam files found in given path")
				settings_expr$bam_path = ""
				temp.DT = NULL
			}
		# compile experiment df with bam paths
			if(!is.null(temp.DT)) {
				colnames(temp.DT)[2] = "bam_file"
				if(is_valid(settings_expr$df.files)) {
			# merge with existing dataframe	
					settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT)
				} else {
			# start anew
					DT = data.table(sample = temp.DT$sample,
						bam_file = "", irf_file = "", cov_file = "", junc_file = "")
					DT[temp.DT, on = "sample", bam_file := i.bam_file] # Update new bam paths
					settings_expr$df.files = as.data.frame(DT)
				}
			}    
    }
		observeEvent(settings_expr$bam_path,{
      req(settings_expr$bam_path)
      Expr_Load_BAMs()
		})
    
		# Run IRFinder
		observeEvent(input$run_irf_expr,{
			req(settings_expr$df.files)
			if(settings_loadref$loadref_path == "") {
				output$txt_run_irf_expr <- renderText("Please load reference")
			} else if(settings_expr$irf_path == "") {
				output$txt_run_irf_expr <- renderText("Please select IRFinder output path")
			} else if(!file.exists(file.path(settings_loadref$loadref_path, "IRFinder.ref.gz"))) {
				output$txt_run_irf_expr <- renderText("IRFinder.ref.gz not found in given reference path")
			} else {				
				df = settings_expr$df.files
				bam_to_run = unname(which(sapply(df$sample, is_valid) & sapply(df$bam_file, is_valid)))
				if("SnowParam" %in% class(BPPARAM)) {
					BPPARAM_mod = BiocParallel::SnowParam(input$expr_Cores)
          message(paste("Using SnowParam", input$expr_Cores, "cores"))
				} else if("MulticoreParam" %in% class(BPPARAM)) {
					BPPARAM_mod = BiocParallel::MulticoreParam(input$expr_Cores)
          message(paste("Using MulticoreParam", input$expr_Cores, "cores"))
				} else {
					BPPARAM_mod = BPPARAM
				}
				BiocParallel::bplapply(bam_to_run, function(i, run_IRF, df, reference_file, output_path) {
					run_IRF(df$bam_file[i], reference_file, file.path(output_path, df$sample[i]))
				}, df = df, run_IRF = run_IRFinder, reference_file = file.path(settings_loadref$loadref_path, "IRFinder.ref.gz"),
					output_path = settings_expr$irf_path, BPPARAM = BPPARAM_mod)
				Expr_Load_IRFs()
			}
		})

    shinyDirChoose(input, "dir_irf_path_load", roots = c(default_volumes, addit_volume), 
      session = session)		
    observeEvent(input$dir_irf_path_load, {
      req(input$dir_irf_path_load)
      settings_expr$irf_path = parseDirPath(c(default_volumes, addit_volume), input$dir_irf_path_load)
    })
    observeEvent(settings_expr$irf_path, {
      output$txt_irf_path_expr <- renderText({
        validate(need(settings_expr$irf_path, "Please select path where IRFinder output should be kept"))
        settings_expr$expr_path = dirname(settings_expr$irf_path)
        settings_expr$irf_path
      })    
    })
    Expr_Load_IRFs = function() {
				# merge irfinder paths
			temp.DT = FindSamples(settings_expr$irf_path, suffix = ".txt.gz", use_subdir = FALSE)
			if(!is.null(temp.DT)) {
        temp.DT = as.data.table(temp.DT)
				if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
					# Assume output names designate sample names
				} else {
					temp.DT = as.data.table(FindSamples(
						settings_expr$irf_path, suffix = ".txt.gz", use_subdir = TRUE))
					if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
				# Else assume subdirectory names designate sample names					
					} else {
						output$txt_irf_path_expr <- renderText("IRFinder file names (or its path names) must be unique")							
						settings_expr$irf_path = ""
						temp.DT = NULL
					}
				}
			} else {
				temp.DT = NULL
			}
			
		# compile experiment df with irfinder paths
			if(!is.null(temp.DT)) {
					colnames(temp.DT)[2] = "irf_file"
				if(is_valid(settings_expr$df.files)) {
			# merge with existing dataframe	
					settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT)
				} else {
			# start anew
					DT = data.table(sample = temp.DT$sample,
						bam_file = "", irf_file = "", cov_file = "", junc_file = "")
					DT[temp.DT, on = "sample", irf_file := i.irf_file] # Update new irf paths
					settings_expr$df.files = as.data.frame(DT)      
				}   
			}
			
			# Attempt to find Coverage files
			temp.DT = FindSamples(settings_expr$irf_path, suffix = ".cov", use_subdir = FALSE)
			if(!is.null(temp.DT)) {
        temp.DT = as.data.table(temp.DT)
				if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
					# Assume output names designate sample names
				} else {
					temp.DT = as.data.table(FindSamples(
						settings_expr$irf_path, suffix = ".cov", use_subdir = TRUE))
					if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
				# Else assume subdirectory names designate sample names					
					} else {
						temp.DT = NULL
					}
				}
			} else {
				temp.DT = NULL
			}
			
		# compile experiment df with irfinder paths
			if(!is.null(temp.DT)) {
				colnames(temp.DT)[2] = "cov_file"
			# merge with existing dataframe	
					settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT)
			}			
    }
		
		observeEvent(settings_expr$irf_path,{
      req(settings_expr$irf_path)
      Expr_Load_IRFs()
		})

		# Add annotation to data frame
    shinyFileChoose(input, "file_expr_path_load", roots = c(default_volumes, addit_volume), 
      session = session)
    observeEvent(input$file_expr_path_load, {
      req(input$file_expr_path_load)
      file_selected<-parseFilePaths(c(default_volumes, addit_volume), 
        input$file_expr_path_load)
      settings_expr$anno_file = as.character(file_selected$datapath)
    })

		observeEvent(settings_expr$anno_file,{
      output$txt_sample_anno_expr <- renderText({
        validate(need(settings_expr$anno_file, "Please select file where sample annotations are kept"))
        settings_expr$anno_file
      })    
      req(settings_expr$anno_file)
			temp.DT = tryCatch(as.data.table(fread(settings_expr$anno_file)),
				error = function(e) NULL)
			req(temp.DT)
			req(nrow(temp.DT) > 0)
			output$txt_sample_anno_expr <- renderText({
				validate(need("sample" %in% colnames(temp.DT), "'sample' must be the header of the first column"))
				""
			})
			req("sample" %in% colnames(temp.DT))
			
      anno_header = names(temp.DT)[!(names(temp.DT) %in% files_header)]
			temp.DT.files = copy(temp.DT)
      if(length(anno_header) > 0) temp.DT.files[, c(anno_header) := NULL]
			if(is_valid(settings_expr$df.files)) {
				settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT.files)
			} else {
				DT = data.table(sample = temp.DT$sample, bam_file = "", irf_file = "",
					cov_file = "", junc_file = "")
				settings_expr$df.files = update_data_frame(DT, temp.DT.files)
			}
			
			temp.DT.anno = copy(temp.DT)
      files_header_exist = intersect(files_header, names(temp.DT))
      if(length(files_header_exist) > 0) temp.DT.anno[, c(files_header_exist):= NULL]
			if(is_valid(settings_expr$df.anno)) {
				settings_expr$df.anno = update_data_frame(settings_expr$df.anno, temp.DT.anno)
			} else {
				settings_expr$df.anno = temp.DT.files
			}
		})

    observe({
      shinyDirChoose(input, "dir_collate_path_load", roots = c(default_volumes, addit_volume), 
        session = session)
      output$txt_collate_path_expr <- renderText({
          validate(need(input$dir_collate_path_load, "Please select path where NxtIRF compiled output should be kept"))
          settings_expr$expr_path = dirname(parseDirPath(c(default_volumes, addit_volume), input$dir_collate_path_load))
          settings_expr$collate_path = parseDirPath(c(default_volumes, addit_volume), 
            input$dir_collate_path_load)
      })      
    })
    Expr_Load_FSTs = function() {
				# merge nxtirf paths
			temp.DT = FindSamples(settings_expr$collate_path, suffix = ".junc.fst", use_subdir = FALSE)
			if(!is.null(temp.DT)) {
        temp.DT = as.data.table(temp.DT)
				if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
					# Assume output names designate sample names
				} else {
					temp.DT = as.data.table(FindSamples(
						settings_expr$collate_path, suffix = ".junc.fst", use_subdir = TRUE))
					if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
				# Else assume subdirectory names designate sample names					
					} else {
						output$txt_collate_path_expr <- renderText("NxtIRF FST file names (or its path names) must be unique")							
						settings_expr$collate_path = ""
						temp.DT = NULL
					}
				}
			} else {
				temp.DT = NULL
			}
			
		# compile experiment df with fst paths
			if(!is.null(temp.DT)) {
					colnames(temp.DT)[2] = "junc_file"
				if(is_valid(settings_expr$df.files)) {
			# merge with existing dataframe	
					settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT)
				} else {
			# start anew
					DT = data.table(sample = temp.DT$sample,
						bam_file = "", irf_file = "", cov_file = "", junc_file = "")
					DT[temp.DT, on = "sample", junc_file := i.junc_file] # Update new fst paths
					settings_expr$df.files = as.data.frame(DT)    
				}
			}    
    }
		observeEvent(settings_expr$collate_path,{
      req(settings_expr$collate_path)
      Expr_Load_FSTs()
		})
	
    output$newcol_expr <- renderUI({
      textInput("newcolumnname_expr", "New Column Name", sprintf("newcol%s", 1+ncol(settings_expr$df.anno)))
    })
 		# Add column
		observeEvent(input$addcolumn_expr, {
      df <- isolate(settings_expr$df.anno)
      newcolumn <- eval(parse(text=sprintf('%s(nrow(df))', isolate(input$type_newcol_expr))))
      settings_expr$df.anno <- setNames(cbind(df, newcolumn, stringsAsFactors=FALSE), 
				c(names(df), isolate(input$newcolumnname_expr)))
    })
 		# Remove column
		observeEvent(input$removecolumn_expr, {
      DT <- as.data.table(isolate(settings_expr$df.anno))
			if(isolate(input$newcolumnname_expr) %in% colnames(DT)) {
        message("removing column")
				DT[, c(input$newcolumnname_expr) := NULL]
				settings_expr$df.anno = as.data.frame(DT)
			}
    })	
    
    # Run CollateData()
    observeEvent(input$run_collate_expr, {
      req(settings_expr$df.files)

      Experiment = na.omit(as.data.table(settings_expr$df.files[, c("sample", "irf_file")]))
      reference_path = settings_loadref$loadref_path
      output_path = settings_expr$collate_path
      # BPPARAM = BPPARAM_mod

      output$txt_run_col_expr <- renderText({
        validate(need(settings_loadref$loadref_path, "Please load a reference before generating NxtIRF FST files"))
        validate(need(settings_expr$collate_path, "Please select a path to store NxtIRF FST files"))
        "running CollateData()"
      })

			cores_to_use = as.numeric(input$expr_Cores)
			if(!is_valid(input$expr_Cores)) cores_to_use = 1
			withProgress(message = 'Collating IRFinder output', value = 0, {
				CollateData(Experiment, reference_path, output_path, n_threads = cores_to_use)#, BPPARAM = BPPARAM)
			})
			Expr_Load_FSTs()
			output$txt_run_col_expr <- renderText("Finished compiling NxtIRF FST files")
    })
		
		shinyFileChoose(input, "loadexpr_expr", roots = c(default_volumes, addit_volume), session = session,
      filetypes = c("csv"))
    observeEvent(input$loadexpr_expr, {
      selectedfile <- parseFilePaths(c(default_volumes, addit_volume), input$loadexpr_expr)
      req(selectedfile$datapath)
      DT = fread(selectedfile$datapath, na.strings = c("", "NA"))
			DT.header = DT[sample == "(Experiment)"]
			DT = DT[sample != "(Experiment)"]

      if(all(c("sample", "bam_file", "irf_file", "cov_file", "junc_file") %in% names(DT)) & nrow(DT.header) == 1) {
        if(all(is.na(DT$bam_file))) DT[, bam_file:=as.character(bam_file)]
        if(all(is.na(DT$irf_file))) DT[, irf_file:=as.character(irf_file)]
        if(all(is.na(DT$cov_file))) DT[, cov_file:=as.character(cov_file)]
        if(all(is.na(DT$junc_file))) DT[, junc_file:=as.character(junc_file)]
        
        DT.files = copy(DT)
        anno_col = names(DT)[!(names(DT) %in% c("sample", files_header))]
        DT.files[, c(anno_col) := NULL]
        DT.anno = copy(DT)
        DT.anno[, c(files_header) := NULL]

        settings_expr$df.files = as.data.frame(DT.files)
        settings_expr$df.anno = as.data.frame(DT.anno)

        settings_expr$bam_path = DT.header$bam_file
        settings_expr$irf_path = DT.header$irf_file
        settings_expr$collate_path = DT.header$junc_file
        
        output$txt_run_save_expr <- renderText({
          paste(selectedfile$datapath, "loaded")
        })
      } else {
        output$txt_run_save_expr <- renderText({
          paste(selectedfile$datapath, "is not a valid NxtIRF experiment file")
        })
      }
		})		

		# Save Experiment
		shinyFileSave(input, "saveexpr_expr", roots = c(default_volumes, addit_volume), session = session,
      filetypes = c("csv"))
    observeEvent(input$saveexpr_expr, {
      req(settings_expr$df.files)
      req(settings_expr$df.anno)
      selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$saveexpr_expr)
      req(selectedfile$datapath)
      
      df = update_data_frame(settings_expr$df.files, settings_expr$df.anno)
      df = rbind(as.data.table(df), list("(Experiment)"), fill = TRUE)
      if(is_valid(settings_expr$bam_path)) df[sample == "(Experiment)", bam_file:=settings_expr$bam_path]
      if(is_valid(settings_expr$irf_file)) df[sample == "(Experiment)", irf_file:=settings_expr$irf_file]
      if(is_valid(settings_expr$collate_path)) df[sample == "(Experiment)", junc_file:=settings_expr$collate_path]
			fwrite(df, selectedfile$datapath)
			output$txt_run_save_expr <- renderText({
        paste("Experiment saved to:", selectedfile$datapath)
      })
		})		
    observeEvent(input$clear_expr, {
			settings_expr$expr_path = ""
			settings_expr$bam_path = ""
			settings_expr$irf_path = ""
			settings_expr$anno_file = ""
			settings_expr$collate_path = ""
			settings_expr$df.files = c()
			settings_expr$df.anno = c()
      output$txt_run_save_expr <- renderText("")
    })
		observeEvent(input$build_expr, {
			output$txt_run_save_expr <- renderText({
				validate(need(settings_expr$collate_path, "Please set path to FST main files first"))
        colData = as.data.table(settings_expr$df.anno)
        # colData = colData[, -c("bam_file", "irf_file", "cov_file", "junc_file")]
        validate(need(ncol(colData) > 1, "Please assign at least 1 column of annotation to the experiment first"))
				se = MakeSE(colData, settings_expr$collate_path)
				settings_SE$se = se
				"SummarizedExperiment Loaded"
			})
		})
		
	# Analyse - Calculate PSIs
		settings_SE <- shiny::reactiveValues(
			se = NULL,
			
			filterSummary = NULL,
      filters = list()
		)
		
    filter1 <- filterModule_server("filter1", conditionList)
    filter2 <- filterModule_server("filter2", conditionList)
    filter3 <- filterModule_server("filter3", conditionList)
    filter4 <- filterModule_server("filter4", conditionList)
    filter5 <- filterModule_server("filter5", conditionList)
    filter6 <- filterModule_server("filter6", conditionList)
    filter7 <- filterModule_server("filter7", conditionList)
    filter8 <- filterModule_server("filter8", conditionList)

    observe({
      settings_SE$filters[[1]] = (reactiveValuesToList(filter1))
      settings_SE$filters[[2]] = (reactiveValuesToList(filter2))
      settings_SE$filters[[3]] = (reactiveValuesToList(filter3))
      settings_SE$filters[[4]] = (reactiveValuesToList(filter4))
      settings_SE$filters[[5]] = (reactiveValuesToList(filter5))
      settings_SE$filters[[6]] = (reactiveValuesToList(filter6))
      settings_SE$filters[[7]] = (reactiveValuesToList(filter7))
      settings_SE$filters[[8]] = (reactiveValuesToList(filter8))
    })
    
    conditionList = reactive({
      if(is(settings_SE$se, "SummarizedExperiment")) {
        colnames(SummarizedExperiment::colData(settings_SE$se))
      } else {
        c("")
      }
    })
    
    se.filterModule <- reactive({
      if(is(settings_SE$se, "SummarizedExperiment")) {
        settings_SE$se
      } else {
        NULL
      }
    })

    processFilters <- function() {
      message("Refreshing filters")
      if(is(settings_SE$se, "SummarizedExperiment")) {
        filterSummary = rep(TRUE, nrow(settings_SE$se))
        for(i in 1:8) {
          print(settings_SE$filters[[i]]$filterVars)
          print(settings_SE$filters[[i]]$trigger)
            if(!is.null(settings_SE$filters[[i]]$trigger)) {
              filterSummary = filterSummary & runFilter(
                settings_SE$filters[[i]]$filterClass,
                settings_SE$filters[[i]]$filterType,
                settings_SE$filters[[i]]$filterVars,
                settings_SE$se)
            } else {
              message(paste("Trigger", i, "is NULL"))
            }
        }
        settings_SE$filterSummary = filterSummary    
      }
    }
    
    observeEvent({
      settings_SE$se
      settings_SE$filters[[1]]$trigger
      settings_SE$filters[[2]]$trigger
      settings_SE$filters[[3]]$trigger
      settings_SE$filters[[4]]$trigger
      settings_SE$filters[[5]]$trigger
      settings_SE$filters[[6]]$trigger
      settings_SE$filters[[7]]$trigger
      settings_SE$filters[[8]]$trigger
    }, {
      processFilters()
    })
    
    observeEvent(input$refresh_filters_Filters, {
      req(input$refresh_filters_Filters)
      processFilters()
    })
    
    observeEvent(settings_SE$filterSummary, {
      req(settings_SE$filterSummary)

      if(is(settings_SE$se, "SummarizedExperiment")) {
        filteredEvents.DT = data.table(EventType = SummarizedExperiment::rowData(settings_SE$se)$EventType,
          keep = settings_SE$filterSummary)
        if(input$graphscale_Filters == "log10") {
          filteredEvents.DT[, Included := log10(sum(keep == TRUE)), by = "EventType"]
          filteredEvents.DT[, Excluded := log10(sum(!is.na(keep))) - log10(sum(keep == TRUE)) , by = "EventType"]
        } else {
          filteredEvents.DT[, Included := sum(keep == TRUE), by = "EventType"]
          filteredEvents.DT[, Excluded := sum(keep != TRUE) , by = "EventType"]        
        }
        filteredEvents.DT = unique(filteredEvents.DT, by = "EventType")
        incl = as.data.frame(filteredEvents.DT[, c("EventType", "Included")]) %>%
          dplyr::mutate(filtered = "Included") %>% dplyr::rename(Events = Included)
        excl = as.data.frame(filteredEvents.DT[, c("EventType", "Excluded")]) %>%
          dplyr::mutate(filtered = "Excluded") %>% dplyr::rename(Events = Excluded)
        # ggplot summary as bar plot
        
        p = ggplot(rbind(incl, excl), aes(x = EventType, y = Events, fill = filtered)) +
          geom_bar(position="stack", stat="identity")
        if(input$graphscale_Filters == "log10") {
          p = p + labs(y = "log10 Events")
        } else {
          p = p + labs(y = "Events")        
        }
        output$plot_filtered_Events <- renderPlotly({
          print(
            ggplotly(p)
          )
        })
      }
    })
    
    # DE
		settings_DE <- shiny::reactiveValues(
			res = NULL,
			res_settings = list(),
			method = NULL,
			batchVar1 = NULL,
			batchVar2 = NULL,
			DE_Var = NULL,
			nom_DE = NULL,
			denom_DE = NULL
		)
    
    observeEvent(input$variable_DE, {
			req(input$variable_DE)
			req(input$variable_DE != "(none)")
			
			colData = SummarizedExperiment::colData(settings_SE$se)
			
			if(class(colData[,input$variable_DE]) != "factor") {
				output$warning_DE = renderText("Contrast must be performed on discrete categories")
        updateSelectInput(session = session, inputId = "variable_DE", 
          choices = c("(none)", colnames(colData)), selected = "(none)")
			} else {
				updateSelectInput(session = session, inputId = "nom_DE", 
          choices = c("(none)", levels(colData[,input$variable_DE])), selected = "(none)")					
        updateSelectInput(session = session, inputId = "denom_DE", 
          choices = c("(none)", levels(colData[,input$variable_DE])), selected = "(none)")
				if(is_valid(settings_DE$nom_DE) && settings_DE$nom_DE %in% levels(colData[,input$variable_DE])) {
					updateSelectInput(session = session, inputId = "denom_DE", selected = settings_DE$nom_DE)
				}
				if(is_valid(settings_DE$nom_DE) && settings_DE$denom_DE %in% levels(colData[,input$variable_DE])) {
					updateSelectInput(session = session, inputId = "denom_DE", selected = settings_DE$denom_DE)
				}
			}
			
			
		})

		observeEvent(settings_DE$method, {
			req(settings_DE$method)
      updateSelectInput(session = session, inputId = "method_DE", selected = settings_DE$method)		
		})    
		observeEvent(settings_DE$DE_Var, {
			req(settings_DE$DE_Var)
      updateSelectInput(session = session, inputId = "variable_DE", selected = settings_DE$DE_Var)		
		})
		observeEvent(settings_DE$batchVar1, {
			req(settings_DE$batchVar1)
			updateSelectInput(session = session, inputId = "batch1_DE", selected = settings_DE$batchVar1)		
		})
		observeEvent(settings_DE$batchVar2, {
			req(settings_DE$batchVar2)
			updateSelectInput(session = session, inputId = "batch2_DE", selected = settings_DE$batchVar2)	
		})
		
    observeEvent(input$perform_DE, {
			req(settings_SE$se)
			output$warning_DE = renderText({
				validate(need(input$variable_DE != "(none)"), "Variable for DE needs to be defined")
				validate(need(input$nom_DE != "(none)"), "Nominator for DE Variable needs to be defined")
				validate(need(input$denom_DE != "(none)"), "Denominator for DE Variable needs to be defined")
				validate(need(input$denom_DE != input$nom_DE), "Denominator and Nominator must be different")
        paste("Running", input$method_DE)
			})
			req(input$variable_DE)
			req(input$nom_DE)
			req(input$denom_DE)
			req(input$variable_DE != "(none)" & input$nom_DE != "(none)" & input$denom_DE != "(none)")
			
			if(length(settings_SE$filterSummary) == nrow(settings_SE$se)) {
				se = settings_SE$se[settings_SE$filterSummary,]
			} else {
				se = settings_SE$se
			}
			rowData = as.data.frame(SummarizedExperiment::rowData(se))
			colData = as.data.frame(SummarizedExperiment::colData(se))
			
			settings_DE$DE_Var = input$variable_DE
			settings_DE$nom_DE = input$nom_DE
			settings_DE$denom_DE = input$denom_DE
			
			if(input$batch1_DE != "(none)" & input$batch1_DE != input$variable_DE) {
				settings_DE$batchVar1 = input$batch1_DE
			} else {
				settings_DE$batchVar1 = NULL
        updateSelectInput(session = session, inputId = "batch1_DE", 
          selected = "(none)")
			}
			if(input$batch2_DE != "(none)" & input$batch2_DE != input$variable_DE & 
					input$batch2_DE != input$batch1_DE) {
				settings_DE$batchVar2 = input$batch2_DE
			} else {
				settings_DE$batchVar2 = NULL
        updateSelectInput(session = session, inputId = "batch2_DE", 
          selected = "(none)")
			}
			req(input$method_DE)
			settings_DE$method = input$method_DE
			
			if(settings_DE$method == "DESeq2") {
				NxtIRF.CheckPackageInstalled("DESeq2", "1.28.0")
				# build design
				if(!is_valid(settings_DE$batchVar2)) {
					dds_formula = paste0("~", paste(
						settings_DE$batchVar1, settings_DE$batchVar2, settings_DE$DE_Var,
						paste0(settings_DE$DE_Var, ":ASE"),
						sep="+"))
				} else if(!is_valid(settings_DE$batchVar1)) {
					dds_formula = paste0("~", paste(
						settings_DE$batchVar1, settings_DE$DE_Var,
						paste0(settings_DE$DE_Var, ":ASE"),
						sep="+"))				
				} else {
					dds_formula = paste0("~", paste(settings_DE$DE_Var,
						paste0(settings_DE$DE_Var, ":ASE"),
						sep="+"))						
				}
				
				# construct dds
				countData = cbind(SummarizedExperiment::assay(se, "Included"), 
					SummarizedExperiment::assay(se, "Excluded"))
				colData_use = rbind(colData, colData)
				rownames(colData_use) = c(
					paste(rownames(colData), "Included", sep="."),
					paste(rownames(colData), "Excluded", sep=".")
				)
        colData_use$ASE = rep(c("Included", "Excluded"), each = nrow(colData))
				colnames(countData) = rownames(colData_use)
				rownames(countData) = rowData$EventName
				countData = round(countData)
				mode(countData) = "integer"

				dds = DESeq2::DESeqDataSetFromMatrix(
					countData = countData,
					colData = colData_use,
					design = as.formula(dds_formula)
				)
				
				DESeq2::sizeFactors(dds) = 1

        cores_to_use = as.numeric(input$expr_Cores)
        if(!is_valid(input$expr_Cores)) cores_to_use = 1
        
        BPPARAM = BiocParallel::bpparam()
        
        if("SnowParam" %in% class(BPPARAM)) {
          BPPARAM_mod = BiocParallel::SnowParam(cores_to_use)
          message(paste("Using SnowParam", BPPARAM_mod$workers, "cores"))
        } else if("MulticoreParam" %in% class(BPPARAM)) {
          BPPARAM_mod = BiocParallel::MulticoreParam(cores_to_use)
          message(paste("Using MulticoreParam", BPPARAM_mod$workers, "cores"))
        } else {
          BPPARAM_mod = BiocParallel::SerialParam()
          message(paste("Using SerialParam mode with", BPPARAM_mod$workers, "cores"))
        }
				
				dds = DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM_mod)
				
        # print(DESeq2::resultsNames(dds))
        message(paste("Factors to contrast are", paste0(settings_DE$DE_Var, settings_DE$nom_DE, ".ASEIncluded"),
          paste0(settings_DE$DE_Var, settings_DE$denom_DE, ".ASEIncluded")))
        
				res = as.data.frame(DESeq2::results(dds,
					list(
						paste0(settings_DE$DE_Var, settings_DE$nom_DE, ".ASEIncluded"),
						paste0(settings_DE$DE_Var, settings_DE$denom_DE, ".ASEIncluded")
					), parallel = TRUE, BPPARAM = BPPARAM_mod)
				)
				res = cbind(as.data.frame(rowData), res)
				settings_DE$res = res

        output$warning_DE = renderText({"Finished"})
        
			} else if(settings_DE$method == "limma") {
				NxtIRF.CheckPackageInstalled("limma", "3.44.0")
        
				countData = cbind(SummarizedExperiment::assay(se, "Included"), 
					SummarizedExperiment::assay(se, "Excluded"))
				colData_use = rbind(colData, colData)
				rownames(colData_use) = c(
					paste(rownames(colData), "Included", sep="."),
					paste(rownames(colData), "Excluded", sep=".")
				)
        colData_use$ASE = rep(c("Included", "Excluded"), each = nrow(colData))
				colnames(countData) = rownames(colData_use)
				rownames(countData) = rowData$EventName
				# countData = round(countData)
				# mode(countData) = "integer"
        
        res.ASE = limma_DT(countData, colData_use, settings_DE$DE_Var, 
          settings_DE$nom_DE, settings_DE$denom_DE, useASE = TRUE)
        res.inc = limma_DT(countData[, seq_len(nrow(colData))], 
          colData_use[seq_len(nrow(colData)),], settings_DE$DE_Var, settings_DE$nom_DE, settings_DE$denom_DE)
        res.exc = limma_DT(countData[, seq(nrow(colData) + 1, nrow(colData) * 2)], 
          colData_use[seq(nrow(colData) + 1, nrow(colData) * 2),], 
            settings_DE$DE_Var, settings_DE$nom_DE, settings_DE$denom_DE)
       
        res.ASE[res.inc, on = "EventName",
          paste("Inc", colnames(res.inc)[1:6], sep=".") := list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
          
        res.ASE[res.exc, on = "EventName",
          paste("Exc", colnames(res.exc)[1:6], sep=".") := list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
          
        setorder(res.ASE, -B)

				settings_DE$res = as.data.frame(res.ASE)

        output$warning_DE = renderText({"Finished"})
      }
			
			req(settings_DE$res)
			# save settings for current res
			settings_DE$res_settings$method = settings_DE$method
			settings_DE$res_settings$DE_Var = settings_DE$DE_Var
			settings_DE$res_settings$nom_DE = settings_DE$nom_DE
			settings_DE$res_settings$denom_DE = settings_DE$denom_DE
			settings_DE$res_settings$batchVar1 = settings_DE$batchVar1
			settings_DE$res_settings$batchVar2 = settings_DE$batchVar2
		})
		
		observeEvent(settings_DE$res, {
			req(settings_DE$res)
			output$DT_DE <- DT::renderDataTable(
				DT::datatable(
					settings_DE$res,
					class = 'cell-border stripe',
					rownames = FALSE,
					filter = 'top'
				)
			)		
		})
		
		shinyFileSave(input, "save_DE", roots = c(default_volumes, addit_volume), session = session,
      filetypes = c("Rds"))
		observeEvent(input$save_DE, {	
			req(settings_DE$res)
			req(length(settings_DE$res_settings) > 0)
      
			selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$save_DE)
			req(selectedfile$datapath)
			
			save_DE = list(res = settings_DE$res, settings = settings_DE$res_settings)
			saveRDS(save_DE,selectedfile$datapath)
		})

    shinyFileChoose(input, "load_DE", roots = c(default_volumes, addit_volume), 
      session = session, filetype = list(RDS = "Rds"))
    observeEvent(input$load_DE, {
      req(input$load_DE)
      file_selected<-parseFilePaths(c(default_volumes, addit_volume), 
        input$load_DE)
      
			load_DE = readRDS(file_selected)
			req(all(c("res", "settings") %in% names(load_DE)))
			
		# check all parameters exist in colData(se)
			req(settings_SE$se)
			colData = SummarizedExperiment::colData(settings_SE$se)
			req(load_DE$settings$DE_Var %in% colnames(colData))
			req(!is_valid(load_DE$settings$batchVar1) || load_DE$settings$batchVar1 %in% colnames(colData))
			req(!is_valid(load_DE$settings$batchVar2) || load_DE$settings$batchVar2 %in% colnames(colData))
			req(any(unlist(colData[,load_DE$settings$DE_Var]) == load_DE$settings$nom_DE))
			req(any(unlist(colData[,load_DE$settings$DE_Var]) == load_DE$settings$denom_DE))
			req(load$settings$method %in% c("DESeq2", "limma", "DSS"))
			
			settings_DE$res = load_DE$res
			settings_DE$res_settings$method = load_DE$settings$method
			settings_DE$res_settings$DE_Var = load_DE$settings$DE_Var
			settings_DE$res_settings$nom_DE = load_DE$settings$nom_DE
			settings_DE$res_settings$denom_DE = load_DE$settings$denom_DE
			settings_DE$res_settings$batchVar1 = load_DE$settings$batchVar1
			settings_DE$res_settings$batchVar2 = load_DE$settings$batchVar2
			
			settings_DE$method = settings_DE$res_settings$method
			settings_DE$DE_Var = settings_DE$res_settings$DE_Var
			settings_DE$nom_DE = settings_DE$res_settings$nom_DE
			settings_DE$denom_DE = settings_DE$res_settings$denom_DE
			settings_DE$batchVar1 = settings_DE$res_settings$batchVar1
			settings_DE$batchVar2 = settings_DE$res_settings$batchVar2			
    })
		
		
	make_matrix <- function(se, event_list, sample_list, method, depth_threshold = 10, logit_max = 5) {

    inc = SummarizedExperiment::assay(se, "Included")[event_list, sample_list]
    exc = SummarizedExperiment::assay(se, "Excluded")[event_list, sample_list]
		
		if(method == "PSI") {
			# essentially M/Cov
			mat = inc/(inc + exc)
			mat[inc + exc < depth_threshold] = NA
			return(mat)
		} else if(method == "logit") {
			mat = inc/(inc + exc)
			mat[inc + exc < depth_threshold] = NA
      mat = boot::logit(mat)
      mat[mat > logit_max] = logit_max
      mat[mat < -logit_max] = -logit_max
      return(mat)
    }
		
	}
  
	make_diagonal <- function(se, event_list, condition, nom_DE, denom_DE, depth_threshold = 10, logit_max = 5) {

    inc = SummarizedExperiment::assay(se, "Included")[event_list, ]
    exc = SummarizedExperiment::assay(se, "Excluded")[event_list, ]
    mat = inc/(inc + exc)
    mat[inc + exc < depth_threshold] = NA

    # use logit method to calculate geometric mean

    mat.nom = boot::logit(mat[, SummarizedExperiment::colData(se)[,condition] == nom_DE])
    mat.denom = boot::logit(mat[, SummarizedExperiment::colData(se)[,condition] == denom_DE])
    
    mat.nom[mat.nom > logit_max] = logit_max
    mat.denom[mat.denom > logit_max] = logit_max
    mat.nom[mat.nom < -logit_max] = -logit_max
    mat.denom[mat.denom < -logit_max] = -logit_max
        
    df = data.frame(EventName = event_list, nom = boot::inv.logit(rowMeans(mat.nom, na.rm = TRUE)),
      denom = boot::inv.logit(rowMeans(mat.denom, na.rm = TRUE)))
		
    return(df)
	}
	
  
# End of server function		
  }

  runApp(shinyApp(ui, server))

}
update_data_frame <- function(existing_df, new_df) {
	# add extra samples to existing df
	DT1 = as.data.table(existing_df)
	DT2 = as.data.table(new_df)
	md1 = melt(DT1, id = "sample")
	md2 = melt(DT2, id = "sample")
	
	res = unique(rbind(md1, md2), by = c("sample", "variable"), fromLast = TRUE)
	return(as.data.frame(dcast(res, sample ~ ...)))
}