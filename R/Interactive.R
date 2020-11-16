

#' @export
startNxtIRF <- function(offline = FALSE, BPPARAM = BiocParallel::bpparam()) {

	assert_that(interactive(),
		msg = "NxtIRF App can only be run in interactive mode (i.e. RStudio).")

  ah = AnnotationHub(localHub = offline)

	filterModule_UI <- function(id, label = "Counter") {
		ns <- NS(id)
		wellPanel(
			h5(label),	# e.g. "Filter #1"
			selectInput(ns("filterClass"), "Filter Class", width = '100%', choices = c("(none)", "Annotation", "Data")),
			selectInput(ns("filterType"), "Filter Type", width = '100%', choices = c("(none)")),
      conditionalPanel(ns = ns,
        condition = "['Transcript_Support_Level'].indexOf(input.filterType) >= 0",
        shinyWidgets::sliderTextInput(ns("slider_TSL_min"), "TSL Threshold", 
          choices = 1:5, selected = 1)
      ),
      conditionalPanel(ns = ns,
        condition = "['Consistency'].indexOf(input.filterType) >= 0",
        shinyWidgets::sliderTextInput(ns("slider_cons_max"), "log-fold maximum", choices = seq(0.2, 5, by = 0.2), selected = 1)
      ),
      conditionalPanel(ns = ns,
        condition = "['Coverage'].indexOf(input.filterType) >= 0",
        sliderInput(ns("slider_cov_min"), "Percent Coverage", min = 0, max = 100, value = 80)
      ),
      conditionalPanel(ns = ns,
        condition = "['Depth'].indexOf(input.filterType) >= 0",
        shinyWidgets::sliderTextInput(ns("slider_depth_min"), 
					"Minimum", choices = c(1,2,3,5,10,20,30,50,100,200,300,500), selected = 20),
      ),
      conditionalPanel(ns = ns,
        condition = "['Depth', 'Coverage'].indexOf(input.filterType) >= 0",
        tagList(
          shinyWidgets::sliderTextInput(ns("slider_mincond"), "Minimum Conditions Satisfy Criteria", 
						choices = c(as.character(1:8), "All"), selected = "All"),
          selectInput(ns("select_conds"), "Condition", width = '100%',
            choices = c("(none)")),
          sliderInput(ns("slider_pcTRUE"), "Percent samples per condition satisfying criteria", min = 0, max = 100, value = 80)
        )
      ),
      conditionalPanel(ns = ns,
        condition = "['Coverage', 'Consistency'].indexOf(input.filterType) >= 0",
        shinyWidgets::sliderTextInput(ns("slider_minDepth"), "Signal Threshold to apply criteria", 
          choices = c(1,2,3,5,10,20,30,50,100,200,300,500), selected = 20),
      ),
      conditionalPanel(ns = ns,
        condition = "['(none)'].indexOf(input.filterClass) < 0",
        selectInput(ns("EventType"), "Splice Type", width = '100%', multiple = TRUE,
          choices = c("IR", "MXE", "SE", "AFE", "ALE", "A5SS", "A3SS"))
      )
		)
	}

	filterModule_server <- function(id, filterdata, conditionList) {
		moduleServer(
			id,
			function(input, output, session) {
				final <- reactiveValues(
          filterClass = "",
          filterType = "",
          filterVars = list()
        )

				# inputs from final -> UI			
				observeEvent(filterdata(), {
					final = filterdata()

          if(is_valid(final$filterClass)) {
            if(final$filterClass == "Annotation") {
              type_choices = c("Protein_Coding", "NMD_Switching", "Transcript_Support_Level")              
            } else if(final$filterClass == "Data") {
              type_choices = c("Depth", "Coverage", "Consistency")
            } else {
              type_choices = c("(none)")
            }
          } else {
            type_choices = c("(none)")
          }
          if(is_valid(final$filterType) && final$filterType %in% type_choices) {
            updateSelectInput(session = session, inputId = "filterType", 
              choices = type_choices, selected = final$filterType)
          } else {
            updateSelectInput(session = session, inputId = "filterType", 
              choices = type_choices)              
          }
          if(is_valid((final$filterVars$minimum))) {
            if(final$filterType == "Depth") {
              shinyWidgets::updateSliderTextInput(session = session, inputId = "slider_depth_min", 
                selected = final$filterVars$minimum)
            } else  if(final$filterType == "Coverage"){
              updateSliderInput(session = session, inputId = "slider_cov_min", 
                value = final$filterVars$minimum)							
            } else  if(final$filterType == "Transcript_Support_Level"){
              shinyWidgets::updateSliderTextInput(session = session, inputId = "slider_TSL_min", 
                selected = final$filterVars$minimum)							
            }
          }
          if(is_valid(final$filterVars$maximum)) {
            shinyWidgets::updateSliderTextInput(session = session, inputId = "slider_cons_max", 
              selected = final$filterVars$maximum)
          }
          if(is_valid(final$filterVars$minDepth)) {
            updateSelectInput(session = session, inputId = "slider_minDepth", 
              selected = final$filterVars$minDepth)
          }
          if(is_valid(final$filterVars$minCond)) {
            shinyWidgets::updateSliderTextInput(session = session, inputId = "slider_mincond", 
              selected = final$filterVars$minCond)
          }
          if(is_valid(final$filterVars$condition)) {
            choices_conds = c("(none)", conditionList())
            if(final$filterVars$condition %in% choices_conds) {
              updateSelectInput(session = session, inputId = "select_conds", 
                choices = choices_conds, selected = final$filterVars$condition)
            }
          }
          if(is_valid(final$filterVars$pcTRUE)){
            updateSliderInput(session = session, inputId = "slider_pcTRUE", 
              value = final$filterVars$pcTRUE)
          }
          if(is_valid(final$filterVars$EventTypes)) {
            updateSelectInput(session = session, inputId = "EventType", 
              selected = final$filterVars$EventTypes)
          } else {
            updateSelectInput(session = session, inputId = "EventType", 
              selected = NULL)          
          }
        })
        
        # outputs from UI -> final
				observeEvent(input$filterClass, {
          final$filterClass = input$filterClass
          if(final$filterClass == "Annotation") {
            type_choices = c("Protein_Coding", "NMD_Switching", "Transcript_Support_Level")
          } else if(final$filterClass == "Data") {
            type_choices = c("Depth", "Coverage", "Consistency")
          } else {
            type_choices = "(none)"
          }
          cur_choice = isolate(input$filterType)
          if(is_valid(cur_choice) && cur_choice %in% type_choices) {
            updateSelectInput(session = session, inputId = "filterType", 
              choices = type_choices, selected = cur_choice)                    
          } else {
            updateSelectInput(session = session, inputId = "filterType", 
              choices = type_choices)          
          }
        })
        observeEvent(input$filterType, {
          final$trigger = NULL
          req(is_valid(input$filterType))
          final$filterType = input$filterType
          final$trigger = runif(1)
          
          if(input$filterType == "Depth") {
            final$filterVars$minimum = input$slider_depth_min
          }
          if(input$filterType == "Coverage"){
						final$filterVars$minimum = input$slider_cov_min
          }
          if(input$filterType == "Transcript_Support_Level"){
            final$filterVars$minimum = as.numeric(input$slider_TSL_min)
          }
        })
        observeEvent(input$slider_depth_min, {
          if(final$filterType == "Depth") {
            final$filterVars$minimum = input$slider_depth_min
          }
        })
        observeEvent(input$slider_cov_min, {
          if(final$filterType == "Coverage"){
						final$filterVars$minimum = input$slider_cov_min
          }
        })
        observeEvent(input$slider_TSL_min,{        
          if(final$filterType == "Transcript_Support_Level"){
            final$filterVars$minimum = as.numeric(input$slider_TSL_min)
          }
        })
        observeEvent(input$slider_cons_max,{        
          final$filterVars$maximum = input$slider_cons_max
        })
        observeEvent(input$slider_minDepth,{        
          final$filterVars$minDepth = input$slider_minDepth
        })
        observeEvent(input$slider_mincond,{        
          final$filterVars$minCond = input$slider_mincond
        })
        observeEvent(input$select_conds,{        
          final$filterVars$condition = input$select_conds
        })
        observeEvent(input$slider_pcTRUE,{        
          final$filterVars$pcTRUE = input$slider_pcTRUE
        })
        observeEvent(input$EventType,{        
          final$filterVars$EventTypes = input$EventType
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
        tags$div(title = paste("Number of threads to run computationally-intensive operations",
          "such as IRFinder, NxtIRF-collate, and DESeq2"),
          numericInput("cores_numeric", "# Threads", min = 1, max = 1, value = 1)
        ),
        tags$div(title = paste("Number of threads to run computationally-intensive operations",
          "such as IRFinder, NxtIRF-collate, and DESeq2"),
          sliderInput('cores_slider', "# Threads", min = 1, max = 1, value = 1)			
        )
			),
		),		
# Reference
		navbarMenu("Reference",
		# New reference
			tabPanel("New", value = "navRef_New",
				fluidRow(
					column(5,
						h4("Select Reference Directory"),
            tags$div(title = "Specify (or create) a directory for NxtIRF to create its IRFinder/NxtIRF reference",
              shinyDirButton("dir_reference_path", label = "Choose reference path", 
                title = "Choose reference path")
            ),
						textOutput("txt_reference_path"),
						br(),
						h4("Select Ensembl Reference from AnnotationHub"),
						selectInput('newrefAH_Species', 'Select Species', width = '100%',
							choices = c("")),
						selectInput('newrefAH_Version_Trans', 'Select Transcriptome Version', width = '100%',
							choices = c("")),
            tags$div(title = paste("Choose the source gtf file to build the reference.",
                "Typically avoid choosing *.abinitio.gtf, '*.chr.gtf' or",
                "'*.chr_patch_hapl_scaff.gtf' assemblies"),
              selectInput('newrefAH_Trans', 'Select Transcriptome Reference', width = '100%',
                choices = c(""))
            ),
						selectInput('newrefAH_Assembly', 'Select Genome Assembly', width = '100%',
							choices = c("")),
						selectInput('newrefAH_Version_Genome', 'Select Genome Version', width = '100%',
							choices = c("")),
            tags$div(title = paste("Choose the source genome sequence to build the reference.",
                "Typically choose the primary assembly, or if none, use",
                "either the 'sm.toplevel' or 'rm.toplevel' assemblies.",
                "Avoid using the cdna or ncrna assemblies"),
              selectInput('newrefAH_Genome', 'Select Genome Reference', width = '100%',
                choices = c(""))
            ),
						br(),
						h4("or select Reference from File"),
						br(),
            tags$div(title = paste("Choose a user-supplied genome fasta file"),
              shinyFilesButton("file_genome", label = "Choose genome FASTA File", title = "Choose genome FASTA File", multiple = FALSE)
						),
            textOutput("txt_genome"),
            tags$div(title = paste("Choose a user-supplied transcript reference gtf file"),
              shinyFilesButton("file_gtf", label = "Choose transcriptome GTF File", title = "Choose transcriptome GTF File", multiple = FALSE)
            ),
						textOutput("txt_gtf")
					),
					column(5,
            tags$div(title = paste("NxtIRF will auto-populate default mappability and non-polyA",
                "reference files for hg38, hg19, mm10 and mm9 genomes"),
              selectInput('newref_genome_type',
                'Select Genome Type to set Mappability and non-PolyA files', 
                c("(not specified)", "hg38", "mm10", "hg19", "mm9"))
            ),
            tags$div(title = paste("Select Mappability Exclusion file. This is typically a 3 columns",
                "of values containing seqnames, start and end coordinates of low-mappability regions"),
              shinyFilesButton("file_mappa", label = "Choose Mappability Exclusion file", 
                title = "Choose Mappability Exclusion file", multiple = FALSE)
            ),
						textOutput("txt_mappa"), actionButton("clear_mappa", "Clear"),
            tags$div(title = paste("Select Non-PolyA reference file. This is used by IRFinder",
                "to calculate reads from known non-polyadenylated transcripts to assess",
                "quality of poly-A enrichment in sample QC"),
              shinyFilesButton("file_NPA", label = "Choose non-PolyA BED file", title = "Choose non-PolyA BED file", multiple = FALSE)
						),
            textOutput("txt_NPA"), actionButton("clear_NPA", "Clear"),
            tags$div(title = paste("Select Blacklist file. This is typically a 3 columns",
                "of values containing seqnames, start and end coordinates of regions",
                "to exclude from IRFinder analysis"),
              shinyFilesButton("file_bl", label = "Choose blacklist BED file", title = "Choose blacklist BED file", multiple = FALSE)
            ),
						textOutput("txt_bl"), actionButton("clear_bl", "Clear"),
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
						div(style="display: inline-block;vertical-align:top; width: 240px;",
              selectizeInput('genes_view', 'Genes', choices = "(none)")
            ),
						div(style="display: inline-block;vertical-align:top; width: 240px;",
              selectizeInput('events_view', 'Events', choices = c("(none)"))
            ),
            textOutput("warning_view_ref")
          ),
					column(8,
						div(style="display: inline-block;vertical-align:top; width: 80px;",
							selectInput("chr_view_ref", label = "Chr", c("(none)"), selected = "(none)")),
						div(style="display: inline-block;vertical-align:top; width: 120px;",
							textInput("start_view_ref", label = "Left", c(""))),
						div(style="display: inline-block;vertical-align:top; width: 120px;",
							textInput("end_view_ref", label = "Right", c(""))),            
            br(),
            shinyWidgets::actionBttn("zoom_out_view_ref", style = "material-circle", color = "danger",
              icon = icon("minus")),
						div(style="display: inline-block;vertical-align:center;width: 60px;padding:20px",
							textOutput("label_zoom_view_ref")),
            shinyWidgets::actionBttn("zoom_in_view_ref", style = "material-circle", color = "danger",
              icon = icon("plus")),
            div(style="display: inline-block;vertical-align:top;padding:10px;",
              shinyWidgets::switchInput("move_labels_view_ref", label = "Movable labels", labelWidth = "100px")
            )
					)
				),
				fluidRow(
					plotlyOutput("plot_view_ref", height = "800px"),
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
			tabPanel("View QC", value = "navQC",
				DT::dataTableOutput("DT_QC")
			),
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
						shinySaveButton("saveAnalysis_Filters", "Save Filters", "Save Filters as...", 
							filetype = list(RDS = "Rds")),
            shinyFilesButton("loadAnalysis_Filters", label = "Load Filters", 
              title = "Load Filters from Rds", multiple = FALSE)
              
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
            actionButton("clear_selected_DE", "Clear Selected Events"),
            DT::dataTableOutput("DT_DE")
          )
        )
      )	# DESeq2 or DSS

		),

		navbarMenu("Display",
			tabPanel("Diagonal Plots",  value = "navDiag",
				fluidRow(
					column(3,	
						shinyWidgets::sliderTextInput(inputId = "number_events_diag", label = "Number of Top Events",
							choices = c(100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000), 
							selected = 10000),
						selectInput("EventType_diag", "Splice Type", width = '100%', multiple = TRUE,
							choices = c("IR", "MXE", "SE", "AFE", "ALE", "A5SS", "A3SS")),
						selectInput('variable_diag', 'Variable', 
							c("(none)")),
            selectInput('nom_diag', 'X-axis condition', 
							c("(none)")),
            selectInput('denom_diag', 'Y-axis condition', 
							c("(none)")),
						actionButton("clear_diag", "Clear settings"),
						textOutput("warning_diag")
					),
					column(9,
						plotlyOutput("plot_diag", height = "800px")
					)
        )
      ),
			tabPanel("Volcano Plots",  value = "navVolcano",
				fluidRow(
					column(3,	
						shinyWidgets::sliderTextInput(inputId = "number_events_volc", label = "Number of Top Events",
							choices = c(100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000), 
							selected = 10000),
						selectInput("EventType_volc", "Splice Type", width = '100%', multiple = TRUE,
							choices = c("IR", "MXE", "SE", "AFE", "ALE", "A5SS", "A3SS")),
            shinyWidgets::switchInput("facet_volc", label = "Facet by Type", labelWidth = "150px"),
						actionButton("clear_volc", "Clear settings"),
						textOutput("warning_volc")
					),
					column(9,
						plotlyOutput("plot_volc", height = "800px")
					)
        )
      ),


			tabPanel("Heatmaps", value = "navHeatmap",
				fluidRow(
					column(3, 
            shinyWidgets::radioGroupButtons("select_events_heat", 
              label = "Select Events from Differential Expression Results", justified = FALSE,
              choices = c("Highlighted", "Top N Filtered Results", "Top N All Results"), 
              checkIcon = list(yes = icon("ok", lib = "glyphicon"))
						),
            shinyWidgets::sliderTextInput("slider_num_events_heat", 
              "Num Events", choices = c(5, 10,25,50,100,200,500), selected = 25),
            shinyWidgets::radioGroupButtons("mode_heat", 
              label = "Mode", justified = FALSE,
              choices = c("PSI", "Logit", "Z-score"), 
              checkIcon = list(yes = icon("ok", lib = "glyphicon"))
						),
						selectInput('color_heat', 'Palette', 
							c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
            )            
					),
					column(9, 
            textOutput("warning_heat"),
						plotlyOutput("plot_heat", height = "800px"),
					)
        )
			),
			
			tabPanel("RNA-seq Coverage", value = "navCoverage",
				fluidRow(style='height:20vh',
					column(6, 
            textOutput("warning_cov"),
						div(style="display: inline-block;vertical-align:top; width: 250px;",
              selectizeInput('genes_cov', 'Genes', choices = "(none)")
            ),
						div(style="display: inline-block;vertical-align:top; width: 350px;",
              selectizeInput('events_cov', 'Events', choices = c("(none)"))
            ),
            shinyWidgets::radioGroupButtons("select_events_cov", 
              label = "Select Events from Differential Expression Results", justified = FALSE,
              choices = c("Highlighted", "Top N Filtered Results", "Top N All Results"), 
              checkIcon = list(yes = icon("ok", lib = "glyphicon"))
						)
          ),
					column(3,
						div(style="display: inline-block;vertical-align:top; width: 80px;",
							selectInput("chr_cov", label = "Chr", c("(none)"), selected = "(none)")),
						div(style="display: inline-block;vertical-align:top; width: 120px;",
							textInput("start_cov", label = "Left", c(""))),
						div(style="display: inline-block;vertical-align:top; width: 120px;",
							textInput("end_cov", label = "Right", c(""))),            
            br(),
            shinyWidgets::actionBttn("zoom_out_cov", style = "material-circle", color = "danger",
              icon = icon("minus")),
						div(style="display: inline-block;vertical-align:center;width: 80px;padding:35px",
							textOutput("label_zoom_cov")),
            shinyWidgets::actionBttn("zoom_in_cov", style = "material-circle", color = "danger",
              icon = icon("plus"))
					),
          column(3,
            div(style="display: inline-block;vertical-align:top;padding:10px;",
              shinyWidgets::radioGroupButtons("graph_mode_cov",label = "Graph Mode", justified = FALSE,
                choices = c("Pan", "Zoom", "Movable Labels"), 
                checkIcon = list(yes = icon("ok", lib = "glyphicon")))
              # shinyWidgets::switchInput("move_labels_cov", label = "Movable labels", labelWidth = "100px")
            ),
            shinyWidgets::sliderTextInput("slider_num_events_cov", 
              "Num Events", choices = c(5, 10,25,50,100,200,500), selected = 25)       
            
          )
				),
				fluidRow(
					column(2, 
						selectInput('event_norm_cov', 'Normalize Event', width = '100%',
							choices = c("(none)")),
						selectInput('mode_cov', 'View', width = '100%',
							choices = c("Individual", "By Condition")),					
						selectInput('condition_cov', 'Condition', width = '100%',
							choices = c("(none)")),							
						selectInput('track1_cov', 'Track 1', width = '100%',
							choices = c("(none)")),
						selectInput('track2_cov', 'Track 2', width = '100%',
							choices = c("(none)")),
						selectInput('track3_cov', 'Track 3', width = '100%',
							choices = c("(none)")),
						selectInput('track4_cov', 'Track 4', width = '100%',
							choices = c("(none)")),
            shinyWidgets::switchInput("stack_tracks_cov", label = "Stack Traces", labelWidth = "150px"),
            shinyWidgets::switchInput("pairwise_t_cov", label = "Pairwise t-test", labelWidth = "150px")
					),
					column(10, 
						plotlyOutput("plot_cov", height = "800px"),
					)
				)
      )
		) # last navbar

	)



	server = function(input, output, session) {

	# nullify NSEs to define as local variable
  gene_display_name <- bam_file <- i.bam_file <- irf_file <- i.irf_file <- junc_file <- i.junc_file <- 
	bam_file <- irf_file <- cov_file <- junc_file <- Included <- keep <- Excluded <- EventType <- 
	Events <- filtered <- padj <- i.logFC <- i.AveExpr <- i.t <- i.P.Value <- i.adj.P.Val <- i.B <- B <- 
	nom <- denom <- EventName <- log2FoldChange <- x <- ci <- track <- t_stat <- NULL

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

    default_volumes <- c("Working Directory" = getwd(), getVolumes()())
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
				
				# seqinfo
				settings = readRDS(file.path(settings_loadref$loadref_path, "settings.Rds"))
				if(settings$ah_genome != "") {
					genome = FetchAH(settings$ah_genome, ah = ah)
				} else {
					genome = rtracklayer::TwoBitFile(file.path(reference_path, "resource", "genome.2bit"))
				}
				settings_ViewRef$seqInfo = seqinfo(genome)
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
						choices = c("(none)"), selected = "(none)")    								
					updateSelectizeInput(session = session, inputId = "genes_view", server = TRUE,
						choices = c("(none)"), selected = "(none)") 
				}

			} else if(input$navSelection == "navThreads") {	
				max_cores = parallel::detectCores() - 2
				updateSliderInput(session = session, inputId = "cores_slider", 
					min = 1, max = max_cores, step = 1, value = min(4, max_cores))
				updateNumericInput(session = session, inputId = "cores_numeric", 
					min = 1, max = max_cores, value = min(4, max_cores))
			} else if(input$navSelection == "navExpr") {
				# Determine IRFinder cores

				output$txt_reference_path <- renderText({
					validate(
						need(settings_loadref$loadref_path, "Please Reference->Load and select reference path")
					)
					settings_loadref$loadref_path
				})
			} else if(input$navSelection == "navQC") {
				output$DT_QC <- DT::renderDataTable({
					validate(need(settings_SE$se, "Load Experiment file first"))
					validate(need(file.exists(
						file.path(settings_expr$collate_path, "se", "stats.fst")
						), "stats.fst does not exist in given NxtIRF output directory"))
					DT = as.data.table(read.fst((file.path(settings_expr$collate_path, "se", "stats.fst"))))
					DT = merge(as.data.table(settings_expr$df.anno), DT, all = TRUE)
					DT::datatable(
						DT,
						class = 'cell-border stripe',
						rownames = DT$sample,
						filter = 'top'
					)
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
			} else if(input$navSelection == "navDiag") {
        output$warning_diag = renderText({
          validate(need(settings_SE$se, "Please load experiment via 'Experiment' tab"))
          "Experiment Loaded"
        })
        req(settings_SE$se)
        colData = SummarizedExperiment::colData(settings_SE$se)
        if(is_valid(input$variable_diag) && input$variable_diag %in% colnames(colData)) {
          selected = isolate(input$variable_diag)
          updateSelectInput(session = session, inputId = "variable_diag", 
            choices = c("(none)", colnames(colData)), selected = selected)
        } else {
          updateSelectInput(session = session, inputId = "variable_diag", 
            choices = c("(none)", colnames(colData)), selected = "(none)")
        }
			} else if(input$navSelection == "navCoverage") {
        if(is_valid(settings_loadref$loadref_path)) {
          load_ref()
        }
        output$warning_cov <- renderText({
            validate(need(settings_loadref$loadref_path, "Please select reference path"))            
            settings_loadref$loadref_path
        })
        req(settings_loadref$loadref_path)
				
				# seqinfo
				settings = readRDS(file.path(settings_loadref$loadref_path, "settings.Rds"))
				if(settings$ah_genome != "") {
					genome = FetchAH(settings$ah_genome, ah = ah)
				} else {
					genome = rtracklayer::TwoBitFile(file.path(reference_path, "resource", "genome.2bit"))
				}
				settings_Cov$seqInfo = seqinfo(genome)
        settings_Cov$gene_list <- getGeneList()
        settings_Cov$elem.DT <- loadViewRef()
        settings_Cov$transcripts.DT <- loadTranscripts()
        
        # Populate events here
        
        if(is_valid(settings_DE$res)) {
          if(input$select_events_cov == "Highlighted") {
            selected = input$DT_DE_rows_selected
          } else if(input$select_events_cov == "Top N Filtered Results") {
            selected = input$DT_DE_rows_all
            if(length(selected) > input$slider_num_events_cov) {
              selected = selected[seq_len(input$slider_num_events_cov)]
            }
          } else {
            selected = seq_len(min(input$slider_num_events_cov, nrow(settings_DE$res)))
          }

          if(length(selected) > 0 & is_valid(settings_DE$res)) {
            updateSelectizeInput(session = session, inputId = "events_view", server = TRUE,
              choices = c("(none)", settings_DE$res$EventName[selected]), selected = "(none)")    								
            updateSelectizeInput(session = session, inputId = "events_cov", server = TRUE,
              choices = c("(none)", settings_DE$res$EventName[selected]), selected = "(none)")    								
          } else {
            updateSelectizeInput(session = session, inputId = "events_view", server = TRUE,
              choices = c("(none)"), selected = "(none)")    								
            updateSelectizeInput(session = session, inputId = "events_cov", server = TRUE,
              choices = c("(none)"), selected = "(none)")    								    
          }        
        }
        
				if(!is.null(settings_Cov$gene_list)) {
          message(paste("Populating drop-down box with", 
            length(unique(settings_Cov$gene_list$gene_display_name)),"genes"))
					updateSelectInput(session = session, inputId = "chr_cov", 
						choices = c("(none)", sort(unique(settings_Cov$gene_list$seqnames))), selected = "(none)")    								          
					updateSelectizeInput(session = session, inputId = "genes_cov", server = TRUE,
						choices = c("(none)", settings_Cov$gene_list$gene_display_name), selected = "(none)")    								
				} else {
					updateSelectInput(session = session, inputId = "chr_cov", 
						choices = c("(none)"), selected = "(none)")    								
					updateSelectizeInput(session = session, inputId = "genes_cov", server = TRUE,
						choices = c("(none)"), selected = "(none)") 
				}
				if(is_valid(settings_SE$se)) {
          # dissect rowData ranges
          rowData = as.data.frame(SummarizedExperiment::rowData(settings_SE$se))
          rowData$seqnames = tstrsplit(rowData$EventRegion, split=":")[[1]]
          temp1 = tstrsplit(rowData$EventRegion, split="/")
          temp2 = tstrsplit(temp1[[1]], split=":")[[2]]
          rowData$start = as.numeric(tstrsplit(temp2, split="-")[[1]])
          rowData$end = as.numeric(tstrsplit(temp2, split="-")[[2]])
          rowData$strand = temp1[[2]]
          
          settings_Cov$event.ranges = as.data.table(rowData)
        
					DT.files = as.data.table(settings_expr$df.files[, c("sample", "cov_file", "junc_file")])
					DT.files = na.omit(DT.files)
          settings_Cov$avail_cov = DT.files$cov_file
          names(settings_Cov$avail_cov) = DT.files$sample
					if(input$mode_cov == "By Condition") {
						colData = SummarizedExperiment::colData(settings_SE$se)
						colData = colData[rownames(colData) %in% DT.files$sample,]
						conditions_avail = colnames(colData)
						updateSelectInput(session = session, inputId = "condition_cov", 
							choices = c("(none)", conditions_avail), selected = "(none)")
					} else if(input$mode_cov == "Individual") {
						avail_samples = DT.files$sample
						updateSelectInput(session = session, inputId = "track1_cov", 
							choices = c("(none)", avail_samples), selected = "(none)")
						updateSelectInput(session = session, inputId = "track2_cov", 
							choices = c("(none)", avail_samples), selected = "(none)")
						updateSelectInput(session = session, inputId = "track3_cov", 
							choices = c("(none)", avail_samples), selected = "(none)")
						updateSelectInput(session = session, inputId = "track4_cov", 
							choices = c("(none)", avail_samples), selected = "(none)")
					}
				}        
			}
		})
    
	# Threads:
		observeEvent(input$cores_numeric, {
			req(input$cores_numeric)
			updateSliderInput(session = session, inputId = "cores_slider", value = input$cores_numeric)
		})
		observeEvent(input$cores_slider, {
			req(input$cores_slider)
			updateNumericInput(session = session, inputId = "cores_numeric", value = input$cores_slider)
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
			output$txt_mappa <- renderText(settings_newref$newref_mappa)
    })
    observeEvent(input$clear_mappa, {
      req(input$clear_mappa)
			settings_newref$newref_mappa = ""
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
			output$txt_NPA <- renderText(settings_newref$newref_NPA)    
    })
    observeEvent(input$clear_NPA, {
      req(input$clear_NPA)
			settings_newref$newref_NPA = ""
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
			output$txt_bl <- renderText(settings_newref$newref_bl)
    })
    observeEvent(input$clear_bl, {
      req(input$clear_bl)
			settings_newref$newref_bl = ""
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
			} else if(input$newref_genome_type == "(not specified)") {
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

				queryfirst = tstrsplit(ah.filtered[1]$sourceurl, split="/", fixed=TRUE)
				query_index = which(sapply(queryfirst, function(x) grepl("release", x)))
				choices = unlist(unique(tstrsplit(ah.filtered$sourceurl, split="/", fixed=TRUE)[query_index]))
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

				queryfirst = tstrsplit(ah.filtered[1]$sourceurl, split="/", fixed=TRUE)
				query_index = which(sapply(queryfirst, function(x) grepl("release", x)))
				choices = unlist(unique(tstrsplit(ah.filtered$sourceurl, split="/", fixed=TRUE)[query_index]))
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
				# fwrite(args.df, file.path(args$reference_path, "settings_newref.csv"), row.names = TRUE)
				if("ah_genome_tmp" %in% names(args)) {
					args$ah_genome = tstrsplit(args$ah_genome_tmp, split=":", fixed=TRUE)[[1]]
					args$ah_genome_tmp = NULL
				}
				if("ah_gtf_tmp" %in% names(args)) {
					args$ah_transcriptome = tstrsplit(args$ah_gtf_tmp, split=":", fixed=TRUE)[[1]]
					args$ah_gtf_tmp = NULL
				}
        withProgress(message = 'Building Reference', value = 0, {
          do.call(BuildReference, args)
        })
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
		settings_loadref <- reactiveValues(
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

    settings_ViewRef <- reactiveValues(
			# data
			seqInfo = NULL,
			gene_list = NULL,
      elem.DT = NULL,
      transcripts.DT = NULL,
			# view settings
      view_chr = "",
      view_start = "",
      view_end = "",
      data_start = 0,
      data_end = 0,

      plotly_relayout = NULL,
      plot_ini = FALSE
		)
		
    loadTranscripts <- function() {
      req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 
			file_path = file.path(settings_loadref$loadref_path, "fst", "Transcripts.fst")
      
      Transcripts.DT = as.data.table(read.fst(file_path))

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
			transcript_id <- NULL
			dir_path = file.path(settings_loadref$loadref_path, "fst")

      exons.DT = as.data.table(read.fst(file.path(dir_path, "Exons.fst"), 
        c("seqnames", "start", "end", "strand", "type", "transcript_id")))
      exons.DT = exons.DT[transcript_id != "protein_coding"]

      protein.DT = as.data.table(read.fst(file.path(dir_path, "Proteins.fst"),
        c("seqnames", "start", "end", "strand", "type", "transcript_id")))
      misc.DT = as.data.table(read.fst(file.path(dir_path, "Misc.fst"),
        c("seqnames", "start", "end", "strand", "type", "transcript_id")))
      # introns.DT = as.data.table(read.fst(file.path(dir_path, "junctions.fst")))
      # introns.DT[, type := "intron"]
      
      total.DT = rbindlist(list(
        exons.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")],
        protein.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")],
        misc.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")]
        # introns.DT[, c("seqnames", "start", "end", "strand", "type", "transcript_id")]
      ))
      return(total.DT)
    }
    
		getGeneList <- function() {
      req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 
			file_path = file.path(settings_loadref$loadref_path, "fst", "Genes.fst")
			if(!file.exists(file_path)) return(NULL)
			
			df = as.data.table(read.fst(file_path))
			return(df)
		}

    observe({
      view_chr = input$chr_view_ref
      view_start = suppressWarnings(as.numeric(input$start_view_ref))
      view_end = suppressWarnings(as.numeric(input$end_view_ref))
			
			has_movable_labels = input$move_labels_view_ref
			
      req(view_chr)
      req(view_start)
      req(view_end)

      message("plot_view_ref_fn")

      if(is.null(settings_ViewRef$elem.DT)) settings_ViewRef$elem.DT <- loadViewRef()
      if(is.null(settings_ViewRef$transcripts.DT)) settings_ViewRef$transcripts.DT <- loadTranscripts()

      req(settings_ViewRef$elem.DT)
      req(settings_ViewRef$transcripts.DT)
      message("stuff loaded")
 
      output$plot_view_ref <- renderPlotly({
        settings_ViewRef$plot_ini = TRUE      
        print(
					plot_view_ref_fn(
						view_chr, view_start, view_end, 
						settings_ViewRef$transcripts.DT, settings_ViewRef$elem.DT
					)
				)
      })
    })

		plot_view_ref_fn <- function(view_chr, view_start, view_end, transcripts, elems, highlight_events, condensed = FALSE) {
      
   transcript_support_level <- transcript_id <- group_id <- type <- gene_id <- plot_level <- 
	 i.gene_name <- i.gene_biotype <- i.transcript_name <- i.transcript_biotype <- display_name <- 
	 group_name <- group_biotype <- disp_x <- i.plot_level <- NULL		
			
			data_start = view_start - (view_end - view_start)
      data_end = view_end + (view_end - view_start)

      transcripts.DT = transcripts[seqnames == view_chr]
      transcripts.DT = transcripts.DT[start <= data_end & end >= data_start]
      setorder(transcripts.DT, transcript_support_level, width)

			# filter transcripts by criteria if applicable

      message(paste(nrow(transcripts.DT), " transcripts"))

      screen.DT = elems[transcript_id %in% transcripts.DT$transcript_id]
      # screen.DT = screen.DT[start <= data_end & end >= data_start]
      if(condensed != TRUE & nrow(transcripts.DT) <= 100) {
        condense_this = FALSE
        transcripts.DT[, group_id := transcript_id]
        screen.DT[, group_id := transcript_id]
        reduced.DT = copy(screen.DT)
        reduced.DT[type %in% c("CDS", "start_codon", "stop_codon"), type := "CDS"]
        reduced.DT[type != "CDS", type := "exon"]
      } else {
        condense_this = TRUE
        transcripts.DT[, group_id := gene_id]     
        screen.DT[transcripts.DT, on = "transcript_id", group_id := gene_id]
        # reduce screen.DT 
        reduced.gr = disjoin(
          makeGRangesFromDataFrame(
            as.data.frame(screen.DT)
          )
        )
        reduced.gr$type = "exon"
        OL = findOverlaps(
          reduced.gr,
          makeGRangesFromDataFrame(
            as.data.frame(screen.DT)
          )
        )
        reduced.gr$group_id[OL@from] = screen.DT$group_id[OL@to]
        OL = findOverlaps(
          reduced.gr,
          makeGRangesFromDataFrame(
            as.data.frame(screen.DT[type %in% c("CDS", "start_codon", "stop_codon")])
          )
        )
        reduced.gr$type[OL@from] = "CDS"
        
        reduced.DT = as.data.table(reduced.gr)
      }

      # add introns to reduced.DT
      introns.DT = as.data.table(grlGaps(
        split(makeGRangesFromDataFrame(as.data.frame(reduced.DT)),
          reduced.DT$group_id)
      ))
      introns.DT[, type := "intron"]
      data.table::setnames(introns.DT, "group_name", "group_id")
      
      reduced.DT = rbind(reduced.DT[, c("seqnames", "start", "end", "strand", "type", "group_id")], 
        introns.DT[, c("seqnames", "start", "end", "strand", "type", "group_id")])
      
      # Highlight events here
      if(!missing(highlight_events)) {
      
      }
      
      group.grl = split(
        makeGRangesFromDataFrame(
          as.data.frame(transcripts.DT)
        ), transcripts.DT$group_id
      )
      group.DT = as.data.table(range(group.grl))
      group.DT$group = NULL
      data.table::setnames(group.DT, "group_name", "group_id")
      # apply plot_order on transcripts.DT
      OL = findOverlaps(
        makeGRangesFromDataFrame(as.data.frame(group.DT)),
        makeGRangesFromDataFrame(as.data.frame(group.DT)),
        ignore.strand = TRUE
      )
      group.DT$plot_level = 1      
      cur_level = 1    
      while(any(group.DT$plot_level == cur_level)) {
        j = match(cur_level, group.DT$plot_level)
        repeat {
          bump_up_trs = unique(OL@to[OL@from == j])
          bump_up_trs = bump_up_trs[bump_up_trs > j]
          bump_up_trs = bump_up_trs[group.DT$plot_level[bump_up_trs] == cur_level]
          if(length(bump_up_trs) > 0) {
            group.DT[bump_up_trs, plot_level := cur_level + 1]
          }
          j = j + match(cur_level, group.DT$plot_level[-seq_len(j)])
          if(is.na(j)) break
        }
        cur_level = cur_level + 1
      }
      
      if(condense_this == TRUE) {
        group.DT[transcripts.DT, on = "group_id", c("group_name", "group_biotype") :=
          list(i.gene_name, i.gene_biotype)]
      } else {
        group.DT[transcripts.DT, on = "group_id", c("group_name", "group_biotype") :=
          list(i.transcript_name, i.transcript_biotype)]      
      }

      group.DT = group.DT[end > view_start & start < view_end]
      group.DT[strand == "+", display_name := paste(group_name, "-", group_biotype, " ->>")]
      group.DT[strand == "-", display_name := paste("<-- ", group_name, "-", group_biotype)]
      group.DT[, disp_x := 0.5 * (start + end)]
      group.DT[start < view_start & end > view_start, disp_x := 0.5 * (view_start + end)]
      group.DT[end > view_end & start < view_end, disp_x := 0.5 * (start + view_end)]
      group.DT[start < view_start & end > view_end, disp_x := 0.5 * (view_start + view_end)]
		
      reduced.DT$group_id = factor(reduced.DT$group_id, unique(group.DT$group_id), ordered = TRUE)
      reduced.DT[group.DT, on = "group_id", 
        plot_level := i.plot_level]
        # c("transcript_name", "transcript_biotype", "transcript_support_level", "display_name", "plot_level") := 
        # list(i.transcript_name, i.transcript_biotype, i.transcript_support_level, i.display_name, i.plot_level)]
      
      p = ggplot(reduced.DT)

      if(nrow(subset(as.data.frame(reduced.DT), type = "intron")) > 0) {
        p = p + geom_segment(data = subset(as.data.frame(reduced.DT), type = "intron"), 
          aes(x = start, xend = end, y = plot_level, yend = plot_level))
      }
      if(nrow(subset(as.data.frame(reduced.DT), type != "intron")) > 0) {
        p = p + 
          geom_rect(data = subset(as.data.frame(reduced.DT), type != "intron"), 
            aes(xmin = start, xmax = end, 
              ymin = plot_level - 0.1 - ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0), 
              ymax = plot_level + 0.1 + ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0)
            )
          )
      }
      
      if(condense_this == TRUE) {
        anno = list(
          x = group.DT$disp_x,
          y = group.DT$plot_level - 0.5 + 0.3 * runif(rep(1, nrow(group.DT))),
          text = group.DT$display_name,
          xref = "x", yref = "y", showarrow = FALSE)
      } else {
        anno = list(
          x = group.DT$disp_x,
          y = group.DT$plot_level - 0.4,
          text = group.DT$display_name,
          xref = "x", yref = "y", showarrow = FALSE)      
      }
      
      if(nrow(group.DT) == 0) {
        max_plot_level = 1
      } else {
        max_plot_level = max(group.DT$plot_level)
      }
      pl = ggplotly(p, source = "plotly_ViewRef", tooltip = "text") %>% layout(
        annotations = anno,
        dragmode = "pan",
        xaxis = list(range = c(view_start, view_end)),
        yaxis = list(range = c(0, 1 + max_plot_level), fixedrange = TRUE)
      )
			return(pl)			
		}

    settings_ViewRef$plotly_relayout = reactive({
      req(settings_ViewRef$plot_ini == TRUE)
      event_data("plotly_relayout", source = "plotly_ViewRef")
    })

    observeEvent(input$zoom_out_view_ref, {
      req(input$zoom_out_view_ref)
			req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 

			view_chr = input$chr_view_ref
			view_start = suppressWarnings(as.numeric(input$start_view_ref))
			view_end = suppressWarnings(as.numeric(input$end_view_ref))

			req(view_chr)
			req(view_start)
			req(view_end)
      req(view_end - view_start > 50)

			seqInfo = settings_ViewRef$seqInfo[input$chr_view_ref]
			seqmax = GenomeInfoDb::seqlengths(seqInfo)
			# get center of current range
			center = round((view_start + view_end) / 2)
			span = view_end - view_start
			# zoom range is 50 * 3^z
			cur_zoom = floor(log(span/50) / log(3))
      
			new_span = round(span * 3)
			if(new_span > seqmax - 1) new_span = seqmax - 1

      new_zoom = floor(log(new_span/50) / log(3))
      
      new_start = max(1, center - round(new_span / 2))
      updateTextInput(session = session, inputId = "start_view_ref", 
        value = new_start)
      updateTextInput(session = session, inputId = "end_view_ref", 
        value = new_start + new_span)
    })

    observeEvent(input$zoom_in_view_ref, {
      req(input$zoom_in_view_ref)
			req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 

			view_start = suppressWarnings(as.numeric(input$start_view_ref))
			view_end = suppressWarnings(as.numeric(input$end_view_ref))

			req(view_start)
			req(view_end)
      req(view_end - view_start > 50)

			# get center of current range
			center = round((view_start + view_end) / 2)
			span = view_end - view_start
			# zoom range is 50 * 3^z
			cur_zoom = floor(log(span/50) / log(3))
      
			new_span = round(span / 3)
			if(new_span < 50) new_span = 50

      new_zoom = floor(log(new_span/50) / log(3))
      
      new_start = max(1, center - round(new_span / 2))
      updateTextInput(session = session, inputId = "start_view_ref", 
        value = new_start)
      updateTextInput(session = session, inputId = "end_view_ref", 
        value = new_start + new_span)
    })


    observeEvent(settings_ViewRef$plotly_relayout(), {
      plotly_relayout = settings_ViewRef$plotly_relayout()
      message(names(plotly_relayout))
      req(length(plotly_relayout) == 2)
      req(all(c("xaxis.range[0]", "xaxis.range[1]") %in% names(plotly_relayout)))

      updateTextInput(session = session, inputId = "start_view_ref", 
        value = max(1, round(plotly_relayout[["xaxis.range[0]"]])))
      updateTextInput(session = session, inputId = "end_view_ref", 
        value = round(plotly_relayout[["xaxis.range[1]"]]))
        
    })

		observeEvent(input$genes_view, {
      req(input$genes_view)
      req(input$genes_view != "(none)")

      gene_id_view = settings_ViewRef$gene_list[gene_display_name == input$genes_view]

      # change settings on input based on this
      updateSelectInput(session = session, inputId = "chr_view_ref", 
        selected = gene_id_view$seqnames[1])
      updateTextInput(session = session, inputId = "start_view_ref", 
        value = gene_id_view$start[1])
      updateTextInput(session = session, inputId = "end_view_ref", 
        value = gene_id_view$end[1])
    })

    observeEvent(input$chr_view_ref, {
      req(input$chr_view_ref != "(none)")
			seqInfo = settings_ViewRef$seqInfo[input$chr_view_ref]
			seqmax = GenomeInfoDb::seqlengths(seqInfo)
      
			req(as.numeric(input$end_view_ref))
      if(as.numeric(input$end_view_ref) > seqmax) {
        updateTextInput(session = session, inputId = "end_view_ref", 
          value = seqmax)      
        req(as.numeric(input$start_view_ref))
        if(seqmax - as.numeric(input$start_view_ref) < 50) {
          updateTextInput(session = session, inputId = "end_view_ref", 
            value = seqmax - 50)        
        }
      }
    })
    
    observeEvent(input$start_view_ref, {
			req(as.numeric(input$start_view_ref))
			req(as.numeric(input$end_view_ref))
			req(as.numeric(input$end_view_ref) - as.numeric(input$start_view_ref) > 50)

			# adjust zoom
			span = as.numeric(input$end_view_ref) - as.numeric(input$start_view_ref)
			cur_zoom = floor(log(span/50) / log(3))
			output$label_zoom_view_ref <- renderText({cur_zoom})
    })
    observeEvent(input$end_view_ref, {
			req(as.numeric(input$start_view_ref))
			req(as.numeric(input$end_view_ref))		
			req(as.numeric(input$end_view_ref) - as.numeric(input$start_view_ref) > 50)

			# adjust zoom
			span = as.numeric(input$end_view_ref) - as.numeric(input$start_view_ref)
			cur_zoom = floor(log(span/50) / log(3))
			output$label_zoom_view_ref <- renderText({cur_zoom})
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
				if(is(BPPARAM, "SnowParam")) {
					BPPARAM_mod = BiocParallel::SnowParam(input$cores_slider)
          message(paste("Using SnowParam", input$cores_slider, "threads"))
				} else if(is(BPPARAM, "MulticoreParam")) {
					BPPARAM_mod = BiocParallel::MulticoreParam(input$cores_slider)
          message(paste("Using MulticoreParam", input$cores_slider, "threads"))
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
      settings_expr$df.anno <- data.table::setnames(cbind(df, newcolumn, stringsAsFactors=FALSE), 
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

			cores_to_use = as.numeric(input$cores_slider)
			if(!is_valid(input$cores_slider)) cores_to_use = 1
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
		
    getFilterData = function(i) {
      if(is_valid(settings_SE$filters)) {
        return(settings_SE$filters[[i]])
      } else {
        return(list())
      }
    }
    reactive_filter1 <- reactive({getFilterData(1)})
    reactive_filter2 <- reactive({getFilterData(2)})
    reactive_filter3 <- reactive({getFilterData(3)})
    reactive_filter4 <- reactive({getFilterData(4)})
    reactive_filter5 <- reactive({getFilterData(5)})
    reactive_filter6 <- reactive({getFilterData(6)})
    reactive_filter7 <- reactive({getFilterData(7)})
    reactive_filter8 <- reactive({getFilterData(8)})
    
    conditionList = reactive({
      req(settings_SE$se)
      if(is(settings_SE$se, "SummarizedExperiment")) {
        colnames(SummarizedExperiment::colData(settings_SE$se))
      } else {
        c("")
      }
    })
    filter1 <- filterModule_server("filter1", reactive_filter1, conditionList)
    filter2 <- filterModule_server("filter2", reactive_filter2, conditionList)
    filter3 <- filterModule_server("filter3", reactive_filter3, conditionList)
    filter4 <- filterModule_server("filter4", reactive_filter4, conditionList)
    filter5 <- filterModule_server("filter5", reactive_filter5, conditionList)
    filter6 <- filterModule_server("filter6", reactive_filter6, conditionList)
    filter7 <- filterModule_server("filter7", reactive_filter7, conditionList)
    filter8 <- filterModule_server("filter8", reactive_filter8, conditionList)
    
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
        if(is_valid(settings_SE$filters)) {
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
        }
        settings_SE$filterSummary = filterSummary
        message(sum(filterSummary == TRUE))
      }
    }
    
    # observeEvent({
      # settings_SE$se
      # settings_SE$filters[[1]]$trigger
      # settings_SE$filters[[2]]$trigger
      # settings_SE$filters[[3]]$trigger
      # settings_SE$filters[[4]]$trigger
      # settings_SE$filters[[5]]$trigger
      # settings_SE$filters[[6]]$trigger
      # settings_SE$filters[[7]]$trigger
      # settings_SE$filters[[8]]$trigger
    # }, {
      # processFilters()
    # })
    
    observeEvent(input$refresh_filters_Filters, {
      req(input$refresh_filters_Filters)
      
      settings_SE$filters[[1]] = (reactiveValuesToList(filter1))
      settings_SE$filters[[2]] = (reactiveValuesToList(filter2))
      settings_SE$filters[[3]] = (reactiveValuesToList(filter3))
      settings_SE$filters[[4]] = (reactiveValuesToList(filter4))
      settings_SE$filters[[5]] = (reactiveValuesToList(filter5))
      settings_SE$filters[[6]] = (reactiveValuesToList(filter6))
      settings_SE$filters[[7]] = (reactiveValuesToList(filter7))
      settings_SE$filters[[8]] = (reactiveValuesToList(filter8))

      processFilters()
    })
    
    observe({
      req(settings_SE$se)
      req(settings_SE$filterSummary)

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
    })
    
    shinyFileSave(input, "saveAnalysis_Filters", roots = c(default_volumes, addit_volume), session = session)
    observeEvent(input$saveAnalysis_Filters, {
      req(settings_SE$filters)
      selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$saveAnalysis_Filters)
      req(selectedfile$datapath)

      settings_SE$filters[[1]] = (reactiveValuesToList(filter1))
      settings_SE$filters[[2]] = (reactiveValuesToList(filter2))
      settings_SE$filters[[3]] = (reactiveValuesToList(filter3))
      settings_SE$filters[[4]] = (reactiveValuesToList(filter4))
      settings_SE$filters[[5]] = (reactiveValuesToList(filter5))
      settings_SE$filters[[6]] = (reactiveValuesToList(filter6))
      settings_SE$filters[[7]] = (reactiveValuesToList(filter7))
      settings_SE$filters[[8]] = (reactiveValuesToList(filter8))

      final = settings_SE$filters
      saveRDS(final, selectedfile$datapath)
    })

		shinyFileChoose(input, "loadAnalysis_Filters", roots = c(default_volumes, addit_volume), session = session,
      filetypes = c("Rds"))
    observeEvent(input$loadAnalysis_Filters, {
      selectedfile <- parseFilePaths(c(default_volumes, addit_volume), input$loadAnalysis_Filters)
      req(selectedfile$datapath)
      settings_SE$filters = readRDS(selectedfile$datapath)
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
			
			if(!is(colData[,input$variable_DE], "factor")) {
				output$warning_DE = renderText("Contrast must be performed on discrete categories")
        updateSelectInput(session = session, inputId = "variable_DE", 
          choices = c("(none)", colnames(colData)), selected = "(none)")
			} else {
				updateSelectInput(session = session, inputId = "nom_DE", 
          choices = c("(none)", levels(colData[,input$variable_DE])), selected = "(none)")					
        updateSelectInput(session = session, inputId = "denom_DE", 
          choices = c("(none)", levels(colData[,input$variable_DE])), selected = "(none)")
				if(is_valid(settings_DE$nom_DE) && settings_DE$nom_DE %in% levels(colData[,input$variable_DE])) {
					updateSelectInput(session = session, inputId = "nom_DE", selected = settings_DE$nom_DE)
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
		observeEvent(settings_DE$nom_DE, {
			req(settings_DE$nom_DE)
      updateSelectInput(session = session, inputId = "nom_DE_DE", selected = settings_DE$nom_DE)		
		})
		observeEvent(settings_DE$denom_DE, {
			req(settings_DE$denom_DE)
      updateSelectInput(session = session, inputId = "denom_DE", selected = settings_DE$denom_DE)		
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

        cores_to_use = as.numeric(input$cores_slider)
        if(!is_valid(input$cores_slider)) cores_to_use = 1
        
        BPPARAM = BiocParallel::bpparam()
        
        if(is(BPPARAM, "SnowParam")) {
          BPPARAM_mod = BiocParallel::SnowParam(cores_to_use)
          message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
        } else if(is(BPPARAM, "MulticoreParam")) {
          BPPARAM_mod = BiocParallel::MulticoreParam(cores_to_use)
          message(paste("Using MulticoreParam", BPPARAM_mod$workers, "threads"))
        } else {
          BPPARAM_mod = BiocParallel::SerialParam()
          message(paste("Using SerialParam mode with", BPPARAM_mod$workers, "threads"))
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
				res = cbind(as.data.frame(rowData[,1:3]), res)
				res = res %>% dplyr::arrange(padj)
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
				
				rowData.DT = as.data.table(rowData)
        rowData.DT = rowData.DT[, 1:3]
				res.ASE = rowData.DT[res.ASE, on = "EventName"]

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
      session = session, filetype = c("Rds"))
    observeEvent(input$load_DE, {
      req(input$load_DE)
      file_selected<-parseFilePaths(c(default_volumes, addit_volume), input$load_DE)
      req(file_selected$datapath)
			load_DE = readRDS(as.character(file_selected$datapath))
			req(all(c("res", "settings") %in% names(load_DE)))
			
		# check all parameters exist in colData(se)
			req(settings_SE$se)
			colData = SummarizedExperiment::colData(settings_SE$se)
			req(load_DE$settings$DE_Var %in% colnames(colData))
			req(!is_valid(load_DE$settings$batchVar1) || load_DE$settings$batchVar1 %in% colnames(colData))
			req(!is_valid(load_DE$settings$batchVar2) || load_DE$settings$batchVar2 %in% colnames(colData))
			req(any(unlist(colData[,load_DE$settings$DE_Var]) == load_DE$settings$nom_DE))
			req(any(unlist(colData[,load_DE$settings$DE_Var]) == load_DE$settings$denom_DE))
			req(load_DE$settings$method %in% c("DESeq2", "limma", "DSS"))
			
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
		
		observeEvent(input$clear_selected_DE, {
			req(settings_DE$res)
			req(input$DT_DE_rows_selected)
			DT::dataTableProxy("DT_DE") %>% DT::selectRows(NULL)
		})
	# Diagonal Plots
	
    settings_Diag = reactiveValues(
      plot_ini = FALSE,
			plotly_click = NULL
    )
    
		output$plot_diag <- renderPlotly({
      # settings_Diag$plot_ini = FALSE
			validate(need(settings_SE$se, "Load Experiment first"))
			validate(need(settings_DE$res, "Load DE Analysis first"))
			validate(need(input$variable_diag, "Select conditions and contrasts"))
			validate(need(input$nom_diag, "Select conditions and contrasts"))
			validate(need(input$denom_diag, "Select conditions and contrasts"))
			validate(need(input$variable_diag != "(none)", "Select conditions and contrasts"))
			validate(need(input$nom_diag != "(none)", "Select conditions and contrasts"))
			validate(need(input$denom_diag != "(none)", "Select conditions and contrasts"))
			
      selected = input$DT_DE_rows_selected      
      
			num_events = input$number_events_diag
			res = as.data.table(settings_DE$res[input$DT_DE_rows_all,])
			if(is_valid(input$EventType_diag)) {
				res = res[EventType %in% input$EventType_diag]
			}
			if(num_events < nrow(res)) {
				res = res[seq_len(num_events)]
			}
			df.diag = make_diagonal(settings_SE$se, res$EventName, input$variable_diag,
				input$nom_diag, input$denom_diag)
			
      if(is_valid(input$DT_DE_rows_selected)) {
        df.diag$selected = (df.diag$EventName %in% settings_DE$res$EventName[selected])
      } else {
        df.diag$selected = FALSE
      }
      
      settings_Diag$plot_ini = TRUE
			print(
				ggplotly(
					ggplot(df.diag, aes(x = nom, y = denom, key = EventName, text = EventName, colour = selected)) + 
            geom_point() + scale_color_manual(values = c("black", "red")),
					tooltip = "text",
          source = "plotly_diagonal"
				) %>% layout(
					yaxis = list(scaleanchor="x", scaleratio=1)
				)
			)
		})
		
    settings_Diag$plotly_click = reactive({
      plot_exist = settings_Diag$plot_ini
      if(plot_exist == TRUE) {
        event_data("plotly_click", source = "plotly_diagonal")
      }
    })
  
    observeEvent(settings_Diag$plotly_click(), {
      req(settings_Diag$plotly_click())
      click = settings_Diag$plotly_click()
      print(click)
      click.id = which(settings_DE$res$EventName == click$key)
      req(click.id)
      
      selected = input$DT_DE_rows_selected
      
      if(click.id %in% selected) {
        selected = selected[-which(selected == click.id)]
      } else {
        selected = c(selected, click.id)
      }
      
      DT::dataTableProxy("DT_DE") %>% DT::selectRows(selected)

    })
    
		observeEvent(input$variable_diag, {
			req(settings_SE$se)
			req(input$variable_diag != "(none)")
			colData = SummarizedExperiment::colData(settings_SE$se)
			req(input$variable_diag %in% colnames(colData))

			if(!is(colData[,input$variable_diag], "factor")) {
				output$warning_diag = renderText("Contrast must be performed on discrete categories")
        updateSelectInput(session = session, inputId = "variable_diag", 
          choices = c("(none)", colnames(colData)), selected = "(none)")
			} else {
				updateSelectInput(session = session, inputId = "nom_diag", 
					 choices = c("(none)", levels(colData[,input$variable_diag])), selected = "(none)")
				updateSelectInput(session = session, inputId = "denom_diag", 
					 choices = c("(none)", levels(colData[,input$variable_diag])), selected = "(none)")
			}
		})

		observeEvent(input$clear_diag, {
			updateSelectInput(session = session, "EventType_diag", selected = NULL)
			shinyWidgets::updateSliderTextInput(session = session, "number_events_diag", selected = 10000)
			if(is_valid(settings_SE$se)) {
				colData = SummarizedExperiment::colData(settings_SE$se)
        updateSelectInput(session = session, inputId = "variable_diag", 
          choices = c("(none)", colnames(colData)), selected = "(none)")
			} else {
        updateSelectInput(session = session, inputId = "variable_diag", 
          choices = c("(none)"), selected = "(none)")			
			}
			updateSelectInput(session = session, inputId = "nom_diag", 
				 choices = c("(none)"), selected = "(none)")			
			updateSelectInput(session = session, inputId = "denom_diag", 
				 choices = c("(none)"), selected = "(none)")			
		})

# Volcano Plots
    settings_Volc = reactiveValues(
      plot_ini = FALSE,
			plotly_click = NULL
    )

    settings_Volc$plotly_click = reactive({
      plot_exist = settings_Volc$plot_ini
      if(plot_exist == TRUE) {
        event_data("plotly_click", source = "plotly_volcano")
      }
    })
  
    observeEvent(settings_Volc$plotly_click(), {
      req(settings_Volc$plotly_click())
      click = settings_Volc$plotly_click()
      print(click)
      click.id = which(settings_DE$res$EventName == click$key)
      req(click.id)
      
      selected = input$DT_DE_rows_selected
      
      if(click.id %in% selected) {
        selected = selected[-which(selected == click.id)]
      } else {
        selected = c(selected, click.id)
      }
      
      DT::dataTableProxy("DT_DE") %>% DT::selectRows(selected)
    })

    settings_Volc$plotly_brush = reactive({
      plot_exist = settings_Volc$plot_ini
      if(plot_exist == TRUE) {
        event_data("plotly_selected", source = "plotly_volcano")
      }
    })
  
    observeEvent(settings_Volc$plotly_brush(), {
      req(settings_Volc$plotly_brush())
      brush = settings_Volc$plotly_brush()
      print(brush)
      brush.id = which(settings_DE$res$EventName %in% brush$key)
      req(brush.id)
      
      selected = input$DT_DE_rows_selected
      selected = unique(c(selected, brush.id))
      
      DT::dataTableProxy("DT_DE") %>% DT::selectRows(selected)
   })

		
		output$plot_volc <- renderPlotly({
      # settings_Diag$plot_ini = FALSE
			validate(need(settings_SE$se, "Load Experiment first"))
			validate(need(settings_DE$res, "Load DE Analysis first"))
			
      selected = input$DT_DE_rows_selected  
      
			num_events = input$number_events_volc
			res = as.data.table(settings_DE$res[input$DT_DE_rows_all,])
			if(is_valid(input$EventType_volc)) {
				res = res[EventType %in% input$EventType_volc]
			}
			if(num_events < nrow(res)) {
				res = res[seq_len(num_events)]
			}
			if(input$method_DE == "DESeq2") {
				df.volc = with(res, data.frame(EventName = EventName, EventType = EventType,
					log2FoldChange = log2FoldChange, pvalue = pvalue, padj = padj))
			} else {
				df.volc = with(res, data.frame(EventName = EventName, EventType = EventType,
					log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val))	
			}
			
      if(is_valid(selected)) {
        df.volc$selected = (df.volc$EventName %in% settings_DE$res$EventName[selected])
      } else {
        df.volc$selected = FALSE
      }
      settings_Volc$plot_ini = TRUE
			
			p = ggplot(df.volc, aes(x = log2FoldChange, y = -log10(padj), 
						key = EventName, text = EventName, colour = selected)) + 
            geom_point() + scale_color_manual(values = c("black", "red"))
			if(input$facet_volc == TRUE) {
				p = p + facet_wrap(vars(EventType))
			}
			print(
				ggplotly(p,
					tooltip = "text",
          source = "plotly_volcano"
				) %>% layout(dragmode = "lasso")
			)
		})
		
		observeEvent(input$clear_volc, {
			updateSelectInput(session = session, "EventType_volc", selected = NULL)
			shinyWidgets::updateSliderTextInput(session = session, "number_events_volc", selected = 10000)
		
		})
		
		# Heatmaps
		
		output$plot_heat <- renderPlotly({
      # settings_Diag$plot_ini = FALSE
			validate(need(settings_SE$se, "Load Experiment first"))
			validate(need(settings_DE$res, "Load DE Analysis first"))
			
			if(input$select_events_heat == "Highlighted") {
				selected = input$DT_DE_rows_selected
			} else if(input$select_events_heat == "Top N Filtered Results") {
				selected = input$DT_DE_rows_all
				if(length(selected) > input$slider_num_events_heat) {
					selected = selected[seq_len(input$slider_num_events_heat)]
				}
			} else {
				selected = seq_len(min(input$slider_num_events_heat, nrow(settings_DE$res)))
			}

      validate(need(length(selected) > 0, "Select some Events first"))
      
			colData = as.data.frame(SummarizedExperiment::colData(settings_SE$se))
      
			if(input$mode_heat == "PSI") {
				mat = make_matrix(settings_SE$se, settings_DE$res$EventName[selected],
					rownames(colData), "PSI")
			} else if(input$mode_heat == "Logit") {
				mat = make_matrix(settings_SE$se, settings_DE$res$EventName[selected],
					rownames(colData), "logit")			
			} else {
				mat = make_matrix(settings_SE$se, settings_DE$res$EventName[selected],
					rownames(colData), "Z-score")
			}
			
      validate(need(nrow(mat) > 0 & ncol(mat) > 0, "No data after filtering results"))
      
      colors.df = RColorBrewer::brewer.pal.info
      color.index = which(rownames(colors.df) == input$color_heat)
      color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = colors.df$maxcolors[color.index],
        name = rownames(colors.df)[color.index])))
        
      na.exclude = (rowSums(!is.na(mat)) == 0)
      if(any(na.exclude == TRUE)) {
        output$warning_heat <- renderText({
            cat("The following events have been excluded due to all NA values:")
            paste(rownames(mat)[which(na.exclude)])
        })
        mat = mat[-which(na.exclude),]
      }
			print(
        # ggplotly(ggplotify::ggplotify(
          # pheatmap::pheatmap(mat, annotation_col = colData, color = color)
        # ))

  			heatmaply::heatmaply(mat, color = color, col_side_colors = colData)        
  		)
		})


		
		# RNA-seq Coverage Plots
    settings_Cov <- shiny::reactiveValues(
			# data
			seqInfo = NULL,
			gene_list = NULL,
      elem.DT = NULL,
      transcripts.DT = NULL,
			# view settings
      view_chr = "",
      view_start = "",
      view_end = "",
      data_start = 0,
      data_end = 0,
			
			view_strand = "*",

      event.ranges = NULL,
			avail_cov = NULL,	# named vector of files

      plotly_relayout = NULL,
      plot_ini = FALSE
		)
		get_track_selection <- function(i) {
			if(i == 1) return(input$track1_cov)
			if(i == 2) return(input$track2_cov)
			if(i == 3) return(input$track3_cov)
			if(i == 4) return(input$track4_cov)
		}
		
    update_norm_events_cov <- function(view_chr, view_start, view_end, default_event = NULL) {
      req(settings_Cov$event.ranges)
      event.ranges.legit = settings_Cov$event.ranges[seqnames == view_chr]
      event.ranges.legit = event.ranges.legit[end > view_start & start < view_end]

      if(!is_valid(default_event) && is_valid(isolate(input$event_norm_cov))) {
        default_event = isolate(input$event_norm_cov)
      }

      if(nrow(event.ranges.legit) == 0) {
        # retain selected event
        if(is_valid(default_event)) {
          updateSelectInput(session = session, inputId = "event_norm_cov", 
            choices = c("(none)", default_event), selected = default_event)
          return(default_event)
        } else {
          updateSelectInput(session = session, inputId = "event_norm_cov", 
            choices = c("(none)"), selected = "(none)")
          return("(none)")
        }
      } else if(is_valid(default_event) && default_event %in% event.ranges.legit$EventName) {
        updateSelectInput(session = session, inputId = "event_norm_cov", 
          choices = c("(none)", event.ranges.legit$EventName), selected = default_event)
        return(default_event)
      } else if(is_valid(default_event)){
        updateSelectInput(session = session, inputId = "event_norm_cov", 
          choices = c("(none)", default_event, event.ranges.legit$EventName), selected = default_event)      
        return(default_event)
      } else {
        updateSelectInput(session = session, inputId = "event_norm_cov", 
          choices = c("(none)", event.ranges.legit$EventName), selected = "(none)")      
        return("(none)")
      }    
    }
    
    observe({
      view_chr = input$chr_cov
      view_start = suppressWarnings(as.numeric(input$start_cov))
      view_end = suppressWarnings(as.numeric(input$end_cov))
			graph_mode = isolate(input$graph_mode_cov)
			
      req(view_chr)
      req(view_start)
      req(view_end)
      req(settings_SE$se)
      
      selected_event = isolate(input$events_cov)
      if(is_valid(selected_event)) {
        cur_event = update_norm_events_cov(view_chr, view_start, view_end, selected_event)
      } else {
        cur_event = NULL
      }
      
      # refresh in-range events here

      conf.int = 0.95

      if(is.null(settings_Cov$elem.DT)) settings_Cov$elem.DT <- loadViewRef()
      if(is.null(settings_Cov$transcripts.DT)) settings_Cov$transcripts.DT <- loadTranscripts()

      req(settings_Cov$elem.DT)
      req(settings_Cov$transcripts.DT)
      message("stuff loaded")
			
			p_ref = plot_view_ref_fn(
				view_chr, view_start, view_end, 
				settings_Cov$transcripts.DT, settings_Cov$elem.DT
			)
			p_track = list()
      
      cur_zoom = floor(log((view_end - view_start)/50) / log(3))

      data.list = list()
      data.t_test = NULL
      fac = NULL
      
      if(is_valid(input$condition_cov) & is_valid(cur_event)) {
        for(i in 1:4) {
          track_samples = get_track_selection(i)
          if(is_valid(track_samples)) {
            colData = SummarizedExperiment::colData(settings_SE$se)
            samples = rownames(colData)[unlist(as.character(colData[, input$condition_cov]) == track_samples)]
            event_norms = SummarizedExperiment::assay(settings_SE$se, "Depth")[cur_event,samples]
            samples = samples[event_norms >= 10]
            event_norms = event_norms[event_norms >= 10]
            
            df = as.data.frame(GetCoverage_DF(samples, settings_Cov$avail_cov[samples],
              view_chr, view_start, view_end, settings_Cov$view_strand))
            # bin anything with cur_zoom > 5
            df = bin_df(df, max(1, 3^(cur_zoom - 5)))
              
            for(todo in seq_len(length(samples))) {
              df[, samples[todo]] = df[, samples[todo]] / event_norms[todo]
            }
            
            if(input$pairwise_t_cov == TRUE) {
              if(is.null(data.t_test)) {
                data.t_test <- as.matrix(df)
                fac = rep(as.character(i), ncol(df) - 1)
              } else {
                data.t_test <- cbind(data.t_test, as.matrix(df[, -1]))
                fac = c(fac, rep(as.character(i), ncol(df) - 1))
              }
            }
            
            df$mean = rowMeans(as.matrix(df[,samples]))
            df$sd = matrixStats::rowSds(as.matrix(df[,samples]))
            n = length(samples)
            df$ci = qt((1 + conf.int)/2,df=	n-1) * df$sd / sqrt(n)

            df$track = as.character(i)
            df = df %>% dplyr::select(x, mean, ci, track)
            data.list[[i]] <- as.data.table(df)
          }
        }
        if(input$stack_tracks_cov == TRUE) {
          df = as.data.frame(rbindlist(data.list))
          if(nrow(df) > 0) {
            p_track[[1]] = ggplotly(
              ggplot(df, aes(x = x)) + 
                geom_ribbon(alpha = 0.2, aes(y = mean, ymin = mean - ci, ymax = mean + ci, fill = track)) +
                geom_line(aes(y = mean, colour = track)) +
                labs(y = "Stacked Tracks Normalized Coverage"),
                tooltip = c("x", "y", "ymin", "ymax", "colour")
            )
            p_track[[1]] = p_track[[1]] %>% layout(
              dragmode = "zoom",
              yaxis = list(range = c(0, 1 + max(df$mean + df$ci)), fixedrange = TRUE)
            )
          }
        } else {
          for(i in 1:4) {
            if(length(data.list) >= i && !is.null(data.list[[i]])) {
              p_track[[i]] = ggplotly(
                ggplot(df, aes(x = x)) + 
                  geom_ribbon(alpha = 0.2, colour = NA, aes(y = mean, ymin = mean - ci, ymax = mean + ci)) +
                  geom_line(aes(y = mean)) +
                  labs(y = paste("Track", i, "Normalized Coverage")),
                tooltip = c("x", "y", "ymin", "ymax")
              )						
              p_track[[i]] = p_track[[i]] %>% layout(
                yaxis = list(range = c(0, 1 + max(df$mean + df$ci)), fixedrange = TRUE)
              )
            }
          }
        }
      } else if(input$mode_cov == "Individual") {
        for(i in 1:4) {
          track_samples = get_track_selection(i)
          if(is_valid(track_samples)) {
            df = GetCoverage_DF(track_samples, settings_Cov$avail_cov[track_samples],
              view_chr, view_start, view_end, settings_Cov$view_strand)
            df = bin_df(df, max(1, 3^(cur_zoom - 5)))
            data.list[[i]] <- as.data.table(df)

            p_track[[i]] = ggplotly(
              ggplot(df, aes_string(x = "x", y = track_samples)) + geom_line() +
                  labs(y = paste(track_samples, " Coverage")),
                tooltip = c("x", "y")
            )
            p_track[[i]] = p_track[[i]] %>% layout(
              yaxis = list(range = c(0, 1 + max(unlist(df[,track_samples]))), fixedrange = TRUE)
            )
          }
        }
      }
            
      if(input$pairwise_t_cov == TRUE && !is.null(fac) && length(unique(fac)) == 2) {
        fac = factor(fac)
        t_test = genefilter::rowttests(data.t_test[, -1], fac)
        DT = data.table(x = data.t_test[, 1])
        DT[, t_stat := -log10(t_test$p.value)]
        
        p_track[[5]] = ggplotly(
          ggplot(as.data.frame(DT), aes_string(x = "x", y = "t_stat")) + geom_line() +
              labs(y = paste("Pairwise T-test -log10(p)")),
            tooltip = c("x", "y")
        )
        p_track[[5]] = p_track[[5]] %>% layout(
          yaxis = list(c(0, 1 + max(DT$t_stat)), fixedrange = TRUE)
        )        
      }
      
      plot_tracks = p_track[unlist(lapply(p_track, function(x) !is.null(x)))]

      plot_tracks[[length(plot_tracks) + 1]] = p_ref

      final_plot = subplot(plot_tracks, nrows = length(plot_tracks), shareX = TRUE)
      
      if(graph_mode == "Pan") {
        final_plot = final_plot %>% layout(
          dragmode = "pan"
        )
      } else if(graph_mode == "Zoom") {
        final_plot = final_plot %>% layout(
          dragmode = "zoom"
        )   
      } else if(graph_mode == "Movable Labels") {
        final_plot = final_plot %>% layout(
          dragmode = FALSE
        ) %>% config(editable = TRUE)
      }
      
      output$plot_cov <- renderPlotly({
        settings_Cov$plot_ini = TRUE      
        print(
					final_plot
				)
      })
    })
		
    observeEvent(input$graph_mode_cov, {
      req(settings_Cov$plot_ini == TRUE)
      if(input$graph_mode_cov == "Pan") {
        plotlyProxy("plot_cov", session) %>% plotlyProxyInvoke("relayout", list(dragmode = "pan")) %>%
          plotlyProxyInvoke("reconfig", editable = FALSE)
      } else if(input$graph_mode_cov == "Zoom") {
        plotlyProxy("plot_cov", session) %>% plotlyProxyInvoke("relayout", list(dragmode = "zoom")) %>%
          plotlyProxyInvoke("reconfig", editable = FALSE)
      } else if(input$graph_mode_cov == "Movable Labels") {
        plotlyProxy("plot_cov", session) %>% plotlyProxyInvoke("relayout", list(dragmode = FALSE)) %>%
          plotlyProxyInvoke("reconfig", editable = TRUE)
      }
    })
    
		observeEvent(input$mode_cov, {
      req(input$mode_cov)
      req(settings_SE$se)
      if(input$mode_cov == "By Condition") {
        colData = SummarizedExperiment::colData(settings_SE$se)
        
        updateSelectInput(session = session, inputId = "condition_cov", 
          choices = c("(none)", colnames(colData)))      
      } else {
        updateSelectInput(session = session, inputId = "condition_cov", 
          choices = c("(none)"))      
      }
    })

		observeEvent(input$condition_cov, {
      req(input$condition_cov)
      req(settings_SE$se)
      if(is_valid(input$condition_cov)) {
        colData = SummarizedExperiment::colData(settings_SE$se)
        updateSelectInput(session = session, inputId = "track1_cov", 
          choices = c("(none)", unique(as.character(unlist(colData[, input$condition_cov])))))    
        updateSelectInput(session = session, inputId = "track2_cov", 
          choices = c("(none)", unique(as.character(unlist(colData[, input$condition_cov])))))    
        updateSelectInput(session = session, inputId = "track3_cov", 
          choices = c("(none)", unique(as.character(unlist(colData[, input$condition_cov])))))    
        updateSelectInput(session = session, inputId = "track4_cov", 
          choices = c("(none)", unique(as.character(unlist(colData[, input$condition_cov])))))    
      } else {
        updateSelectInput(session = session, inputId = "track1_cov", 
          choices = c("(none)"))     
        updateSelectInput(session = session, inputId = "track2_cov", 
          choices = c("(none)"))  
        updateSelectInput(session = session, inputId = "track3_cov", 
          choices = c("(none)"))    
        updateSelectInput(session = session, inputId = "track4_cov", 
          choices = c("(none)"))    
      }
    })

    
    settings_Cov$plotly_relayout = reactive({
      req(settings_Cov$plot_ini == TRUE)
      event_data("plotly_relayout", source = "plotly_ViewRef")
    })
    observeEvent(settings_Cov$plotly_relayout(), {
      plotly_relayout = settings_Cov$plotly_relayout()
      message(names(plotly_relayout))
      req(length(plotly_relayout) == 2)
      req(all(c("xaxis.range[0]", "xaxis.range[1]") %in% names(plotly_relayout)))

      updateTextInput(session = session, inputId = "start_cov", 
        value = max(1, round(plotly_relayout[["xaxis.range[0]"]])))
      updateTextInput(session = session, inputId = "end_cov", 
        value = round(plotly_relayout[["xaxis.range[1]"]]))
    })
    observeEvent(input$zoom_out_cov, {
      req(input$zoom_out_cov)
			req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 

			view_chr = input$chr_cov
			view_start = suppressWarnings(as.numeric(input$start_cov))
			view_end = suppressWarnings(as.numeric(input$end_cov))

			req(view_chr)
			req(view_start)
			req(view_end)
      req(view_end - view_start > 50)

			seqInfo = settings_Cov$seqInfo[input$chr_cov]
			seqmax = GenomeInfoDb::seqlengths(seqInfo)
			# get center of current range
			center = round((view_start + view_end) / 2)
			span = view_end - view_start
			# zoom range is 50 * 3^z
			cur_zoom = floor(log(span/50) / log(3))
      
			new_span = round(span * 3)
			if(new_span > seqmax - 1) new_span = seqmax - 1

      new_zoom = floor(log(new_span/50) / log(3))
      
      new_start = max(1, center - round(new_span / 2))
      updateTextInput(session = session, inputId = "start_cov", 
        value = new_start)
      updateTextInput(session = session, inputId = "end_cov", 
        value = new_start + new_span)
    })		
    observeEvent(input$zoom_in_cov, {
      req(input$zoom_in_cov)
			req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 

			view_start = suppressWarnings(as.numeric(input$start_cov))
			view_end = suppressWarnings(as.numeric(input$end_cov))

			req(view_start)
			req(view_end)
      req(view_end - view_start > 50)

			# get center of current range
			center = round((view_start + view_end) / 2)
			span = view_end - view_start
			# zoom range is 50 * 3^z
			cur_zoom = floor(log(span/50) / log(3))
      
			new_span = round(span / 3)
			if(new_span < 50) new_span = 50

      new_zoom = floor(log(new_span/50) / log(3))
      
      new_start = max(1, center - round(new_span / 2))
      updateTextInput(session = session, inputId = "start_cov", 
        value = new_start)
      updateTextInput(session = session, inputId = "end_cov", 
        value = new_start + new_span)
    })
		observeEvent(input$events_cov, {
      req(input$events_cov)
      req(input$events_cov != "(none)")

      events_id_view = settings_Cov$event.ranges[EventName == input$events_cov]

      # change settings on input based on this
      updateSelectInput(session = session, inputId = "chr_cov", 
        selected = events_id_view$seqnames[1])
      
      # default is zoom out by factor of 1
      span = events_id_view$end[1] - events_id_view$start[1]
      view_start = max(1, events_id_view$start[1] - span)
      view_end = view_start + 3 * span
      updateTextInput(session = session, inputId = "start_cov", 
        value = view_start)
      updateTextInput(session = session, inputId = "end_cov", 
        value = view_end)
    })		
		observeEvent(input$genes_cov, {
      req(input$genes_cov)
      req(input$genes_cov != "(none)")

      gene_id_view = settings_Cov$gene_list[gene_display_name == input$genes_cov]

      # change settings on input based on this
      updateSelectInput(session = session, inputId = "chr_cov", 
        selected = gene_id_view$seqnames[1])
      updateTextInput(session = session, inputId = "start_cov", 
        value = gene_id_view$start[1])
      updateTextInput(session = session, inputId = "end_cov", 
        value = gene_id_view$end[1])
    })		

    observeEvent(input$chr_cov, {
      req(input$chr_cov != "(none)")
			seqInfo = settings_Cov$seqInfo[input$chr_cov]
			seqmax = GenomeInfoDb::seqlengths(seqInfo)
      
			req(as.numeric(input$end_cov))
      if(as.numeric(input$end_cov) > seqmax) {
        updateTextInput(session = session, inputId = "end_cov", 
          value = seqmax)      
        req(as.numeric(input$start_cov))
        if(seqmax - as.numeric(input$start_cov) < 50) {
          updateTextInput(session = session, inputId = "end_cov", 
            value = seqmax - 50)        
        }
      }
    })
   observeEvent(input$start_cov, {
			req(as.numeric(input$start_cov))
			req(as.numeric(input$end_cov))
			req(as.numeric(input$end_cov) - as.numeric(input$start_cov) > 50)

			# adjust zoom
			span = as.numeric(input$end_cov) - as.numeric(input$start_cov)
			cur_zoom = floor(log(span/50) / log(3))
			output$label_zoom_cov <- renderText({16 - cur_zoom})
    })
    observeEvent(input$end_cov, {
			req(as.numeric(input$start_cov))
			req(as.numeric(input$end_cov))		
			req(as.numeric(input$end_cov) - as.numeric(input$start_cov) > 50)

			# adjust zoom
			span = as.numeric(input$end_cov) - as.numeric(input$start_cov)
			cur_zoom = floor(log(span/50) / log(3))
			output$label_zoom_cov <- renderText({16 - cur_zoom})
    })
    
  # DE events row selection:
  observeEvent(input$select_events_cov, {
    # Populate events
    req(input$DT_DE_rows_all)
    req(settings_DE$res)
    
    if(input$select_events_cov == "Highlighted") {
      selected = input$DT_DE_rows_selected
    } else if(input$select_events_cov == "Top N Filtered Results") {
      selected = input$DT_DE_rows_all
      if(length(selected) > input$slider_num_events_cov) {
        selected = selected[seq_len(input$slider_num_events_cov)]
      }
    } else {
      selected = seq_len(min(input$slider_num_events_cov, nrow(settings_DE$res)))
    }

    if(length(selected) > 0 & is_valid(settings_DE$res)) {
      updateSelectizeInput(session = session, inputId = "events_view", server = TRUE,
        choices = c("(none)", settings_DE$res$EventName[selected]), selected = "(none)")    								
      updateSelectizeInput(session = session, inputId = "events_cov", server = TRUE,
        choices = c("(none)", settings_DE$res$EventName[selected]), selected = "(none)")    								
    } else {
      updateSelectizeInput(session = session, inputId = "events_view", server = TRUE,
        choices = c("(none)"), selected = "(none)")    								
      updateSelectizeInput(session = session, inputId = "events_cov", server = TRUE,
        choices = c("(none)"), selected = "(none)")    								    
    }
  })
    
    
# End of server function		
  }

  runApp(shinyApp(ui, server))

}

# Temp Utilities

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
	} else if(method == "Z-score") {
		mat = inc/(inc + exc)
		mat[inc + exc < depth_threshold] = NA
    mat = mat - rowMeans(mat)
    mat = mat / matrixStats::rowSds(mat)
		return(mat)
	}
	
}

make_diagonal <- function(se, event_list, condition, nom_DE, denom_DE, depth_threshold = 10, logit_max = 5) {

	inc = SummarizedExperiment::assay(se, "Included")[event_list, ]
	exc = SummarizedExperiment::assay(se, "Excluded")[event_list, ]
	mat = inc/(inc + exc)
	mat[inc + exc < depth_threshold] = NA

	# use logit method to calculate geometric mean

	mat.nom = logit(mat[, SummarizedExperiment::colData(se)[,condition] == nom_DE])
	mat.denom = logit(mat[, SummarizedExperiment::colData(se)[,condition] == denom_DE])
	
	mat.nom[mat.nom > logit_max] = logit_max
	mat.denom[mat.denom > logit_max] = logit_max
	mat.nom[mat.nom < -logit_max] = -logit_max
	mat.denom[mat.denom < -logit_max] = -logit_max
			
	df = data.frame(EventName = event_list, nom = inv.logit(rowMeans(mat.nom, na.rm = TRUE)),
		denom = inv.logit(rowMeans(mat.denom, na.rm = TRUE)))
	
	return(df)
}

update_data_frame <- function(existing_df, new_df) {
	# add extra samples to existing df
	DT1 = as.data.table(existing_df)
	DT2 = as.data.table(new_df)

  common_cols = intersect(names(DT1)[-1], names(DT2)[-1])
  new_cols = names(DT2)[!(names(DT2) %in% names(DT1))]

  if(!all(DT2$sample %in% DT1$sample)) {
    DT_add = DT2[!(sample %in% DT1$sample)]
    if(length(new_cols) > 0) DT_add = DT_add[, c(new_cols) := NULL]
    newDT = rbind(DT1, DT_add, fill = TRUE)
  } else {
    newDT = copy(DT1)
  }

  if(length(new_cols) > 0) {
    DT_tomerge = copy(DT2)
    if(length(common_cols) > 0) {
      DT_tomerge[, c(common_cols) := NULL]
    }
    newDT = merge(newDT, DT_tomerge, all = TRUE, by = "sample")
  }

  # now update conflicting values
  if(length(common_cols) > 0 & any(DT2$sample %in% DT1$sample)) {
    DT_toupdate = DT2[(sample %in% DT1$sample)]
    if(length(new_cols) > 0) DT_toupdate = DT_toupdate[, c(new_cols) := NULL]

    newDT[DT_toupdate, on=.(sample), (common_cols) := mget(paste0("i.", common_cols))]
  }
  return(as.data.frame(newDT))
}







