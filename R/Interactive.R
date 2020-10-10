
#' @export
startNxtIRF <- function(offline = FALSE) {

	assertthat::assert_that(interactive(),
		msg = "NxtIRF App can only be run in interactive mode (i.e. RStudio).")

  ah = AnnotationHub::AnnotationHub(localHub = offline)

	ui <- navbarPage("NxtIRF", id = "navSelection",
# Title Page
		tabPanel("About", value = "navTitle",
			img(src="https://pbs.twimg.com/profile_images/1310789966293655553/7HawCItY_400x400.jpg")
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
							c("", "hg38", "mm10", "hg19", "mm9", "other")),
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
			tabPanel("Load", value = "navRef_Load",
				fluidRow(
					column(9,
						h4("Select Reference Directory"),
						shinyDirButton("dir_reference_path_load", label = "Choose reference path", title = "Choose reference path"),
							textOutput("txt_reference_path_load"),
						br(),
						actionButton("loadRef", "Load Reference"),
						actionButton("clearLoadRef", "Clear settings"),
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
			)
		),
		
		tabPanel("Experiment", value = "navExpr",
			fluidRow(
				column(4,
					textOutput("txt_reference_path_expr"),
					br(),
					
					shinyDirButton("dir_bam_path_load", 
						label = "Choose BAM path", title = "Choose BAM path"),
					textOutput("txt_bam_path_expr"),
					br(),
					
					shinyDirButton("dir_irf_path_load", 
						label = "Choose IRFinder output path", title = "Choose IRFinder output path"),					
					textOutput("txt_irf_path_expr"),
					br(),
					
					actionButton("run_irf_expr", "Run IRFinder on selected bam files"),
					br(),
					
					shinyFilesButton("file_expr_path_load", label = "Choose Sample Annotation Table", 
						title = "Choose Sample Annotation Table", multiple = FALSE),
					textOutput("txt_sample_anno_expr"),
					br(),						
					
					wellPanel(
						h5("Add annotation column"),
						uiOutput("newcol_expr"),
						radioButtons("type_newcol_expr", "Type", c("integer", "double", "character")),
						actionButton("addcolumn_expr", "Add"), actionButton("removecolumn_expr", "Remove")
					),
					
					actionButton("run_collate_expr", "Compile Experiment"),
					br(),					
				),
				column(8,
					rHandsontableOutput("hot_expr")
				)	# last column
			),
		),
		
		navbarMenu("Analyze"

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

    default_volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
    addit_volume <- reactive({
      req(input$navSelection)
      if(input$navSelection == "navRef_New") {
        req(settings_newref$newref_path)
        return(c(Ref = settings_newref$newref_path))        
      } else if(input$navSelection == "navExpr") {
        req(settings_expr$expr_path)
        return(c(Ref = settings_expr$expr_path))              
      }
    })
    
		observe({  
			shinyDirChoose(input, "dir_reference_path", roots = c(addit_volume, default_volumes), session = session)
			output$txt_reference_path <- renderText({
					validate(need(input$dir_reference_path, "Please select reference path"))
					settings_newref$newref_path = parseDirPath(c(addit_volume, default_volumes), input$dir_reference_path)
			})
		})
		observe({
			shinyFileChoose(input, "file_genome", roots = c(addit_volume, default_volumes), 
				session = session, filetypes = c("fa", "fasta", "gz"))
			if(!is.null(input$file_genome)){
			 file_selected<-parseFilePaths(c(addit_volume, default_volumes), input$file_genome)
			 settings_newref$newref_fasta = as.character(file_selected$datapath)
			 output$txt_genome <- renderText(as.character(file_selected$datapath))
			}
		})
		observe({  
			shinyFileChoose(input, "file_gtf", roots = c(addit_volume, default_volumes), 
				session = session, filetypes = c("gtf", "gz"))
			if(!is.null(input$file_gtf)){
			 file_selected<-parseFilePaths(c(addit_volume, default_volumes), input$file_gtf)
			 settings_newref$newref_gtf = as.character(file_selected$datapath)
			 output$txt_gtf <- renderText(as.character(file_selected$datapath))
			}
		})
		observe({  
		shinyFileChoose(input, "file_mappa", roots = c(addit_volume, default_volumes), 
			session = session, filetypes = c("txt", "gz"))
			if(!is.null(input$file_mappa)){
			 file_selected<-parseFilePaths(c(addit_volume, default_volumes), input$file_mappa)
			 settings_newref$newref_mappa = as.character(file_selected$datapath)
			 output$txt_mappa <- renderText(as.character(file_selected$datapath))
			}
		})
		observe({  
			shinyFileChoose(input, "file_NPA", roots = c(addit_volume, default_volumes), 
				session = session, filetypes = c("bed", "txt", "gz"))
			if(!is.null(input$file_NPA)){
				file_selected<-parseFilePaths(c(addit_volume, default_volumes), input$file_NPA)
				settings_newref$newref_NPA = as.character(file_selected$datapath)
				output$txt_NPA <- renderText(as.character(file_selected$datapath))
			}
		})
		observe({  
			shinyFileChoose(input, "file_bl", roots = c(addit_volume, default_volumes), 
				session = session, filetypes = c("bed", "txt", "gz"))
			if(!is.null(input$file_bl)){
			 file_selected<-parseFilePaths(c(addit_volume, default_volumes), input$file_bl)
			 settings_newref$newref_bl = as.character(file_selected$datapath)
			 output$txt_bl <- renderText(as.character(file_selected$datapath))
			}
		})
		observeEvent(input$newref_genome_type, {
			if(input$newref_genome_type == "hg38") {
				settings_newref$newref_NPA = system.file("extra-input-files/Human_hg38_nonPolyA_ROI.bed", package = "NxtIRF")
			} else if(input$newref_genome_type == "hg19")  {
				settings_newref$newref_NPA = system.file("extra-input-files/Human_hg19_nonPolyA_ROI.bed", package = "NxtIRF")
			} else if(input$newref_genome_type == "mm10")  {
				settings_newref$newref_NPA = system.file("extra-input-files/Mouse_mm10_nonPolyA_ROI.bed", package = "NxtIRF")
			} else if(input$newref_genome_type == "mm9")  {
				settings_newref$newref_NPA = system.file("extra-input-files/Mouse_mm9_nonPolyA_ROI.bed", package = "NxtIRF")
			} else if(input$newref_genome_type == "") {
				settings_newref$newref_NPA = ""
			} else {

			}
			output$txt_NPA <- renderText(settings_newref$newref_NPA)
		})

		observeEvent(input$navSelection, {
			if(input$navSelection == "navRef_New") {
				ah.filtered = ah[ah$dataprovider == "Ensembl"]
				ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
				ah.filtered = ah.filtered[ah.filtered$sourcetype == "GTF"]
				updateSelectInput(session = session, inputId = "newrefAH_Species", 
					choices = c("", sort(unique(ah.filtered$species))))
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
				genome_type = "other", nonPolyARef = settings_newref$newref_NPA, MappabilityRef = settings_newref$newref_mappa,
				BlacklistRef = settings_newref$newref_bl)

			is_valid <- function(.) !is.null(.) && . != ""
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
				data.table::fwrite(args.df, paste(args$reference_path, "settings_newref.csv", sep="/"), row.names = TRUE)
				if("ah_genome_tmp" %in% names(args)) {
					args$ah_genome = data.table::tstrsplit(args$ah_genome_tmp, split=":", fixed=TRUE)[[1]]
					args$ah_genome_tmp = NULL
				}
				if("ah_gtf_tmp" %in% names(args)) {
					args$ah_transcriptome = data.table::tstrsplit(args$ah_gtf_tmp, split=":", fixed=TRUE)[[1]]
					args$ah_gtf_tmp = NULL
				}
				do.call(BuildReference, args)
			}
		})
		
		# clearNewRef Button
		observeEvent(input$clearNewRef, {
			newref_path = getwd()
			newref_fasta = ""
			newref_gtf = ""
			newref_AH_fasta = ""
			newref_AH_gtf = ""
			newref_mappa = ""
			newref_NPA = ""
			newref_bl = ""
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
		shinyDirChoose(input, "dir_reference_path_load", roots = c(addit_volume, default_volumes), session = session)
		observeEvent(input$dir_reference_path_load,{  
			output$txt_reference_path_load <- renderText({
					validate(need(input$dir_reference_path_load, "Please select reference path"))
					settings_loadref$loadref_path = parseDirPath(c(addit_volume, default_volumes), input$dir_reference_path_load)
			})
    })
		observeEvent(settings_loadref$loadref_path,{ 
      req(settings_loadref$loadref_path)
			if(file.exists(paste(settings_loadref$loadref_path, "settings.csv", sep="/"))) {
        settings_loadref$settings = as.data.frame(fread(
          paste(settings_loadref$loadref_path, "settings.csv", sep="/"), header = TRUE
          ))
        colnames(settings_loadref$settings) = c("fields", "value")
        if("reference_path" %in% settings_loadref$settings$fields) {
          if("ah_genome_tmp" %in% settings_loadref$settings$fields) {
            output$loadRef_field1 <- renderText({
              paste("AnnotationHub genome:\n",
              settings_loadref$settings$value[which(settings_loadref$settings$fields == "ah_genome_tmp")]
              )
            })
          }
          if("ah_gtf_tmp" %in% settings_loadref$settings$fields) {
            output$loadRef_field2 <- renderText({
              paste("AnnotationHub gene annotations:\n",
              settings_loadref$settings$value[which(settings_loadref$settings$fields == "ah_gtf_tmp")]
              )
            })
          }
          if("fasta" %in% settings_loadref$settings$fields) {
            output$loadRef_field3 <- renderText({
              paste("Genome FASTA file (user):\n",
              settings_loadref$settings$value[which(settings_loadref$settings$fields == "fasta")]
              )
            })					
          }
          if("gtf" %in% settings_loadref$settings$fields) {
            output$loadRef_field4 <- renderText({
              paste("Annotation GTF file (user):\n",
              settings_loadref$settings$value[which(settings_loadref$settings$fields == "gtf")]
              )
            })					
          }
          if("MappabilityRef" %in% settings_loadref$settings$fields) {
            output$loadRef_field5 <- renderText({
              paste("Mappability Exclusions:\n",
              settings_loadref$settings$value[which(settings_loadref$settings$fields == "MappabilityRef")]
              )
            })					
          }
          if("nonPolyARef" %in% settings_loadref$settings$fields) {
            output$loadRef_field6 <- renderText({
              paste("Non-polyA Ref:\n",
              settings_loadref$settings$value[which(settings_loadref$settings$fields == "nonPolyARef")]
              )
            })					
          }
          if("BlacklistRef" %in% settings_loadref$settings$fields) {
            output$loadRef_field7 <- renderText({
              paste("Blacklist file:\n",
              settings_loadref$settings$value[which(settings_loadref$settings$fields == "BlacklistRef")]
              )
            })					
          }
        }
			} else {
          output$loadRef_field1 <- renderText({
            paste(paste(settings_loadref$loadref_path, "settings.csv", sep="/"), "not found")
          })
      }
		})
		
# Design Experiment page
		settings_expr <- shiny::reactiveValues(
			expr_path = "",
			bam_path = "",
			irf_path = "",
			collate_path = "",
			DT = c()
		)
		
		## Handsontable auto-updates settings_expr$DT on user edit
    observe({
      if (!is.null(input$hot_expr)) {
        DT = as.data.table(hot_to_r(input$hot_expr))
      } else {
        if (is.null(settings_expr$DT))
          DT <- DT
        else
          DT <- settings_expr$DT
      }
      settings_expr$DT <- DT
    })
		output$hot_expr <- renderRHandsontable({
			DF <- as.data.frame(settings_expr$DT)
			if (!is.null(DF)) {
				rhandsontable(DF, useTypes = TRUE, stretchH = "all")
			}
		})		
    observe({  
      shinyDirChoose(input, "dir_bam_path_load", roots = c(addit_volume, default_volumes), session = session)
			output$txt_bam_path_expr <- renderText({
					validate(need(input$dir_bam_path_load, "Please select path where BAMs are kept"))
          settings_expr$expr_path = dirname(parseDirPath(c(addit_volume, default_volumes), input$dir_bam_path_load))
					settings_expr$bam_path = parseDirPath(c(addit_volume, default_volumes), input$dir_bam_path_load)
			})
    })
		observeEvent(settings_expr$bam_path,{
      req(settings_expr$bam_path)
		# First assume bams are named by subdirectory names
			temp.df = as.data.table(FindSamples(
				settings_expr$bam_path, suffix = ".bam", use_subdir = TRUE))
			if(nrow(temp.df) > 0) {
				if(length(unique(temp.df$sample)) == nrow(temp.df)) {
					# Assume subdirectory names designate sample names
				} else {
					temp.df = as.data.table(FindSamples(
						settings_expr$bam_path, suffix = ".bam", use_subdir = FALSE))
					if(length(unique(temp.df$sample)) == nrow(temp.df)) {
				# Else assume bam names designate sample names					
					} else {
						output$txt_bam_path_expr <- renderText("BAM file names (or its path names) must be unique")							
						settings_expr$bam_path = ""
						temp.df = NULL
					}
				}
			} else {
				output$txt_bam_path_expr <- renderText("No bam files found in given path")
				settings_expr$bam_path = ""
				temp.df = NULL
			}
			
		# compile experiment df with bam paths
			if(!is.null(temp.df)) {
					colnames(temp.df)[2] = "bam_file"
					temp.DT = as.data.table(temp.df)
				if(!is.null(settings_expr$DT)) {
			# merge with existing dataframe	
					settings_expr$DT = rbind(settings_expr$DT, temp.DT[!(sample %in% settings_expr$DT$sample)],
							fill = TRUE) # Add samples not in original DT
					settings_expr$DT[temp.DT, on = "sample", bam_file := i.bam_file] # Update new bam paths
				} else {
			# start anew
					settings_expr$DT = data.table(sample = temp.df$sample,
						bam_file = "", irf_file = "", collated_file = "")
					settings_expr$DT[temp.DT, on = "sample", bam_file := i.bam_file] # Update new bam paths
				}			
			}
			# output$hot_expr <- renderRHandsontable({
				# DF <- as.data.frame(settings_expr$DT)
				# if (!is.null(DF)) {
					# rhandsontable(DF, useTypes = TRUE, stretchH = "all")
				# }
			# })
		})
    observe({
      shinyDirChoose(input, "dir_irf_path_load", roots = c(addit_volume, default_volumes), 
        session = session)
      output$txt_irf_path_expr <- renderText({
          validate(need(input$dir_irf_path_load, "Please select path where IRFinder output should be kept"))
          settings_expr$irf_path = parseDirPath(c(addit_volume, default_volumes), 
            input$dir_irf_path_load)
      })        
    })
		observeEvent(settings_expr$irf_path,{
      req(settings_expr$irf_path)
				# merge irfinder paths
			temp.df = as.data.table(FindSamples(
				settings_expr$irf_path, suffix = ".txt.gz", use_subdir = FALSE))
			if(nrow(temp.df) > 0) {
				if(length(unique(temp.df$sample)) == nrow(temp.df)) {
					# Assume output names designate sample names
				} else {
					temp.df = as.data.table(FindSamples(
						settings_expr$irf_path, suffix = ".txt.gz", use_subdir = TRUE))
					if(length(unique(temp.df$sample)) == nrow(temp.df)) {
				# Else assume subdirectory names designate sample names					
					} else {
						output$txt_irf_path_expr <- renderText("IRFinder file names (or its path names) must be unique")							
						settings_expr$irf_path = ""
						temp.df = NULL
					}
				}
			} else {
				output$txt_irf_path_expr <- renderText("No IRFinder files found in given path")
				settings_expr$irf_path = ""
				temp.df = NULL
			}
			
		# compile experiment df with irfinder paths
			if(!is.null(temp.df)) {
					colnames(temp.df)[2] = "irf_file"
					temp.DT = as.data.table(temp.df)
				if(!is.null(settings_expr$DT)) {
			# merge with existing dataframe	
					settings_expr$DT = rbind(settings_expr$DT, temp.DT[!(sample %in% settings_expr$DT$sample)],
							fill = TRUE) # Add samples not in original DT
					settings_expr$DT[temp.DT, on = "sample", irf_file := i.irf_file] # Update new bam paths
				} else {
			# start anew
					settings_expr$DT = data.table(sample = temp.df$sample,
						bam_file = "", irf_file = "", collated_file = "")
					settings_expr$DT[temp.DT, on = "sample", irf_file := i.irf_file] # Update new bam paths
				}			
			}
		})
		

    output$newcol_expr <- renderUI({
      textInput("newcolumnname_expr", "Name", sprintf("newcol%s", 1+ncol(settings_expr$DT)))
    })
 		# Add column
		observeEvent(input$addcolumn_expr, {
      DT <- isolate(settings_expr$DT)
      newcolumn <- eval(parse(text=sprintf('%s(nrow(DF))', isolate(input$type_newcol_expr))))
      settings_expr$DT <- setNames(cbind(DF, newcolumn, stringsAsFactors=FALSE), 
				c(names(DT), isolate(input$newcolumnname_expr)))
    })
 		# Remove column
		observeEvent(input$removecolumn_expr, {
      DT <- isolate(settings_expr$DT)
			if(isolate(input$newcolumnname_expr) %in% colnames(DT)) {
				DT[, c(input$newcolumnname_expr) := NULL]
				settings_expr$DT = DT
			}
    })	
# End of server function		
  }

  runApp(shinyApp(ui, server))

}

