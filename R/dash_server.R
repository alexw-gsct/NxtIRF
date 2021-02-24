dash_server = function(input, output, session) {

    default_volumes <- c("Working Directory" = getwd(), "Home" = "~", getVolumes()())
    addit_volume = c()

    settings_system <- setreactive_system()
    settings_newref <- setreactive_newref()
    settings_loadref <- setreactive_loadref()
    settings_expr <- setreactive_expr()
    settings_SE <- setreactive_SE()
    settings_DE <- setreactive_DE()
    settings_Diag = setreactive_Diag()
    settings_Volc = setreactive_Diag()
    settings_Heat = setreactive_Diag()
    settings_Cov <- setreactive_Cov()
        
        # settings_SaveObj <- reactiveValues(
            # obj = NULL
        # )

        # settings_to_file <- function(filename) {
            # final = list(
                # settings_system <- isolate(reactiveValuesToList(settings_system)),
                # settings_newref <- isolate(reactiveValuesToList(settings_newref)),
                # settings_loadref <- isolate(reactiveValuesToList(settings_loadref)),
                # settings_expr <- isolate(reactiveValuesToList(settings_expr)),
                # settings_SE <- isolate(reactiveValuesToList(settings_SE)),
                # settings_DE <- isolate(reactiveValuesToList(settings_DE)),
                # settings_Diag <- isolate(reactiveValuesToList(settings_Diag)),
                # settings_Volc <- isolate(reactiveValuesToList(settings_Volc)),
                # settings_Cov <- isolate(reactiveValuesToList(settings_Cov))
            # )
            # names(final) = c("settings_system","settings_newref","settings_loadref",
                # "settings_expr","settings_SE","settings_DE",
                # "settings_Diag","settings_Volc","settings_Cov")
            # saveRDS(final, filename)
        # }

        # settings_from_file <- function(filename) {
            # if(file.exists(filename)) {
                # final <- readRDS(filename)
                # if("settings_system" %in% names(final)) settings_system <<- do.call(reactiveValues, final$settings_system)
                # if("settings_newref" %in% names(final)) settings_newref <<- do.call(reactiveValues, final$settings_newref)
                # if("settings_loadref" %in% names(final)) settings_loadref <<- do.call(reactiveValues, final$settings_loadref)
                # if("settings_expr" %in% names(final)) settings_expr <<- do.call(reactiveValues, final$settings_expr)
                # if("settings_SE" %in% names(final)) settings_SE <<- do.call(reactiveValues, final$settings_SE)
                # if("settings_DE" %in% names(final)) settings_DE <<- do.call(reactiveValues, final$settings_DE)
                # if("settings_Diag" %in% names(final)) settings_Diag <<- do.call(reactiveValues, final$settings_Diag)
                # if("settings_Volc" %in% names(final)) settings_Volc <<- do.call(reactiveValues, final$settings_Volc)
                # if("settings_Cov" %in% names(final)) settings_Cov <<- do.call(reactiveValues, final$settings_Cov)
            # }
        # }

        # shinyFileSave(input, "file_savestate", 
            # roots = c(default_volumes, addit_volume), session = session,
            # filetype = list(RDS = "Rds"))
        # observeEvent(input$file_savestate, {
            # selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$file_savestate)
            # req(selectedfile$datapath)
            # settings_to_file(selectedfile$datapath)
        # })

        # shinyFileChoose(input, "file_loadstate", 
            # roots = c(default_volumes, addit_volume), session = session)
		# observeEvent(input$file_loadstate, {
            # req(input$file_loadstate)
            # file_selected<-parseFilePaths(c(default_volumes, addit_volume), input$file_loadstate)
            # req(file_selected$datapath)
            # settings_from_file(as.character(file_selected$datapath))
		# })
        
        # observeEvent(settings_system$n_threads, {
            # updateSliderInput(session = session, 
                # inputId = "cores_slider", value = settings_system$n_threads)
            # updateNumericInput(session = session, 
                # inputId = "cores_numeric", value = settings_system$n_threads)
        # })
        
    settings_app <- reactiveValues(
        ah = NULL
    )
    initialize_ah <- function() {
        if(!is_valid(settings_app$ah)) {
            settings_app$ah = .ah
        }        
    }
# tabEvent Observer
    observeEvent(input$navSelection, {
        if(input$navSelection == "navTitle") {
            # Initialize App
            initialize_ah()
        } else if(input$navSelection == "navRef_New") {
            initialize_ah()
            ah = settings_app$ah
            # autopopulate if previous settings detected
            
            ah.filtered = ah[ah$dataprovider == "Ensembl"]
            ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
            ah.filtered = ah.filtered[ah.filtered$sourcetype == "GTF"]
            choices = c("", sort(unique(ah.filtered$species)))

            if(is_valid(settings_newref$ui_newrefAH_Species) &&
                settings_newref$ui_newrefAH_Species %in% choices) {
                updateSelectInput(session = session, inputId = "newrefAH_Species",
                    choices = choices, selected = settings_newref$ui_newrefAH_Species)				
            } else {
                updateSelectInput(session = session, inputId = "newrefAH_Species",
                    choices = choices)	        
            }
            
            choices = c("(not specified)", "hg38", "mm10", "hg19", "mm9")
            if(is_valid(settings_newref$ui_newref_genome_type) &&
                settings_newref$ui_newref_genome_type %in% choices) {
                updateSelectInput(session = session, inputId = "newref_genome_type",
                    choices = choices, selected = settings_newref$ui_newref_genome_type)				
            } else {
                updateSelectInput(session = session, inputId = "newref_genome_type",
                    choices = choices)	        
            }
            
            
            if(settings_newref$newref_path != ""){
                output$txt_reference_path <- renderText(settings_newref$newref_path)
            }
        } else if(input$navSelection == "navSystem") {

        } else if(input$navSelection == "navExpr") {
            output$se_expr_infobox <- renderUI({
                ui_infobox_expr(ifelse(
                    is_valid(settings_SE$se), 2,
                    ifelse(
                        is_valid(settings_expr$collate_path) &&
                        file.exists(file.path(
                            settings_expr$collate_path,
                            "colData.Rds"
                        ))
                    ,1,0)
                ))
            })
        } else if(input$navSelection == "navQC") {
            if(file.exists(file.path(settings_expr$collate_path, "stats.fst"))) {
                settings_SE$QC = 
                    as.data.table(read.fst((file.path(settings_expr$collate_path, "stats.fst"))))                
                settings_SE$QC = 
                    merge(as.data.table(settings_expr$df.anno), settings_SE$QC, all = TRUE)
            }
            output$DT_QC <- DT::renderDataTable({
                validate(need(settings_SE$se, "Load Experiment file first"))
                validate(need(file.exists(
                    file.path(settings_expr$collate_path, "stats.fst")
                    ), "stats.fst does not exist in given NxtIRF output directory"))
                DT::datatable(
                    as.data.frame(settings_SE$QC),
                    class = 'cell-border stripe',
                    rownames = settings_SE$QC$sample,
                    filter = 'top'
                )
            })
            if(is_valid(settings_SE$QC)) {
                choices = colnames(settings_SE$QC)
                choices = choices[!(choices %in% colnames(settings_expr$df.anno))]
                choices = choices[!(choices %in% 
                    c("paired", "strand", "path")
                )]
                choices = c("(none)", choices)
                updateSelectInput(session = session, inputId = "QC_xaxis",
                    choices = choices)	        
                updateSelectInput(session = session, inputId = "QC_yaxis",
                    choices = choices)	        
            } else {
                 updateSelectInput(session = session, inputId = "QC_xaxis",
                    choices = "(none)")	        
                updateSelectInput(session = session, inputId = "QC_yaxis",
                    choices = "(none)")	               
            }
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
        } else if(input$navSelection == "navHeatmap") {
            if(is_valid(settings_SE$se)) {
                updateSelectInput(session = session, inputId = "anno_col_heat", 
                    choices = colnames(SummarizedExperiment::colData(settings_SE$se)), selected = NULL)
            }
        } else if(input$navSelection == "navCoverage") {
            initialize_ah()
            ah = settings_app$ah
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
                genome = rtracklayer::TwoBitFile(
                    file.path(settings_loadref$loadref_path, "resource", "genome.2bit"))
            }
            settings_Cov$seqInfo = seqinfo(genome)
            settings_Cov$gene_list <- getGeneList(settings_loadref$loadref_path)
            settings_Cov$elem.DT <- loadViewRef(settings_loadref$loadref_path)
            settings_Cov$transcripts.DT <- loadTranscripts(settings_loadref$loadref_path)
            
            # Populate events here

            if(!is.null(settings_Cov$gene_list)) {
                # Only refresh if new reference added
                if(settings_Cov$loaded_reference != settings_loadref$loadref_path) {
                    message(paste("Populating drop-down box with", 
                        length(unique(settings_Cov$gene_list$gene_display_name)),"genes"))
                    updateSelectInput(session = session, inputId = "chr_cov", 
                        choices = c("(none)", sort(unique(settings_Cov$gene_list$seqnames))), selected = "(none)")    								          
                    updateSelectizeInput(session = session, inputId = "genes_cov", server = TRUE,
                        choices = c("(none)", settings_Cov$gene_list$gene_display_name), selected = "(none)")
                    settings_Cov$loaded_reference = settings_loadref$loadref_path
                }
            } else {
                updateSelectInput(session = session, inputId = "chr_cov", 
                    choices = c("(none)"), selected = "(none)")    								
                updateSelectizeInput(session = session, inputId = "genes_cov", server = TRUE,
                    choices = c("(none)"), selected = "(none)") 
            }
            
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
                    if(is_valid(input$events_cov)) {
                        selected_event = input$events_cov
                        if(!(selected_event %in% settings_DE$res$EventName[selected])) {
                            selected_event = "(none)"
                        }
                    } else {
                        selected_event = "(none)"
                    }
                    updateSelectizeInput(session = session, inputId = "events_cov", server = TRUE,
                        choices = c("(none)", settings_DE$res$EventName[selected]), selected = selected_event)    								
                } else {
                    updateSelectizeInput(session = session, inputId = "events_cov", server = TRUE,
                        choices = c("(none)"), selected = "(none)")    								    
                }
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
            
                DT.files = as.data.table(settings_expr$df.files[, c("sample", "cov_file")])
                DT.files = na.omit(DT.files)
                settings_Cov$avail_cov = DT.files$cov_file
                names(settings_Cov$avail_cov) = DT.files$sample
                settings_Cov$avail_cov = settings_Cov$avail_cov[file.exists(settings_Cov$avail_cov)]

                refresh_tracks_cov()
            }        
        }
    })

# Function that returns the number of threads to use
    get_threads <- function() {
        if(input$thread_option == "Single-Thread"){
            n_threads = 1
        } else if(input$thread_option == "Multi-Thread (Low)") {
            n_threads = floor( (parallel::detectCores() - 2) / 2)
            if(n_threads < 1) n_threads = 1
        } else if(input$thread_option == "Multi-Thread (High)") {
            n_threads = parallel::detectCores() - 2
            if(n_threads < 1) n_threads = 1
        } else if(input$thread_option == "Custom") {
            n_threads = input$cores_numeric
        }
        n_threads
    }
		
# TAB: Build Reference
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
        req(input$newref_genome_type)
        resource_path = "https://raw.github.com/alexw-gsct/NxtIRF_resources/main/data"
        
        if(input$newref_genome_type == "hg38") {
            settings_newref$newref_NPA = system.file(
                "extra-input-files/Human_hg38_nonPolyA_ROI.bed", package = "NxtIRF")
            settings_newref$newref_mappa = 
                paste(resource_path, "Mappability_Regions_hg38_v94.txt.gz", sep="/")
        
        } else if(input$newref_genome_type == "hg19")  {
            settings_newref$newref_NPA = system.file(
                "extra-input-files/Human_hg19_nonPolyA_ROI.bed", package = "NxtIRF")
            settings_newref$newref_mappa = 
                paste(resource_path, "Mappability_Regions_hg19_v75.txt.gz", sep="/")
        
        } else if(input$newref_genome_type == "mm10")  {
            settings_newref$newref_NPA = system.file(
                "extra-input-files/Mouse_mm10_nonPolyA_ROI.bed", package = "NxtIRF")
            settings_newref$newref_mappa = 
                paste(resource_path, "Mappability_Regions_mm10_v94.txt.gz", sep="/")
        
        } else if(input$newref_genome_type == "mm9")  {
            settings_newref$newref_NPA = system.file(
                "extra-input-files/Mouse_mm9_nonPolyA_ROI.bed", package = "NxtIRF")
            settings_newref$newref_mappa = 
                paste(resource_path, "Mappability_Regions_mm9_v67.txt.gz", sep="/")
        
        } else if(input$newref_genome_type == "(not specified)") {
    # do nothing. This allows user to first select the default and then change to user-defined files
        } else {
            settings_newref$newref_NPA = ""
            settings_newref$newref_mappa = ""
        }
        settings_newref$ui_newref_genome_type = input$newref_genome_type
    })

    observeEvent(input$newrefAH_Species, {
        initialize_ah()
        ah = settings_app$ah
        if(input$newrefAH_Species != "") {
            ah.filtered = ah[ah$dataprovider == "Ensembl"]
            ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
            ah.filtered = ah.filtered[ah.filtered$sourcetype == "GTF"]
            ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]

            queryfirst = tstrsplit(ah.filtered[1]$sourceurl, split="/", fixed=TRUE)
            query_index = which(vapply(queryfirst, function(x) grepl("release", x), logical(1)))
            choices = unlist(unique(tstrsplit(ah.filtered$sourceurl, split="/", fixed=TRUE)[query_index]))
            choices = c("", sort(choices))
            # updateSelectInput(session = session, inputId = "newrefAH_Version_Trans", 
                # choices = c("", sort(choices)))
            update_select_without_clearing(session, "newrefAH_Version_Trans", choices, input)
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
        settings_newref$ui_newrefAH_Species = input$newrefAH_Species
    })
    observeEvent(input$newrefAH_Version_Trans, {
        initialize_ah()
        ah = settings_app$ah
        if(input$newrefAH_Version_Trans != "") {
            ah.filtered = ah[ah$dataprovider == "Ensembl"]
            ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
            ah.filtered = ah.filtered[ah.filtered$sourcetype == "GTF"]
            ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]
            ah.filtered = ah.filtered[grepl(input$newrefAH_Version_Trans, ah.filtered$sourceurl)]
            choices = c("", paste(names(ah.filtered), basename(ah.filtered$sourceurl), sep=": "))
            # updateSelectInput(session = session, inputId = "newrefAH_Trans", 
                # choices = c("", paste(names(ah.filtered), basename(ah.filtered$sourceurl), sep=": ")))
            update_select_without_clearing(session, "newrefAH_Trans", choices, input)
            
            # Also search for compatible genome
            genomes_avail = ah.filtered$genome
            ah.filtered = ah[ah$dataprovider == "Ensembl"]
            ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
            ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]
            ah.filtered = ah.filtered[ah.filtered$sourcetype == "FASTA"]

            if(any(ah.filtered$genome %in% genomes_avail)) {
                ah.filtered = ah.filtered[ah.filtered$genome %in% genomes_avail]    
                choices = c("", unique(ah.filtered$genome))
                # updateSelectInput(session = session, inputId = "newrefAH_Assembly", 
                        # choices = c("", unique(ah.filtered$genome)))
                update_select_without_clearing(session, "newrefAH_Assembly", choices, input)
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
        settings_newref$ui_newrefAH_Version_Trans = input$newrefAH_Version_Trans
    })
    observeEvent(input$newrefAH_Assembly, {
        initialize_ah()
        ah = settings_app$ah
        if(input$newrefAH_Assembly != "") {
            ah.filtered = ah[ah$dataprovider == "Ensembl"]
            ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
            ah.filtered = ah.filtered[ah.filtered$sourcetype == "FASTA"]
            ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]
            ah.filtered = ah.filtered[ah.filtered$genome == input$newrefAH_Assembly]

            queryfirst = tstrsplit(ah.filtered[1]$sourceurl, split="/", fixed=TRUE)
            query_index = which(vapply(queryfirst, function(x) grepl("release", x), logical(1)))
            choices = unlist(unique(tstrsplit(ah.filtered$sourceurl, split="/", fixed=TRUE)[query_index]))
            choices = c("", sort(choices))
            update_select_without_clearing(session, "newrefAH_Version_Genome", choices, input)
        } else {
            updateSelectInput(session = session, inputId = "newrefAH_Version_Genome", 
                choices = c(""))
            updateSelectInput(session = session, inputId = "newrefAH_Genome", 
                choices = c(""))        
        }
        settings_newref$ui_newrefAH_Assembly = input$newrefAH_Assembly
    })
    observeEvent(input$newrefAH_Version_Genome, {
        initialize_ah()
        ah = settings_app$ah
        if(input$newrefAH_Version_Genome != "") {
            ah.filtered = ah[ah$dataprovider == "Ensembl"]
            ah.filtered = ah.filtered[grepl("release", ah.filtered$sourceurl)]
            ah.filtered = ah.filtered[ah.filtered$sourcetype == "FASTA"]
            ah.filtered = ah.filtered[ah.filtered$species == input$newrefAH_Species]
            ah.filtered = ah.filtered[ah.filtered$genome == input$newrefAH_Assembly]
            ah.filtered = ah.filtered[grepl(input$newrefAH_Version_Genome, ah.filtered$sourceurl)]
            choices = c("", paste(names(ah.filtered), basename(ah.filtered$sourceurl), sep=": "))
            update_select_without_clearing(session, "newrefAH_Genome", choices, input)
        } else {
            updateSelectInput(session = session, inputId = "newrefAH_Genome", 
                choices = c(""))        
        }
        settings_newref$ui_newrefAH_Version_Genome = input$newrefAH_Version_Genome
    })
    observeEvent(input$newrefAH_Genome, {
        if(input$newrefAH_Genome != "") {
            settings_newref$newref_AH_fasta = input$newrefAH_Genome
        } else {
            settings_newref$newref_AH_fasta = ""
        }
        settings_newref$ui_newrefAH_Genome = input$newrefAH_Genome
    })
    observeEvent(input$newrefAH_Trans, {
        if(input$newrefAH_Trans != "") {
            settings_newref$newref_AH_gtf = input$newrefAH_Trans
        } else {
            settings_newref$newref_AH_gtf = ""
        }
        settings_newref$ui_newrefAH_Trans = input$newrefAH_Trans
    })
		
    match_genome_type <- function() {
        if(!is_valid(settings_newref$newref_NPA) | !is_valid(settings_newref$newref_mappa)) {
            return("Interactive")            
        }
        resource_path = "https://raw.github.com/alexw-gsct/NxtIRF_resources/main/data"
        if(
            settings_newref$newref_NPA == system.file(
                "extra-input-files/Human_hg38_nonPolyA_ROI.bed", package = "NxtIRF") &
            settings_newref$newref_mappa == 
                paste(resource_path, "Mappability_Regions_hg38_v94.txt.gz", sep="/")
        ) {
            return("hg38")
        } else if(
            settings_newref$newref_NPA == system.file(
                "extra-input-files/Human_hg19_nonPolyA_ROI.bed", package = "NxtIRF") &
            settings_newref$newref_mappa == 
                paste(resource_path, "Mappability_Regions_hg19_v75.txt.gz", sep="/")
        ) {
            return("hg19")
        } else if(
            settings_newref$newref_NPA == system.file(
                "extra-input-files/Mouse_mm10_nonPolyA_ROI.bed", package = "NxtIRF") &
            settings_newref$newref_mappa == 
                paste(resource_path, "Mappability_Regions_mm10_v94.txt.gz", sep="/")
        ) {
            return("mm10")
        } else if(
            settings_newref$newref_NPA == system.file(
                "extra-input-files/Mouse_mm9_nonPolyA_ROI.bed", package = "NxtIRF") &
            settings_newref$newref_mappa == 
                paste(resource_path, "Mappability_Regions_mm9_v67.txt.gz", sep="/")
        ) {
            return("mm9")
        } else {
            return("Interactive")            
        }
    }
        
	# buildRef Button
    observeEvent(input$buildRef, {
        args <- list(reference_path = settings_newref$newref_path,
            ah_genome_tmp = settings_newref$newref_AH_fasta, ah_gtf_tmp = settings_newref$newref_AH_gtf, 
            fasta = settings_newref$newref_fasta, gtf = settings_newref$newref_gtf,
            genome_type = match_genome_type(), nonPolyARef = settings_newref$newref_NPA, 
            MappabilityRef = settings_newref$newref_mappa,
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
                sendSweetAlert(
                    session = session,
                    title = "Reference Build complete!",
                    type = "success"
                )           
            } else {
                sendSweetAlert(
                    session = session,
                    title = "Reference Build failed. An error must have occurred",
                    type = "error"
                )               
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

# TAB: Load Reference Page

	shinyDirChoose(input, "dir_reference_path_load", roots = c(default_volumes, addit_volume), session = session)
		observeEvent(input$dir_reference_path_load,{  
        req(input$dir_reference_path_load)
        settings_loadref$loadref_path = parseDirPath(c(default_volumes, addit_volume), 
            input$dir_reference_path_load)
    })
    load_ref = function() {
        req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds")))
        initialize_ah()
        ah = settings_app$ah
        
        settings_loadref$settings = readRDS(file.path(settings_loadref$loadref_path, "settings.Rds"))
        if("reference_path" %in% names(settings_loadref$settings)) {
            if("ah_genome" %in% names(settings_loadref$settings)) {
                output$fasta_source_infobox <- renderInfoBox({
                    infoBox(
                        "Genome - AnnotationHub", "",
                        basename(
                            ah$sourceurl[
                                which(names(ah) == settings_loadref$settings$ah_genome)
                            ]
                        ),
                        icon = icon("dna", lib = "font-awesome"),
                        color = "green"
                    )
                })      
            } else if("fasta_file" %in% names(settings_loadref$settings)) {
                output$fasta_source_infobox <- renderInfoBox({
                    infoBox(
                        "Genome - User FASTA",  "",
                        basename(settings_loadref$settings$fasta_file), 
                        icon = icon("dna", lib = "font-awesome"),
                        color = "green"
                    )
                })          
            } else {
                output$fasta_source_infobox <- renderInfoBox(NULL)
            }
            if("ah_transcriptome" %in% names(settings_loadref$settings)) {
                output$gtf_source_infobox <- renderInfoBox({
                    infoBox(
                        "Gene Annotation - AnnotationHub",  "",
                        basename(
                            ah$sourceurl[
                                which(names(ah) == 
                                    settings_loadref$settings$ah_transcriptome)
                            ]
                        ),
                        icon = icon("book-medical", lib = "font-awesome"),
                        color = "orange"
                    )
                })               
            } else if("gtf_file" %in% names(settings_loadref$settings)) {
                output$gtf_source_infobox <- renderInfoBox({
                    infoBox(
                        "Gene Annotation - User GTF",  "",
                        basename(settings_loadref$settings$gtf_file), 
                        icon = icon("book-medical", lib = "font-awesome"),
                        color = "orange"
                    )
                })
            } else {
                output$gtf_source_infobox <- renderInfoBox(NULL)
            }
            if("MappabilityRef" %in% names(settings_loadref$settings)) {
                output$mappa_source_infobox <- renderInfoBox({
                    infoBox(
                        "Mappability",  "",
                        basename(settings_loadref$settings$MappabilityRef), 
                        icon = icon("map", lib = "font-awesome"),
                        color = "blue"
                    )
                })  			
            } else {
                output$mappa_source_infobox <- renderInfoBox(NULL)
            }
            if("nonPolyARef" %in% names(settings_loadref$settings)) {
                output$NPA_source_infobox <- renderInfoBox({
                    infoBox(
                        "Non-PolyA",  "",
                        basename(settings_loadref$settings$nonPolyARef), 
                        icon = icon("font", lib = "font-awesome"),
                        color = "purple"
                    )
                })  			
            } else {
                output$NPA_source_infobox <- renderInfoBox(NULL)
            }
            if("BlacklistRef" %in% names(settings_loadref$settings)) {
                output$BL_source_infobox <- renderInfoBox({
                    infoBox(
                        "BlackList",  "",
                        basename(settings_loadref$settings$BlacklistRef), 
                        icon = icon("list-alt", lib = "font-awesome"),
                        color = "red"
                    )
                })  			
            } else {
                output$BL_source_infobox <- renderInfoBox(NULL)
            }

        } else {
            settings_loadref$loadref_path = ""
            output$fasta_source_infobox <- renderInfoBox(NULL)
            output$gtf_source_infobox <- renderInfoBox(NULL)
            output$mappa_source_infobox <- renderInfoBox(NULL)
            output$NPA_source_infobox <- renderInfoBox(NULL)
            output$BL_source_infobox <- renderInfoBox(NULL)
        }
    }
        
    observeEvent(input$clearLoadRef,{
        settings_loadref$loadref_path = ""
        output$fasta_source_infobox <- renderInfoBox(NULL)
        output$gtf_source_infobox <- renderInfoBox(NULL)
        output$mappa_source_infobox <- renderInfoBox(NULL)
        output$NPA_source_infobox <- renderInfoBox(NULL)
        output$BL_source_infobox <- renderInfoBox(NULL)            
    })

    observeEvent(settings_loadref$loadref_path,{
        req(settings_loadref$loadref_path)
        if(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) {
            load_ref()
        }
        output$txt_reference_path_load <- renderText({
            if(!is_valid(settings_loadref$loadref_path)) {
                ""
            } else {
                settings_loadref$loadref_path
            }
        })
    })
		
# TAB: Design Experiment page
    files_header = c("bam_file", "irf_file", "cov_file")

    # Make sure df.anno exists when df.files exist:
    observeEvent(settings_expr$df.files, {
        req(settings_expr$df.files)
        if(!is_valid(settings_expr$df.anno)) {
            DT = as.data.table(settings_expr$df.files)
            settings_expr$df.anno = DT[, "sample"]
        } else {
            # merge with existing samples
            DT = as.data.table(settings_expr$df.files)
            settings_expr$df.anno = update_data_frame(settings_expr$df.anno, DT[, "sample"])
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
            rhandsontable(settings_expr$df.files, useTypes = TRUE, 
            selectCallback = TRUE, stretchH = "all")
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
		
    observeEvent(settings_loadref$loadref_path, {
        output$ref_expr_infobox <- renderUI({
            if(!is_valid(settings_loadref$loadref_path)) {
                ui_infobox_ref("")            
            } else {
                ui_infobox_ref(file.path(settings_loadref$loadref_path, "settings.Rds"))
            }
        }) 
    })
    
    observeEvent(input$ref_expr_goto_loadref, {
        req(input$ref_expr_goto_loadref)
        updateTabItems(session, "navSelection", "navRef_Load")
    })

    observe({  
        shinyDirChoose(input, "dir_bam_path_load", roots = c(default_volumes, addit_volume), session = session)

        req(input$dir_bam_path_load)
        settings_expr$expr_path = 
            dirname(parseDirPath(c(default_volumes, addit_volume), input$dir_bam_path_load))
        settings_expr$bam_path = 
            parseDirPath(c(default_volumes, addit_volume), input$dir_bam_path_load)
    })
    Expr_Load_BAMs = function() {
    # First assume bams are named by subdirectory names
        if(!is_valid(settings_expr$bam_path)) return(-1)
        temp.DT = FindSamples(settings_expr$bam_path, suffix = ".bam", use_subdir = TRUE)
        if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
            temp.DT = as.data.table(temp.DT)
            if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
                # Assume subdirectory names designate sample names
            } else {
                temp.DT = as.data.table(FindSamples(
                    settings_expr$bam_path, suffix = ".bam", use_subdir = FALSE))
                if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Else assume bam names designate sample names					
                } else {
                    # output$txt_bam_path_expr <- renderText("BAM file names (or its path names) must be unique")
                    sendSweetAlert(
                        session = session,
                        title = "Incompatible BAM file names",
                        text = paste("Could not determine sample names.",
                            "Please ensure either BAMs are uniquely named by sample name,",
                            "or its parent directories are uniquely named."
                        ),
                        type = "error"
                    )
                    settings_expr$bam_path = ""
                    temp.DT = NULL
                }
            }
        } else {
            # output$txt_bam_path_expr <- renderText("No bam files found in given path")
            sendSweetAlert(
                session = session,
                title = "No BAM files found",
                text = "No BAM files found",
                type = "error"
            )            
            settings_expr$bam_path = ""
            temp.DT = NULL
        }
    # compile experiment df with bam paths
        if(!is.null(temp.DT) && nrow(temp.DT) > 0)  {
            colnames(temp.DT)[2] = "bam_file"
            if(is_valid(settings_expr$df.files)) {
        # merge with existing dataframe	
                settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT)
            } else {
        # start anew
                DT = data.table(sample = temp.DT$sample,
                    bam_file = "", irf_file = "", cov_file = "")
                DT[temp.DT, on = "sample", bam_file := i.bam_file] # Update new bam paths
                settings_expr$df.files = as.data.frame(DT)
            }
            return(0)
        } else {
            return(-1)
        }
    }
    observeEvent(settings_expr$bam_path,{
        Expr_Load_BAMs()
        if(is_valid(settings_expr$df.files)) {
            if(is_valid(settings_expr$bam_path) &&
                    "bam_file" %in% colnames(settings_expr$df.files) && 
                    all(file.exists(settings_expr$df.files$bam_file))) {
                    
                output$bam_expr_infobox <- renderUI({
                    ui_infobox_bam(settings_expr$bam_path, settings_expr$df.files$bam_file)
                })
                
            } else if("irf_file" %in% colnames(settings_expr$df.files) && 
                    all(file.exists(settings_expr$df.files$irf_file))) {
                    
                output$bam_expr_infobox <- renderUI({
                    ui_infobox_bam(settings_expr$bam_path, escape = TRUE)
                })
                
            } else if(is_valid(settings_expr$collate_path) && 
                    file.exists(file.path(
                    settings_expr$collate_path, "colData.Rds"
                    ))){
                    
                output$bam_expr_infobox <- renderUI({
                    ui_infobox_bam(settings_expr$bam_path, escape = TRUE)
                })
                
            } else if("bam_file" %in% colnames(settings_expr$df.files)) {
                output$bam_expr_infobox <- renderUI({
                    ui_infobox_bam(settings_expr$bam_path, settings_expr$df.files$bam_file)
                })
            }
        } else {
            output$bam_expr_infobox <- renderUI({
                ui_infobox_bam(settings_expr$bam_path)
            })                 
        } 
    })

    # Run IRFinder
    observeEvent(input$run_irf_expr,{
        req(settings_expr$df.files)
        req("bam_file" %in% colnames(settings_expr$df.files))
        if(!is_valid(input$hot_files_expr_select$select$r)) {
             sendSweetAlert(
                session = session,
                title = "No BAM files selected",
                text = "Please select bam files to run IRFinder",
                type = "error"
            )           
        }
        req(input$hot_files_expr_select$select$r)
        selected_rows = seq(
            input$hot_files_expr_select$select$r,
            input$hot_files_expr_select$select$r2                
        )
        selected_cols = seq(
            input$hot_files_expr_select$select$c,
            input$hot_files_expr_select$select$c2                
        )
        bam_col = which(colnames(settings_expr$df.files) == "bam_file")
        bam_files = settings_expr$df.files$bam_file[selected_rows]
        if(settings_loadref$loadref_path == "") {
            sendSweetAlert(
                session = session,
                title = "Missing Reference",
                text = "Please load Reference before running IRFinder",
                type = "error"
            )
        } else if(!(bam_col %in% selected_cols)) {
            sendSweetAlert(
                session = session,
                title = "No BAM files selected",
                text = "Please select bam files to run IRFinder",
                type = "error"
            )                
        } else if(!all(file.exists(bam_files))) {
            sendSweetAlert(
                session = session,
                title = "Missing BAMs",
                text = "Please check all selected bam files exist",
                type = "error"
            )                
        } else if(!file.exists(file.path(settings_loadref$loadref_path, "IRFinder.ref.gz"))) {
            sendSweetAlert(
                session = session,
                title = "Missing Reference",
                text = "IRFinder.ref.gz is missing",
                type = "error"
            )
        } else if(!is_valid(settings_expr$irf_path) || !dir.exists(settings_expr$irf_path)) {
            sendSweetAlert(
                session = session,
                title = "Missing IRFinder output path",
                text = "Please set IRFinder output path",
                type = "error"
            )
        } else {
            settings_expr$selected_rows = selected_rows
            n_threads = min(get_threads(), length(selected_rows))
            if(n_threads < length(selected_rows)) {
                n_rounds = ceiling(length(selected_rows) / n_threads)
                n_threads = ceiling(length(selected_rows) / n_rounds)
            }
            msg = paste("Run IRFinder on", length(selected_rows), "samples?",
                "Estimated runtime", 10 * 
                    ceiling(length(selected_rows) / n_threads),
                "minutes using", n_threads, 
                "threads (10min per BAM @ 100 million reads per sample)"
            )
            ask_confirmation(
                inputId = "irf_confirm",
                type = "warning",
                title = msg,
                btn_labels = c("Cancel", "Run IRFinder"),
                btn_colors = c("#00BFFF", "#FE2E2E")
            )
        }
    })

    # Function that ACTUALLY runs IRFinder as it requires the sweetalert confirmation to proceed
    observeEvent(input$irf_confirm, {
        if(input$irf_confirm == FALSE) {
            # do nothing
        } else {
            selected_rows = settings_expr$selected_rows

            n_threads = min(get_threads(), length(selected_rows))
            if(n_threads == 1) {
                # run IRFinder using single thread
                withProgress(message = 'Running IRFinder', value = 0, {
                    i_done = 0
                    incProgress(0.001, 
                        message = paste('Running IRFinder',
                            i_done, "of", length(selected_rows), "done")
                    )
                    for(i in selected_rows) {
                        # run_IRFinder(
                            # settings_expr$df.files$bam_file[i], 
                            # file.path(settings_loadref$loadref_path, "IRFinder.ref.gz"), 
                            # file.path(output_path, settings_expr$df.files$sample[i])
                        # )
                        IRFinder(
                            bamfiles = settings_expr$df.files$bam_file[i],
                            sample_names = settings_expr$df.files$sample[i],
                            reference_path = settings_loadref$loadref_path,
                            output_path = settings_expr$irf_path,
                            n_threads = 1,
                            run_featureCounts = FALSE,
                            localHub = FALSE,
                            ah = AnnotationHub(localHub = localHub)                            
                        )
                        i_done = i_done + 1
                        incProgress(1 / length(selected_rows), 
                            message = paste(i_done, "of", length(selected_rows), "done")
                        )
                    }
                })
            } else if(n_threads <= length(selected_rows)) {
                n_rounds = ceiling(length(selected_rows) / n_threads)
                n_threads = ceiling(length(selected_rows) / n_rounds)

                BPPARAM = BiocParallel::bpparam()
                if(n_threads == 1) {
                    BPPARAM_mod = BiocParallel::SerialParam()
                    message(paste("Using SerialParam", BPPARAM_mod$workers, "threads"))
                } else if(Sys.info()["sysname"] == "Windows") {
                    BPPARAM_mod = BiocParallel::SnowParam(n_threads)
                    message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
                } else {
                    BPPARAM_mod = BiocParallel::MulticoreParam(n_threads)
                    message(paste("Using MulticoreParam", BPPARAM_mod$workers, "threads"))
                }

                # extract subset to run in parallel
                row_starts = seq(selected_rows[1], by = n_threads,
                    length.out = n_rounds)
                withProgress(message = 'Running IRFinder - Multi-threaded', value = 0, {
                    i_done = 0
                    incProgress(0.001, 
                        message = paste('Running IRFinder - Multi-threaded,',
                        i_done, "of", length(selected_rows), "done")
                    )
                    for(i in seq_len(n_rounds)) {
                        selected_rows_subset = seq(row_starts[i], 
                            min(length(selected_rows), row_starts[i] + n_threads - 1)
                        )
                        # BiocParallel::bplapply(selected_rows_subset, 
                            # function(i, run_IRF, df, reference_file, output_path) {
                                # run_IRF(df$bam_file[i], reference_file, 
                                    # file.path(output_path, df$sample[i]))
                            # }, 

                            # df = settings_expr$df.files, 
                            # run_IRF = run_IRFinder, 
                            # reference_file = file.path(settings_loadref$loadref_path, "IRFinder.ref.gz"),
                            # output_path = settings_expr$irf_path, BPPARAM = BPPARAM_mod
                        # )
                        IRFinder(
                            bamfiles = settings_expr$df.files$bam_file[selected_rows_subset],
                            sample_names = settings_expr$df.files$sample[selected_rows_subset],
                            reference_path = settings_loadref$loadref_path,
                            output_path = settings_expr$irf_path,
                            n_threads = n_threads,
                            run_featureCounts = FALSE,
                            localHub = FALSE,
                            ah = AnnotationHub(localHub = localHub)                            
                        )                        
                        i_done = i_done + n_threads
                        incProgress(n_threads / length(selected_rows), 
                            message = paste(i_done, "of", length(selected_rows), "done")
                        )
                    }
                })
            }
            sendSweetAlert(
                session = session,
                title = "IRFinder run completed",
                type = "success"
            )                
        }
        settings_expr$selected_rows = c()
        Expr_Load_IRFs()
        if(is_valid(settings_expr$df.files) && "irf_file" %in% colnames(settings_expr$df.files)) {
            irf_files = settings_expr$df.files$irf_file
        } else {
            irf_files = NULL
        }
        
        output$irf_expr_infobox <- renderUI({
            ui_infobox_irf(settings_expr$irf_path, irf_files)
        })      
    })

    shinyDirChoose(input, "dir_irf_path_load", 
        roots = c(default_volumes, addit_volume), session = session)		
    observeEvent(input$dir_irf_path_load, {
        req(input$dir_irf_path_load)
        settings_expr$irf_path = parseDirPath(c(default_volumes, addit_volume), input$dir_irf_path_load)
    })

    Expr_Load_IRFs = function() {
        if(!is_valid(settings_expr$irf_path)) return(-1)
		# merge irfinder paths
        temp.DT = FindSamples(settings_expr$irf_path, suffix = ".txt.gz", use_subdir = FALSE)
        if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
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
                    # settings_expr$irf_path = ""
                    temp.DT = NULL
                }
            }
        } else {
            temp.DT = NULL
        }
        
    # compile experiment df with irfinder paths
        if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
            colnames(temp.DT)[2] = "irf_file"
            if(is_valid(settings_expr$df.files)) {
        # merge with existing dataframe	
                settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT)
            } else {
        # start anew
                DT = data.table(sample = temp.DT$sample,
                    bam_file = "", irf_file = "", cov_file = "")
                DT[temp.DT, on = "sample", irf_file := i.irf_file] # Update new irf paths
                settings_expr$df.files = as.data.frame(DT)      
            }   
        }
        
        # Attempt to find Coverage files
        temp.DT2 = FindSamples(settings_expr$irf_path, suffix = ".cov", use_subdir = FALSE)
        if(!is.null(temp.DT2) && nrow(temp.DT2) > 0) {
            temp.DT2 = as.data.table(temp.DT2)
            if(length(unique(temp.DT2$sample)) == nrow(temp.DT2)) {
                # Assume output names designate sample names
            } else {
                temp.DT2 = as.data.table(FindSamples(
                    settings_expr$irf_path, suffix = ".cov", use_subdir = TRUE))
                if(length(unique(temp.DT2$sample)) == nrow(temp.DT2)) {
            # Else assume subdirectory names designate sample names					
                } else {
                    temp.DT2 = NULL
                }
            }
        } else {
            temp.DT2 = NULL
        }
        
    # compile experiment df with irfinder paths
        if(!is.null(temp.DT2) && nrow(temp.DT2) > 0) {
            colnames(temp.DT2)[2] = "cov_file"
        # merge with existing dataframe	
            settings_expr$df.files = update_data_frame(settings_expr$df.files, temp.DT2)
        }			
        
        if(!is.null(temp.DT) && nrow(temp.DT) > 0)  {
            return(0)
        } else {
            return(-1)
        }  
    }
		
    observeEvent(settings_expr$irf_path,{
        ret = Expr_Load_IRFs()
        if(is_valid(settings_expr$df.files) && "irf_file" %in% colnames(settings_expr$df.files)) {
            irf_files = settings_expr$df.files$irf_file
        } else {
            irf_files = NULL
        }
        
        output$irf_expr_infobox <- renderUI({
            ui_infobox_irf(settings_expr$irf_path, irf_files)
        })
        # ret = is_valid(irf_files) && all(file.exists(irf_files))
        # if(ret == TRUE && !is_valid(settings_expr$bam_path)) {
            # output$bam_expr_infobox <- renderUI({
                # ui_infobox_bam(escape = TRUE)
            # })
        # }
    })

	# Add annotation to data frame
    shinyFileChoose(input, "file_expr_path_load", 
        roots = c(default_volumes, addit_volume), session = session)
    observeEvent(input$file_expr_path_load, {
        req(input$file_expr_path_load)
        file_selected<-parseFilePaths(c(default_volumes, addit_volume), 
            input$file_expr_path_load)
        settings_expr$anno_file = as.character(file_selected$datapath)
    })

    observeEvent(settings_expr$anno_file,{
        req(settings_expr$anno_file)
        temp.DT = tryCatch(as.data.table(fread(settings_expr$anno_file)),
            error = function(e) NULL)
        req(temp.DT)
        req(nrow(temp.DT) > 0)
        output$txt_sample_anno_expr <- renderText({
            validate(need("sample" %in% colnames(temp.DT), 
                "'sample' must be the header of the first column"))
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
            cov_file = "")
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
        req(input$dir_collate_path_load)
        settings_expr$expr_path = 
            dirname(parseDirPath(c(default_volumes, addit_volume), input$dir_collate_path_load))
        settings_expr$collate_path = parseDirPath(c(default_volumes, addit_volume), 
        input$dir_collate_path_load)
    })

    observeEvent(settings_expr$collate_path,{
        # ret = Expr_Load_FSTs()
        if(is_valid(settings_expr$collate_path) &&
            file.exists(file.path(settings_expr$collate_path, "colData.Rds"))) {
            # Load Experiment from "colData.Rds"
            colData.Rds = readRDS(file.path(settings_expr$collate_path, "colData.Rds"))
            req_columns = c("df.anno", "df.files")
            if(all(req_columns %in% names(colData.Rds))) {
                settings_expr$df.files = colData.Rds$df.files
                settings_expr$df.anno = colData.Rds$df.anno
                if("bam_path" %in% names(colData.Rds)) {
                    settings_expr$bam_path = colData.Rds$bam_path
                }
                if("irf_path" %in% names(colData.Rds)) {
                    settings_expr$irf_path = colData.Rds$irf_path
                }
            }
            output$se_expr_infobox <- renderUI({
                ui_infobox_expr(ifelse(is_valid(settings_SE$se),2,1),
                    "Ready to Build Experiment")
            })
        } else if(is_valid(settings_expr$collate_path) &&
            is_valid(settings_expr$df.files) &&
            all(file.exists(settings_expr$df.files$irf_file))) {
            # Ready to run NxtIRF-Collate
            output$se_expr_infobox <- renderUI({
                ui_infobox_expr(1, "Ready to run NxtIRF-Collate")
            })
        } else if(is_valid(settings_expr$collate_path)) {
            # Flag that some IRFinder files need to be built
            output$se_expr_infobox <- renderUI({
                ui_infobox_expr(1, "IRFinder files incomplete")
            })
        } else {
            output$se_expr_infobox <- renderUI({
                ui_infobox_expr(0)
            })        
        }
    })

    observeEvent(input$save_expr,{
        req(input$save_expr)
        if(is_valid(settings_expr$collate_path) &&
            file.exists(file.path(settings_expr$collate_path, "colData.Rds"))) {
            colData.Rds = list(
                df.anno = settings_expr$df.anno,
                df.files = settings_expr$df.files,
                bam_path = settings_expr$bam_path,
                irf_path = settings_expr$irf_path
            )
            saveRDS(colData.Rds, file.path(settings_expr$collate_path, "colData.Rds"))
        } else {
        
        }
    })

    observeEvent(input$load_expr,{
        req(input$load_expr)
        if(is_valid(settings_expr$collate_path) &&
            file.exists(file.path(settings_expr$collate_path, "colData.Rds"))) {
            # Load Experiment from "colData.Rds"
            colData.Rds = readRDS(file.path(settings_expr$collate_path, "colData.Rds"))
            req_columns = c("df.anno", "df.files")
            if(all(req_columns %in% names(colData.Rds))) {
                settings_expr$df.files = colData.Rds$df.files
                settings_expr$df.anno = colData.Rds$df.anno
                if("bam_path" %in% names(colData.Rds)) {
                    settings_expr$bam_path = colData.Rds$bam_path
                }
                if("irf_path" %in% names(colData.Rds)) {
                    settings_expr$irf_path = colData.Rds$irf_path
                }
                output$se_expr_infobox <- renderUI({
                    ui_infobox_expr(ifelse(is_valid(settings_SE$se),2,1),
                        "Ready to Build Experiment")
                })                
            }
        } else {
        
        }
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

        Experiment = na.omit(as.data.table(settings_expr$df.files[, c("sample", "irf_file", "cov_file")]))
        reference_path = settings_loadref$loadref_path
        output_path = settings_expr$collate_path
        # BPPARAM = BPPARAM_mod
        if(!is_valid(settings_loadref$loadref_path)) {
            sendSweetAlert(
                session = session,
                title = "Missing Reference",
                text = "Please load Reference before running NxtIRF::CollateData",
                type = "error"
            )
        } else if(!is_valid(settings_expr$collate_path)) {
            sendSweetAlert(
                session = session,
                title = "Missing NxtIRF Path",
                text = "Please select NxtIRF path before running NxtIRF::CollateData",
                type = "error"
            )
        }
        req(settings_expr$collate_path)

        withProgress(message = 'Collating IRFinder output', value = 0, {
            CollateData(Experiment, reference_path, output_path, n_threads = get_threads())
        })
        # Update colData with df.files and df.anno
        if(file.exists(file.path(
                settings_expr$collate_path, "colData.Rds"))) {
        
            colData.Rds = readRDS(file.path(settings_expr$collate_path, "colData.Rds"))
            if(all(colData.Rds$df.anno$sample %in% settings_expr$df.anno$sample)) {
                colData.Rds$df.anno = settings_expr$df.anno
                colData.Rds$df.files = settings_expr$df.files
                colData.Rds$bam_path = settings_expr$bam_path
                colData.Rds$irf_path = settings_expr$irf_path
                saveRDS(colData.Rds, file.path(settings_expr$collate_path, "colData.Rds"))
                sendSweetAlert(
                    session = session,
                    title = "NxtIRF-Collate run completed",
                    type = "success"
                )
            } else {
                sendSweetAlert(
                    session = session,
                    title = "NxtIRF-Collate did not collate all samples",
                    type = "warning"
                )
            }
        
        } else {

        }
        # Expr_Load_FSTs()
    })
		
    observeEvent(input$clear_expr, {
        settings_expr$expr_path = ""
        settings_expr$bam_path = ""
        settings_expr$irf_path = ""
        settings_expr$anno_file = ""
        settings_expr$collate_path = ""
        settings_expr$df.files = c()
        settings_expr$df.anno = c()
        settings_SE$se = NULL
    })
    observeEvent(input$build_expr, {
        if(is_valid(settings_expr$collate_path) &&
                file.exists(file.path(
                    settings_expr$collate_path, "colData.Rds"))) {
                    
            colData = as.data.table(settings_expr$df.anno)
            settings_SE$se = MakeSE(settings_expr$collate_path, colData)
            
        }
    })
		
# TAB: QC plots
    observeEvent(input$QCmode, {
        req(settings_SE$QC)
        df = as.data.frame(settings_SE$QC)
        rownames(df) = df$sample
        choices = colnames(df)
        choices = choices[!(choices %in% colnames(settings_expr$df.anno))]
        choices = choices[!(choices %in% 
            c("paired", "strand", "path")
        )]
        if(input$QCmode == "PCA") {
            mat = as.matrix(df[, choices])
            rownames(mat) = df$sample
            colVar = matrixStats::colVars(mat)
            mat = mat[,colVar > 0]
            PCA = prcomp(mat, scale. = TRUE)
            output$QC_plot <- renderPlotly({
                print(
                    ggplotly(
                        ggplot(as.data.frame(PCA$x), 
                            aes(x = PC1, y = PC2, text = rownames(PCA$x))
                        ) + geom_point() + geom_text(aes(label = rownames(PCA$x))),
                        tooltip = "text"
                    )
                )
            })
        } else if(input$QCmode == "Graphs") {
            output$QC_plot <- renderPlotly({
                validate(need(is_valid(input$QC_xaxis) | is_valid(input$QC_yaxis),
                    "Specify X or Y axis"))
                if(is_valid(input$QC_xaxis) & is_valid(input$QC_yaxis)) {
                    df.plot = data.frame(sample = df$sample,
                        xaxis = unname(unlist(df[,input$QC_xaxis])),
                        yaxis = unname(unlist(df[,input$QC_yaxis]))
                    )
                    colnames(df.plot)[2:3] = c(
                        input$QC_xaxis, input$QC_yaxis
                    )
                    print(ggplotly(
                        ggplot(df.plot, aes_string(
                            x = input$QC_xaxis, y = input$QC_yaxis,
                            text = "sample")) +
                        geom_point() + geom_text(aes(label = sample)),
                        tooltip = "text"
                    )) 
                } else if(is_valid(input$QC_xaxis)) {
                     df.plot = data.frame(sample = df$sample,
                        xaxis = unname(unlist(df[,input$QC_xaxis]))
                    )
                    colnames(df.plot)[2] = input$QC_xaxis
                    print(ggplotly(
                        ggplot(df.plot, aes_string(
                            x = input$QC_xaxis, y = "sample",
                            text = "sample")) +
                        geom_bar(stat="identity"),
                        tooltip = "text"
                    ))          
                } else if(is_valid(input$QC_yaxis)) {
                     df.plot = data.frame(sample = df$sample,
                        yaxis = unname(unlist(df[,input$QC_yaxis]))
                    )
                    colnames(df.plot)[2] = input$QC_yaxis
                    print(ggplotly(
                        ggplot(df.plot, aes_string(
                            y = input$QC_yaxis, x = "sample",
                            text = "sample")) +
                        geom_bar(stat="identity"),
                        tooltip = "text"
                    ))
                }
            })
        }
    })
        
# TAB: Analyse - Calculate PSIs

    observeEvent(settings_SE$se, {
        req(settings_SE$se)
        req(is(settings_SE$se, "SummarizedExperiment"))
        output$se_expr_infobox <- renderUI({
            ui_infobox_expr(2)
        })

        if(!is_valid(settings_expr$irf_path)) {
            output$irf_expr_infobox <- renderUI({
                ui_infobox_irf(escape = TRUE)
            })
        }            
    })

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
                            settings_SE$se
                        )
                    } else {
                        message(paste("Trigger", i, "is NULL"))
                    }
                }
            }
            settings_SE$filterSummary = filterSummary
            message(sum(filterSummary == TRUE))
        }
    }
        
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
    
    shinyFileSave(input, "saveAnalysis_Filters", 
        roots = c(default_volumes, addit_volume), session = session)
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

    shinyFileChoose(input, "loadAnalysis_Filters", 
        roots = c(default_volumes, addit_volume), session = session,
        filetypes = c("Rds"))
    observeEvent(input$loadAnalysis_Filters, {
        selectedfile <- parseFilePaths(c(default_volumes, addit_volume), input$loadAnalysis_Filters)
        req(selectedfile$datapath)
        settings_SE$filters = readRDS(selectedfile$datapath)
    })

    observeEvent(input$loadDefault_Filters, {
        settings_SE$filters = get_default_filters()
    })

    
# TAB: DE
    
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
                choices = c("(none)", levels(colData[,input$variable_DE])), 
                selected = "(none)")					
            updateSelectInput(session = session, inputId = "denom_DE", 
                choices = c("(none)", levels(colData[,input$variable_DE])), 
                selected = "(none)")
            if(is_valid(settings_DE$nom_DE) && 
                    settings_DE$nom_DE %in% levels(colData[,input$variable_DE])) {
            
                updateSelectInput(session = session, inputId = "nom_DE", selected = settings_DE$nom_DE)
            }
            if(is_valid(settings_DE$nom_DE) && 
                    settings_DE$denom_DE %in% levels(colData[,input$variable_DE])) {
                    
                updateSelectInput(session = session, 
                    inputId = "denom_DE", selected = settings_DE$denom_DE)
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
            settings_DE$batchVar1 = ""
            updateSelectInput(session = session, inputId = "batch1_DE", 
            selected = "(none)")
        }
        if(input$batch2_DE != "(none)" & input$batch2_DE != input$variable_DE & 
            input$batch2_DE != input$batch1_DE) {
            settings_DE$batchVar2 = input$batch2_DE
        } else {
            settings_DE$batchVar2 = ""
            updateSelectInput(session = session, inputId = "batch2_DE", 
            selected = "(none)")
        }
        req(input$method_DE)
        settings_DE$method = input$method_DE

        if(settings_DE$method == "DESeq2") {
        
            res.ASE = DESeq_ASE(se, 
                settings_DE$DE_Var, settings_DE$nom_DE, settings_DE$denom_DE,
                settings_DE$batchVar1, settings_DE$batchVar2,
                n_threads = get_threads()
            )
            if(!input$adjP_DE) {
                setorder(res.ASE, pvalue)
            } else {
                setorder(res.ASE, padj)            
            }
            settings_DE$res = as.data.frame(res.ASE)
            
        } else if(settings_DE$method == "limma") {

            res.ASE = limma_ASE(se, 
                settings_DE$DE_Var, settings_DE$nom_DE, settings_DE$denom_DE,
                settings_DE$batchVar1, settings_DE$batchVar2,
            )
            if(!input$adjP_DE) {
                setorder(res.ASE, P.Value)
            } else {
                setorder(res.ASE, adj.P.Val)
            }
        }
        
        settings_DE$res = as.data.frame(res.ASE)
        output$warning_DE = renderText({"Finished"})

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
		
    shinyFileSave(input, "save_DE", 
        roots = c(default_volumes, addit_volume), session = session,
        filetypes = c("Rds"))
    observeEvent(input$save_DE, {	
        req(settings_DE$res)
        req(length(settings_DE$res_settings) > 0)

        selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$save_DE)
        req(selectedfile$datapath)

        save_DE = list(
            res = settings_DE$res, 
            settings = settings_DE$res_settings, 
            filters = settings_SE$filters)
        saveRDS(save_DE,selectedfile$datapath)
    })

    shinyFileChoose(input, "load_DE", 
        roots = c(default_volumes, addit_volume), 
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
        req(load_DE$settings$method %in% c("DESeq2", "limma"))

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

        if("filters" %in% names(load_DE)) {
            settings_SE$filters = load_DE$filters
        }
    })
		
    observeEvent(input$clear_selected_DE, {
        req(settings_DE$res)
        req(input$DT_DE_rows_selected)
        DT::dataTableProxy("DT_DE") %>% DT::selectRows(NULL)
    })
    
# TAB: Diagonal Plots


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
        df.diag$NMD_direction = settings_DE$res$NMD_direction[
            match(df.diag$EventName, settings_DE$res$EventName)]
        
        settings_Diag$plot_ini = TRUE
        if(input$NMD_diag == TRUE) {
            df.diag = df.diag[df.diag$NMD_direction != 0, ]
            df.diag$nom_NMD = ifelse(df.diag$NMD_direction == 1, df.diag$nom, df.diag$denom)
            df.diag$denom_NMD = ifelse(df.diag$NMD_direction == -1, df.diag$nom, df.diag$denom)
            p = ggplot(df.diag, 
                    aes(x = nom_NMD, y = denom_NMD, key = EventName, 
                        text = EventName, colour = selected)
                ) + geom_point() + scale_color_manual(values = c("black", "red")) +
                labs(
                    x = paste(input$nom_diag, "NMD substrate"),
                    y = paste(input$denom_diag, "NMD substrate")
                )
        } else {
            p = ggplot(df.diag, aes(x = nom, y = denom, key = EventName, 
                    text = EventName, colour = selected)
                ) + geom_point() + scale_color_manual(values = c("black", "red")) +
            labs(
                x = paste(input$nom_diag),
                y = paste(input$denom_diag)
            )         
        }   
        
        settings_Diag$final_plot = ggplotly(p, 
            tooltip = "text",
            source = "plotly_diagonal") %>% 
            layout(
                dragmode = "lasso",
                yaxis = list(scaleanchor="x", scaleratio=1)
            )      
        print(
            settings_Diag$final_plot
        )
    })

    shinyFileSave(input, "saveplot_diag", 
        roots = c(default_volumes, addit_volume), session = session,
        filetypes = c("pdf"))
    observeEvent(input$saveplot_diag, {	
        req(settings_Diag$final_plot)
        selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$saveplot_diag)
        req(selectedfile$datapath)

        obj = isolate(settings_Diag$final_plot)
        plotly::orca(obj, make.path.relative(getwd(), selectedfile$datapath))
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
        # print(click)
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

    settings_Diag$plotly_brush = reactive({
        plot_exist = settings_Diag$plot_ini
        if(plot_exist == TRUE) {
            event_data("plotly_selected", source = "plotly_diagonal")
        }
    })
  
    observeEvent(settings_Diag$plotly_brush(), {
        req(settings_Diag$plotly_brush())
        brush = settings_Diag$plotly_brush()
        # print(brush)
        brush.id = which(settings_DE$res$EventName %in% brush$key)
        req(brush.id)

        selected = input$DT_DE_rows_selected
        selected = unique(c(selected, brush.id))

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

    settings_Volc$plotly_click = reactive({
        plot_exist = settings_Volc$plot_ini
        if(plot_exist == TRUE) {
            event_data("plotly_click", source = "plotly_volcano")
        }
    })
  
    observeEvent(settings_Volc$plotly_click(), {
        req(settings_Volc$plotly_click())
        click = settings_Volc$plotly_click()
        # print(click)
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
        # print(brush)
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
            NMD_direction = NMD_direction,
            log2FoldChange = log2FoldChange, pvalue = pvalue, padj = padj))
        } else {
            df.volc = with(res, data.frame(EventName = EventName, EventType = EventType,
            NMD_direction = NMD_direction,
            log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val))	
        }

        if(is_valid(selected)) {
            df.volc$selected = (df.volc$EventName %in% settings_DE$res$EventName[selected])
        } else {
            df.volc$selected = FALSE
        }
        if(input$NMD_volc == TRUE) {
            df.volc = df.volc[df.volc$NMD_direction != 0, ]
            df.volc$log2FoldChange = df.volc$log2FoldChange * df.volc$NMD_direction
        }

        settings_Volc$plot_ini = TRUE
        if(input$adjP_volc == TRUE) {
            p = ggplot(df.volc, aes(x = log2FoldChange, y = -log10(padj),
                key = EventName, text = EventName, colour = selected))           
        } else {
            p = ggplot(df.volc, aes(x = log2FoldChange, y = -log10(pvalue),
                key = EventName, text = EventName, colour = selected))               
        }

        p = p + geom_point() + scale_color_manual(values = c("black", "red"))

        if(input$facet_volc == TRUE) {
            p = p + facet_wrap(vars(EventType))
        }
        if(input$NMD_volc == TRUE) {
            p = p + labs(x = "Log2 Fold Change NMD substrate")
        } else {
            p = p + labs(x = "Log2 Fold Change")            
        }
        if(input$adjP_volc == TRUE) {
            p = p + labs(y = "Adjusted P Value (-log10)")
        } else {
            p = p + labs(x = "Nominal P Value (-log10)")            
        }
        settings_Volc$final_plot = ggplotly(p, 
            tooltip = "text",
            source = "plotly_volcano") %>% layout(dragmode = "lasso")
        print(
        settings_Volc$final_plot
        )
    })

    shinyFileSave(input, "saveplot_volc", 
        roots = c(default_volumes, addit_volume), session = session,
        filetypes = c("pdf"))
    observeEvent(input$saveplot_volc, {	
        req(settings_Volc$final_plot)
        selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$saveplot_volc)
        req(selectedfile$datapath)

        obj = isolate(settings_Volc$final_plot)
        plotly::orca(obj, make.path.relative(getwd(), selectedfile$datapath))
    })
    
    observeEvent(input$clear_volc, {
        updateSelectInput(session = session, "EventType_volc", selected = NULL)
        shinyWidgets::updateSliderTextInput(session = session, "number_events_volc", selected = 10000)
    })
		
# TAB: Heatmaps
		
    output$plot_heat <- renderPlotly({
        
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
        color = grDevices::colorRampPalette(
            rev(RColorBrewer::brewer.pal(n = colors.df$maxcolors[color.index],
            name = rownames(colors.df)[color.index]))
        )

        # Hopefully the fixed filtering in limma pipeline will also fix the NA issues here:
        na.exclude = (rowSums(!is.na(mat)) == 0)
        if(any(na.exclude == TRUE)) {
        output$warning_heat <- renderText({
            cat("The following events have been excluded due to all NA values:")
            paste(rownames(mat)[which(na.exclude)])
        })
        mat = mat[-which(na.exclude),]
        }

        if(is_valid(input$anno_col_heat) && all(input$anno_col_heat %in% colnames(colData))) {
            settings_Heat$final_plot = heatmaply::heatmaply(
                mat, color = color, 
                col_side_colors = colData[, input$anno_col_heat, drop=FALSE]
            )
        } else {
            settings_Heat$final_plot = heatmaply::heatmaply(mat, color = color)            
        }      
        print(
            settings_Heat$final_plot
        )
    })

    shinyFileSave(input, "saveplot_heat", 
        roots = c(default_volumes, addit_volume), session = session,
        filetypes = c("pdf"))
    observeEvent(input$saveplot_heat, {	
        req(settings_Heat$final_plot)
        selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$saveplot_heat)
        req(selectedfile$datapath)
        
        obj = isolate(settings_Heat$final_plot)
        plotly::orca(obj, make.path.relative(getwd(), selectedfile$datapath))
    })
		
# TAB: RNA-seq Coverage Plots
    
    get_track_selection <- function(i) {
        if(i == 1) return(input$track1_cov)
        if(i == 2) return(input$track2_cov)
        if(i == 3) return(input$track3_cov)
        if(i == 4) return(input$track4_cov)
    }

    # getGeneList <- function() {
        # req(file.exists(file.path(settings_loadref$loadref_path, "settings.Rds"))) 
        # file_path = file.path(settings_loadref$loadref_path, "fst", "Genes.fst")
        # if(!file.exists(file_path)) return(NULL)

        # df = as.data.table(read.fst(file_path))
        # return(df)
    # }
	
    get_inrange_events <- function(view_chr, view_start, view_end) {
        req(settings_Cov$event.ranges)
        event.ranges.legit = settings_Cov$event.ranges[seqnames == view_chr]
        event.ranges.legit = event.ranges.legit[end > view_start & start < view_end]
        return(event.ranges.legit$EventName)
    }
    
    observe({
        view_chr = input$chr_cov
        view_start = suppressWarnings(as.numeric(input$start_cov))
        view_end = suppressWarnings(as.numeric(input$end_cov))

        req(view_chr)
        req(view_start)
        req(view_end)
        
        selected_event = isolate(input$events_cov)
        cur_event = isolate(input$event_norm_cov)

        event_choices = c("(none)")
        if(is_valid(selected_event)) {
            event_choices = c(event_choices, selected_event)
        } else if(is_valid(cur_event)) {
            event_choices = c(event_choices, cur_event)        
        }
        event_choices = unique(c(event_choices, 
            get_inrange_events(view_chr, view_start, view_end)))

        if(is_valid(selected_event)) {
            updateSelectInput(session = session, inputId = "event_norm_cov", 
                choices = event_choices, selected = selected_event)
        } else if(is_valid(cur_event)) {
            updateSelectInput(session = session, inputId = "event_norm_cov", 
                choices = event_choices, selected = cur_event)
        } else {
            updateSelectInput(session = session, inputId = "event_norm_cov", 
                choices = event_choices, selected = "(none)")        
        }
    })
    
    observeEvent(input$refresh_coverage, {
        view_chr = input$chr_cov
        view_start = suppressWarnings(as.numeric(input$start_cov))
        view_end = suppressWarnings(as.numeric(input$end_cov))
        graph_mode = input$graph_mode_cov
			
        req(view_chr)
        req(view_start)
        req(view_end)
        req(settings_SE$se)
        
        req(view_end - view_start > 0)
        
      # refresh in-range events here
        if(is_valid(input$event_norm_cov)) {
            norm_event = input$event_norm_cov
        } else {
            norm_event = NULL
        }

        rowData = SummarizedExperiment::rowData(settings_SE$se)
        events_to_highlight = list()
        if(!is.null(norm_event) && norm_event %in% rowData$EventName) {
            if(rowData$EventType[match(norm_event, rowData$EventName)] 
                %in% c("MXE", "SE")) {
                events_to_highlight[[1]] = c(rowData$Event1a[match(norm_event, rowData$EventName)],
                    rowData$Event2a[match(norm_event, rowData$EventName)])
            } else {
                events_to_highlight[[1]] = rowData$Event1a[match(norm_event, rowData$EventName)]
            }
            if(rowData$EventType[match(norm_event, rowData$EventName)] 
                %in% c("MXE")) {
                events_to_highlight[[2]] = c(rowData$Event1b[match(norm_event, rowData$EventName)],
                    rowData$Event2b[match(norm_event, rowData$EventName)])
            } else if(rowData$EventType[match(norm_event, rowData$EventName)] 
                %in% c("SE", "A3SS", "A5SS", "ALE", "AFE")){
                events_to_highlight[[2]] = rowData$Event1b[match(norm_event, rowData$EventName)]
            }           
        }
        conf.int = 0.95

        if(is.null(settings_Cov$elem.DT)) settings_Cov$elem.DT <- 
            loadViewRef(settings_loadref$loadref_path)
        if(is.null(settings_Cov$transcripts.DT)) settings_Cov$transcripts.DT <- 
            loadTranscripts(settings_loadref$loadref_path)

        req(settings_Cov$elem.DT)
        req(settings_Cov$transcripts.DT)
     
        tracks = list()
        for(i in seq_len(4)) {
            tracks[[i]] = get_track_selection(i)       
        }
     
        obj = plot_cov_fn(
            view_chr, view_start, view_end, input$strand_cov, 
            norm_event, input$condition_cov, tracks = tracks, 
            track_names = "",
            se = settings_SE$se, settings_Cov$avail_cov,
            settings_Cov$transcripts.DT, settings_Cov$elem.DT,
            events_to_highlight,
            graph_mode = graph_mode,
            stack_tracks = input$stack_tracks_cov,
            t_test = input$pairwise_t_cov,
            condensed = input$condense_cov
        )
        
        settings_Cov$final_plot = obj$final_plot
        # settings_SaveObj$obj = obj$ggplot
        
        output$plot_cov <- renderPlotly({
            settings_Cov$plot_ini = TRUE      
            print(
                settings_Cov$final_plot
            )
        })
    })

    observeEvent(input$graph_mode_cov, {
        req(settings_Cov$plot_ini == TRUE)
        if(input$graph_mode_cov == "Pan") {
            plotlyProxy("plot_cov", session) %>% 
                plotlyProxyInvoke("relayout", list(dragmode = "pan")) %>%
                plotlyProxyInvoke("reconfig", editable = FALSE)
        } else if(input$graph_mode_cov == "Zoom") {
            plotlyProxy("plot_cov", session) %>% 
                plotlyProxyInvoke("relayout", list(dragmode = "zoom")) %>%
                plotlyProxyInvoke("reconfig", editable = FALSE)
        } else if(input$graph_mode_cov == "Movable Labels") {
            plotlyProxy("plot_cov", session) %>% 
                plotlyProxyInvoke("relayout", list(dragmode = FALSE)) %>%
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

            updateSelectInput(session = session, inputId = "track1_cov", 
                choices = c("(none)"), selected = "(none)")     
            updateSelectInput(session = session, inputId = "track2_cov", 
                choices = c("(none)"), selected = "(none)")  
            updateSelectInput(session = session, inputId = "track3_cov", 
                choices = c("(none)"), selected = "(none)")    
            updateSelectInput(session = session, inputId = "track4_cov", 
                choices = c("(none)"), selected = "(none)")             
        } else {
            updateSelectInput(session = session, inputId = "condition_cov", 
                choices = c("(none)"))
        }        
        refresh_tracks_cov()
    })

    observeEvent(input$condition_cov, {
        req(input$condition_cov)
        req(settings_SE$se)
        refresh_tracks_cov()
    })

    refresh_tracks_cov <- function() {
        if(input$mode_cov == "By Condition") {
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
            }  else {
                updateSelectInput(session = session, inputId = "track1_cov", 
                    choices = "(none)")    
                updateSelectInput(session = session, inputId = "track2_cov", 
                    choices = "(none)")    
                updateSelectInput(session = session, inputId = "track3_cov", 
                    choices = "(none)")    
                updateSelectInput(session = session, inputId = "track4_cov", 
                    choices = "(none)")    
            }
        } else {
            avail_samples = names(settings_Cov$avail_cov)
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
  
    shinyFileSave(input, "saveplot_cov", roots = c(default_volumes, addit_volume), session = session,
        filetypes = c("pdf"))
    observeEvent(input$saveplot_cov, {	
        req(settings_Cov$final_plot)
        selectedfile <- parseSavePath(c(default_volumes, addit_volume), input$saveplot_cov)
        req(selectedfile$datapath)
        
        obj = isolate(settings_Cov$final_plot)
        plotly::orca(settings_Cov$final_plot, make.path.relative(getwd(), selectedfile$datapath),
            width = 1920, height = 1080)
    })
 
# End of server function		
 }