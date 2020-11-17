
ui_sidebar <- function() {
    dashboardSidebar(
        sidebarMenu(id = "navSelection",
            menuItem("About", tabName = "navTitle"),
            menuItem("System", tabName = "navSystem"),
            menuItem("Reference", tabName = "navRef_New"),
                 # menuSubItem("New Reference", tabName = "navRef_New"),
                 # menuSubItem("Load Reference", tabName = "navRef_Load"),
                 # menuSubItem("View Reference", tabName = "navRef_View")
            menuItem("Experiment", tabName = "navExpr"),
            menuItem("Analysis",
                 menuSubItem("Experiment QC", tabName = "navQC"),
                 menuSubItem("Filters", tabName = "navFilter"),
                 menuSubItem("Differential Expression Analysis", tabName = "navAnalyse")
            ),
            menuItem("Display",
                 menuSubItem("Diagonal", tabName = "navDiag"),
                 menuSubItem("Volcano", tabName = "navVolcano"),
                 menuSubItem("Heatmap", tabName = "navHeatmap"),
                 menuSubItem("Coverage", tabName = "navCoverage")
            )
        )
    )
}

ui_tab_title <- function() {
    tabItem(tabName = "navTitle",
        img(src="https://pbs.twimg.com/profile_images/1310789966293655553/7HawCItY_400x400.jpg")		
    )
}

ui_tab_system <- function() {
    tabItem(tabName = "navSystem",
        # box(
            # tags$div(title = paste("Load NxtIRF Session from File"),
                # shinyFilesButton("file_loadstate", label = "Choose Save State to Load", 
                    # title = "Choose Save State to Load", multiple = FALSE,
                    # filetype = list(RDS = "Rds"))
            # ),
            # tags$div(title = paste("Saves NxtIRF Session to File"),
                # shinySaveButton("file_savestate", "Choose file to Save Session", "Choose file to Save Session", 
                    # filetype = list(RDS = "Rds")),
            # )
        # ),
        box(
            tags$div(title = paste("Number of threads to run computationally-intensive operations",
                "such as IRFinder, NxtIRF-collate, and DESeq2"),
                numericInput("cores_numeric", "# Threads", min = 1, max = 1, value = 1)
            ),
            tags$div(title = paste("Number of threads to run computationally-intensive operations",
                "such as IRFinder, NxtIRF-collate, and DESeq2"),
                sliderInput('cores_slider', "# Threads", min = 1, max = 1, value = 1)			
            )
        )
    )
}

ui_tab_ref_new <- function() {
    tabItem(tabName = "navRef_New",
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
    )
}

ui_tab_ref_load <- function() {
    tabItem(tabName = "navRef_Load",
        fluidRow(
            box(
                h4("Select Reference Directory"),
                shinyDirButton("dir_reference_path_load", label = "Choose reference path", title = "Choose reference path"),
                textOutput("txt_reference_path_load"),                        
            ),
            box(
                actionButton("clearLoadRef", "Clear settings"), # TODO                        
            )

        ),
        conditionalPanel(
            condition = "output.txt_reference_path_load != ''",
            fluidRow(
                infoBoxOutput("fasta_source_infobox"),
                infoBoxOutput("gtf_source_infobox")
            ),
            fluidRow(
                infoBoxOutput("mappa_source_infobox"),
                infoBoxOutput("NPA_source_infobox"),
                infoBoxOutput("BL_source_infobox")
            )
        )
    )
}

ui_ddb_ref_load <- function() {
    shinyWidgets::dropdownButton(
        # box(
            h4("Select Reference Directory"),
            shinyDirButton("dir_reference_path_load", label = "Choose reference path", title = "Choose reference path"),
            textOutput("txt_reference_path_load"),                        
        # ),
        # box(
            actionButton("clearLoadRef", "Clear settings"), # TODO                        
        # ),
        
        circle = TRUE, status = "danger",
        icon = icon("folder-open", lib = "font-awesome"), 
        width = "300px",
        inputId = "ref_ddb"
    )
}

ui_infobox_ref <- function(settings_file) {
    box1 = infoBox(
        title = "Reference", 
        value = ifelse(file.exists(settings_file),
            "LOADED", "MISSING"),
        subtitle = ifelse(file.exists(settings_file),
            dirname(settings_file), "click here to set reference"),
        icon = icon("dna", lib = "font-awesome"),
        color = ifelse(file.exists(settings_file),
            "green", "red")
    )
    ddb = ui_ddb_ref_load()
    
    box1$children[[1]]$attribs$id = ddb$children[[1]]$attribs$id
    box1$children[[1]]$attribs$`data-toggle` = "dropdown"
    box1$children[[2]] = ddb$children[[2]]
    return(box1)
}

ui_ddb_bam_path <- function() {
    shinyWidgets::dropdownButton(
        tags$h4("Set BAM Path"),
        shinyDirButton("dir_bam_path_load", 
            label = "Choose BAM path", 
            title = "Choose BAM path"), # done
        textOutput("txt_bam_path_expr"),
        
        circle = TRUE, status = "danger",
        icon = icon("folder-open", lib = "font-awesome"), 
        width = "300px",
        inputId = "bam_ddb"
    )
}

ui_infobox_bam <- function(bam_path, bam_files, escape = FALSE) {
    if(escape == TRUE) {
        box1 = infoBox(
            title = "bam path", 
            value = "NOT REQUIRED",
            icon = icon("folder-open", lib = "font-awesome"),
            color = "green"
        )
    } else {
        ret = !missing(bam_files) &&is_valid(bam_files) && all(file.exists(bam_files))
        box1 = infoBox(
            title = "bam path", 
            value = ifelse(!is_valid(bam_path),
                "MISSING", ifelse(ret == TRUE, "LOADED", "No BAMs found")),
            subtitle = ifelse(is_valid(bam_path),
                bam_path, ""),
            icon = icon("folder-open", lib = "font-awesome"),
            color = ifelse(!is_valid(bam_path),
                "red", ifelse(ret == TRUE, "green", "yellow"))
        )
    }
    ddb = ui_ddb_bam_path()
    
    box1$children[[1]]$attribs$id = ddb$children[[1]]$attribs$id
    box1$children[[1]]$attribs$`data-toggle` = "dropdown"
    box1$children[[2]] = ddb$children[[2]]
    return(box1)
}

ui_ddb_irf_path <- function() {
    shinyWidgets::dropdownButton(
        tags$h4("Set IRFinder Path"),
        shinyDirButton("dir_irf_path_load", 
        
            label = "Choose IRFinder output path", 
            title = "Choose IRFinder output path"), # done					
        textOutput("txt_irf_path_expr"),
        
        circle = TRUE, status = "danger",
        icon = icon("folder-open", lib = "font-awesome"), 
        width = "300px",
        inputId = "irf_ddb"
    )
}

ui_infobox_irf <- function(irf_path, irf_files, escape = FALSE) {
    if(escape == TRUE) {
        box1 = infoBox(
            title = "irfinder output", 
            value = "NOT REQUIRED",
            icon = icon("align-center", lib = "font-awesome"),
            color = "green"
        )
    } else {
        ret = is_valid(irf_files) && all(file.exists(irf_files))
        box1 =  infoBox(
            title = "irfinder output", 
            value = ifelse(!is_valid(irf_path),
                "MISSING", ifelse(ret == TRUE, "LOADED", "Some IRF files missing")),
            subtitle = ifelse(is_valid(irf_path),
                irf_path, ""),
            icon = icon("align-center", lib = "font-awesome"),
            color = ifelse(!is_valid(irf_path),
            "red", ifelse(ret == TRUE, "green", "yellow"))
        )
    }
    ddb = ui_ddb_irf_path()
    
    box1$children[[1]]$attribs$id = ddb$children[[1]]$attribs$id
    box1$children[[1]]$attribs$`data-toggle` = "dropdown"
    box1$children[[2]] = ddb$children[[2]]
    return(box1)
}

ui_ddb_nxt_path <- function() {
    shinyWidgets::dropdownButton(
        tags$h4("Set NxtIRF Path"),
        shinyDirButton("dir_collate_path_load", 
            label = "Choose NxtIRF FST output path", 
            title = "Choose NxtIRF FST output path"), # done
        textOutput("txt_collate_path_expr"), # done
        
        circle = TRUE, status = "danger",
        icon = icon("folder-open", lib = "font-awesome"), 
        width = "300px",
        inputId = "nxt_ddb"
    )
}

ui_infobox_nxt <- function(nxt_path, nxt_files, escape = FALSE) {
    if(escape == TRUE) {
        box1 = infoBox(
            title = "NxtIRF output", 
            value = "NOT REQUIRED",
            icon = icon("align-center", lib = "font-awesome"),
            color = "green"
        )
    } else {
        ret = is_valid(nxt_files) && all(file.exists(nxt_files))
        box1 =  infoBox(
            title = "NxtIRF output", 
            value = ifelse(!is_valid(nxt_path),
                "MISSING", ifelse(ret == TRUE, "LOADED", "Some IRF files missing")),
            subtitle = ifelse(is_valid(nxt_path),
                nxt_path, ""),
            icon = icon("layer-group", lib = "font-awesome"),
            color = ifelse(!is_valid(nxt_path),
            "red", ifelse(ret == TRUE, "green", "yellow"))
        )
    }
    ddb = ui_ddb_nxt_path()
    
    box1$children[[1]]$attribs$id = ddb$children[[1]]$attribs$id
    box1$children[[1]]$attribs$`data-toggle` = "dropdown"
    box1$children[[2]] = ddb$children[[2]]
    return(box1)
}

ui_ddb_build_expr <- function() {
    shinyWidgets::dropdownButton(
        tags$h4("Construct Experiment"),
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
        
        circle = TRUE, status = "danger",
        icon = icon("folder-open", lib = "font-awesome"), 
        width = "300px",
        inputId = "nxt_ddb"
    )
}

ui_infobox_expr <- function(status = 0) {

    box1 =  infoBox(
        title = "SummarizedExperiment Object", 
        value = ifelse(status == 0,
            "MISSING", ifelse(status == 2, "BUILT / LOADED", "Please BUILD Experiment")),
        icon = icon("flask", lib = "font-awesome"),
        color = ifelse(status == 0,
        "red", ifelse(status == 2, "green", "yellow"))
    )

    ddb = ui_ddb_build_expr()
    
    box1$children[[1]]$attribs$id = ddb$children[[1]]$attribs$id
    box1$children[[1]]$attribs$`data-toggle` = "dropdown"
    box1$children[[2]] = ddb$children[[2]]
    return(box1)
}


ui_tab_expr <- function() {
    tabItem(tabName = "navExpr",
        conditionalPanel(
            condition = "output.txt_reference_path_load != '' && input.ref_ddb_state == true",
            fluidRow(
                infoBoxOutput("fasta_source_infobox"),
                infoBoxOutput("gtf_source_infobox")
            ),
            fluidRow(
                infoBoxOutput("mappa_source_infobox"),
                infoBoxOutput("NPA_source_infobox"),
                infoBoxOutput("BL_source_infobox")
            )
        ),    
        fluidRow(
            uiOutput("ref_expr_infobox"),
            uiOutput("bam_expr_infobox"),
            uiOutput("irf_expr_infobox"),
            uiOutput("nxt_expr_infobox"),
            uiOutput("se_expr_infobox"),
        ),
        fluidRow(
            column(4,

                shinyWidgets::dropdownButton(
                    tags$h4("Run IRFinder on Selected BAMs"),
                    actionButton("run_irf_expr", 
                        "Run IRFinder"), # TODO
                    textOutput("txt_run_irf_expr"),
                
                    circle = TRUE, status = "danger",
                    icon = icon("align-center", lib = "font-awesome"), 
                    width = "300px"                               
                ), br(),


                shinyWidgets::dropdownButton(
                    tags$h4("Run NxtIRF CollateData on IRFinder output"),
                    actionButton("run_collate_expr", 
                        "Compile NxtIRF FST files"),
                    textOutput("txt_run_col_expr"),
                
                    circle = TRUE, status = "danger",
                    icon = icon("layer-group", lib = "font-awesome"), 
                    width = "300px"                           
                ), br(),
                
                conditionalPanel(
                    condition = "['Annotations'].indexOf(input.hot_switch_expr) >= 0",
                    shinyWidgets::dropdownButton(
                        tags$h4("Annotation Columns"),
                        uiOutput("newcol_expr"), # done
                        div(class='row',
                            div(class= "col-sm-6",
                                radioButtons("type_newcol_expr", "Type", c("character", "integer", "double"))
                            ),
                            div(class = "col-sm-6", 
                                actionButton("addcolumn_expr", "Add"), br(),  # done
                                actionButton("removecolumn_expr", "Remove") # done
                            )
                        ),
                        circle = TRUE, status = "danger",
                        icon = icon("columns", lib = "font-awesome"), 
                        width = "300px"                                      
                    )
                ),
            ),
            column(8,
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
        )
    )
}

ui_tab_filter <- function() {
    tabItem(tabName = "navFilter",
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
                    title = "Load Filters from Rds", multiple = FALSE),
                actionButton("loadDefault_Filters", "Load Default Filters"),
                    
            )
        )
    )
}

ui_tab_analyse <- function() {
    tabItem(tabName = "navAnalyse",
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
                div(style = 'overflow-x: scroll',  DT::dataTableOutput('DT_DE'))
            )
        )
    )
}

ui_tab_diag <- function() {
    tabItem(tabName = "navDiag",
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
    )
}

ui_tab_volcano <- function() {
    tabItem(tabName = "navVolcano",
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
    )
}

ui_tab_heatmap <- function() {
    tabItem(tabName = "navHeatmap",
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
                    c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "
                        RdGy", "RdYlBu", "RdYlGn", "Spectral")
                )            
            ),
            column(9, 
                textOutput("warning_heat"),
                plotlyOutput("plot_heat", height = "800px"),
            )
        )
    )
}

ui_tab_coverage <- function() {
    tabItem(tabName = "navCoverage",
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
                    choices = c("Top N All Results", "Top N Filtered Results", "Highlighted"), 
                    checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                )
            ),
            column(4,
                div(style="display: inline-block;vertical-align:top; width: 80px;",
                    selectInput("chr_cov", label = "Chr", c("(none)"), selected = "(none)")),
                div(style="display: inline-block;vertical-align:top; width: 120px;",
                    textInput("start_cov", label = "Left", c(""))),
                div(style="display: inline-block;vertical-align:top; width: 120px;",
                    textInput("end_cov", label = "Right", c(""))),            
                br(),
                shinyWidgets::actionBttn("zoom_out_cov", style = "material-circle", color = "danger",
                    icon = icon("minus")),
                div(style="display: inline-block;vertical-align:center;width: 50px;padding:25px",
                    textOutput("label_zoom_cov")),
                shinyWidgets::actionBttn("zoom_in_cov", style = "material-circle", color = "danger",
                    icon = icon("plus")),
                div(style="display: inline-block;vertical-align:center;padding:15px",
                    shinyWidgets::radioGroupButtons("strand_cov", 
                        label = "Strand", justified = FALSE,
                        choices = c("*", "+", "-"), 
                        checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                    )
                )
            ),
            column(2,
                div(style="display: inline-block;vertical-align:top;padding:10px;",
                    shinyWidgets::radioGroupButtons("graph_mode_cov",label = "Graph Mode", justified = FALSE,
                        choices = c("Pan", "Zoom", "Movable Labels"), 
                        checkIcon = list(yes = icon("ok", lib = "glyphicon")))
                ),
                shinyWidgets::sliderTextInput("slider_num_events_cov", 
                    "Num Events", choices = c(5, 10,25,50,100,200,500), selected = 25)       
            )
        ),
        fluidRow(
            column(2, 
                selectInput('mode_cov', 'View', width = '100%',
                    choices = c("Individual", "By Condition")),					
                conditionalPanel(
                    condition = "['By Condition'].indexOf(input.mode_cov) >= 0",
                    selectInput('event_norm_cov', 'Normalize Event', width = '100%',
                        choices = c("(none)")),
                    selectInput('condition_cov', 'Condition', width = '100%',
                        choices = c("(none)"))					
                ),
                selectInput('track1_cov', 'Track 1', width = '100%',
                    choices = c("(none)")),
                selectInput('track2_cov', 'Track 2', width = '100%',
                    choices = c("(none)")),
                selectInput('track3_cov', 'Track 3', width = '100%',
                    choices = c("(none)")),
                selectInput('track4_cov', 'Track 4', width = '100%',
                    choices = c("(none)")),
                shinyWidgets::switchInput("stack_tracks_cov", label = "Stack Traces", labelWidth = "150px"),
                shinyWidgets::switchInput("pairwise_t_cov", label = "Pairwise t-test", labelWidth = "150px"),
                shinySaveButton("saveplot_cov", "Save Plot as PDF", "Save Plot as PDF...", 
                    filetype = list(PDF = "pdf")),
            ),
            column(10, 
                plotlyOutput("plot_cov", height = "800px"),
            )
        )
    )
}