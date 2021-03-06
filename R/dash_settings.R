setreactive_system <- function() {
    reactiveValues(
        threads_initialized = FALSE,
        n_threads = 1
    )
}

setreactive_newref <- function() {
    reactiveValues(
        newref_path = "",
        newref_fasta = "",
        newref_gtf = "",
        newref_AH_fasta = "",
        newref_AH_gtf = "",
        newref_mappa = "",
        newref_NPA = "",
        newref_bl = "",
        
        ui_newrefAH_Species = "",
        ui_newrefAH_Version_Trans = "",
        ui_newrefAH_Trans = "",
        ui_newrefAH_Assembly = "",
        ui_newrefAH_Version_Genome = "",
        ui_newrefAH_Genome = "",
        
        ui_newref_genome_type = ""
    )
}

setreactive_loadref <- function() {
    reactiveValues(
        loadref_path = "",
        settings = c()
    )
}

setreactive_expr <- function() {
    reactiveValues(
        expr_path = "",
        bam_path = "",
        irf_path = "",
        anno_file = "",
        collate_path = "",
        selected_rows = c(),
        df = c(),
        df.files = c(),
        df.anno = c()
    )
}

setreactive_SE <- function() {
    reactiveValues(
        se = NULL,
        filterSummary = NULL,
        filters = list(),
        QC = c()
    )
}

setreactive_DE <- function() {
    reactiveValues(
        res = NULL,
        res_settings = list(),
        method = NULL,
        batchVar1 = NULL,
        batchVar2 = NULL,
        DE_Var = NULL,
        nom_DE = NULL,
        denom_DE = NULL
    )
}

setreactive_Diag <- function() {
    reactiveValues(
        plot_ini = FALSE,
		plotly_click = NULL,
        final_plot = NULL
    )
}

setreactive_Heat <- function() {
    reactiveValues(
        final_plot = NULL
    )
}

setreactive_Cov <- function() {
    reactiveValues(
    # Memory
        loaded_reference = "",
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
        plot_ini = FALSE,
        
        final_plot = NULL
    )
}