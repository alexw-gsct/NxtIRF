GetCoverage_DF <- function(samples, files, seqname, start, end, strand) {
  covData = list()
  for(i in seq_len(length(files))) {
    cov = GetCoverage(files[i], seqname, start - 1, end, ifelse(strand == "+", 0, ifelse(strand == "-", 1, 2)))
    view = IRanges::Views(cov, start, end)
    view.df = as.data.frame(view[[1]])
    covData[[i]] = view.df
  }
  df = do.call(cbind, covData)
  colnames(df) = samples
  x = seq(start,end)
  df = cbind(x, df)
  return(df)
}

bin_df <- function(df, binwidth = 3) {
    DT = as.data.table(df)
    brks = seq(1, nrow(DT), length = round(nrow(DT) / binwidth))
    bin <- NULL
    DT[, bin := findInterval(seq_len(nrow(DT)), brks)]
    DT2 <- DT[, lapply(.SD, mean, na.rm = TRUE), by = bin]
    DT2[, bin := NULL]
    return(as.data.frame(DT2))
}

update_select_without_clearing <- function(session, inputId, choices, input) {
    req(inputId %in% names(input))
    selected = input[[inputId]]
    if(selected %in% choices) {
        updateSelectInput(session = session, inputId = inputId,
            choices = choices, selected = selected)
    } else {
        updateSelectInput(session = session, inputId = inputId,
            choices = choices, selected = selected)    
    }
}

get_psi <- function(se_path, 
        view_chr, view_start, view_end, view_strand = "*"
){
    junc_fst = file.path(se_path, "se", "junc_PSI.fst")
    junc_fst_index = file.path(se_path, "se", "junc_PSI_index.fst")
    assert_that(all(file.exists(c(junc_fst, junc_fst_index))),
        msg = "Some junction fst files do not exist")
        
    junc_index = as.data.table(fst::read.fst(junc_fst_index)) 
    junc_index$index = seq_len(nrow(junc_index))
    junc_index = junc_index[
        seqnames == view_chr &
        start < view_end & end > view_start
    ]
    index_start = min(junc_index$index)
    index_end = max(junc_index$index)
    
    junc_index = as.data.table(fst::read.fst(junc_fst_index, 
        from = index_start, to = index_end))

    junc_data = as.matrix(fst::read.fst(junc_fst, 
        from = index_start, to = index_end))
    junc_data[is.na(junc_data)] = 0
    junc_data = as.data.table(junc_data)
    
    junc = cbind(junc_index, junc_data)
    junc = junc[
        seqnames == view_chr &
        start < view_end & end > view_start
    ]
    junc
}

determine_compatible_events <- function(reduced.DT, highlight_events) {

    introns = reduced.DT[get("type") == "intron"]
    introns[, c("highlight") := "0"]
    exons = reduced.DT[get("type") == "exon"]
    exons[, c("highlight") := "0"]
    misc = reduced.DT[get("type") == "CDS"]
    misc[, c("highlight") := "0"]

    tr_filter = c()
    if(length(highlight_events) == 1) {
        # This is IR
        gr = NxtIRF.CoordToGR(highlight_events[[1]])
        introns.gr = makeGRangesFromDataFrame(as.data.frame(introns))
        OL = findOverlaps(gr, introns.gr)
        introns[OL@to, c("highlight") := 1]
    } else if(length(highlight_events) == 2) {
        # This is AS
        AS_count = 1;
        for(event in highlight_events) {
            gr = NxtIRF.CoordToGR(event)
            introns.gr = makeGRangesFromDataFrame(as.data.frame(introns))
            OL = findOverlaps(gr, introns.gr, type = "equal")
            introns[OL@to, c("highlight") := as.character(AS_count)]

            OL_s1 = findOverlaps(gr[1], introns.gr, type = "equal")
            tr1 = unique(introns$transcript_id[OL_s1@to])
            if(length(gr) == 2) {
                OL_s2 = findOverlaps(gr[2], introns.gr, type = "equal")
                tr1 = unique(intersect(tr1, introns$transcript_id[OL_s2@to]))
            }
            tr_filter = c(tr_filter, tr1)
            coord_keys = c(BiocGenerics::start(gr[1]) - 1, BiocGenerics::end(gr[1]) + 1)
            if(length(gr) == 2) {
                coord_keys = c(coord_keys,
                    BiocGenerics::start(gr[2]) - 1, BiocGenerics::end(gr[2]) + 1)
            }
            exons[get("transcript_id") %in% tr1 & 
                (get("start") %in% coord_keys | get("end")%in% coord_keys),
                c("highlight") := as.character(AS_count)]
            AS_count = AS_count + 1
        }  
    }
    # if(any(exons$highlight != "0")) {
        # OL_cds = findOverlaps(
            # makeGRangesFromDataFrame(as.data.frame(exons[get("highlight") != "0"])),
            # makeGRangesFromDataFrame(as.data.frame(misc))
        # )
        # misc[OL_cds@to, c("highlight") := TRUE]
        # misc[!(get("transcript_id") %in% tr_filter), c("highlight") := FALSE]           
    # }

    return(rbind(introns, exons, misc))
}

plot_view_ref_fn <- function(view_chr, view_start, view_end, 
    transcripts, elems, highlight_events, condensed = FALSE) {
			
    data_start = view_start - (view_end - view_start)
    data_end = view_end + (view_end - view_start)

    transcripts.DT = transcripts[get("seqnames") == view_chr]
    transcripts.DT = transcripts.DT[
        get("start") <= data_end & 
        get("end") >= data_start]
    setorderv(transcripts.DT, c("transcript_support_level", "width"))
			# filter transcripts by criteria if applicable
    message(paste(nrow(transcripts.DT), " transcripts"))

    screen.DT = elems[
        get("transcript_id") %in% transcripts.DT$transcript_id &
        get("type") %in% c("CDS", "start_codon", "stop_codon", "exon")
    ]
    if(condensed != TRUE & nrow(transcripts.DT) <= 100) {
        condense_this = FALSE
        transcripts.DT[, c("group_id") := get("transcript_id")]
        screen.DT[, c("group_id") := get("transcript_id")]        
    } else {
        condense_this = TRUE
        transcripts.DT[, c("group_id") := get("gene_id")]     
        screen.DT[transcripts.DT, on = "transcript_id", 
            c("group_id") := get("gene_id")]
    }

    reduced.DT = copy(screen.DT)
    reduced.DT[get("type") %in% c("CDS", "start_codon", "stop_codon"), c("type") := "CDS"]
    reduced.DT[get("type") != "CDS", c("type") := "exon"]
    
    # add introns to reduced.DT
    introns.DT = as.data.table(grlGaps(
        split(makeGRangesFromDataFrame(as.data.frame(reduced.DT)),
            reduced.DT$transcript_id)
    ))
    introns.DT[, c("type") := "intron"]
    setnames(introns.DT, "group_name", "transcript_id")
    introns.DT[reduced.DT, on = "transcript_id",
        "group_id" := get("i.group_id")]

    reduced.DT = rbind(reduced.DT[, 
        c("seqnames", "start", "end", "strand", "type", "group_id", "transcript_id")], 
    introns.DT[, 
        c("seqnames", "start", "end", "strand", "type", "group_id", "transcript_id")])

    # Highlight events here
    # highlight_events is of syntax chrX:10000-11000/-
    # if(!missing(highlight_events) & condense_this == FALSE) {
    if(!missing(highlight_events)) {    
        reduced.DT = determine_compatible_events(reduced.DT, highlight_events)
    }

    group.grl = split(
        makeGRangesFromDataFrame(
            as.data.frame(transcripts.DT)
        ), 
        transcripts.DT$group_id
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
                group.DT[bump_up_trs, 
                    c("plot_level") := cur_level + 1]
            }
            j = j + match(cur_level, group.DT$plot_level[-seq_len(j)])
            if(is.na(j)) break
        }
        cur_level = cur_level + 1
    }
      
    if(condense_this == TRUE) {
        group.DT[transcripts.DT, on = "group_id", 
            c("group_name", "group_biotype") :=
            list(get("i.gene_name"), get("i.gene_biotype"))]
    } else {
        group.DT[transcripts.DT, on = "group_id", 
            c("group_name", "group_biotype") :=
            list(get("i.transcript_name"), get("i.transcript_biotype"))]      
    }

    group.DT = group.DT[get("end") > view_start & get("start") < view_end]
    group.DT[get("strand") == "+", c("display_name") := 
        paste(get("group_name"), "-", get("group_biotype"), " ->>")]
    group.DT[get("strand") == "-", c("display_name") := 
        paste("<-- ", get("group_name"), "-", get("group_biotype"))]
    group.DT[, c("disp_x") := 0.5 * (get("start") + get("end"))]
    group.DT[get("start") < view_start & get("end") > view_start, 
        c("disp_x") := 0.5 * (view_start + get("end"))]
    group.DT[get("end") > view_end & get("start") < view_end, 
        c("disp_x") := 0.5 * (get("start") + view_end)]
    group.DT[get("start") < view_start & get("end") > view_end, 
        c("disp_x") := 0.5 * (view_start + view_end)]
		
    reduced.DT$group_id = factor(reduced.DT$group_id, unique(group.DT$group_id), ordered = TRUE)
    reduced.DT[group.DT, on = "group_id", 
        c("plot_level") := get("i.plot_level")]
    
    if(missing(highlight_events)) {
        reduced.DT[, c("highlight") := FALSE]
    } else {
        setorderv(reduced.DT, "highlight")
    }
    p = ggplot(reduced.DT)

    if(nrow(subset(as.data.frame(reduced.DT), type = "intron")) > 0) {
        p = p + geom_segment(data = subset(as.data.frame(reduced.DT), type = "intron"), 
            aes(x = start, xend = end, y = plot_level, yend = plot_level,
            color = highlight))
    }
    if(nrow(subset(as.data.frame(reduced.DT), type != "intron")) > 0) {
        p = p + 
            geom_rect(data = subset(as.data.frame(reduced.DT), type != "intron"), 
                aes(xmin = start, xmax = end, 
                ymin = plot_level - 0.1 - 
                    ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0), 
                ymax = plot_level + 0.1 + 
                    ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0),
                fill = highlight
            )
        )
    }

    if(!missing(highlight_events)) {
        p = p + scale_color_manual(values = c("black", "red", "blue")) +
            scale_fill_manual(values = c("black", "red", "blue"))
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
    gp = p + geom_text(data = data.frame(x = anno[["x"]], y = anno[["y"]], text = anno[["text"]]), 
        aes(x = x, y = y, label = text)) + theme_white
    pl = ggplotly(p, source = "plotly_ViewRef", tooltip = "text") %>% 
    layout(
        annotations = anno,
        dragmode = "pan",
        xaxis = list(range = c(view_start, view_end),
            title = paste("Chromosome/Scaffold", view_chr)),
        yaxis = list(range = c(0, 1 + max_plot_level), 
            fixedrange = TRUE)
    )
    
    return(list(gp = gp, pl = pl))		
}

plot_cov_fn <- function(view_chr, view_start, view_end, view_strand,
    norm_event, condition, tracks = list(), track_names = "", se, avail_files,
    transcripts, elems, highlight_events, stack_tracks, graph_mode,
    conf.int = 0.95,
    t_test = FALSE, condensed = FALSE) {

    p_ref = plot_view_ref_fn(
        view_chr, view_start, view_end, 
        transcripts, elems, highlight_events,
        condensed = condensed
    )
    gp_track = list()
    pl_track = list()
  
    cur_zoom = floor(log((view_end - view_start)/50) / log(3))

    data.list = list()
    data.t_test = NULL
    fac = NULL
    
    print(tracks)
    
    if(is_valid(condition) & is_valid(norm_event)) {
        for(i in 1:4) {
            if(length(tracks) >= i && is_valid(tracks[[i]])) {
                track_samples = tracks[[i]]
                colData = SummarizedExperiment::colData(se)
                samples = rownames(colData)[
                    unlist(as.character(colData[, condition]) 
                        == track_samples)]
                event_norms = SummarizedExperiment::assay(
                    se, "Depth")[norm_event,samples]
                samples = samples[event_norms >= 10]
                event_norms = event_norms[event_norms >= 10]

                if(length(avail_files[samples]) > 0 &&
                        all(file.exists(avail_files[samples]))) {

                    df = as.data.frame(GetCoverage_DF(
                        samples, avail_files[samples],
                        view_chr, view_start, view_end, view_strand))
                    # bin anything with cur_zoom > 5
                    df = bin_df(df, max(1, 3^(cur_zoom - 5)))

                    for(todo in seq_len(length(samples))) {
                        df[, samples[todo]] = 
                            df[, samples[todo]] / event_norms[todo]
                    }

                    if(t_test == TRUE) {
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
                    DT = as.data.table(df)
                    DT = DT[, c("x", "mean", "ci", "track")]
                    data.list[[i]] <- DT 
                }
            }
        }
        if(stack_tracks == TRUE) {
            df = as.data.frame(rbindlist(data.list))
            if(nrow(df) > 0) {
                gp_track[[1]] = ggplot() + 
                    geom_ribbon(data = df, alpha = 0.2, 
                        aes(x = x, y = mean, ymin = mean - ci, ymax = mean + ci, fill = track)) +
                    geom_line(data = df, aes(x = x, y = mean, colour = track)) +
                    labs(y = "Stacked Tracks Normalized Coverage") +
                    theme_white_legend
                pl_track[[1]] = ggplotly(gp_track[[1]],
                    tooltip = c("x", "y", "ymin", "ymax", "colour")
                )
                pl_track[[1]] = pl_track[[1]] %>% layout(
                    dragmode = "zoom",
                    yaxis = list(range = c(0, 1 + max(df$mean + df$ci)), fixedrange = TRUE)
                )
                max_tracks = length(unique(df$track))
                for(j in seq_len(max_tracks)) {
                    pl_track[[1]]$x$data[[j]]$showlegend = FALSE
                    pl_track[[1]]$x$data[[j + max_tracks]]$showlegend = TRUE
                    if(!missing(track_names) && length(track_names) >= max_tracks) {
                        pl_track[[1]]$x$data[[j]]$name = track_names[j]
                        pl_track[[1]]$x$data[[j + max_tracks]]$name = track_names[j]                
                    } else {
                        pl_track[[1]]$x$data[[j]]$name = paste(condition, tracks[[j]])
                        pl_track[[1]]$x$data[[j + max_tracks]]$name = paste(condition, tracks[[j]])
                    }
                }
            }
        } else {
            for(i in 1:4) {
                if(length(data.list) >= i && !is.null(data.list[[i]])) {
                    df = as.data.frame(data.list[[i]])
                    gp_track[[i]] = ggplot() + 
                        geom_ribbon(data = df, alpha = 0.2, colour = NA, 
                            aes(x = x, y = mean, ymin = mean - ci, ymax = mean + ci)) +
                        geom_line(data = df, aes(x = x, y = mean)) +
                        labs(y = paste("Track", i, "Normalized Coverage"))
                    pl_track[[i]] = ggplotly(gp_track[[i]],
                        tooltip = c("x", "y", "ymin", "ymax")
                    )						
                    pl_track[[i]] = pl_track[[i]] %>% layout(
                        yaxis = list(
                            range = c(0, 1 + max(df$mean + df$ci)), 
                            fixedrange = TRUE)
                    )
                    pl_track[[i]]$x$data[[1]]$showlegend = FALSE
                    pl_track[[i]]$x$data[[2]]$showlegend = FALSE
                    if(!missing(track_names) && length(track_names) >= i) {
                        pl_track[[i]]$x$data[[1]]$name = track_names[i]
                        pl_track[[i]]$x$data[[2]]$name = track_names[i]
                    } else {
                        pl_track[[i]]$x$data[[1]]$name = paste(condition, tracks[[i]]) # paste("Track",i)
                        pl_track[[i]]$x$data[[2]]$name = paste(condition, tracks[[i]])
                    }
                }
            }
        }
    } else if(!is_valid(condition)){
        for(i in 1:4) {
            if(length(tracks) >= i && is_valid(tracks[[i]])) {
                track_samples = tracks[[i]]
                filename = avail_files[which(names(avail_files) == track_samples)]
                if(length(filename) == 1 && file.exists(filename)) {
                    df = GetCoverage_DF("sample", filename,
                        view_chr, view_start, view_end, view_strand)
                    df = bin_df(df, max(1, 3^(cur_zoom - 5)))
                    data.list[[i]] <- as.data.table(df)
                    if("sample" %in% colnames(df)) {
                        gp_track[[i]] = ggplot() + geom_line(data = df, aes_string(x = "x", y = "sample"))
                        pl_track[[i]] = ggplotly(gp_track[[i]],
                            tooltip = c("x", "y")
                        )
                        pl_track[[i]] = pl_track[[i]] %>% layout(
                            yaxis = list(
                                range = c(0, 1 + max(unlist(df[,"sample"]))), 
                                fixedrange = TRUE,
                                title = paste(track_samples, " Coverage")
                            )
                        )
                        pl_track[[i]]$x$data[[1]]$showlegend = FALSE
                        # pl_track[[i]]$x$data[[2]]$showlegend = TRUE
                        if(!missing(track_names) && length(track_names) >= i) {
                            pl_track[[i]]$x$data[[1]]$name = track_names[i]
                            # pl_track[[i]]$x$data[[2]]$name = track_names[i]
                        } else {
                            pl_track[[i]]$x$data[[1]]$name = track_samples # paste("Track",i)
                            # pl_track[[i]]$x$data[[2]]$name = track_samples
                        }
                    }
                }
            }
        }
    }

    if(t_test == TRUE && !is.null(fac) && length(unique(fac)) == 2) {
        fac = factor(fac)
        t_test = genefilter::rowttests(data.t_test[, -1], fac)
        DT = data.table(x = data.t_test[, 1])
        DT[, c("t_stat") := -log10(t_test$p.value)]
        gp_track[[5]] = ggplot() + geom_line(data = as.data.frame(DT), 
                mapping = aes_string(x = "x", y = "t_stat"))
        pl_track[[5]] = ggplotly(gp_track[[5]],
            # labs(y = paste("Pairwise T-test -log10(p)")),
            tooltip = c("x", "y")
        )
        pl_track[[5]] = pl_track[[5]] %>% layout(
            yaxis = list(c(0, 1 + max(DT$t_stat)), 
                fixedrange = TRUE,
                title = paste("T-test -log10(p)")
            )
        )
        pl_track[[5]]$x$data[[1]]$showlegend = FALSE
    }

    plot_tracks = pl_track[unlist(lapply(pl_track, function(x) !is.null(x)))]

    for(i in seq_len(length(p_ref$pl$x$data))) {
        p_ref$pl$x$data[[i]]$showlegend = FALSE
    }
    
    plot_tracks[[length(plot_tracks) + 1]] = p_ref$pl
    
    gp_track[[6]] = p_ref$gp
    
    # Work out which x axis ticks to use, based on zoom level
    # cur_zoom = 0 {50, 150}    10
    # cur_zoom = 1 {151, 450}   50
    # cur_zoom = 2 {451, 1350}  100
    # cur_zoom = 3 {1351, 1350}  100
    view_range = view_end - view_start
    min_tick_size = view_range / 15
    # round up tick size to nearest 1, 2, 5
    if(min_tick_size / 10 ^ floor(log10(min_tick_size)) > 5) {
        tick_size = (10 ^ floor(log10(a))) * (min_tick_size) * 10
    } else if(min_tick_size / 10 ^ floor(log10(min_tick_size)) > 2) {
        tick_size = (10 ^ floor(log10(a))) * (min_tick_size) * 5
    } else {
        tick_size = (10 ^ floor(log10(a))) * (min_tick_size) * 2
    }
    first_tick = ceiling(view_start / tick_size) * tick_size
    
    final_plot = subplot(plot_tracks, nrows = length(plot_tracks), 
        shareX = TRUE, titleY = TRUE) %>%
      layout(
        xaxis = list(
          dtick = tick_size, 
          tick0 = first_tick, 
          tickmode = "linear"
      ))

    if(graph_mode == "Pan") {
        final_plot = final_plot %>% 
            layout(dragmode = "pan")
    } else if(graph_mode == "Zoom") {
        final_plot = final_plot %>% 
            layout(dragmode = "zoom")   
    } else if(graph_mode == "Movable Labels") {
        final_plot = final_plot %>% 
            layout(dragmode = FALSE) %>%
            config(editable = TRUE)
    }

    # ggplot equivalent: list of ggplots. Allows advanced end-users to apply final edits to ggplots
    
    return(list(ggplot = gp_track, final_plot = final_plot))

}

get_default_filters <- function() {
    filterUnit <- list()
    filterUnit$filterVars = list(
        1, 20, "All", "(none)", 80
    )
    names(filterUnit$filterVars) = 
        c("maximum", "minDepth", "minCond", "condition", "pcTRUE")
    filterUnit$filterClass = "(none)"
    filterUnit$filterType = "(none)"
    
    filters = list()
    for(i in seq_len(8)) {
        filters[[i]] = filterUnit
    }

    filters[[1]]$filterClass = "Data"    
    filters[[1]]$filterType = "Depth"    
    filters[[1]]$filterVars$minimum =  20

    filters[[2]]$filterClass = "Data"    
    filters[[2]]$filterType = "Coverage"
    filters[[2]]$filterVars$minimum =  90
    filters[[2]]$filterVars$minDepth =  5
    filters[[2]]$filterVars$EventTypes =  "IR"

    filters[[3]]$filterClass = "Data"    
    filters[[3]]$filterType = "Coverage"
    filters[[3]]$filterVars$minimum =  60
    filters[[3]]$filterVars$minDepth =  20
    filters[[3]]$filterVars$EventTypes =
        c("MXE", "SE", "AFE", "ALE", "A5SS", "A3SS")

    filters[[4]]$filterClass = "Data"    
    filters[[4]]$filterType = "Consistency"
    filters[[4]]$filterVars$maximum =  2
    filters[[4]]$filterVars$minDepth =  20
    filters[[4]]$filterVars$EventTypes = c("MXE", "SE")

    return(filters)
}

make_matrix <- function(se, event_list, sample_list, method, depth_threshold = 10, logit_max = 5, na.percent.max = 0.1) {

	inc = SummarizedExperiment::assay(se, "Included")[event_list, sample_list]
	exc = SummarizedExperiment::assay(se, "Excluded")[event_list, sample_list]
    mat = inc/(inc + exc)
    mat[inc + exc < depth_threshold] = NA
    mat = mat[rowSums(is.na(mat)) < na.percent.max * ncol(mat),]	
	if(method == "PSI") {
		# essentially M/Cov
		return(mat)
	} else if(method == "logit") {
		mat = boot::logit(mat)
		mat[mat > logit_max] = logit_max
		mat[mat < -logit_max] = -logit_max
		return(mat)
	} else if(method == "Z-score") {
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

