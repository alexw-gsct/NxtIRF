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
  x = start:end
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


plot_view_ref_fn <- function(view_chr, view_start, view_end, 
    transcripts, elems, highlight_events, condensed = FALSE) {
      
   # transcript_support_level <- transcript_id <- group_id <- type <- gene_id <- plot_level <- 
	 # i.gene_name <- i.gene_biotype <- i.transcript_name <- i.transcript_biotype <- display_name <- 
	 # group_name <- group_biotype <- disp_x <- i.plot_level <- NULL		
			
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
        get("transcript_id") %in% transcripts.DT$transcript_id
    ]
    if(condensed != TRUE & nrow(transcripts.DT) <= 100) {
        condense_this = FALSE
        transcripts.DT[, c("group_id") := get("transcript_id")]
        screen.DT[, c("group_id") := get("transcript_id")]
        reduced.DT = copy(screen.DT)
        reduced.DT[get("type") %in% c("CDS", "start_codon", "stop_codon"), c("type") := "CDS"]
        reduced.DT[get("type") != "CDS", c("type") := "exon"]
    } else {
        condense_this = TRUE
        transcripts.DT[, c("group_id") := get("gene_id")]     
        screen.DT[transcripts.DT, on = "transcript_id", 
            c("group_id") := get("gene_id")]
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
                as.data.frame(screen.DT[
                        get("type") %in% c("CDS", "start_codon", "stop_codon")
                    ]
                )
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
    introns.DT[, c("type") := "intron"]
    setnames(introns.DT, "group_name", "group_id")

    reduced.DT = rbind(reduced.DT[, 
        c("seqnames", "start", "end", "strand", "type", "group_id")], 
    introns.DT[, 
        c("seqnames", "start", "end", "strand", "type", "group_id")])

    # Highlight events here
    # highlight_events is of syntax chrX:10000-11000/-
    if(!missing(highlight_events)) {
        event.df = as.data.table(NxtIRF.CoordToInt(highlight_events))
        event.df = event.df[get("seqnames") == view_chr]
        
        
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
      
    p = ggplot(reduced.DT)

    if(nrow(subset(as.data.frame(reduced.DT), type = "intron")) > 0) {
        p = p + geom_segment(data = subset(as.data.frame(reduced.DT), type = "intron"), 
            aes(x = start, xend = end, y = plot_level, yend = plot_level))
    }
    if(nrow(subset(as.data.frame(reduced.DT), type != "intron")) > 0) {
        p = p + 
            geom_rect(data = subset(as.data.frame(reduced.DT), type != "intron"), 
                aes(xmin = start, xmax = end, 
                ymin = plot_level - 0.1 - 
                    ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0), 
                ymax = plot_level + 0.1 + 
                    ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 0.1, 0)
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
    
    pl = ggplotly(p, source = "plotly_ViewRef", tooltip = "text") %>% 
    layout(
        annotations = anno,
        dragmode = "pan",
        xaxis = list(range = c(view_start, view_end)),
        yaxis = list(range = c(0, 1 + max_plot_level), fixedrange = TRUE)
    )
    
    return(pl)			
}

plot_cov_fn <- function(view_chr, view_start, view_end, view_strand,
    norm_event, condition, tracks = list(), se, avail_files,
    transcripts, elems, highlight_events, stack_tracks, graph_mode,
    conf.int = 0.95,
    t_test = FALSE, condensed = FALSE) {

    p_ref = plot_view_ref_fn(
        view_chr, view_start, view_end, 
        transcripts, elems
    )
    p_track = list()
  
    cur_zoom = floor(log((view_end - view_start)/50) / log(3))

    data.list = list()
    data.t_test = NULL
    fac = NULL
    
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
                p_track[[1]] = ggplotly(
                    ggplot(df, aes(x = x)) + 
                    geom_ribbon(alpha = 0.2, 
                        aes(y = mean, ymin = mean - ci, ymax = mean + ci, fill = track)) +
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
                    df = as.data.frame(data.list[[i]])
                    p_track[[i]] = ggplotly(
                        ggplot(df, aes(x = x)) + 
                        geom_ribbon(alpha = 0.2, colour = NA, 
                            aes(y = mean, ymin = mean - ci, ymax = mean + ci)) +
                        geom_line(aes(y = mean)) +
                        labs(y = paste("Track", i, "Normalized Coverage")),
                        tooltip = c("x", "y", "ymin", "ymax")
                    )						
                    p_track[[i]] = p_track[[i]] %>% layout(
                        yaxis = list(
                            range = c(0, 1 + max(df$mean + df$ci)), 
                            fixedrange = TRUE)
                    )
                }
            }
        }
    } else if(!is_valid(condition)){
        for(i in 1:4) {
            if(length(tracks) >= i && is_valid(tracks[[i]])) {
                track_samples = tracks[[i]]
                if(length(avail_files[track_samples]) > 0 &&
                        all(file.exists(avail_files[track_samples]))) {
                    df = GetCoverage_DF(track_samples, avail_files[track_samples],
                        view_chr, view_start, view_end, view_strand)
                    df = bin_df(df, max(1, 3^(cur_zoom - 5)))
                    data.list[[i]] <- as.data.table(df)

                    p_track[[i]] = ggplotly(
                        ggplot(df, aes_string(x = "x", y = track_samples)) + geom_line() +
                        labs(y = paste(track_samples, " Coverage")),
                        tooltip = c("x", "y")
                    )
                    p_track[[i]] = p_track[[i]] %>% layout(
                        yaxis = list(
                            range = c(0, 1 + max(unlist(df[,track_samples]))), 
                            fixedrange = TRUE
                        )
                    )
                }
            }
        }
    }

    if(t_test == TRUE && !is.null(fac) && length(unique(fac)) == 2) {
        fac = factor(fac)
        t_test = genefilter::rowttests(data.t_test[, -1], fac)
        DT = data.table(x = data.t_test[, 1])
        DT[, c("t_stat") := -log10(t_test$p.value)]

        p_track[[5]] = ggplotly(
            ggplot(as.data.frame(DT), 
                aes_string(x = "x", y = "t_stat")) + geom_line() +
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
    
    return(final_plot)

}