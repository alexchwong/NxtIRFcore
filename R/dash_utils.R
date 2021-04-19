GetCoverage_DF <- function(samples, files, seqname, start, end, strand) {
    covData = list()
    for(i in seq_len(length(files))) {
        cov = GetCoverage(files[i], seqname, start - 1, end, 
            ifelse(strand == "+", 0, ifelse(strand == "-", 1, 2)))
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
    if(!all(file.exists(c(junc_fst, junc_fst_index)))) {
        stop(paste("In get_psi(),",
            "Some junction fst files do not exist"
        ), call. = FALSE)
    }
        
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
        OL2 = findOverlaps(gr, introns.gr, type = "equal")
        introns[OL2@to, c("highlight") := 2]

    } else if(length(highlight_events) == 2) {
        # This is AS
        AS_count = 1;
        for(event in highlight_events) {
            gr = NxtIRF.CoordToGR(event)
            introns.gr = .grDT(introns)
            OL = findOverlaps(gr, introns.gr, type = "equal")
            introns[OL@to, c("highlight") := as.character(AS_count)]

            OL_s1 = findOverlaps(gr[1], introns.gr, type = "equal")
            tr1 = unique(introns$transcript_id[OL_s1@to])
            if(length(gr) == 2) {
                OL_s2 = findOverlaps(gr[2], introns.gr, type = "equal")
                tr1 = unique(intersect(tr1, introns$transcript_id[OL_s2@to]))
            }
            tr_filter = c(tr_filter, tr1)
            coord_keys = c(start(gr[1]) - 1, end(gr[1]) + 1)
            if(length(gr) == 2) {
                coord_keys = c(coord_keys,
                    start(gr[2]) - 1, end(gr[2]) + 1)
            }
            exons[get("transcript_id") %in% tr1 & 
                (get("start") %in% coord_keys | get("end")%in% coord_keys),
                c("highlight") := as.character(AS_count)]
            AS_count = AS_count + 1
        }  
    }
    return(rbind(introns, exons, misc))
}

plot_view_ref_fn <- function(view_chr, view_start, view_end, 
    transcripts, elems, highlight_events, condensed = FALSE,
    selected_transcripts) {

    data_start = view_start - (view_end - view_start)
    data_end = view_end + (view_end - view_start)

    transcripts.DT = transcripts[get("seqnames") == view_chr]
    transcripts.DT = transcripts.DT[
        get("start") <= data_end & 
        get("end") >= data_start]
    setorderv(transcripts.DT, c("transcript_support_level", "width"))
    # filter transcripts if applicable
    if(!missing(selected_transcripts)) {
        transcripts.DT = transcripts.DT[
            get("transcript_id") %in% selected_transcripts |
            get("transcript_name") %in% selected_transcripts]
    }
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
    reduced.DT[get("type") %in% c("CDS", "start_codon", "stop_codon"), 
        c("type") := "CDS"]
    reduced.DT[get("type") != "CDS", c("type") := "exon"]
    
    # add introns to reduced.DT
    introns.DT = as.data.table(.grlGaps(
        split(makeGRangesFromDataFrame(as.data.frame(reduced.DT)),
            reduced.DT$transcript_id)
    ))
    introns.DT[, c("type") := "intron"]
    setnames(introns.DT, "group_name", "transcript_id")
    introns.DT[reduced.DT, on = "transcript_id",
        "group_id" := get("i.group_id")]

    filter_cols = c("seqnames", "start", "end", "strand", 
        "type", "group_id", "transcript_id")
    reduced.DT = rbind(reduced.DT[, filter_cols, with = FALSE], 
        introns.DT[, filter_cols, with = FALSE])

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
            bump_up_trs = bump_up_trs[
                group.DT$plot_level[bump_up_trs] == cur_level]
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

    reduced.DT$group_id = factor(reduced.DT$group_id, 
        unique(group.DT$group_id), ordered = TRUE)
    reduced.DT[group.DT, on = "group_id", 
        c("plot_level") := get("i.plot_level")]
    
    if(missing(highlight_events)) {
        reduced.DT[, c("highlight") := FALSE]
    } else {
        setorderv(reduced.DT, "highlight")
    }
    p = ggplot(reduced.DT)

    if(nrow(subset(as.data.frame(reduced.DT), type = "intron")) > 0) {
        p = p + geom_segment(data = subset(as.data.frame(reduced.DT), 
                type = "intron"), 
            aes(x = get("start"), xend = get("end"), 
                y = get("plot_level"), yend = get("plot_level"),
            color = get("highlight")))
    }
    if(nrow(subset(as.data.frame(reduced.DT), type != "intron")) > 0) {
        p = p + 
            geom_rect(data = subset(as.data.frame(reduced.DT), 
                    type != "intron"), 
                aes(xmin = get("start"), xmax = get("end"), 
                ymin = get("plot_level") - 0.1 - 
                    ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 
                        0.1, 0), 
                ymax = get("plot_level") + 0.1 + 
                    ifelse(type %in% c("CDS", "start_codon", "stop_codon"), 
                        0.1, 0),
                fill = get("highlight")
            )
        )
    }

    if(!missing(highlight_events)) {
        p = p + scale_color_manual(values = c("black", "blue", "red")) +
            scale_fill_manual(values = c("black", "blue", "red"))
    }

    p = p + theme_white_legend +
        theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
        labs(fill = "")

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
    gp = p + geom_text(data = data.frame(x = anno[["x"]], y = anno[["y"]], 
            text = anno[["text"]]), 
        aes(x = get("x"), y = get("y"), label = get("text"))) + 
        coord_cartesian(xlim = c(view_start, view_end))
    pl = ggplotly(p, 
        # source = "plotly_ViewRef", 
        tooltip = "text") %>% 
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
    transcripts, elems, highlight_events, selected_transcripts, 
    stack_tracks, graph_mode,
    conf.int = 0.95,
    t_test = FALSE, condensed = FALSE) {

    if(!missing(selected_transcripts)) {
        p_ref = plot_view_ref_fn(
            view_chr, view_start, view_end, 
            transcripts, elems, highlight_events,
            condensed = condensed,
            selected_transcripts = selected_transcripts
        )   
    } else {
        p_ref = plot_view_ref_fn(
            view_chr, view_start, view_end, 
            transcripts, elems, highlight_events,
            condensed = condensed
        )    
    }
    gp_track = list()
    pl_track = list()
    
    cur_zoom = floor(log((view_end - view_start)/50) / log(3))

    data.list = list()
    data.t_test = NULL
    fac = NULL
    
    if(is_valid(condition) & is_valid(norm_event)) {
        max_tracks = 0
        for(i in seq_len(4)) {
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
                    message(paste("Group GetCoverage performed for", condition))
                    for(todo in seq_len(length(samples))) {
                        df[, samples[todo]] = 
                            df[, samples[todo]] / event_norms[todo]
                    }

                    if(t_test == TRUE) {
                        if(is.null(data.t_test)) {
                            data.t_test <- as.matrix(df)
                            fac = rep(as.character(i), ncol(df) - 1)
                        } else {
                            data.t_test <- cbind(
                                data.t_test, as.matrix(df[, -1]))
                            fac = c(fac, rep(as.character(i), ncol(df) - 1))
                        }
                    }

                    df$mean = rowMeans(as.matrix(df[,samples]))
                    df$sd = rowSds(as.matrix(df[,samples]))
                    n = length(samples)
                    df$ci = qt((1 + conf.int)/2,df = n-1) * df$sd / sqrt(n)

                    if(length(track_names) == length(tracks)) {
                        df$track = track_names[i]
                    } else {
                        df$track = as.character(i)
                    }
                    DT = as.data.table(df)
                    DT = DT[, c("x", "mean", "ci", "track")]
                    data.list[[i]] <- DT 
                    max_tracks = max_tracks + 1
                }
            }
        }
        if(stack_tracks == TRUE) {
            df = as.data.frame(rbindlist(data.list))
            if(nrow(df) > 0) {
                if(length(track_names) == length(tracks)) {
                    df$track = factor(df$track, track_names)
                }
                gp_track[[1]] = ggplot() + 
                    geom_hline(yintercept = 0) +
                    geom_ribbon(data = df, alpha = 0.2, 
                        aes(x = get("x"), y = get("mean"), 
                        ymin = get("mean") - get("ci"), 
                        ymax = get("mean") + get("ci"), 
                        fill = get("track"))) +
                    geom_line(data = df, aes(x = get("x"), 
                        y = get("mean"), colour = get("track"))) +
                    labs(y = "Normalized Coverage") +
                    theme_white_legend
                pl_track[[1]] = ggplotly(gp_track[[1]],
                    tooltip = c("x", "y", "ymin", "ymax", "colour")
                )
                pl_track[[1]] = pl_track[[1]] %>% layout(
                    dragmode = "zoom",
                    yaxis = list(rangemode = "tozero")
                )
                for(j in seq_len(max_tracks)) {
                    pl_track[[1]]$x$data[[1 + j]]$showlegend = FALSE
                    pl_track[[1]]$x$data[[1 + j + max_tracks]]$showlegend = TRUE
                    if(!missing(track_names) && 
                            length(track_names) >= max_tracks) {
                        pl_track[[1]]$x$data[[1 + j]]$name = track_names[j]
                        pl_track[[1]]$x$data[[1 + j + max_tracks]]$name = 
                            track_names[j]                
                    } else {
                        pl_track[[1]]$x$data[[1 + j]]$name = 
                            paste(condition, tracks[[j]])
                        pl_track[[1]]$x$data[[1 + j + max_tracks]]$name = 
                            paste(condition, tracks[[j]])
                    }
                }
            }
        } else {
            for(i in seq_len(4)) {
                if(length(data.list) >= i && !is.null(data.list[[i]])) {
                    df = as.data.frame(data.list[[i]])
                    gp_track[[i]] = ggplot() + 
                        geom_hline(yintercept = 0) +
                        geom_ribbon(data = df, alpha = 0.2, colour = NA, 
                            aes(x = get("x"), y = get("mean"), 
                                ymin = get("mean") - get("ci"), 
                                ymax = get("mean") + get("ci"))) +
                        geom_line(data = df, 
                            aes(x = get("x"), y = get("mean"))) +
                        labs(y = paste(condition, tracks[[i]])) + 
                        theme_white_legend
                    pl_track[[i]] = ggplotly(gp_track[[i]],
                        tooltip = c("x", "y", "ymin", "ymax")
                    )
                    pl_track[[i]] = pl_track[[i]] %>% layout(
                        yaxis = list(rangemode = "tozero")
                    )
                    pl_track[[i]]$x$data[[2]]$showlegend = FALSE
                    pl_track[[i]]$x$data[[3]]$showlegend = FALSE
                    if(!missing(track_names) && length(track_names) >= i) {
                        pl_track[[i]]$x$data[[2]]$name = track_names[i]
                        pl_track[[i]]$x$data[[3]]$name = track_names[i]
                    } else {
                        pl_track[[i]]$x$data[[2]]$name = paste(
                            condition, tracks[[i]])
                        pl_track[[i]]$x$data[[3]]$name = 
                            paste(condition, tracks[[i]])
                    }
                }
            }
        }
    } else if(!is_valid(condition)){
        for(i in seq_len(4)) {
            if(length(tracks) >= i && is_valid(tracks[[i]])) {
                track_samples = tracks[[i]]
                filename = avail_files[which(
                    names(avail_files) == track_samples)]
                if(length(filename) == 1 && file.exists(filename)) {
                    df = GetCoverage_DF("sample", filename,
                        view_chr, view_start, view_end, view_strand)
                    df = bin_df(df, max(1, 3^(cur_zoom - 5)))
                    data.list[[i]] <- as.data.table(df)
                    if("sample" %in% colnames(df)) {
                        gp_track[[i]] = ggplot() + 
                        geom_hline(yintercept = 0) +
                        geom_line(data = df, 
                            aes(x = get("x"), y = get("sample"))) +
                        theme_white_legend
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
                        pl_track[[i]]$x$data[[2]]$showlegend = FALSE
                        if(!missing(track_names) && length(track_names) >= i) {
                            pl_track[[i]]$x$data[[2]]$name = track_names[i]
                        } else {
                            pl_track[[i]]$x$data[[2]]$name = track_samples
                        }
                    }
                }
            }
        }
    }

    if(t_test == TRUE && !is.null(fac) && length(unique(fac)) == 2) {
        fac = factor(fac)
        t_test = genefilter::rowttests(data.t_test[, -1], fac)
        message(paste("rowttests performed for", condition))

        DT = data.table(x = data.t_test[, 1])
        DT[, c("t_stat") := -log10(t_test$p.value)]
        gp_track[[5]] = ggplot() + 
            geom_hline(yintercept = 0) +
            geom_line(data = as.data.frame(DT), 
                mapping = aes(x = get("x"), y = get("t_stat"))) +
                theme_white_legend
        pl_track[[5]] = ggplotly(gp_track[[5]],
            # labs(y = paste("Pairwise T-test -log10(p)")),
            tooltip = c("x", "y")
        )
        pl_track[[5]] = pl_track[[5]] %>% layout(
            yaxis = list(
                # c(0, 1 + max(DT$t_stat)), 
                # fixedrange = TRUE,
                rangemode = "tozero",
                title = paste("T-test -log10(p)")
            )
        )
        pl_track[[5]]$x$data[[2]]$showlegend = FALSE
    }

    plot_tracks = pl_track[unlist(lapply(pl_track, function(x) !is.null(x)))]

    for(i in seq_len(length(p_ref$pl$x$data))) {
        p_ref$pl$x$data[[i]]$showlegend = FALSE
    }
    
    plot_tracks[[length(plot_tracks) + 1]] = p_ref$pl
    
    gp_track[[6]] = p_ref$gp
    
    # Work out which x axis ticks to use, based on zoom level
    view_range = view_end - view_start
    min_tick_size = view_range / 15
    tick_order_magn = 10 ^ floor(log10(min_tick_size))
    # round up tick size to nearest 1, 2, 5
    if(min_tick_size / tick_order_magn > 5) {
        tick_size = tick_order_magn * 20
    } else if(min_tick_size / tick_order_magn > 2) {
        tick_size = tick_order_magn * 10
    } else {
        tick_size = tick_order_magn * 5
    }
    first_tick = ceiling(view_start / tick_size) * tick_size
    
    final_plot = subplot(plot_tracks, nrows = length(plot_tracks), 
        shareX = TRUE, titleY = TRUE) %>%
        layout(
            xaxis = list(
                dtick = tick_size, 
                tick0 = first_tick, 
                tickmode = "linear"
            )
        )

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
    
    # Zero all but last subplot
    n_plots = length(plot_tracks) - 1
    
    final_plot = final_plot %>% layout(yaxis = list(
        rangemode = "tozero", tick0 = 0))
    if(n_plots > 1) final_plot = final_plot %>% layout(yaxis2 = 
        list(rangemode = "tozero", tick0 = 0))
    if(n_plots > 2) final_plot = final_plot %>% layout(yaxis3 = 
        list(rangemode = "tozero", tick0 = 0))
    if(n_plots > 3) final_plot = final_plot %>% layout(yaxis4 = 
        list(rangemode = "tozero", tick0 = 0))
    if(n_plots > 4) final_plot = final_plot %>% layout(yaxis5 = 
        list(rangemode = "tozero", tick0 = 0))

    # ggplot equivalent: list of ggplots. 
    # Allows advanced end-users to apply final edits to ggplots
    
    message("Cov Plot finished")
    return(list(ggplot = gp_track, final_plot = final_plot))

}

# Exported wrapper functions to plot coverage plots from command line:

#' Generate plotly / ggplot RNA-seq coverage plots from command line
#'
#' This function generates a coverage plot illustrating differential expression
#' of intron retention or alternative splicing (differential exon coverage). 
#' It can normalise coverage of samples per condition using NxtIRF's SpliceOver
#' parameter that normalises the flanking exon boundaries to 1.
#' 
#' @param se A SummarizedExperiment object. It must contain the Event to be 
#'   displayed
#' @param Event The `EventName` or the IR / alternative splicing event to be 
#'   normalised. Valid names are all entries within rownames(se)
#' @param cov_data A list containing all the data required by this function. 
#'   Generate this data using prepare_covplot_data()
#' @param strand Whether to show coverage of both strands "*" (default), or
#'   from the "+" or "-" strand only.
#' @param tracks The names of individual samples (if condition is not set),
#'   or the names of the different conditions to be plotted.
#' @param track_names The names of the tracks to be displayed. Defaults to 
#'   `tracks`
#' @param condition The name of the column (of `colData(se)` containing the 
#'   conditions to be contrasted. If this is not set, `tracks` are assumed to be
#'   names of individual samples
#' @param selected_transcripts Transcript ID or transcript names of transcripts
#'   to be displayed on the gene annotation track. Useful to remove overlapping
#'   transcripts that are not relevant to the samples being displayed.
#' @param condense_tracks Whether to collapse the transcript tracks by gene.
#' @param stack_tracks Whether to graph all the conditions on a single coverage
#'   track. If set to true, each condition will be displayed in a different
#'   colour
#' @param t_test Whether to perform a pair-wise T-test. Only used if there are
#'   TWO condition tracks.
#' @param norm_event Whether to normalise by an event different to that given
#'   in "Event". The difference between this and Event is that the genomic
#'   coordinates are only centered around Event. norm_event overrides Event
#'   if both are given.
#' @param stack_tracks Whether to graph all the conditions on a single coverage
#'   track. If set to true, each condition
#' @param zoom_factor Zoom out from event. Each level of zoom zooms out by a
#'   factor of 3.
#' @param bases_flanking How many bases flanking the zoomed window. Useful when
#'    used in conjunction with zoom_factor == 0
#' @param Gene Whether to use the range for the given Gene. If given and valid,
#'   overrides Event (but `Event` or `norm_event` will be used to normalise by
#'  condition). Valid Gene entries include gene_id (Ensembl ID) or gene_name
#'  (Gene Symbol)
#' @param seqnames Use chromosome. Only used if `seqnames`, `start` and `end`
#'   are also given and valid.
#'   If these coordinates are set, they will override Event and Gene, although
#'   Event / norm_event is still required for normalising by condition
#' @param start Overrides with given start coordinate. See `seqnames`
#' @param end Overrides with given end coordinate. See `seqnames`
#' 
#' @return A list containing two objects. final_plot is the plotly object. 
#'   ggplot is a list of ggplot tracks.\cr\cr
#'   `ggplot[[n]]` is the nth track (this function supports up to 4 tracks).
#'   \cr\cr
#'   `ggplot[[5]]` contains the T-test track if valid.\cr\cr
#'   `ggplot[[6]]` always contains the genome track.\cr\cr
#' @examples
#' se = NxtIRF_example_NxtSE()
#' 
#' colData(se)$treatment = rep(c("A", "B"), each = 3)
#'
#' Plot_Coverage(se, Event = rowData(se)$EventName[1], 
#'   cov_data = ref(se), tracks = colnames(se)[1:2])
#' @md
#' @export
Plot_Coverage <- function(se, Event, cov_data,
    strand = c("*", "+", "-"),
    tracks,                 
    track_names = tracks,
    condition,              
    selected_transcripts,   
    condense_tracks = FALSE,
    stack_tracks = TRUE,  
    t_test = TRUE,  
    norm_event,
    zoom_factor = 1,        
    bases_flanking = 100,   
    Gene,
    seqnames, start, end   # Optional
    ) {

# Assertions
    args = as.list(match.call())
    do.call(.plot_cov_validate_args, args)

    strand = match.arg(strand)
    if(strand == "") strand <- "*"

    # Prepare zoom window
    if((!missing(seqnames) & !missing(start) & !missing(end))) {
        view_chr = as.character(seqnames)
        view_start = start
        view_end = end
    } else if(!missing(Gene)) {        
        if(Gene %in% cov_data$gene_list$gene_id) {
            gene.df = as.data.frame(
                cov_data$gene_list[get("gene_id") == get("Gene")])
        } else {
            gene.df = as.data.frame(
                cov_data$gene_list[get("gene_name") == get("Gene")])
        }
        view_chr = as.character(gene.df$seqnames)
        view_start = gene.df$start
        view_end = gene.df$end        
    } else {
        rowData = as.data.frame(rowData(se))
        rowData = rowData[Event,]
        view_chr = tstrsplit(rowData$EventRegion, split=":")[[1]]
        temp1 = tstrsplit(rowData$EventRegion, split="/")
        temp2 = tstrsplit(temp1[[1]], split=":")[[2]]
        view_start = as.numeric(tstrsplit(temp2, split="-")[[1]])
        view_end = as.numeric(tstrsplit(temp2, split="-")[[2]])
    }

    view_center = (view_start + view_end) / 2
    view_length = view_end - view_start

    # Apply zoom 
    new_view_length = view_length * 3 ^ zoom_factor + 2 * bases_flanking
    view_start = round(view_center - new_view_length / 2)
    view_end = round(view_center + new_view_length / 2)
    
    # Validate genomic window and shift if invalid
    if(view_start < 1) view_start = 1
    seqInfo = cov_data$seqInfo[view_chr]
    seqmax = GenomeInfoDb::seqlengths(seqInfo)
    if(view_end > seqmax) view_end = seqmax - 1
    
    if(missing(norm_event)) {
        if(!missing(Event)) {
            norm_event = Event
        } else {
            norm_event = ""       
        }
    }
    if(missing(condition)) {
        condition = ""
    }
    
    # Last sanity check view_chr, view_start, view_end here:
    
    args = list(
        view_chr = view_chr,
        view_start = view_start,
        view_end = view_end,
        view_strand = strand,
        norm_event = norm_event,
        condition = condition,
        tracks = as.list(tracks),
        track_names = track_names,
        se = se,
        avail_files = covfile(se),
        transcripts = cov_data$transcripts.DT,
        elems = cov_data$elem.DT,
        stack_tracks = stack_tracks,
        graph_mode = "Pan",
        conf.int = 0.95,
        t_test = t_test,
        condensed = condense_tracks
    )
    if(norm_event != "") {
        events_to_highlight = list()
        rowData = as.data.frame(SummarizedExperiment::rowData(se))

        if(rowData$EventType[match(norm_event, rowData$EventName)] 
            %in% c("MXE", "SE")) {
            events_to_highlight[[1]] = c(
                rowData$Event1a[match(norm_event, rowData$EventName)],
                rowData$Event2a[match(norm_event, rowData$EventName)])
        } else {
            events_to_highlight[[1]] = rowData$Event1a[
                match(norm_event, rowData$EventName)]
        }
        if(rowData$EventType[match(norm_event, rowData$EventName)] 
            %in% c("MXE")) {
            events_to_highlight[[2]] = c(
                rowData$Event1b[match(norm_event, rowData$EventName)],
                rowData$Event2b[match(norm_event, rowData$EventName)])
        } else if(rowData$EventType[match(norm_event, rowData$EventName)] 
            %in% c("SE", "A3SS", "A5SS", "ALE", "AFE")){
            events_to_highlight[[2]] = rowData$Event1b[
                match(norm_event, rowData$EventName)]
        }
        args$highlight_events = events_to_highlight
    }
    if(!missing(selected_transcripts)) {
        args$selected_transcripts = selected_transcripts
    }
    return(
        do.call(plot_cov_fn, args)
    )
}

.plot_cov_validate_args <- function(se, Event, cov_data,
    strand = c("*", "+", "-"),
    tracks,                 
    track_names = tracks,
    condition,              
    selected_transcripts,   
    condense_tracks = FALSE,
    stack_tracks = TRUE,  
    t_test = TRUE,  
    norm_event,
    zoom_factor = 1,        
    bases_flanking = 100,   
    Gene,
    seqnames, start, end
) {
    # Requires all cov files in the SE to be present
    if(!all(colnames(se) %in% names(covfile(se)))) {
        stop(paste("In Plot_Coverage,",
            "Some samples do not have COV files.",
            "Make sure metadata(se)$cov_file is a named vector",
            "containing COV file paths,",
            "and names(metadata(se)$cov_file) correspond to sample names"
        ), call. = FALSE)
    }
    if(!all(c("seqInfo", "gene_list", "elem.DT", "transcripts.DT") %in% 
            names(cov_data))) {
        stop(paste("In Plot_Coverage,",
            "cov_data must be a valid object created by prepare_covplot_data()"
        ), call. = FALSE)
    }
    # Check condition and tracks
    if(!missing(condition)) {
        if(!(condition %in% names(colData(se)))) {
            stop(paste("In Plot_Coverage,",
                "condition must be a valid column name in colData(se)"
            ), call. = FALSE)
        }
        condition_options = unique(colData(se)[, condition])
        if(!all(tracks %in% condition_options)) {
            stop(paste("In Plot_Coverage,",
                "some tracks do not match valid condition names in", condition
            ), call. = FALSE)
        }   
    } else {
        if(!all(tracks %in% colnames(se))) {
            stop(paste("In Plot_Coverage,",
                "some tracks do not match valid sample names in se"
            ), call. = FALSE)
        }
    }
    # Check we know where to plot
    if(missing(Event) & missing(Gene) & 
            (missing(seqnames) | missing(start) | missing(end))
    ) {
        stop(paste("In Plot_Coverage,",
            "Event or Gene cannot be empty, unless coordinates are provided"
        ), call. = FALSE)
    } else if((!missing(seqnames) & !missing(start) & !missing(end))) {
        view_chr = as.character(seqnames)
        view_start = start
        view_end = end
    } else if(!missing(Gene)) {
        if(!(Gene %in% cov_data$gene_list$gene_id) & 
                !(Gene %in% cov_data$gene_list$gene_name)) {
            stop(paste("In Plot_Coverage,",
                Gene, "is not a valid gene symbol or Ensembl gene id"
            ), call. = FALSE)
        }       
        if(!(Gene %in% cov_data$gene_list$gene_id)) {
            gene.df = as.data.frame(
                cov_data$gene_list[get("gene_name") == get("Gene")])
            if(nrow(gene.df) != 1) {
                stop(paste("In Plot_Coverage,",
                    Gene, "is an ambiguous name referring to 2 or more genes.",
                    "Please provide its gene_id instead"
                ), call. = FALSE)
            }
        } else {
            gene.df = as.data.frame(
                cov_data$gene_list[get("gene_id") == get("Gene")])        
        }
        view_chr = as.character(gene.df$seqnames)
        view_start = gene.df$start
        view_end = gene.df$end             
    } else {
        rowData = as.data.frame(rowData(se))
        if(!(Event %in% rownames(rowData))) {
            stop(paste("In Plot_Coverage,",
                Event, "is not a valid IR or alternate splicing event",
                "in rowData(se)"
            ), call. = FALSE)
        }
        rowData = rowData[Event,]
        view_chr = tstrsplit(rowData$EventRegion, split=":")[[1]]
        temp1 = tstrsplit(rowData$EventRegion, split="/")
        temp2 = tstrsplit(temp1[[1]], split=":")[[2]]
        view_start = as.numeric(tstrsplit(temp2, split="-")[[1]])
        view_end = as.numeric(tstrsplit(temp2, split="-")[[2]])
    }
    view_center = (view_start + view_end) / 2
    view_length = view_end - view_start    
    
    if(!(view_chr %in% names(cov_data$seqInfo))) {
        stop(paste("In Plot_Coverage,", view_chr, 
            "is not a valid chromosome reference name in the given genome"
        ), call. = FALSE)
    }
    if(!is.numeric(zoom_factor) || zoom_factor < 0) {
        stop(paste("In Plot_Coverage,",
            "zoom_factor must be a non-negative number"
        ), call. = FALSE)
    }
    if(!is.numeric(bases_flanking) || bases_flanking < 0) {
        stop(paste("In Plot_Coverage,",
            "bases_flanking must be a non-negative number"
        ), call. = FALSE)
    }
    if(!is.numeric(view_length) || view_length < 0) {
        stop(paste("In Plot_Coverage,",
            "view_length must be a non-negative number"
        ), call. = FALSE)
    }
}

#' Retrieves a list of default recommended filters
#'
#' This function returns the recommended filters. These are:\cr\cr
#' (1) Depth filter of 20,\cr\cr
#' (2) Coverage filter requiring 90% coverage in IR events\cr\cr
#' (3) Coverage filter requiring 60% coverage in AS events
#'   (i.e. Included + Excluded isoforms must cover at least 60% of all junction
#'   events across the given region)\cr\cr
#' (4) Consistency filter requring log difference of 2 (for skipped exon and
#'  mutually exclusive exon events, each junction must comprise at least 1/(2^2)
#'  = 1/4 of all reads associated with each isoform)\cr\cr
#' In all filters, we require at least 80% samples `pcTRUE = 80` from the entire
#'   dataset `minCond = "All"`.
#' Events with read depth (reads supporting either included or excluded isoforms)
#'   lower than 20 `minDepth = 20` are excluded from filters 2,3,4.
#' @return A list of filters to be used in `apply_filters()`. Alternatively,
#'   individual filters can be run using `runFilter()`
#' @examples
#' filters = get_default_filters()
#' @md
#' @export
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


#' Constructs a matrix containing PSI values of the given ASE events
#'
#' This function takes an input SummarizedExperiment `se`, a list of 
#'   alternative splicing events `event_list`, and (optionally) a list of
#'   sample names `sample_list`. 
#' It returns a matrix containing PSI values with columns as samples and rows
#'   as ASE events.
#' 
#' @param se A NxtIRF SummarizedExperiment
#' @param event_list A character vector containing the row names of ASE events
#'   (as given by the `EventName` column of differential ASE results table 
#'   using `limma_ASE()` or `DESeq_ASE()`)
#' @param sample_list (default = `colnames(se)`) A list of sample names
#'   referring to the subset of samples in the given experiment to be included
#'   in the returned matrix
#' @param method The values to be returned (default = "PSI"). It can
#'   alternately be "logit" which returns logit-transformed PSI values, or 
#'   "Z-score" which returns Z-score-transformed PSI values
#' @param depth_threshold (default = 10) If any PSI is derived from raw values
#'   with both isoforms represented by reads below this value, it is given as
#'   `NA` as the uncertainty of PSI would be deemed too highlight
#' @param logit_max The max or min logit values to be capped at, because 
#'   `logit(0) == -Inf` and `logit(1) = Inf`,
#'   this function caps logit values using logit_max. 
#'   (Only used if `method = "logit"`)
#' @param na.percent.max (default = 0.1) The maximum number of NA values in 
#'   an event for the PSI values for each event to be returned. Most heatmap
#'   functions will spring an error if there are too many NA values in any
#'   given row. This option caps the number of NA values to avoid returning
#'   this error.
#' @return A matrix of PSI (or alternate) values, with columns as samples and
#'   rows as ASE events.
#' @md
#' @examples
#' se = NxtIRF_example_NxtSE()
#' 
#' colData(se)$treatment = rep(c("A", "B"), each = 3)
#' 
#' event_list = rowData(se)$EventName
#'
#' mat = make_matrix(se, event_list[1:10])
#' @export
make_matrix <- function(se, event_list, sample_list = colnames(se), 
    method = c("PSI", "logit", "Z-score"), 
    depth_threshold = 10, logit_max = 5, na.percent.max = 0.1) {

    method = match.arg(method)
    inc = assay(se, "Included")[event_list, sample_list, drop = FALSE]
    exc = assay(se, "Excluded")[event_list, sample_list, drop = FALSE]
    mat = inc/(inc + exc)
    mat[inc + exc < depth_threshold] = NA
    mat = mat[rowSums(is.na(mat)) < na.percent.max * ncol(mat),, drop = FALSE]
    if(method == "PSI") {
        # essentially M/Cov
        return(mat)
    } else if(method == "logit") {
        mat = qlogis(mat)
        mat[mat > logit_max] = logit_max
        mat[mat < -logit_max] = -logit_max
        return(mat)
    } else if(method == "Z-score") {
        mat = mat - rowMeans(mat)
        mat = mat / rowSds(mat)
        return(mat)
    }
    
}

#' Constructs a data frame containing average PSI values of the two contrasted 
#'   conditions
#'
#' This function takes an input SummarizedExperiment `se`, a list of 
#'   alternative splicing events `event_list`, the `condition` column 
#'   containing the experimental annotation, and the nominator `nom_DE` and
#'   denominator `denom_DE`.\cr\cr
#' Note that this function takes the geometric mean of PSI, by first converting
#'   all values to logit(PSI), taking the average logit(PSI) values of each
#'   condition, and then converting back to PSI using inv.logit.
#' 
#' @param se A NxtIRF SummarizedExperiment
#' @param event_list A character vector containing the row names of ASE events
#'   (as given by the `EventName` column of differential ASE results table
#'   using `limma_ASE()` or `DESeq_ASE()`)
#' @param condition The name of the column containing the condition values in 
#'   `colData(se)`
#' @param nom_DE The condition to be contrasted, e.g. `nom_DE = "treatment"`
#' @param denom_DE The condition to be contrasted against, e.g. 
#'   `denom_DE = "control"`
#' @param depth_threshold (default = 10) Samples with the number of reads 
#'   supporting either included or excluded isoforms below this values are 
#'   excluded
#' @param logit_max (default = 5). This function works by converting all PSI 
#'   values as logit(PSI), taking the average, then converting back to PSI 
#'   using inverse logit. Because `logit(0) == -Inf` and `logit(1) = Inf`,
#'   this function caps logit values using logit_max
#' @return A a 3 column data frame, with the first column containing 
#'   `event_list` list of ASE events, and the last 2 columns containing the 
#'   average PSI values of the nominator and denominator conditions.
#' @examples
#' se = NxtIRF_example_NxtSE()
#'
#' event_list = rowData(se)$EventName
#'
#' colData(se)$treatment = rep(c("A", "B"), each = 3)
#'
#' diag_values <- make_diagonal(se, event_list,
#'   condition = "treatment", nom_DE = "A", denom_DE = "B"
#' )
#' @md
#' @export
make_diagonal <- function(se, event_list = rownames(se), 
        condition, nom_DE, 
        denom_DE, depth_threshold = 10, logit_max = 5) {

    inc = assay(se, "Included")[event_list, ]
    exc = assay(se, "Excluded")[event_list, ]
    mat = inc/(inc + exc)
    mat[inc + exc < depth_threshold] = NA

    # use logit method to calculate geometric mean

    mat.nom = qlogis(mat[, colData(se)[,condition] == nom_DE])
    mat.denom = qlogis(mat[, colData(se)[,condition] == denom_DE])
    
    mat.nom[mat.nom > logit_max] = logit_max
    mat.denom[mat.denom > logit_max] = logit_max
    mat.nom[mat.nom < -logit_max] = -logit_max
    mat.denom[mat.denom < -logit_max] = -logit_max

    df = data.frame(EventName = event_list, 
        nom = plogis(rowMeans(mat.nom, na.rm = TRUE)),
        denom = plogis(rowMeans(mat.denom, na.rm = TRUE)))
    
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
        if(length(new_cols) > 0) {
            DT_toupdate = DT_toupdate[, c(new_cols) := NULL]
        }
        newDT[DT_toupdate, on = .(sample), 
            (common_cols) := mget(paste0("i.", common_cols))]
    }
    return(as.data.frame(newDT))
}

NxtIRF.SpliceCurve = function(xstart,xend,ystart,yend,y_height,info) {
    source.df = data.frame(
        xstart = xstart, xend = xend,
        ystart = ystart, yend = yend,
        y_height = y_height, info = info,
        stringsAsFactors = FALSE
    )
    final = c()
    for(i in seq_len(nrow(source.df))) {
        temp = with(source.df, 
            data.frame(info = info[i],
                x = seq(xstart[i],xend[i], length.out = 20),
                stringsAsFactors = FALSE)
            )
        temp$y = with(source.df, seq(ystart[i],yend[i], length.out = 20))
        temp$y = with(source.df, temp$y + 
            y_height[i] * sinpi((temp$x - xstart[i]) / (xend[i] - xstart[i]))
            )
        final = rbind(final, temp)
    }
    return(final)
}

#' Draws an arc plot representing the spliced reads of a given sample
#'
#' NB Experimental / WIP. Internal function only
#' @param fst_path The path to the output generated by [CollateData()]
#' @param seqnames,start,end,strand The region to display. `seqnames` and
#'   `strand` should be characters, whereas start and end should be integers
#' @param sample_name The sample name.
#' @return A ggplot object containing the rendered arc plot.
Plot_Junctions <- function(fst_path, 
        seqnames, start, end, strand, sample_name) {
    if(!file.exists(file.path(fst_path, "junc_counts.fst"))) {
        stop(paste("In Plot_Junctions(),",
            "The file", file.path(fst_path, "junc_counts.fst"), "was not found"
        ), call. = FALSE)
    }
    data = fst::read.fst(file.path(fst_path, "junc_counts.fst"))
    if(!(sample_name %in% colnames(data))) {
        stop(paste("In Plot_Junctions(),",
            sample_name, "was not a sample in the given data set"
        ), call. = FALSE)
    }
    
    rownames(data) = data$rownames
    data$rownames = NULL
    
    index.gr = NxtIRF.CoordToGR(rownames(data))
    region.gr = GenomicRanges::GRanges(seqnames = seqnames, 
        ranges = IRanges::IRanges(start = start, end = end),
        strand = strand)
        
    OL = GenomicRanges::findOverlaps(index.gr, region.gr)
    data = data[sort(unique(OL@from)),sample_name, drop = FALSE]
    index = as.data.frame(index.gr[sort(unique(OL@from))])
    data = cbind(data, index)
    colnames(data)[1] = "sample"
    data = data[data$sample > 0,]

    df = NxtIRF.SpliceCurve(data$start, data$end, 0, 0,
        data$sample, rownames(data))
    y_range = max(df$y)
    
    data_mod = data
    data_mod$info = rownames(data)
    data_mod$x = ifelse(data_mod$start < start, start,
        ifelse(data_mod$end > end, end, (data_mod$start + data_mod$end) / 2))
    data_mod$start = ifelse(data_mod$start < start, start, data_mod$start)
    data_mod$end = ifelse(data_mod$end < end, end, data_mod$end)
    
    pl = ggplot() + 
        geom_line(data = df, 
            mapping = aes(
                x = get("x"), y = get("y"), 
                color = get("info"), group = get("info")
            )
        ) + geom_text(data = data_mod, mapping = aes(
            x = get("x"), y = get("sample"), label = get("sample")
        ), nudge_y = 0.05 * y_range) + theme_white
    return(pl)
}