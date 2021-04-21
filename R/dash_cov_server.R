server_cov <- function(id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected) {
    moduleServer(id, function(input, output, session) {
        settings_Cov <- setreactive_Cov()
        get_ref <- function(){
            se = get_se()
            ref(se)
        }
        observeEvent(refresh_tab(), {
            req(refresh_tab())
            output$warning_cov <- renderText({
                validate(need(get_se(), "Please build experiment first"))
            })
            req(get_se())
            .server_cov_refresh(session, get_ref()$gene_list,
                get_de(), rows_all(), rows_selected(),
                input$slider_num_events_cov, input$events_cov,
                input$select_events_cov)
            se = get_se()
            settings_Cov$event.ranges = as.data.table(
                NxtIRF.CoordToGR(rowData(se)$EventRegion))
            settings_Cov$event.ranges$EventName <- rowData(se)$EventName
            .server_cov_refresh_tracks_cov(session, input$mode_cov, 
                input$condition_cov, se)
        })
    # Delayed (debounced) reactives
    chr_r <- reactive({
        req(is_valid(input$chr_cov), 
            input$chr_cov %in% names(get_ref()$seqInfo))
        input$chr_cov
    })
    start_r <- reactive({
        req(suppressWarnings(as.numeric(input$start_cov)))
        as.numeric(input$start_cov)
    })
    end_r <- reactive({
        req(suppressWarnings(as.numeric(input$end_cov)))
        as.numeric(input$end_cov)
    })
    tracks_r <- reactive({
        server_cov_get_all_tracks(input)
    })
    trigger_r <- reactive({
        settings_Cov$trigger
    })
    chr_rd <- chr_r %>% debounce(1000)
    start_rd <- start_r %>% debounce(1000)
    end_rd <- end_r %>% debounce(1000)
    tracks_rd <- tracks_r %>% debounce(3000)
    trigger_rd <- trigger_r %>% debounce(1000)
    
    observeEvent(list(chr_rd(), start_rd(), end_rd()), {
        .server_cov_update_norm_event(input, session, settings_Cov$event.ranges)
    })
    observeEvent(list(input$refresh_coverage, input$stack_tracks_cov,
        input$graph_mode_cov, input$pairwise_t_cov, input$condense_cov), {
        settings_Cov$trigger <- runif(1)
    })
    observeEvent(list(trigger_rd(), tracks_rd()), {
        tracks <- tracks_r()
        settings_Cov$plot_params <- .server_cov_refresh_plot_args(get_se(), get_ref(), 
            input$event_norm_cov, 
            input$chr_cov, suppressWarnings(as.numeric(input$start_cov)), 
            suppressWarnings(as.numeric(input$end_cov)), tracks, 
            settings_Cov$plot_params, input)
        obj <- do.call(plot_cov_fn, settings_Cov$plot_params)
        req(obj)
        settings_Cov$final_plot = obj$final_plot
        settings_Cov$final_plot$x$source = "plotly_ViewRef"
        output$plot_cov <- renderPlotly({
            settings_Cov$plot_ini = TRUE      
            print(
                settings_Cov$final_plot
            )
        })
    })
    observeEvent(input$graph_mode_cov, {
        req(settings_Cov$plot_ini == TRUE)
        .server_cov_plot_change_mode(session, input$graph_mode_cov)
    })
    observeEvent(input$mode_cov, {
        .server_cov_refresh_track_condition(session, input$mode_cov, get_se())
        .server_cov_refresh_tracks_cov(session, input$mode_cov, 
            input$condition_cov, get_se())
    })
    observeEvent(input$condition_cov, {
        .server_cov_refresh_tracks_cov(session, input$mode_cov, 
            input$condition_cov, get_se())
    })
    settings_Cov$plotly_relayout = reactive({
        req(settings_Cov$plot_ini == TRUE)
        event_data("plotly_relayout", source = "plotly_ViewRef")
    })
    observeEvent(settings_Cov$plotly_relayout(), {
        print(settings_Cov$plotly_relayout())
        req(length(settings_Cov$plotly_relayout()) == 2)
        req(all(c("xaxis.range[0]", "xaxis.range[1]") %in% 
            names(settings_Cov$plotly_relayout())))
        updateTextInput(session = session, inputId = "start_cov", 
            value = max(1, 
            round(settings_Cov$plotly_relayout()[["xaxis.range[0]"]])))
        updateTextInput(session = session, inputId = "end_cov", 
            value = round(settings_Cov$plotly_relayout()[["xaxis.range[1]"]]))
    })
    observeEvent(input$zoom_out_cov, {
        req(input$zoom_out_cov, input$chr_cov, 
            input$chr_cov %in% names(get_ref()$seqInfo))
        seqInfo = get_ref()$seqInfo[input$chr_cov]
        .server_cov_zoom_out(input, session, seqInfo)
    })
    observeEvent(input$zoom_in_cov, {
        req(input$zoom_in_cov)
        .server_cov_zoom_in(input, session)
    })
    observeEvent(input$events_cov, {
        req(input$events_cov)
        req(input$events_cov != "(none)")
        events_id_view = settings_Cov$event.ranges[
            get("EventName") == input$events_cov]
        .server_cov_locate_events(input, session, events_id_view)
    })
    observeEvent(input$genes_cov, {
        req(input$genes_cov)
        req(input$genes_cov != "(none)")
        gene_id_view = get_ref()$gene_list[
            get("gene_display_name") == input$genes_cov]
        .server_cov_locate_genes(input, session, gene_id_view)
    })        
    observeEvent(chr_rd(), {
        seqInfo = get_ref()$seqInfo[chr_rd()]
        seqmax = as.numeric(GenomeInfoDb::seqlengths(seqInfo))
        .server_cov_change_chr(input, session, seqmax)
        settings_Cov$trigger <- runif(1)
    })
    observeEvent(list(start_rd(), end_rd()), {
        req(input$chr_cov, input$chr_cov %in% names(get_ref()$seqInfo))
        seqInfo = get_ref()$seqInfo[input$chr_cov]
        seqmax = as.numeric(GenomeInfoDb::seqlengths(seqInfo))
        req(seqmax > 50)
        output <- .server_cov_change_start_end(input, session, output, seqmax)
        settings_Cov$trigger <- runif(1)
    })
    # Populate events
    observeEvent(input$select_events_cov, {
        req(rows_all())
        req(get_de())
        .server_cov_change_event_list(session, input$select_events_cov, 
                input$slider_num_events_cov,
                get_de(), rows_all(), rows_selected())
    })
    observe({
        shinyFileSave(input, "saveplot_cov", roots = volumes(), 
            session = session, filetypes = c("pdf"))    
    })
    observeEvent(input$saveplot_cov, {    
        req(settings_Cov$final_plot)
        selectedfile <- parseSavePath(volumes(), input$saveplot_cov)
        req(selectedfile$datapath)
        plotly::orca(settings_Cov$final_plot, 
            make.path.relative(getwd(), selectedfile$datapath),
            width = 1920, height = 1080)
    })
    
    })
}

.server_cov_get_track_selection <- function(input, i) {
    return(input[[paste0("track", as.character(i), "_cov")]])
}

server_cov_get_all_tracks <- function(input) {
    tracks = list()
    for(i in seq_len(4)) {
        tracks[[i]] = .server_cov_get_track_selection(input, i)       
    }
    return(tracks)
}

.server_cov_refresh <- function(session, gene_list, DE,
        rows_all, rows_selected, num_events, selected_event,
        mode) {
    if(!is.null(gene_list)) {
        message(paste("Populating drop-down box with", 
            length(unique(gene_list$gene_display_name)), "genes"))
        updateSelectInput(session = session, inputId = "chr_cov", 
            choices = c("(none)", 
                as.character(sort(unique(gene_list$seqnames)))),
            selected = "(none)")
        updateSelectizeInput(session = session, inputId = "genes_cov",
            server = TRUE, choices = c("(none)", gene_list$gene_display_name), 
            selected = "(none)")
    } else {
        updateSelectInput(session = session, inputId = "chr_cov", 
            choices = c("(none)"), selected = "(none)")
        updateSelectizeInput(session = session, inputId = "genes_cov", 
            server = TRUE, choices = c("(none)"), selected = "(none)") 
    }
    if(is_valid(DE)) {
        if(mode == "Highlighted") {
            selected = rows_selected
        } else if(mode == "Top N Filtered Results") {
            selected = rows_all
        } else {
            selected = seq_len(nrow(DE))
        }
        if(length(selected) > num_events) {
            selected = selected[seq_len(num_events)]
        }
        if(length(selected) > 0 & is_valid(DE)) {
            if(is_valid(selected_event)) {
                if(!(selected_event %in% DE$EventName[selected])) {
                    selected_event = "(none)"
                }
            } else {
                selected_event = "(none)"
            }
            updateSelectizeInput(session = session, 
                inputId = "events_cov", server = TRUE,
                choices = c("(none)", 
                    DE$EventName[selected]), 
                    selected = selected_event)
        } else {
            updateSelectizeInput(session = session, 
                inputId = "events_cov", server = TRUE,
                choices = c("(none)"), selected = "(none)")
        }
    }
}

.server_cov_get_inrange_events <- function(view_chr, view_start, view_end,
        event.ranges) {
    req(event.ranges)
    DT = event.ranges[get("seqnames") == view_chr &
        get("end") > view_start & get("start") < view_end]
    DT$EventName
}

.server_cov_update_norm_event <- function(input, session, event.ranges) {
    view_chr = isolate(input$chr_cov)
    view_start = suppressWarnings(as.numeric(isolate(input$start_cov)))
    view_end = suppressWarnings(as.numeric(isolate(input$end_cov)))
    selected_event = isolate(input$events_cov)
    cur_event = isolate(input$event_norm_cov)
    req(view_chr)
    req(view_start)
    req(view_end)
    
    event_choices = c("(none)")
    if(is_valid(selected_event)) {
        event_choices = c(event_choices, selected_event)
    } else if(is_valid(cur_event)) {
        event_choices = c(event_choices, cur_event)        
    }
    event_choices = unique(c(event_choices, 
        .server_cov_get_inrange_events(view_chr, view_start, view_end,
            event.ranges)))
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
}

.server_cov_refresh_plot_args <- function(se, ref, norm_event, 
        view_chr, view_start, view_end, tracks, plot_params, input) {
    req(view_chr, view_start, view_end, se)

    rowData = rowData(se)
    events_to_highlight = list()
    if(is_valid(norm_event) && norm_event %in% rowData$EventName) {
        if(rowData$EventType[match(norm_event, rowData$EventName)] 
            %in% c("MXE", "SE")) {
            events_to_highlight[[1]] = c(
                rowData$Event1a[match(norm_event, rowData$EventName)],
                rowData$Event2a[match(norm_event, rowData$EventName)]
            )
        } else {
            events_to_highlight[[1]] = rowData$Event1a[
                match(norm_event, rowData$EventName)
            ]
        }
        if(rowData$EventType[match(norm_event, rowData$EventName)] 
                %in% c("MXE")) {
            events_to_highlight[[2]] = c(
                rowData$Event1b[match(norm_event, rowData$EventName)],
                rowData$Event2b[match(norm_event, rowData$EventName)]
            )
        } else if(rowData$EventType[match(norm_event, rowData$EventName)] 
                %in% c("SE", "A3SS", "A5SS", "ALE", "AFE")){
            events_to_highlight[[2]] = rowData$Event1b[
                match(norm_event, rowData$EventName)]
        }           
    }
    conf.int = 0.95
    req(ref$elem.DT)
    req(ref$transcripts.DT)
    args = list(
        view_chr = view_chr, view_start = view_start, view_end = view_end, 
        view_strand = input$strand_cov, norm_event = norm_event,
        condition = input$condition_cov, tracks = tracks, 
        track_names = "", se = se, 
        avail_files = covfile(se)[file.exists(covfile(se))],
        transcripts = ref$transcripts.DT, 
        elems = ref$elem.DT,
        highlight_events = events_to_highlight,
        stack_tracks = input$stack_tracks_cov,
        graph_mode = input$graph_mode_cov,
        t_test = input$pairwise_t_cov,
        condensed = input$condense_cov
    )
    req(is.null(plot_params) || !identical(plot_params, args))

    return(args)
}

.server_cov_plot_change_mode <- function(session, mode) {
    if(mode == "Pan") {
        plotlyProxy("plot_cov", session) %>% 
            plotlyProxyInvoke("relayout", list(dragmode = "pan")) %>%
            plotlyProxyInvoke("reconfig", editable = FALSE)
    } else if(mode == "Zoom") {
        plotlyProxy("plot_cov", session) %>% 
            plotlyProxyInvoke("relayout", list(dragmode = "zoom")) %>%
            plotlyProxyInvoke("reconfig", editable = FALSE)
    } else if(mode == "Movable Labels") {
        plotlyProxy("plot_cov", session) %>% 
            plotlyProxyInvoke("relayout", list(dragmode = FALSE)) %>%
            plotlyProxyInvoke("reconfig", editable = TRUE)
    }
}

.server_cov_refresh_track_condition <- function(session, mode, se) {
    req(se)
    if(mode == "By Condition") {
        colData = colData(se)
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
}

.server_cov_refresh_tracks_cov <- function(session, mode, condition, se) {
    req(se)
    if(mode == "By Condition") {
        if(is_valid(condition)) {
            colData = colData(se)
            updateSelectInput(session = session, inputId = "track1_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
            updateSelectInput(session = session, inputId = "track2_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
            updateSelectInput(session = session, inputId = "track3_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
            updateSelectInput(session = session, inputId = "track4_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
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
        avail_samples = names(covfile(se)[file.exists(covfile(se))])
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

.server_cov_zoom_out <- function(input, session, seqInfo) {
    view_start = suppressWarnings(as.numeric(input$start_cov))
    view_end = suppressWarnings(as.numeric(input$end_cov))
    req(view_start, view_end, view_end - view_start >= 50)
    seqmax = as.numeric(GenomeInfoDb::seqlengths(seqInfo))
    # get center of current range
    center = round((view_start + view_end) / 2)
    span = view_end - view_start
    # zoom range is 50 * 3^z
    cur_zoom = floor(log(span/50) / log(3))

    new_span = round(span * 3)
    # if(new_span > seqmax - 1) new_span = seqmax - 1
    new_start = max(1, center - round(new_span / 2))
    updateTextInput(session = session, inputId = "start_cov", 
        value = new_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = new_start + new_span)
}

.server_cov_zoom_in <- function(input, session) {
    view_start = suppressWarnings(as.numeric(input$start_cov))
    view_end = suppressWarnings(as.numeric(input$end_cov))
    # get center of current range
    req(view_start, view_end, view_end - view_start >= 50)
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
}

.server_cov_locate_events <- function(input, session, events_id_view) {
    updateSelectInput(session = session, inputId = "chr_cov", 
        selected = events_id_view$seqnames[1])
    span = events_id_view$end[1] - events_id_view$start[1]
    view_start = max(1, events_id_view$start[1] - span)
    view_end = view_start + 3 * span
    updateTextInput(session = session, inputId = "start_cov", 
        value = view_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = view_end)
}

.server_cov_locate_genes <- function(input, session, gene_id_view) {
    updateSelectInput(session = session, inputId = "chr_cov", 
        selected = gene_id_view$seqnames[1])
    updateTextInput(session = session, inputId = "start_cov", 
        value = gene_id_view$start[1])
    updateTextInput(session = session, inputId = "end_cov", 
        value = gene_id_view$end[1])
}

.server_cov_change_chr <- function(input, session, seqmax) {
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
}

.server_cov_change_start_end <- function(input, session, output, seqmax) {
    if(as.numeric(input$end_cov) > seqmax) {
        updateTextInput(session = session, inputId = "end_cov", 
            value = seqmax)
    }
    if(as.numeric(input$end_cov) - as.numeric(input$start_cov) < 50) {
        if(as.numeric(input$end_cov) > 50) {
            updateTextInput(session = session, inputId = "start_cov", 
                value = as.numeric(input$end_cov) - 50)
        } else {
            updateTextInput(session = session, inputId = "start_cov", 
                value = 1)
            updateTextInput(session = session, inputId = "end_cov", 
                value = 51)
        }
    }
    span = as.numeric(input$end_cov) - as.numeric(input$start_cov)
    cur_zoom = floor(log(span/50) / log(3))
    output$label_zoom_cov <- renderText({16 - cur_zoom})
    return(output)
}

.server_cov_change_event_list <- function(session, mode, num_events,
        DE, rows_all, rows_selected) {
    if(mode == "Highlighted") {
        selected = rows_selected
    } else if(mode == "Top N Filtered Results") {
        selected = rows_all
    if(length(selected) > num_events) {
        selected = selected[seq_len(num_events)]
    }
    } else {
        selected = seq_len(min(num_events, nrow(DE)))
    }
    if(length(selected) > 0 & is_valid(DE)) {
        updateSelectizeInput(session = session, inputId = "events_view", 
            server = TRUE, choices = c("(none)", 
                DE$EventName[selected]), selected = "(none)")
        updateSelectizeInput(session = session, inputId = "events_cov", 
            server = TRUE, choices = c("(none)", 
                DE$EventName[selected]), selected = "(none)")
    } else {
        updateSelectizeInput(session = session, inputId = "events_view", 
            server = TRUE, choices = c("(none)"), selected = "(none)")
        updateSelectizeInput(session = session, inputId = "events_cov", 
            server = TRUE, choices = c("(none)"), selected = "(none)")
    }
}