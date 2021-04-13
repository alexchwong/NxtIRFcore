ui_filters <- function(id) {
    ns <- NS(id)
    wellPanel(
        fluidRow(
            column(4,
                wellPanel(style = "overflow-y:scroll; max-height: 800px",
                    filterModule_UI(ns("filter1"), "Filter #1"),
                    filterModule_UI(ns("filter2"), "Filter #2"),
                    filterModule_UI(ns("filter3"), "Filter #3"),
                    filterModule_UI(ns("filter4"), "Filter #4")
                )
            ),
            column(4,
                wellPanel(style = "overflow-y:scroll; max-height: 800px",
                    filterModule_UI(ns("filter5"), "Filter #5"),
                    filterModule_UI(ns("filter6"), "Filter #6"),
                    filterModule_UI(ns("filter7"), "Filter #7"),
                    filterModule_UI(ns("filter8"), "Filter #8")
                )
            ),
            column(4,
                textOutput(ns("current_expr_Filters")), br(),
                textOutput(ns("current_ref_Filters")), br(),
                # actionButton("load_filterdata_Filters", "Load Data"),
                actionButton(ns("refresh_filters_Filters"), "Refresh Filters"),
                plotlyOutput(ns("plot_filtered_Events")),
                selectInput(ns('graphscale_Filters'), 'Y-axis Scale', width = '100%',
                    choices = c("linear", "log10")), 
                shinySaveButton(ns("saveAnalysis_Filters"), "Save Filters", "Save Filters as...", 
                    filetype = list(RDS = "Rds")),
                shinyFilesButton(ns("loadAnalysis_Filters"), label = "Load Filters", 
                    title = "Load Filters from Rds", multiple = FALSE),
                actionButton(ns("loadDefault_Filters"), "Load Default Filters"),
            )
        )
    )
}

server_filters <- function(id, refresh_tab, volumes, 
        get_se, get_filters_from_DE) {
    moduleServer(id, function(input, output, session) {
        settings_filter <- setreactive_filtered_SE()

        observeEvent(refresh_tab(), {
            if(is_valid(get_se())) {
                output$current_expr_Filters = 
                    renderText("SummarizedExperiment loaded")
            } else {
                output$current_expr_Filters = 
                    renderText("Please load SummarizedExperiment first")
            }
            if(is_valid(get_se())) {
                processFilters()
            }
        })
        getFilterData = function(i) {
            if(is_valid(settings_filter$filters)) {
                return(settings_filter$filters[[i]])
            } else {
                return(list())
            }
        }
        r_filter1 <- reactive({getFilterData(1)})
        r_filter2 <- reactive({getFilterData(2)})
        r_filter3 <- reactive({getFilterData(3)})
        r_filter4 <- reactive({getFilterData(4)})
        r_filter5 <- reactive({getFilterData(5)})
        r_filter6 <- reactive({getFilterData(6)})
        r_filter7 <- reactive({getFilterData(7)})
        r_filter8 <- reactive({getFilterData(8)})
        conditionList = reactive({
            req(get_se())
            if(is(get_se(), "NxtSE")) {
                colnames(colData(get_se()))
            } else {
                c("")
            }
        })
        filter1 <- filterModule_server("filter1", r_filter1, conditionList)
        filter2 <- filterModule_server("filter2", r_filter2, conditionList)
        filter3 <- filterModule_server("filter3", r_filter3, conditionList)
        filter4 <- filterModule_server("filter4", r_filter4, conditionList)
        filter5 <- filterModule_server("filter5", r_filter5, conditionList)
        filter6 <- filterModule_server("filter6", r_filter6, conditionList)
        filter7 <- filterModule_server("filter7", r_filter7, conditionList)
        filter8 <- filterModule_server("filter8", r_filter8, conditionList)

        processFilters <- function() {
            message("Refreshing filters")
            if(is(get_se(), "NxtSE")) {
                filterSummary = rep(TRUE, nrow(get_se()))
                if(is_valid(settings_filter$filters)) {
                    for(i in seq_len(8)) {
                        print(settings_filter$filters[[i]]$filterVars)
                        print(settings_filter$filters[[i]]$trigger)
                        if(!is.null(settings_filter$filters[[i]]$trigger)) {
                            filterSummary = filterSummary & runFilter(
                                settings_filter$filters[[i]]$filterClass,
                                settings_filter$filters[[i]]$filterType,
                                settings_filter$filters[[i]]$filterVars,
                                get_se()
                            )
                        } else {
                            message(paste("Trigger", i, "is NULL"))
                        }
                    }
                }
                settings_filter$filterSummary = filterSummary
                message(sum(filterSummary == TRUE))
            }
        }
            
        observeEvent({list(
            input$refresh_filters_Filters
        )}, {
            req(input$refresh_filters_Filters)
            settings_filter$filters[[1]] = (reactiveValuesToList(filter1))
            settings_filter$filters[[2]] = (reactiveValuesToList(filter2))
            settings_filter$filters[[3]] = (reactiveValuesToList(filter3))
            settings_filter$filters[[4]] = (reactiveValuesToList(filter4))
            settings_filter$filters[[5]] = (reactiveValuesToList(filter5))
            settings_filter$filters[[6]] = (reactiveValuesToList(filter6))
            settings_filter$filters[[7]] = (reactiveValuesToList(filter7))
            settings_filter$filters[[8]] = (reactiveValuesToList(filter8))
            processFilters()
        })
        
        observeEvent({list(
            settings_filter$filterSummary,
            input$graphscale_Filters
        )}, {
            req(get_se())
            req(settings_filter$filterSummary)
            DT = data.table(
                EventType = rowData(get_se())$EventType,
                keep = settings_filter$filterSummary
            )
            output$plot_filtered_Events <- renderPlotly({
                print(
                    Filters_Plot_Summary(DT, input$graphscale_Filters)
                )
            })
        })
        observe({
            shinyFileSave(input, "saveAnalysis_Filters", 
                roots = volumes(), session = session)        
        })
        observeEvent(input$saveAnalysis_Filters, {
            req(settings_filter$filters)
            selectedfile <- parseSavePath(volumes(), input$saveAnalysis_Filters)
            req(selectedfile$datapath)

            settings_filter$filters[[1]] = (reactiveValuesToList(filter1))
            settings_filter$filters[[2]] = (reactiveValuesToList(filter2))
            settings_filter$filters[[3]] = (reactiveValuesToList(filter3))
            settings_filter$filters[[4]] = (reactiveValuesToList(filter4))
            settings_filter$filters[[5]] = (reactiveValuesToList(filter5))
            settings_filter$filters[[6]] = (reactiveValuesToList(filter6))
            settings_filter$filters[[7]] = (reactiveValuesToList(filter7))
            settings_filter$filters[[8]] = (reactiveValuesToList(filter8))

            final = settings_filter$filters
            saveRDS(final, selectedfile$datapath)
        })

        observe({
            shinyFileChoose(input, "loadAnalysis_Filters", 
                roots = volumes(), session = session,
                filetypes = c("Rds"))        
        })
        observeEvent(input$loadAnalysis_Filters, {
            selectedfile <- parseFilePaths(volumes(), 
                input$loadAnalysis_Filters)
            req(selectedfile$datapath)
            settings_filter$filters = readRDS(selectedfile$datapath)
        })
        # Import filters from loading DE object
        observeEvent(get_filters_from_DE(), {
            req(get_filters_from_DE())
            settings_filter$filters = get_filters_from_DE()
        })
        observeEvent(input$loadDefault_Filters, {
            settings_filter$filters = get_default_filters()
        })
        return(settings_filter)
    })
}

Filters_Plot_Summary <- function(DT, scale) {
    if(scale == "log10") {
        DT[, c("Included", "Excluded") := .(
            log10(sum(get("keep") == TRUE)),
            log10(sum(!is.na(get("keep")))) - 
                log10(sum(get("keep") == TRUE))
        ), by = "EventType"]
    } else {
        DT[, c("Included", "Excluded") := .(
            sum(get("keep") == TRUE),
            sum(get("keep") != TRUE)
        ), by = "EventType"] 
    }
    DT = unique(DT, by = "EventType")
    incl <- as.data.frame(DT[, c("EventType", "Included")])
    incl$filtered = "Included"
    colnames(incl)[colnames(incl) == "Included"] <- "Events"
    excl <- as.data.frame(DT[, c("EventType", "Excluded")])
    excl$filtered = "Excluded"
    colnames(excl)[colnames(excl) == "Excluded"] <- "Events"    

    # ggplot summary as bar plot
    p = ggplot(rbind(incl, excl), 
        aes(x = get("EventType"), y = get("Events"), fill = get("filtered"))
    ) + geom_bar(position="stack", stat="identity")
    if(scale == "log10") {
        p = p + labs(y = "log10 Events")
    } else {
        p = p + labs(y = "Events")        
    }
    p = p + labs(fill = "Filtered")
    return(p)
}