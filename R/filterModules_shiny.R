filterModule_UI <- function(id, label = "Counter") {
    ns <- NS(id)
    wellPanel(
        h5(label),  # e.g. "Filter #1"
        selectInput(ns("filterClass"), "Filter Class", 
            width = '100%', choices = c("(none)", "Annotation", "Data")),
        selectInput(ns("filterType"), "Filter Type", 
            width = '100%', choices = c("(none)")),
        conditionalPanel(ns = ns,
        condition = paste0("['Transcript_Support_Level'].",
            "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_TSL_min"), 
                "TSL Threshold", 
                choices = seq_len(5), selected = 1)
        ),
        conditionalPanel(ns = ns,
            condition = paste0("['Consistency'].",
                "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_cons_max"), 
                "log-fold maximum", choices = seq(0.2, 5, by = 0.2), 
                selected = 1)
        ),
        conditionalPanel(ns = ns,
            condition = "['Coverage'].indexOf(input.filterType) >= 0",
            sliderInput(ns("slider_cov_min"), "Percent Coverage", 
                min = 0, max = 100, value = 80)
        ),
        conditionalPanel(ns = ns,
            condition = "['Depth'].indexOf(input.filterType) >= 0",
            shinyWidgets::sliderTextInput(ns("slider_depth_min"), 
                "Minimum", choices = c(1,2,3,5,10,20,30,50,100,200,300,500), 
                selected = 20),
        ),
        conditionalPanel(ns = ns,
            condition = "['Depth', 'Coverage'].indexOf(input.filterType) >= 0",
            tagList(
                shinyWidgets::sliderTextInput(ns("slider_mincond"), 
                    "Minimum Conditions Satisfy Criteria", 
                    choices = c(as.character(seq_len(8)), "All"), 
                    selected = "All"),
                selectInput(ns("select_conds"), "Condition", width = '100%',
                    choices = c("(none)")),
                sliderInput(ns("slider_pcTRUE"), 
                    "Percent samples per condition satisfying criteria", 
                    min = 0, max = 100, value = 80)
            )
        ),
        conditionalPanel(ns = ns,
            condition = paste0("['Coverage', 'Consistency'].",
                "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_minDepth"), 
                "Signal Threshold to apply criteria", 
                choices = c(1,2,3,5,10,20,30,50,100,200,300,500), 
                selected = 20),
        ),
        conditionalPanel(ns = ns,
            condition = "['(none)'].indexOf(input.filterClass) < 0",
            selectInput(ns("EventType"), "Splice Type", width = '100%', 
                multiple = TRUE,
                choices = c("IR", "MXE", "SE", "AFE", "ALE", "A5SS", "A3SS"))
        )
    )
}

filterModule_server <- function(id, filterdata, conditionList) {
    moduleServer(id, function(input, output, session) {
        final <- reactiveValues(
            filterClass = "",
            filterType = "",
            filterVars = list()
        )

        observeEvent(conditionList(), {
            choices_conds = c("(none)", conditionList())
            if(is_valid(final$filterVars$condition) && 
                    final$filterVars$condition %in% choices_conds) {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = final$filterVars$condition)
            } else {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)")            
            }
        })

        # inputs from final -> UI
        observeEvent(filterdata(), {
            final = filterdata()

            if(is_valid(final$filterClass)) {
                if(final$filterClass == "Annotation") {
                    type_choices = c("Protein_Coding", 
                        "NMD_Switching", "Transcript_Support_Level")              
                } else if(final$filterClass == "Data") {
                    type_choices = c("Depth", "Coverage", "Consistency")
                } else {
                    type_choices = c("(none)")
                }
            } else {
                type_choices = c("(none)")
            }
            if(is_valid(final$filterType) && 
                    final$filterType %in% type_choices) {
                updateSelectInput(session = session, inputId = "filterType", 
                    choices = type_choices, selected = final$filterType)
                updateSelectInput(session = session, inputId = "filterClass", 
                    choices = c("(none)", "Annotation", "Data"), 
                    selected = final$filterClass)
            } else {
                # Invalid filter; destroy this record
                final$filterClass = "(none)"
                final$filterType = "(none)"
                updateSelectInput(session = session, inputId = "filterClass", 
                    choices = c("(none)", "Annotation", "Data"))
                updateSelectInput(session = session, inputId = "filterType", 
                    choices = c("(none)"))

                return()
            }
         
            if(is_valid((final$filterVars$minimum))) {
                if(final$filterType == "Depth") {
                    shinyWidgets::updateSliderTextInput(
                        session = session, inputId = "slider_depth_min", 
                        selected = final$filterVars$minimum)
                } else  if(final$filterType == "Coverage"){
                    updateSliderInput(session = session, 
                        inputId = "slider_cov_min", 
                        value = final$filterVars$minimum)
                } else  if(final$filterType == "Transcript_Support_Level"){
                    shinyWidgets::updateSliderTextInput(
                        session = session, inputId = "slider_TSL_min", 
                        selected = final$filterVars$minimum)
                }
            }
            if(is_valid(final$filterVars$maximum)) {
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_cons_max", 
                    selected = final$filterVars$maximum)
            }
            if(is_valid(final$filterVars$minDepth)) {
                updateSelectInput(session = session, 
                    inputId = "slider_minDepth", 
                    selected = final$filterVars$minDepth)
            }
            if(is_valid(final$filterVars$minCond)) {
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_mincond", 
                    selected = final$filterVars$minCond)
            }
            choices_conds = c("(none)", conditionList())
            if(is_valid(final$filterVars$condition) && 
                    final$filterVars$condition %in% choices_conds) {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = final$filterVars$condition)
            } else {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)")            
            }
            if(is_valid(final$filterVars$pcTRUE)){
                updateSliderInput(session = session, 
                    inputId = "slider_pcTRUE", 
                    value = final$filterVars$pcTRUE)
            }
            if(is_valid(final$filterVars$EventTypes)) {
                updateSelectInput(session = session, 
                    inputId = "EventType", 
                    selected = final$filterVars$EventTypes)
            } else {
                updateSelectInput(session = session, 
                    inputId = "EventType", selected = NULL)          
            }
        })

        # outputs from UI -> final
        observeEvent(input$filterClass, {
            final$filterClass = input$filterClass
            if(final$filterClass == "Annotation") {
                type_choices = c("Protein_Coding", 
                    "NMD_Switching", "Transcript_Support_Level")
            } else if(final$filterClass == "Data") {
                type_choices = c("Depth", "Coverage", "Consistency")
            } else {
                type_choices = "(none)"
            }
            cur_choice = isolate(final$filterType)
            if(is_valid(cur_choice) && cur_choice %in% type_choices) {
                updateSelectInput(session = session, 
                    inputId = "filterType", 
                    choices = type_choices, selected = cur_choice)                    
            } else {
                updateSelectInput(session = session, 
                    inputId = "filterType", 
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
            } else if(input$filterType == "Coverage"){
                final$filterVars$minimum = input$slider_cov_min
            } else if(input$filterType == "Transcript_Support_Level"){
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
    })
}