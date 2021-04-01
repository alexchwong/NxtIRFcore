ui_cov <- function(id) {
    ns <- NS(id)
    wellPanel(
        fluidRow(style='height:20vh',
            column(6, 
                textOutput(ns("warning_cov")),
                div(style=paste0(
                    "display: inline-block;",
                    "vertical-align:top;",
                    "width: 250px;"), 
                    selectizeInput(ns('genes_cov'), 'Genes', 
                        choices = "(none)")
                ),
                div(style=paste0("display: inline-block;",
                    "vertical-align:top;",
                    "width: 350px;"), 
                    selectizeInput(ns('events_cov'), 'Events', 
                        choices = c("(none)"))
                ),
                shinyWidgets::radioGroupButtons(ns("select_events_cov"), 
                    label = 
                        "Select Events from Differential Expression Results",
                        justified = FALSE,
                        choices = c(
                            "Top N All Results", 
                            "Top N Filtered Results", 
                            "Highlighted"
                        ),
                    checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                )
            ),
            column(4,
                div(style=paste0("display: inline-block;",
                    "vertical-align:top;",
                    "width: 80px;"),
                    selectInput(ns("chr_cov"), label = "Chr", 
                        choices = c("(none)"), selected = "(none)")),
                div(style=paste0("display: inline-block;",
                    "vertical-align:top; width: 120px;"),
                    textInput(ns("start_cov"), label = "Left", 
                        value = c(""))),
                div(style=paste0("display: inline-block;",
                    "vertical-align:top; width: 120px;"),
                    textInput(ns("end_cov"), label = "Right", 
                        value = c(""))),           
                br(),
                shinyWidgets::actionBttn(ns("zoom_out_cov"), 
                    style = "material-circle", 
                    color = "danger",icon = icon("minus")),
                div(style=paste0("display: inline-block;",
                    "vertical-align:center;width: 50px;padding:25px"),
                    textOutput(ns("label_zoom_cov"))),
                shinyWidgets::actionBttn(ns("zoom_in_cov"), 
                    style = "material-circle", 
                    color = "danger", icon = icon("plus")),
                div(style=paste0("display: inline-block;",
                    "vertical-align:center;padding:15px"),
                    shinyWidgets::radioGroupButtons(ns("strand_cov"), 
                        label = "Strand", justified = FALSE,
                        choices = c("*", "+", "-"), 
                        checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                    )
                )
            ),
            column(2,
                div(style=paste0("display: inline-block;",
                    "vertical-align:top;padding:10px;"),
                    shinyWidgets::radioGroupButtons(ns("graph_mode_cov"), 
                        label = "Graph Mode", justified = FALSE,
                        choices = c("Pan", "Zoom", "Movable Labels"), 
                        checkIcon = list(yes = icon("ok", lib = "glyphicon")))
                ),
                shinyWidgets::sliderTextInput(ns("slider_num_events_cov"), 
                    "Num Events", choices = c(5, 10, 25, 50, 100, 200, 500),
                    selected = 25) 
            )
        ),
        fluidRow(
            column(2, 
                actionButton(ns("refresh_coverage"), "Refresh Plot"),
                selectInput(ns('mode_cov'), 'View', width = '100%',
                    choices = c("Individual", "By Condition")),					
                selectInput(ns('event_norm_cov'), 'Normalize Event', 
                    width = '100%', choices = c("(none)")),
                conditionalPanel(
                    ns = ns,
                    condition = "['By Condition'].indexOf(input.mode_cov) >= 0",
                    selectInput(ns('condition_cov'), 'Condition', 
                        width = '100%', choices = c("(none)"))					
                ),
                selectInput(ns('track1_cov'), 'Track 1', width = '100%',
                    choices = c("(none)")),
                selectInput(ns('track2_cov'), 'Track 2', width = '100%',
                    choices = c("(none)")),
                selectInput(ns('track3_cov'), 'Track 3', width = '100%',
                    choices = c("(none)")),
                selectInput(ns('track4_cov'), 'Track 4', width = '100%',
                    choices = c("(none)")),
                shinyWidgets::switchInput(ns("stack_tracks_cov"), 
                    label = "Stack Traces", labelWidth = "150px"),
                shinyWidgets::switchInput(ns("pairwise_t_cov"), 
                    label = "Pairwise t-test", labelWidth = "150px"),
                shinyWidgets::switchInput(ns("condense_cov"), 
                    label = "Condensed Tracks", labelWidth = "150px"),
                shinySaveButton(ns("saveplot_cov"), "Save Plot as PDF", "Save Plot as PDF...", 
                    filetype = list(PDF = "pdf")),
            ),
            column(10, 
                plotlyOutput(ns("plot_cov"), height = "800px"),
            )
        )    
    )
}
