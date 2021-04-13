#' Launches the NxtIRF Graphics User Interface using Shiny Dashboard
#' 
#' This function launches the NxtIRF interactive app using Shiny Dashboard
#' @return None
#' @examples
#' # nxtIRF() # Launches interactive ShinyDashboard NxtIRF app
#' @export
nxtIRF <- function() {
    if(!interactive()) {
        stop(paste("In nxtIRF(),",
            "NxtIRF App can only be run in interactive mode (i.e. RStudio)."
        ), call. = FALSE)
    }

    ui_dash <- dashboardPage(
        dashboardHeader(title = "NxtIRF"),
        ui_sidebar(),
        dashboardBody(
            tabItems(
                ui_tab_title(),
                
                ui_tab_system(),
                
                ui_tab_ref_new(),        

                ui_tab_expr(),
                
                ui_tab_expr_load(),
                ui_tab_qc(),
                ui_tab_filter(),
                ui_tab_analyse(),

                ui_tab_diag(),
                ui_tab_volcano(),
                ui_tab_heatmap(),
                ui_tab_coverage()
            )
        )
    )
    runApp(shinyApp(ui_dash, dash_server))
}