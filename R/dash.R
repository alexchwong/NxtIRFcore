#' Launches the NxtIRF Graphics User Interface using Shiny Dashboard
#' 
#' This function launches the NxtIRF interactive app using Shiny Dashboard
#' @param mode (default `"dialog"`) `"dialog"` displays NxtIRF in a dialog box 
#'   with specified width and height. `"browser"` opens NxtIRF in a browser-
#'   like resizable window.
#' @param width,height If `mode` is set to `"dialog"`, the specified width
#'   and height of the NxtIRF app.
#' @return An interactive shinydashboard NxtIRF app runs.
#' @examples
#' \donttest{
#' # Launches interactive ShinyDashboard NxtIRF app as fixed-size dialog box
#' nxtIRF(mode = "dialog", width = 1600, height = 900) 
#'
#' # Launches interactive ShinyDashboard NxtIRF app as browser window
#' nxtIRF(mode = "browser") 

#' }
#' @md
#' @export
nxtIRF <- function(mode = c("dialog", "browser"), 
        width = 1600, height = 900) {
    if(!interactive()) {
        stop(paste("In nxtIRF(),",
            "NxtIRF App can only be run in interactive mode (i.e. RStudio)."
        ), call. = FALSE)
    }
    mode = match.arg(mode)
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
    if(mode == "dialog") {
        runGadget(
            shinyApp(ui_dash, dash_server),
            viewer = dialogViewer('NxtIRF', width = width, height = height)
        )
    } else {
        runApp(shinyApp(ui_dash, dash_server))
    }
}