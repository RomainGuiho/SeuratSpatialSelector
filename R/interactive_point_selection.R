#' Interactive selection of cells from a ggplot (Seurat spatial plot)
#'
#' @importFrom ggplot2 ggplot
#'
#' @param ggplot_object A ggplot object (e.g. from ImageDimPlot)
#' @return A data.frame of selected points (x, y, cell, etc.)
#' @export
interactive_point_selection <- function(ggplot_object) {
  ui <- shiny::fluidPage(
    shiny::titlePanel("Interactive Cell Selection"),
    shiny::fluidRow(
      shiny::column(
        width = 9,
        plotly::plotlyOutput("scatterPlot", height = "800px"),
        shiny::verbatimTextOutput("selectedPoints")
      ),
      shiny::column(
        width = 3,
        shiny::actionButton("return_points", "Validate selection", class = "btn-primary")
      )
    )
  )

  server <- function(input, output, session) {
    output$scatterPlot <- plotly::renderPlotly({
      ggplot2::ggplot_build(ggplot_object)
      ggplot_object$data$text <- ggplot_object$data$cell
      plotly::ggplotly(ggplot_object, tooltip = "text", source = "scatter")
    })

    selectedData <- shiny::reactive({
      plotly::event_data("plotly_selected", source = "scatter")
    })

    output$selectedPoints <- shiny::renderPrint({
      if (is.null(selectedData())) "No selection yet" else utils::head(selectedData())
    })

    shiny::observeEvent(input$return_points, {
      selected <- selectedData()
      if (!is.null(selected)) {
        temp_x <- selected$x
        selected$x <- selected$y
        selected$y <- temp_x
        shiny::stopApp(selected)
      }
    })
  }

  shiny::runApp(shiny::shinyApp(ui = ui, server = server))
}
