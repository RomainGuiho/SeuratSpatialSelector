#' Interactive selection of cells from a ggplot (Seurat spatial plot)
#'
#' @description
#' Launches an interactive Shiny app that allows users to select spatial points
#' directly on a Seurat spatial plot (converted via Plotly). The user should use
#' the **lasso selection tool** to draw around desired cells, then click
#' **"Validate selection"** to return the selected coordinates.
#'
#' @details
#' Some Seurat spatial objects may display inverted
#' coordinate axes relative to their internal representation.
#' If your selection does not map correctly to the expected cells, try setting
#' `invert_coordinates = TRUE` to swap the `x` and `y` coordinates before returning.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom shiny fluidPage titlePanel fluidRow column
#' @importFrom shiny actionButton verbatimTextOutput renderPrint observeEvent stopApp runApp shinyApp
#' @importFrom plotly ggplotly plotlyOutput renderPlotly event_data config
#' @importFrom utils head
#'
#' @param ggplot_object A ggplot object (e.g. from ImageDimPlot)
#' @param invert_coordinates Logical. If TRUE, swaps x and y coordinates before returning
#'                           (useful if selections appear mirrored or mismatched). Default: FALSE.
#' @param point_size Numeric. Size of points displayed in the interactive plot (default: 3).
#'
#'
#' @return A data.frame of selected points (x, y, cell, etc.)
#' @export
interactive_point_selection <- function(ggplot_object,
                                        invert_coordinates = FALSE,
                                        point_size = 3) {

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
        shiny::tags$div(
          class = "alert alert-info",
          style = "font-size:14px; margin-bottom:15px;",
          shiny::tags$b("Tips for selection:"),
          shiny::tags$ul(
            shiny::tags$li("Use the lasso tool (icon in toolbar) to select cells."),
            shiny::tags$li("Click and drag around the region of interest."),
            shiny::tags$li("Hold down Shift to add areas."),
            shiny::tags$li("Click 'Validate selection' when done."),
            shiny::tags$li(
              shiny::em("If selection looks misaligned, try ",
                        shiny::code("invert_coordinates = TRUE"))
            )
          )
        ),
        shiny::actionButton(
          "return_points",
          "Validate selection",
          class = "btn btn-success",
          style = "width:100%; margin-top:10px; font-size:16px;"
        ),
        shiny::actionButton(
          "cancel",
          "Cancel",
          class = "btn btn-secondary",
          style = "width:100%; margin-top:10px; font-size:16px;"
        )
      )
    )
  )

  server <- function(input, output, session) {

    # Convert ggplot to plotly and display
    output$scatterPlot <- plotly::renderPlotly({
      ggplot2::ggplot_build(ggplot_object)
      ggplot_object$data$text <- ggplot_object$data$cell
      # Overwrite default ggplot2 point size for consistency
      ggplot_object$layers[[1]]$aes_params$size <- point_size

      p <- plotly::ggplotly(ggplot_object, tooltip = "text", source = "scatter")

      # Default to lasso mode
      p <- plotly::config(p, modeBarButtonsToAdd = list("lasso2d"))
      p$x$layout$dragmode <- "lasso"
      p
    })

    # Reactive: selected data
    selectedData <- shiny::reactive({
      plotly::event_data("plotly_selected", source = "scatter")
    })

    # Display preview of selection
    output$selectedPoints <- shiny::renderPrint({
      if (is.null(selectedData())) "No selection yet" else utils::head(selectedData())
    })

    # Handle validation
    shiny::observeEvent(input$return_points, {
      selected <- selectedData()
      if (!is.null(selected) && nrow(selected) > 0) {
        # Optionally invert coordinates for Seurat alignment
        if (invert_coordinates) {
          selected$x <- selected$x
          selected$y <- selected$y
        } else {
          temp_x <- selected$x
          selected$x <- selected$y
          selected$y <- temp_x
        }
        shiny::stopApp(selected)
      } else {
        shiny::showNotification("No cells selected", type = "warning")
      }
    })

    # Cancel button: safely close app without returning anything
    shiny::observeEvent(input$cancel, {
      message("Selection cancelled by user.")
      shiny::stopApp(NULL)
    })
  }

  shiny::runApp(shiny::shinyApp(ui = ui, server = server))
}
