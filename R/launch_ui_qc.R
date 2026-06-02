launch_ui_qc <- function(res)
{
    ## =========================
    ## 3D QC Explorer
    ## =========================
    try({
        print("Preparing data frame")
        qc_df <- res$QC_metrics
        qc_df$bam_id <- rownames(qc_df)
        numeric_cols <- names(qc_df)[sapply(qc_df, is.numeric)]
        qc_df <- qc_df[, colSums(!is.na(qc_df)) > 0]
    })
    ## =========================
    library(shiny)
    library(plotly)
    library(dplyr)
    library(DT)
    ## =========================
    if (!exists("qc_df"))
    {
        set.seed(1)
        qc_df <- data.frame(
            bam_id = paste0("bam_", 1:200),
            metric1 = rnorm(200, 50, 10),
            metric2 = rnorm(200, 100, 20),
            metric3 = rnorm(200, 0.5, 0.1),
            metric4 = rnorm(200, 1000, 200)
        )
    }
    ## =========================
    qc_df <- qc_df %>%
        dplyr::select(where(~ !all(is.na(.))))
    if (!"bam_id" %in% colnames(qc_df))
    {
        qc_df$bam_id <- rownames(qc_df)
    }
    numeric_cols <- names(qc_df)[sapply(qc_df, is.numeric)]
    ## ---- UI ----
    ui <- fluidPage(
        titlePanel("3D BAM QC Explorer (Multidimensional Gating)"),
        sidebarLayout(
            sidebarPanel(
                h4("3D View"),
                selectInput("xvar", "X-axis", numeric_cols, selected = numeric_cols[1]),
                selectInput("yvar", "Y-axis", numeric_cols, selected = numeric_cols[2]),
                selectInput("zvar", "Z-axis", numeric_cols, selected = numeric_cols[3]),
                hr(),
                h4("Adjust bounds (current axes)"),
                uiOutput("x_ui"),
                uiOutput("y_ui"),
                uiOutput("z_ui"),
                hr(),
                downloadButton("download_filtered", "Download filtered BAMs"),
                downloadButton("download_bounds", "Download all bounds")
            ),
            mainPanel(
                plotlyOutput("plot3d", height = "600px"),
                hr(),
                DTOutput("table")
            )
        )
    )
    ## =========================
    ## ---- SERVER ----
    server <- function(input, output, session)
    {
        ## ---- GLOBAL BOUNDS STORE ----
        bounds <- reactiveValues()
        ## Initialize bounds once
        observe({
            for (col in numeric_cols) {
                if (is.null(bounds[[col]])) {
                    rng <- range(qc_df[[col]], na.rm = TRUE)
                    bounds[[col]] <- rng
                }
            }
        })
        ## ---- SLIDER UI ----
        make_slider <- function(varname, id) {
            req(bounds[[varname]])
            rng <- range(qc_df[[varname]], na.rm = TRUE)
            sliderInput(
                id,
                varname,
                min = rng[1],
                max = rng[2],
                value = bounds[[varname]]
            )
        }
        output$x_ui <- renderUI({ make_slider(input$xvar, "xrange") })
        output$y_ui <- renderUI({ make_slider(input$yvar, "yrange") })
        output$z_ui <- renderUI({ make_slider(input$zvar, "zrange") })
        ## ---- UPDATE GLOBAL BOUNDS ----
        observeEvent(input$xrange, {
            bounds[[input$xvar]] <- input$xrange
        })
        observeEvent(input$yrange, {
            bounds[[input$yvar]] <- input$yrange
        })
        observeEvent(input$zrange, {
            bounds[[input$zvar]] <- input$zrange
        })
        ## ---- CLASSIFICATION (ALL DIMENSIONS) ----
        classified_data <- reactive({
            df <- qc_df
            df$pass <- TRUE
            for (col in numeric_cols) {
                b <- bounds[[col]]
                if (!is.null(b)) {
                    df$pass <- df$pass &
                        df[[col]] >= b[1] &
                        df[[col]] <= b[2]
                }
            }
            df
        })
        ## ---- 3D PLOT ----
        output$plot3d <- renderPlotly({
            df <- classified_data()
            plot_ly(
                df,
                x = ~.data[[input$xvar]],
                y = ~.data[[input$yvar]],
                z = ~.data[[input$zvar]],
                color = ~pass,
                colors = c("red", "blue"),
                text = ~bam_id,
                type = "scatter3d",
                mode = "markers"
            ) %>%
                layout(
                    scene = list(
                        xaxis = list(title = input$xvar),
                        yaxis = list(title = input$yvar),
                        zaxis = list(title = input$zvar)
                    )
                )
        })
        ## ---- TABLE ----
        output$table <- renderDT({
            classified_data()
        })
        ## ---- DOWNLOAD FILTERED ----
        output$download_filtered <- downloadHandler(
            filename = function() "filtered_bams.csv",
            content = function(file) {
                write.csv(
                    classified_data() %>% filter(pass),
                    file,
                    row.names = FALSE
                )
            }
        )
        ## ---- DOWNLOAD BOUNDS ----
        output$download_bounds <- downloadHandler(
            filename = function() "bounds.csv",
            content = function(file) {
                bounds_df <- data.frame(
                    metric = numeric_cols,
                    min = sapply(numeric_cols, function(x) bounds[[x]][1]),
                    max = sapply(numeric_cols, function(x) bounds[[x]][2])
                )
                write.csv(bounds_df, file, row.names = FALSE)
            }
        )
    }
    ## =========================
    shinyApp(ui, server)
}
