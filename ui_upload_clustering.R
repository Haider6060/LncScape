lustering_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Clustering & UMAP"),
    sliderInput(ns("resolutionInput"), "Clustering Resolution:",
                min = 0.2, max = 1.2, value = 0.5, step = 0.05),
    actionButton(ns("runClustering"), "Run Clustering and UMAP", icon = icon("play")),
    br(), br(),
    downloadButton(ns("downloadPlot"), "Download UMAP (PNG)", icon = icon("download")),
    br(), br(),
    plotOutput(ns("umapPlot"), height = "700px")
  )
}