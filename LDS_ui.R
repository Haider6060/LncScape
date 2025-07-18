lds_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),  # Initialize shinyjs
    
    h2("LDS: lncRNA Dynamics Score"),
    fileInput(ns("lds_csv"), "Upload All lncRNAs CSV:", accept = ".csv"),
    h4("Uploaded lncRNA Data"),
    DT::DTOutput(ns("lnc_table")),
    
    br(),
    fluidRow(
      column(6, actionButton(ns("runLDS"), "▶ Run LDS Analysis")),
      column(6, actionButton(ns("plotLDS"), "📈 Generate LDS Plot"))
    ),
    br(),
    fluidRow(
      column(4, downloadButton(ns("downloadAllLDS"), "⬇ Download All LDS (CSV)")),
      column(4, downloadButton(ns("downloadTop50LDS"), "⬇ Download Top 50 LDS (CSV)")),
      column(4, downloadButton(ns("downloadPlotLDS"), "⬇ Download LDS Plot (PNG)"))
    ),
    br(),
    h4("LDS Expression Plot (Top 5 Genes)"),
    plotOutput(ns("lds_plot"))
  )
}
