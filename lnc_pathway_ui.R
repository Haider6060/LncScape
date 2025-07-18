lnc_pathway_ui <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = "lncPathway",
    fluidPage(
      h3("lncRNA-Centric Pathway Enrichment and Visualization"),
      
      tags$b("Upload lncRNA CSV (must have a 'gene_name' column):"),
      fileInput(
        ns("lnc_pathway_file"),
        label = NULL,
        accept = ".csv",
        width = "50%"
      ),
      br(),
      
      uiOutput(ns("dataset_selector_ui")),
      br(),
      
      h4("Uploaded lncRNAs Preview"),
      DT::dataTableOutput(ns("lnc_pathway_table")),
      br(),
      
      fluidRow(
        column(
          width = 6,
          actionButton(
            ns("run_lnc_pathway"),
            "Run Pathway Enrichment",
            icon = icon("play"),
            style = "width: 100%;"
          )
        ),
        column(
          width = 6,
          actionButton(
            ns("plot_lnc_pathway"),
            "Generate Plot",
            icon = icon("chart-bar"),
            style = "width: 100%;"
          )
        )
      ),
      br(),
      
      fluidRow(
        column(
          width = 6,
          downloadButton(
            ns("download_lnc_pathway"),
            "Download Enrichment Results (CSV)",
            style = "width: 100%;"
          )
        ),
        column(
          width = 6,
          downloadButton(
            ns("download_lnc_pathway_plot"),
            "Download Plot (PNG)",
            style = "width: 100%;"
          )
        )
      ),
      br(),
      
      plotOutput(ns("lnc_pathway_plot"), height = "700px")
    )
  )
}
