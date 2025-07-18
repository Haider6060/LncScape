gsva_pathway_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("GSVA Pathway Enrichment"),
    
    # Dataset Selector
    selectInput(
      ns("gsva_dataset"),
      "Select GSVA Dataset:",
      choices = c(
        "Cancer Hallmark" = "GSVA_Cancer_Hallmark_ssGSEA_scores.rds",
        "Disease Type" = "GSVA_Disease_Type_scores.rds",
        "Drug" = "gsva_scores_drug.rds",
        "Survival" = "gsva_scores_survival.rds",
        "lncRNADisease" = "gsva_lncrnadisease_scores.rds",
        "NPInter" = "gsva_scores_npinter.rds"
      ),
      selected = "GSVA_Cancer_Hallmark_ssGSEA_scores.rds"
    ),
    
    # File upload
    fileInput(
      ns("gsva_file"),
      "Upload lncRNA CSV File (All_lncRNAs.csv):",
      accept = c(".csv")
    ),
    
    # Action buttons in 2x2 layout
    fluidRow(
      column(6, actionButton(ns("run_gsva"), "Run GSVA Enrichment", icon = icon("play"), width = "100%")),
      column(6, actionButton(ns("plot_gsva"), "Generate GSVA Plot", icon = icon("chart-area"), width = "100%"))
    ),
    br(),
    fluidRow(
      column(6, downloadButton(ns("download_gsva_csv"), "Download Enrichment CSV", icon = icon("download"), width = "100%")),
      column(6, downloadButton(ns("download_gsva_plot"), "Download GSVA Plot (PNG)", icon = icon("download"), width = "100%"))
    ),
    br(),
    
    # Outputs
    h4("GSVA Enrichment Table"),
    DT::dataTableOutput(ns("gsva_table")),
    br(),
    h4("GSVA Pathway Activity Plot"),
    plotOutput(ns("gsva_plot"), height = "700px")
  )
}
