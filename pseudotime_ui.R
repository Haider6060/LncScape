pseudotime_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Pseudotime Analysis"),
    
    # Upload Top 20 CSV
    fileInput(
      ns("lnc_file"),
      "Upload Top 20 lncRNA CSV:",
      accept = ".csv"
    ),
    
    h4("Uploaded Top 20 lncRNAs"),
    DT::dataTableOutput(ns("lnc_table")),  # Displays uploaded file contents
    
    br(),
    
    # Select lncRNA dynamically
    selectInput(
      ns("selected_gene"),
      "Select lncRNA for Pseudotime Plot:",
      choices = NULL  # Populated dynamically in server
    ),
    
    br(),
    
    # Generate plot and download buttons
    actionButton(
      ns("run_plot"),
      "Generate Pseudotime Plot",
      icon = icon("chart-line")
    ),
    
    br(), br(),
    
    downloadButton(
      ns("download_plot"),
      "Download Plot (PNG)",
      icon = icon("download")
    ),
    
    br(), br(),
    
    # Pseudotime plot display
    plotOutput(ns("pseudo_plot"), height = "700px")
  )
}
