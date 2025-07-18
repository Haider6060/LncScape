tld_ui <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = "tld_analysis",
    
    fluidPage(
      titlePanel("TF–lncRNA Dynamics (TLD) Analysis"),
      
      fluidRow(
        column(12,
               fileInput(ns("lds_file"), "Upload lncRNA LDS CSV File:", accept = c(".csv"))
        )
      ),
      
      fluidRow(
        column(6, actionButton(ns("run_tld"), "Run TLD Analysis")),
        column(6, actionButton(ns("generate_plot"), "Generate Plot"))
      ),
      
      br(),
      
      fluidRow(
        column(6, downloadButton(ns("download_tld_results"), "Download CSV")),
        column(6, downloadButton(ns("download_tld_barplot"), "Download Plot (PNG)"))
      ),
      
      br(), br(),
      
      fluidRow(
        column(12,
               h4("Barplot of Top 20 TFs for Top Dynamic lncRNA"),
               plotOutput(ns("tld_barplot"), height = "600px")
        )
      )
    )
  )
}
