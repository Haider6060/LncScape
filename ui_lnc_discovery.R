lnc_discovery_ui <- tabItem(
  tabName = "lncDiscovery",
  
  fluidPage(
    h3("lncRNA Discovery (GENCODE v44)"),
    
    actionButton("runLncDiscovery", "Run lncRNA Discovery", icon = icon("search")),
    br(), br(),
    
    ## lncRNAs Table & Downloads
    h4("Discovered lncRNAs"),
    dataTableOutput("knownLncTable"),
    br(),
    
    fluidRow(
      column(6, downloadButton("downloadKnownLnc", "Download All lncRNAs", icon = icon("download"))),
      column(6, downloadButton("downloadTopKnownLnc", "Download Top 20 lncRNAs", icon = icon("download")))
    ),
    br(),
    
    ## Heatmap & Download
    actionButton("plotKnownHeatmap", "Generate Top 20 lncRNA Heatmap", icon = icon("chart-area")),
    br(), br(),
    downloadButton("downloadKnownHeatmap", "Download Heatmap (PNG)", icon = icon("download")),
    br(), br(),
    
    plotOutput("knownHeatmapPlot", height = "700px")
  )
)