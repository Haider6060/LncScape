lncscape_ui <- dashboardPage(
  
  dashboardHeader(title = "lncScape"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Input & Detection", tabName = "dataInput", icon = icon("upload")),
      menuItem("Clustering & UMAP", tabName = "clustering", icon = icon("project-diagram")),
      menuItem("lncRNA Discovery", tabName = "lncDiscovery", icon = icon("search")),
      menuItem("Pseudotime Analysis", tabName = "pseudotime", icon = icon("chart-line")),
      menuItem("Pathway Enrichment", tabName = "lncPathway", icon = icon("network-wired")),
      menuItem("Pathway Activity (GSVA)", tabName = "gsvaPathway", icon = icon("dna")),
      menuItem("LDS Method", tabName = "ldsMethod", icon = icon("chart-bar")),
      menuItem("TF–lncRNA Dynamics (TLD)", tabName = "tld_analysis", icon = icon("exchange-alt"))
    )
  ),
  
  dashboardBody(
    tabItems(
      
      tabItem(
        tabName = "dataInput",
        h3("Data Input & Detection"),
        fileInput("fileUpload", "Upload Seurat Object (.rds):", accept = ".rds"),
        h4("Seurat Object Information"),
        verbatimTextOutput("fileInfo"),
        h4("Dataset Type"),
        textOutput("datasetType"),
        h4("Metadata Preview"),
        verbatimTextOutput("metaDataSummary")
      ),
      
      tabItem(
        tabName = "clustering",
        h3("Clustering & UMAP"),
        sliderInput("resolutionInput", "Clustering Resolution:",
                    min = 0.2, max = 1.2, value = 0.5, step = 0.05),
        actionButton("runClustering", "Run Clustering and UMAP", icon = icon("play")),
        br(), br(),
        downloadButton("downloadPlot", "Download UMAP (PNG)", icon = icon("download")),
        br(), br(),
        plotOutput("umapPlot", height = "700px")
      ),
      
      tabItem(
        tabName = "lncDiscovery",
        h3("lncRNA Discovery (GENCODE v44)"),
        actionButton("runLncDiscovery", "Run lncRNA Discovery", icon = icon("search")),
        br(), br(),
        h4("Detected lncRNAs"),
        DT::dataTableOutput("knownLncTable"),
        br(),
        fluidRow(
          column(6, downloadButton("downloadKnownLnc", "Download All lncRNAs", icon = icon("download"))),
          column(6, downloadButton("downloadTopKnownLnc", "Download Top 20 lncRNAs", icon = icon("download")))
        ),
        br(),
        actionButton("plotKnownHeatmap", "Generate Top 20 lncRNA Heatmap", icon = icon("chart-area")),
        br(), br(),
        downloadButton("downloadKnownHeatmap", "Download Heatmap (PNG)", icon = icon("download")),
        br(), br(),
        plotOutput("knownHeatmapPlot", height = "600px")
      ),
      
      tabItem(tabName = "pseudotime", pseudotime_ui("pseudotime")),
      tabItem(tabName = "lncPathway", lnc_pathway_ui("lncPathway")),
      tabItem(tabName = "gsvaPathway", gsva_pathway_ui("gsvaPathway")),
      tabItem(tabName = "ldsMethod", lds_ui("ldsMethod")),
      tabItem(tabName = "tld_analysis", tld_ui("tldModule"))
    )
  )
)
