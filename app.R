library(shiny)
library(shinydashboard)

# ==== Load Globals and Modules ====
source("global.R")

# Load UI Modules FIRST
source("pseudotime_ui.R")      # ✅ Load this BEFORE main_ui.R
source("lnc_pathway_ui.R")
source("gsva_pathway_ui.R")
source("LDS_ui.R")
source("ui_tld_module.R")
source("ui_upload_clustering.R")

# Load Main UI LAST (because it calls the modules)
source("main_ui.R")             # ✅ Load this AFTER all ui modules

# Load Server Modules
source("server_upload_clustering.R")
source("server_lnc_discovery.R")
source("pseudotime_server.R")
source("lnc_pathway_server.R")
source("gsva_pathway_server.R")
source("LDS_server.R")
source("server_tld_module.R")

# Reactive Containers
seurat_obj_reactive <- reactiveVal(NULL)
dataset_type <- reactiveVal(NULL)
cds_reactive <- reactiveVal(NULL)

# Server
server <- function(input, output, session) {
  upload_clustering_server(input, output, session, seurat_obj_reactive, dataset_type, cds_reactive)
  lnc_discovery_server(input, output, session, seurat_obj_reactive)
  pseudotime_server("pseudotime", cds_reactive)
  lnc_pathway_server("lncPathway")
  gsva_pathway_server("gsvaPathway")
  lds_server("ldsMethod", cds_reactive)
  tld_server("tldModule", seurat_obj_reactive)
}

# Launch App
shinyApp(ui = lncscape_ui, server = server)
