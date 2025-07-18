upload_clustering_server <- function(input, output, session,
                                     seurat_obj_reactive,
                                     dataset_type,
                                     cds_reactive) {
  
  observeEvent(input$fileUpload, {
    req(input$fileUpload)
    
    # Load Seurat object safely
    seurat_obj <- tryCatch({
      readRDS(input$fileUpload$datapath)
    }, error = function(e) {
      showNotification("❌ Failed to load Seurat object.", type = "error")
      message("❌ Failed to load Seurat object: ", e$message)
      return(NULL)
    })
    
    if (is.null(seurat_obj)) return(NULL)
    
    seurat_obj_reactive(seurat_obj)
    showNotification("✅ Seurat object loaded successfully.", type = "message")
    message("✅ Seurat object loaded successfully.")
    
    # Dataset type detection
    dataset_type_value <- if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
      if (any(grepl("^GSM", seurat_obj$orig.ident))) {
        "GEO Dataset"
      } else {
        "10x Dataset"
      }
    } else {
      "Unknown Format"
    }
    dataset_type(dataset_type_value)
    
    output$fileInfo <- renderPrint({ seurat_obj })
    output$datasetType <- renderText({ paste("Dataset Type:", dataset_type_value) })
    output$metaDataSummary <- renderPrint({ head(seurat_obj@meta.data, 10) })
    
    # === CDS preparation ===
    tryCatch({
      showNotification("⏳ Preparing CDS for pseudotime...", type = "message")
      message("🟩 Starting CDS preparation pipeline...")
      
      # Check available layers safely
      available_layers <- tryCatch({
        names(seurat_obj[["RNA"]]@layers)
      }, error = function(e) {
        NULL
      })
      
      counts_layers <- if (!is.null(available_layers)) {
        grep("^counts", available_layers, value = TRUE)
      } else {
        NULL
      }
      
      # Extract counts matrix robustly
      if (!is.null(counts_layers) && length(counts_layers) > 0) {
        counts_list <- lapply(counts_layers, function(layer) {
          Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = layer)
        })
        counts_matrix <- do.call(cbind, counts_list)
      } else {
        counts_matrix <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
      }
      
      # Align metadata
      cell_metadata <- seurat_obj@meta.data[colnames(counts_matrix), , drop = FALSE]
      gene_metadata <- data.frame(
        gene_short_name = rownames(counts_matrix),
        row.names = rownames(counts_matrix)
      )
      
      # Create CDS object
      cds <- monocle3::new_cell_data_set(
        expression_data = as.matrix(counts_matrix),
        cell_metadata = cell_metadata,
        gene_metadata = gene_metadata
      )
      
      cds <- monocle3::preprocess_cds(cds, num_dim = 50)
      cds <- monocle3::reduce_dimension(cds)
      cds <- monocle3::cluster_cells(cds)
      cds <- monocle3::learn_graph(cds)
      
      # 🚩 Avoid GUI by explicitly providing root_cells
      root_cells <- colnames(cds)[1]
      cds <- monocle3::order_cells(cds, root_cells = root_cells)
      
      # ✅ Attach the Seurat object for LDS module
      attr(cds, "seurat_obj") <- seurat_obj
      
      cds_reactive(cds)
      showNotification("✅ Monocle3 pseudotime preparation completed.", type = "message")
      message("✅ Monocle3 pseudotime preparation completed.")
      
    }, error = function(e) {
      showNotification(paste("❌ Error during CDS preparation:", e$message), type = "error")
      message("❌ Error during CDS preparation: ", e$message)
    })
  })
  
  # === Clustering & UMAP on button click only ===
  observeEvent(input$runClustering, {
    req(seurat_obj_reactive())
    
    seurat_obj <- seurat_obj_reactive()
    showNotification("⏳ Running clustering and UMAP...", type = "message")
    message("🟩 Running clustering and UMAP...")
    
    tryCatch({
      seurat_obj <- Seurat::NormalizeData(seurat_obj)
      seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
      seurat_obj <- Seurat::ScaleData(seurat_obj)
      seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 50)
      seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
      seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = input$resolutionInput)
      
      # === FIX: CLEAN CLUSTER LABELS ===
      clusters <- seurat_obj$seurat_clusters
      if (is.factor(clusters)) {
        seurat_obj$seurat_clusters <- as.integer(as.character(clusters))
      } else if (is.character(clusters)) {
        seurat_obj$seurat_clusters <- as.integer(clusters)
      }
      
      seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:30)
      
      seurat_obj_reactive(seurat_obj)
      
      output$umapPlot <- renderPlot({
        Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
          ggplot2::ggtitle("Seurat UMAP Clustering")
      })
      
      output$downloadPlot <- downloadHandler(
        filename = function() { paste0("UMAP_Plot_", Sys.Date(), ".png") },
        content = function(file) {
          png(file, width = 2000, height = 1600, res = 1000)
          print(Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
                  ggplot2::ggtitle("Seurat UMAP Clustering"))
          dev.off()
        }
      )
      
      showNotification("✅ Clustering and UMAP completed.", type = "message")
      message("✅ Clustering and UMAP completed.")
      
    }, error = function(e) {
      showNotification(paste("❌ Error during clustering:", e$message), type = "error")
      message("❌ Error during clustering: ", e$message)
    })
  })
}
