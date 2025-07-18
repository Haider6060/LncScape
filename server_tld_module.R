tld_server <- function(id, seurat_obj_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run_tld, {
      
      showModal(modalDialog(
        "Running TF–lncRNA Dynamics Analysis...",
        footer = NULL
      ))
      
      seurat_obj <- seurat_obj_reactive()
      all_layers <- Layers(seurat_obj[["RNA"]])
      counts_layer <- grep("counts", all_layers, value = TRUE)[1]
      counts <- GetAssayData(seurat_obj, assay = "RNA", layer = counts_layer)
      cells_in_counts <- colnames(counts)
      meta <- seurat_obj@meta.data
      meta <- meta[cells_in_counts, ]
      
      library(SingleCellExperiment)
      library(slingshot)
      
      sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = meta
      )
      
      # UMAP check and run if missing
      if (!"umap" %in% names(seurat_obj@reductions)) {
        seurat_obj <- RunPCA(seurat_obj)
        seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
      }
      
      umap <- Embeddings(seurat_obj, "umap")
      umap <- umap[cells_in_counts, ]
      reducedDims(sce) <- SimpleList(UMAP = umap)
      
      # Clustering check and run if missing
      if (is.null(seurat_obj$seurat_clusters)) {
        seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
        seurat_obj <- FindClusters(seurat_obj)
      }
      
      sce <- slingshot(
        sce,
        clusterLabels = seurat_obj$seurat_clusters[cells_in_counts],
        reducedDim = 'UMAP'
      )
      
      pseudo <- slingPseudotime(sce)
      pseudo_df <- as.data.frame(pseudo)
      pseudo_df$Cell_ID <- rownames(pseudo_df)
      
      # --- Load LDS and TF Pairs ---
      TF_pairs <- readRDS("E:/New_Tool_Split_App/Data/TF_target.rds")
      lds_scores <- read.csv(input$lds_file$datapath)
      
      top_lncRNAs <- lds_scores[order(-lds_scores$LDS_score), "gene_name"][1:100]
      filtered_pairs <- TF_pairs[TF_pairs$lncRNA %in% top_lncRNAs, ]
      
      data_layer <- grep("data", all_layers, value = TRUE)[1]
      expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = data_layer)
      
      TF_list <- unique(filtered_pairs$TF)
      lncRNA_list <- unique(filtered_pairs$lncRNA)
      
      common_TFs <- intersect(TF_list, rownames(expr_matrix))
      common_lncRNAs <- intersect(lncRNA_list, rownames(expr_matrix))
      
      TF_expr <- expr_matrix[common_TFs, ]
      lncRNA_expr <- expr_matrix[common_lncRNAs, ]
      
      pseudo_vector <- pseudo_df[, 1]
      names(pseudo_vector) <- pseudo_df$Cell_ID
      common_cells <- intersect(colnames(expr_matrix), names(pseudo_vector))
      
      pseudo_vector <- pseudo_vector[common_cells]
      TF_expr <- TF_expr[, common_cells]
      lncRNA_expr <- lncRNA_expr[, common_cells]
      
      top_pairs <- filtered_pairs[1:100, ]
      cor_results <- list()
      
      for (i in 1:nrow(top_pairs)) {
        tf <- top_pairs$TF[i]
        lnc <- top_pairs$lncRNA[i]
        if (tf %in% rownames(TF_expr) && lnc %in% rownames(lncRNA_expr)) {
          tf_values <- as.numeric(TF_expr[tf, ])
          lnc_values <- as.numeric(lncRNA_expr[lnc, ])
          if (sd(tf_values) > 0 && sd(lnc_values) > 0) {
            cor_val <- cor(tf_values, lnc_values, method = "pearson")
          } else {
            cor_val <- NA
          }
          cor_results[[i]] <- data.frame(TF = tf, lncRNA = lnc, Correlation = cor_val)
        }
      }
      
      cor_table <- do.call(rbind, cor_results)
      merged_table <- merge(cor_table, lds_scores, by.x = "lncRNA", by.y = "gene_name", all.x = TRUE)
      merged_table$TLD_Score <- merged_table$Correlation * merged_table$LDS_score
      
      # Save merged_table to reactive value for plotting later
      tld_results <- reactiveVal(merged_table)
      
      # --- Download Handler for TLD Results ---
      output$download_tld_results <- downloadHandler(
        filename = function() { "TLD_results.csv" },
        content = function(file) { write.csv(tld_results(), file, row.names = FALSE) }
      )
      
      # --- Plot Generation only when button is clicked ---
      observeEvent(input$generate_plot, {
        merged_table <- tld_results()
        
        top20 <- merged_table %>%
          filter(lncRNA == "MALAT1" & !is.na(TLD_Score)) %>%
          arrange(-abs(TLD_Score)) %>%
          head(20)
        
        top20 <- top20 %>%
          mutate(Regulation = ifelse(TLD_Score > 0, "Activation (Blue)", "Repression (Red)"))
        
        top20$TF <- factor(top20$TF, levels = unique(top20$TF[order(top20$TLD_Score)]))
        
        p1 <- ggplot(top20, aes(x = TF, y = TLD_Score, fill = Regulation)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          labs(title = "Top 20 TFs Interacting with MALAT1",
               x = "TF", y = "TLD-Score", fill = "Interaction Type") +
          scale_fill_manual(values = c("Activation (Blue)" = "steelblue", "Repression (Red)" = "red")) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 10),
                legend.position = "right",
                plot.title = element_text(size = 14, face = "bold"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white"))
        
        output$tld_barplot <- renderPlot({ p1 })
        
        output$download_tld_barplot <- downloadHandler(
          filename = function() { "TLD_barplot.png" },
          content = function(file) { ggsave(file, p1, width = 8, height = 6, dpi = 1000) }
        )
        
        showNotification("✅ TLD Barplot Generated!", type = "message", duration = 4)
      })
      
      # --- Final Notifications ---
      removeModal()
      showNotification("✅ TLD Analysis Completed! Click 'Generate Plot' to visualize.", type = "message", duration = 5)
      
    })
  })
}
