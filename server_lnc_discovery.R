## server_lnc_discovery.R

lnc_discovery_server <- function(input, output, session, seurat_obj_reactive) {
  
  lnc_long_reactive <- reactiveVal(NULL)
  
  observeEvent(input$runLncDiscovery, {
    req(seurat_obj_reactive())
    
    showNotification("⏳ Running lncRNA detection...", type = "message", id = "lncNote")
    
    seu <- seurat_obj_reactive()
    expressed_genes <- rownames(seu)
    
    if (!exists("lncRNA_list")) {
      showNotification("❌ lncRNA_list not loaded. Check your global.R.", type = "error")
      removeNotification("lncNote")
      return(NULL)
    }
    
    detected <- lncRNA_list %>%
      dplyr::filter(gene_name %in% expressed_genes | gene_id %in% expressed_genes) %>%
      dplyr::distinct(gene_id, gene_name, gene_type, chromosome, start, end, strand)
    
    if (nrow(detected) > 0) {
      avg_expr <- Seurat::AverageExpression(
        seu,
        features = unique(detected$gene_name),
        assays = "RNA",
        group.by = "seurat_clusters",
        return.seurat = FALSE
      )$RNA
      
      avg_df <- avg_expr %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_name") %>%
        tidyr::pivot_longer(-gene_name, names_to = "Cluster", values_to = "Avg_Expression") %>%
        dplyr::mutate(
          Cluster = gsub("^g", "", Cluster)  # remove 'g' prefix for clarity
        )
      
      lnc_final <- dplyr::left_join(avg_df, detected, by = "gene_name") %>%
        dplyr::filter(Avg_Expression > 0)
    } else {
      lnc_final <- data.frame(Message = "No lncRNAs detected in this dataset.")
    }
    
    lnc_long_reactive(lnc_final)
    
    output$knownLncTable <- DT::renderDataTable({
      req(lnc_long_reactive())
      DT::datatable(lnc_long_reactive(), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    removeNotification("lncNote")
    showNotification("✅ lncRNA detection completed.", type = "message")
  })
  
  ## ==== Download All Detected lncRNAs ====
  output$downloadKnownLnc <- downloadHandler(
    filename = function() paste0("All_lncRNAs_", Sys.Date(), ".csv"),
    content = function(file) {
      req(lnc_long_reactive())
      write.csv(lnc_long_reactive(), file, row.names = FALSE)
    }
  )
  
  ## ==== Download Top 20 lncRNAs ====
  output$downloadTopKnownLnc <- downloadHandler(
    filename = function() paste0("Top20_lncRNAs_", Sys.Date(), ".csv"),
    content = function(file) {
      req(lnc_long_reactive())
      top_df <- lnc_long_reactive() %>%
        dplyr::group_by(gene_name) %>%
        dplyr::slice_max(order_by = Avg_Expression, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(desc(Avg_Expression)) %>%
        dplyr::slice_head(n = 20)
      write.csv(top_df, file, row.names = FALSE)
    }
  )
  
  ## ==== Generate and Download Top 20 lncRNA Heatmap ====
  observeEvent(input$plotKnownHeatmap, {
    req(lnc_long_reactive())
    
    top20 <- lnc_long_reactive() %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarise(mean_exp = mean(Avg_Expression), .groups = "drop") %>%
      dplyr::arrange(desc(mean_exp)) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::pull(gene_name)
    
    plot_data <- lnc_long_reactive() %>%
      dplyr::filter(gene_name %in% top20) %>%
      dplyr::select(gene_name, Cluster, Avg_Expression) %>%
      dplyr::mutate(Cluster = gsub("^g", "", Cluster)) %>%
      tidyr::pivot_wider(names_from = Cluster, values_from = Avg_Expression, values_fill = 0) %>%
      tibble::column_to_rownames("gene_name") %>%
      as.matrix()
    
    output$knownHeatmapPlot <- renderPlot({
      ComplexHeatmap::Heatmap(
        plot_data,
        name = "Avg_Expression",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        column_title = "Top 20 lncRNAs Heatmap",
        heatmap_legend_param = list(title = "Expression")
      )
    })
    
    output$downloadKnownHeatmap <- downloadHandler(
      filename = function() paste0("Top20_lncRNAs_Heatmap_", Sys.Date(), ".png"),
      content = function(file) {
        png(file, width = 2000, height = 1600, res = 1000)
        ComplexHeatmap::draw(
          ComplexHeatmap::Heatmap(
            plot_data,
            name = "Avg_Expression",
            row_names_gp = grid::gpar(fontsize = 10),
            column_names_gp = grid::gpar(fontsize = 10),
            column_title = "Top 20 lncRNAs Heatmap",
            heatmap_legend_param = list(title = "Expression")
          )
        )
        dev.off()
      }
    )
  })
}
