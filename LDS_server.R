lds_server <- function(id, cds) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    lnc_data <- reactiveVal(NULL)
    lds_result <- reactiveVal(NULL)
    lds_plot_obj <- reactiveVal(NULL)
    
    # Disable Top 50 button initially
    shinyjs::disable(ns("downloadTop50LDS"))
    
    # Upload CSV
    observeEvent(input$lds_csv, {
      req(input$lds_csv)
      df <- readr::read_csv(input$lds_csv$datapath, show_col_types = FALSE)
      lnc_data(df)
      output$lnc_table <- DT::renderDataTable({
        DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
      })
      showNotification("✅ lncRNA file uploaded.", type = "message")
    })
    
    # Run LDS
    observeEvent(input$runLDS, {
      if (is.null(cds())) {
        showNotification("❌ Pseudotime not available yet.", type = "error")
        return(NULL)
      }
      
      if (is.null(attr(cds(), "seurat_obj"))) {
        showNotification("❌ Seurat object not attached to CDS.", type = "error")
        return(NULL)
      }
      
      ptime <- monocle3::pseudotime(cds())
      valid_cells <- names(ptime[is.finite(ptime)])
      
      seurat_obj <- attr(cds(), "seurat_obj")
      available_layers <- names(seurat_obj[["RNA"]]@layers)
      counts_layers <- grep("^counts", available_layers, value = TRUE)
      
      counts_list <- lapply(counts_layers, function(layer) {
        Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = layer)
      })
      combined_counts <- do.call(cbind, counts_list)
      
      lnc_genes <- unique(lnc_data()$gene_name)
      lnc_expr <- combined_counts[rownames(combined_counts) %in% lnc_genes, valid_cells]
      
      ptime_vec <- ptime[colnames(lnc_expr)]
      lnc_names <- rownames(lnc_expr)
      
      lds_scores <- sapply(lnc_names, function(gene) {
        expr <- as.numeric(lnc_expr[gene, ])
        if (all(expr == 0)) return(0)
        var(expr) * abs(cor(expr, ptime_vec, method = "spearman"))
      })
      
      lds_df <- data.frame(gene_name = lnc_names, LDS_score = lds_scores)
      lds_df <- lds_df[order(-lds_df$LDS_score), ]
      
      lds_result(lds_df)
      
      # Enable top 50 button after calculation
      shinyjs::enable(ns("downloadTop50LDS"))
      
      showNotification("✅ LDS calculation completed.", type = "message")
    })
    
    # Plot LDS
    observeEvent(input$plotLDS, {
      req(lds_result())
      req(cds())
      
      ptime <- monocle3::pseudotime(cds())
      valid_cells <- names(ptime[is.finite(ptime)])
      
      seurat_obj <- attr(cds(), "seurat_obj")
      available_layers <- names(seurat_obj[["RNA"]]@layers)
      counts_layers <- grep("^counts", available_layers, value = TRUE)
      
      counts_list <- lapply(counts_layers, function(layer) {
        Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = layer)
      })
      combined_counts <- do.call(cbind, counts_list)
      
      top_genes <- head(lds_result()$gene_name, 5)
      
      plot_df <- do.call(rbind, lapply(top_genes, function(gene) {
        data.frame(
          gene = gene,
          pseudotime = ptime[valid_cells],
          expression = as.numeric(combined_counts[gene, valid_cells])
        )
      }))
      
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pseudotime, y = expression)) +
        ggplot2::geom_point(alpha = 0.3, size = 0.5, color = "blue") +
        ggplot2::geom_smooth(method = "loess", color = "red", se = FALSE) +
        ggplot2::facet_wrap(~gene, scales = "free_y") +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "Top 5 lncRNAs by LDS Score", x = "Pseudotime", y = "lncRNA Expression")
      
      lds_plot_obj(p)
      output$lds_plot <- renderPlot({ p })
    })
    
    # Downloads
    output$downloadAllLDS <- downloadHandler(
      filename = function() { "All_LDS_Scores.csv" },
      content = function(file) {
        write.csv(lds_result(), file, row.names = FALSE)
      }
    )
    
    output$downloadTop50LDS <- downloadHandler(
      filename = function() { "Top50_LDS_Scores.csv" },
      content = function(file) {
        write.csv(head(lds_result(), 50), file, row.names = FALSE)
      }
    )
    
    output$downloadPlotLDS <- downloadHandler(
      filename = function() { "LDS_Top5_Plot.png" },
      content = function(file) {
        png(file, width = 10, height = 8, units = "in", res = 1000)
        print(lds_plot_obj())
        dev.off()
      }
    )
  })
}
