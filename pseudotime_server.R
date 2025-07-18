pseudotime_server <- function(id, cds) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    top20_data <- reactiveVal(NULL)
    
    # === Upload CSV and update dropdown ===
    observeEvent(input$lnc_file, {
      req(input$lnc_file)
      df <- tryCatch({
        readr::read_csv(input$lnc_file$datapath, show_col_types = FALSE)
      }, error = function(e) {
        showNotification("❌ Failed to read CSV.", type = "error")
        return(NULL)
      })
      if (!"gene_name" %in% names(df)) {
        showNotification("❌ 'gene_name' column missing in CSV.", type = "error")
        return(NULL)
      }
      top20_data(df)
      output$lnc_table <- DT::renderDataTable({
        DT::datatable(df, options = list(pageLength = 5, scrollX = TRUE))
      })
      
      # ✅ Only allow gene names that exist in CDS
      cds_obj <- cds()
      if (!is.null(cds_obj)) {
        available_genes <- rownames(SingleCellExperiment::counts(cds_obj))
        valid_genes <- intersect(df$gene_name, available_genes)
        updateSelectInput(session, "selected_gene", choices = sort(valid_genes), selected = NULL)
      } else {
        updateSelectInput(session, "selected_gene", choices = NULL, selected = NULL)
      }
      
      showNotification("✅ Top 20 lncRNAs loaded.", type = "message")
    })
    
    # === Generate publication-style violin plot ===
    observeEvent(input$run_plot, {
      req(cds())
      req(input$selected_gene)
      cds_obj <- cds()
      gene <- input$selected_gene
      
      if (!(gene %in% rownames(SingleCellExperiment::counts(cds_obj)))) {
        showNotification(paste0("❌ Gene ", gene, " not found in CDS."), type = "error")
        return(NULL)
      }
      
      output$pseudo_plot <- renderPlot({
        suppressWarnings({
          pseudotime_vals <- monocle3::pseudotime(cds_obj)
          expr_vals <- SingleCellExperiment::counts(cds_obj)[gene, ]
          expr_log <- log1p(expr_vals)
          
          df_plot <- data.frame(
            Pseudotime = pseudotime_vals,
            Expression = expr_log
          )
          df_plot <- df_plot[is.finite(df_plot$Pseudotime), ]
          df_plot$PseudotimeBin <- cut(df_plot$Pseudotime, breaks = 10, labels = FALSE)
          
          ggplot(df_plot, aes(x = factor(PseudotimeBin), y = Expression)) +
            geom_violin(fill = "#56B4E9", alpha = 0.6, scale = "width") +
            geom_jitter(width = 0.2, size = 0.3, alpha = 0.2) +
            stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", size = 1) +
            scale_y_continuous(trans = "log1p") +
            labs(
              x = "Pseudotime Bins",
              y = "Normalized Expression (log1p)",
              title = paste("Pseudotime Expression of", gene)
            ) +
            theme_classic(base_size = 14) +
            theme(
              plot.title = element_text(hjust = 0.5),
              axis.text = element_text(color = "black"),
              axis.title = element_text(face = "bold")
            )
        })
      })
    })
    
    # === Download publication-ready plot ===
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("Pseudotime_", input$selected_gene, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(cds())
        req(input$selected_gene)
        cds_obj <- cds()
        gene <- input$selected_gene
        
        if (!(gene %in% rownames(SingleCellExperiment::counts(cds_obj)))) {
          showNotification(paste0("❌ Gene ", gene, " not found in CDS."), type = "error")
          return(NULL)
        }
        
        pseudotime_vals <- monocle3::pseudotime(cds_obj)
        expr_vals <- SingleCellExperiment::counts(cds_obj)[gene, ]
        expr_log <- log1p(expr_vals)
        
        df_plot <- data.frame(
          Pseudotime = pseudotime_vals,
          Expression = expr_log
        )
        df_plot <- df_plot[is.finite(df_plot$Pseudotime), ]
        df_plot$PseudotimeBin <- cut(df_plot$Pseudotime, breaks = 10, labels = FALSE)
        
        p <- ggplot(df_plot, aes(x = factor(PseudotimeBin), y = Expression)) +
          geom_violin(fill = "#56B4E9", alpha = 0.6, scale = "width") +
          geom_jitter(width = 0.2, size = 0.3, alpha = 0.2) +
          stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", size = 1) +
          scale_y_continuous(trans = "log1p") +
          labs(
            x = "Pseudotime Bins",
            y = "Normalized Expression (log1p)",
            title = paste("Pseudotime Expression of", gene)
          ) +
          theme_classic(base_size = 14) +
          theme(
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(color = "black"),
            axis.title = element_text(face = "bold")
          )
        
        ggsave(file, plot = p, width = 8, height = 6, dpi = 1000)
      }
    )
  })
}
