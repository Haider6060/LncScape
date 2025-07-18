gsva_pathway_server <- function(id) { 
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    uploaded_lncrna <- reactiveVal(NULL)
    gsva_scores <- reactiveVal(NULL)
    plot_ready <- reactiveVal(FALSE)
    
    # Reset plot when dataset changes
    observeEvent(input$gsva_dataset, {
      plot_ready(FALSE)
      gsva_scores(NULL)
    })
    
    # Load uploaded lncRNA CSV
    observeEvent(input$gsva_file, {
      req(input$gsva_file)
      df <- readr::read_csv(input$gsva_file$datapath, show_col_types = FALSE)
      uploaded_lncrna(df)
      
      output$gsva_table <- DT::renderDataTable({
        DT::datatable(df, options = list(pageLength = 10))
      })
    })
    
    # Run GSVA enrichment
    observeEvent(input$run_gsva, {
      req(uploaded_lncrna())
      tryCatch({
        selected_file <- input$gsva_dataset
        gsva_file_path <- file.path("E:/New_Tool_Split_App/prepare_data", selected_file)
        gsva_data <- readRDS(gsva_file_path)
        
        # Only for NPINTER, extract using Biobase
        if (selected_file == "gsva_scores_npinter.rds" && inherits(gsva_data, "ExpressionSet")) {
          gsva_data <- Biobase::exprs(gsva_data)
        }
        
        if (!is.null(rownames(gsva_data))) {
          pathway_names <- sapply(strsplit(rownames(gsva_data), ";"), `[`, 1)
          rownames(gsva_data) <- pathway_names
        }
        
        df_gsva <- data.frame(
          Pathway = rownames(gsva_data),
          gsva_data,
          check.names = FALSE
        )
        
        gsva_scores(df_gsva)
        plot_ready(FALSE)
        showNotification("✅ GSVA enrichment completed.", type = "message")
      }, error = function(e) {
        showNotification(paste("❌ GSVA error:", e$message), type = "error")
      })
    })
    
    # Enable plot generation explicitly
    observeEvent(input$plot_gsva, {
      req(gsva_scores())
      plot_ready(TRUE)
    })
    
    # Download CSV
    output$download_gsva_csv <- downloadHandler(
      filename = function() {
        paste0("GSVA_Pathway_Enrichment_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(gsva_scores())
        write.csv(gsva_scores(), file, row.names = FALSE)
      }
    )
    
    # Reactive plot object
    gsva_plot_obj <- reactive({
      req(plot_ready())
      req(gsva_scores())
      
      selected_file <- input$gsva_dataset
      
      df <- gsva_scores()
      df_long <- tidyr::pivot_longer(
        df,
        cols = -Pathway,
        names_to = "Cluster",
        values_to = "EnrichmentScore"
      )
      
      df_long$Cluster <- factor(df_long$Cluster, levels = gtools::mixedsort(unique(df_long$Cluster)))
      
      # Apply NPINTER-specific top 20 bar logic ONLY for NPINTER
      if (selected_file == "gsva_scores_npinter.rds") {
        pathway_vars <- df_long %>%
          group_by(Pathway) %>%
          summarise(Var = var(EnrichmentScore, na.rm = TRUE)) %>%
          arrange(desc(Var)) %>%
          slice_head(n = 20)
        selected_paths <- pathway_vars$Pathway
        df_plot <- df_long %>% filter(Pathway %in% selected_paths)
        plot_type <- "bar"
      } else {
        # Your original logic untouched for all other files
        total_paths <- length(unique(df_long$Pathway))
        if (total_paths <= 30) {
          selected_paths <- unique(df_long$Pathway)
        } else {
          pathway_vars <- df_long %>%
            group_by(Pathway) %>%
            summarise(Var = var(EnrichmentScore, na.rm = TRUE)) %>%
            arrange(desc(Var)) %>%
            slice_head(n = 30)
          selected_paths <- pathway_vars$Pathway
        }
        df_plot <- df_long %>% filter(Pathway %in% selected_paths)
        
        plot_type <- switch(
          selected_file,
          "GSVA_Cancer_Hallmark_ssGSEA_scores.rds" = "boxplot",
          "GSVA_Disease_Type_scores.rds" = "violin",
          "gsva_scores_drug.rds" = "bar",
          "gsva_scores_survival.rds" = "dot",
          "gsva_lncrnadisease_scores.rds" = "bar_horizontal",
          "boxplot"
        )
      }
      
      # Use your exact plotting logic unchanged
      if (plot_type == "boxplot") {
        p <- ggplot(df_plot, aes(x = Cluster, y = EnrichmentScore, fill = Pathway)) +
          geom_boxplot(outlier.size = 0.4, alpha = 0.8, position = position_dodge(width = 0.8))
      } else if (plot_type == "violin") {
        p <- ggplot(df_plot, aes(x = Cluster, y = EnrichmentScore, fill = Pathway)) +
          geom_violin(alpha = 0.7, scale = "width", trim = TRUE)
      } else if (plot_type == "bar") {
        p <- ggplot(df_plot, aes(x = Cluster, y = EnrichmentScore, fill = Pathway)) +
          geom_bar(stat = "identity", position = position_dodge(width = 0.8))
      } else if (plot_type == "dot") {
        p <- ggplot(df_plot, aes(x = EnrichmentScore, y = reorder(Pathway, EnrichmentScore), color = Cluster)) +
          geom_point(size = 3, alpha = 0.8) +
          facet_wrap(~Cluster, scales = "free_x")
      } else if (plot_type == "bar_horizontal") {
        p <- ggplot(df_plot, aes(x = EnrichmentScore, y = reorder(Pathway, EnrichmentScore), fill = Cluster)) +
          geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
          coord_flip()
      }
      
      p <- p +
        labs(
          title = paste0("GSVA Enrichment Scores (", tools::toTitleCase(gsub("_", " ", plot_type)), ")"),
          x = ifelse(plot_type %in% c("dot", "bar_horizontal"), "Enrichment Score", "Cluster"),
          y = ifelse(plot_type %in% c("dot", "bar_horizontal"), "Pathway", "Enrichment Score")
        ) +
        theme_classic(base_size = 16) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.title = element_blank(),
          legend.position = "bottom",
          panel.grid = element_blank(),
          panel.border = element_blank()
        ) +
        scale_fill_viridis_d(option = "plasma") +
        scale_color_viridis_d(option = "plasma")
      
      return(p)
    })
    
    # Render plot
    output$gsva_plot <- renderPlot({
      gsva_plot_obj()
    })
    
    # Download high-res PNG
    output$download_gsva_plot <- downloadHandler(
      filename = function() {
        paste0("GSVA_Pathway_Enrichment_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(plot_ready())
        ggsave(
          filename = file,
          plot = gsva_plot_obj(),
          device = "png",
          dpi = 1000,
          width = 12,
          height = 8,
          units = "in"
        )
      }
    )
  })
}
