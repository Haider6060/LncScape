lnc_pathway_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    uploaded_data <- reactiveVal(NULL)
    enrichment_results <- reactiveVal(NULL)
    plot_ready <- reactiveVal(FALSE)
    
    output$dataset_selector_ui <- renderUI({
      selectInput(
        ns("selected_dataset"),
        "Select Enrichment Dataset:",
        choices = c(
          "Cancer Hallmark" = "Cancer_Hallmark_LncSEA.rds",
          "Disease Phenotype" = "Disease_Type_LncSEA.rds",
          "Drug Associations" = "Drug_LncSEA.rds",
          "Survival Associations" = "Survival_LncSEA.rds",
          "lncRNADisease Clean" = "lncRNADisease_clean.rds",
          "NPInter lncRNA Interaction" = "NPInter_lncRNA_interaction_clean.rds"
        ),
        selected = "Cancer_Hallmark_LncSEA.rds"
      )
    })
    
    observeEvent(input$selected_dataset, {
      plot_ready(FALSE) # Prevent auto plot generation
    })
    
    observeEvent(input$lnc_pathway_file, {
      req(input$lnc_pathway_file)
      df <- readr::read_csv(input$lnc_pathway_file$datapath, show_col_types = FALSE)
      uploaded_data(df)
      output$lnc_pathway_table <- DT::renderDataTable({
        DT::datatable(df, options = list(pageLength = 5))
      })
    })
    
    observeEvent(input$run_lnc_pathway, {
      req(uploaded_data())
      tryCatch({
        selected_file <- input$selected_dataset
        gsc <- readRDS(file.path("E:/New_Tool_Split_App/prepare_data", selected_file))
        
        if (selected_file == "lncRNADisease_clean.rds") {
          gene_sets <- split(gsc$lncRNA, gsc$disease)
          names(gene_sets) <- names(gene_sets)
        } else if (selected_file == "NPInter_lncRNA_interaction_clean.rds") {
          gene_sets <- split(as.character(gsc$ncName), gsc$tarType)
          names(gene_sets) <- names(gene_sets)
        } else {
          pathway_names <- sapply(gsc, GSEABase::setName)
          gene_sets <- lapply(gsc, GSEABase::geneIds)
          names(gene_sets) <- sapply(strsplit(pathway_names, ";"), `[`, 1)
        }
        
        all_lncrnas <- uploaded_data()
        lnc_list <- unique(all_lncrnas$gene_name)
        
        pathway_hits <- sapply(gene_sets, function(genes) {
          length(intersect(toupper(genes), toupper(lnc_list)))
        })
        pathway_hits <- sort(pathway_hits, decreasing = TRUE)
        
        all_matched_labels <- sapply(names(pathway_hits), function(pw) {
          genes <- gene_sets[[pw]]
          matched_lncs <- intersect(toupper(genes), toupper(lnc_list))
          matched_lncs <- unique(matched_lncs)
          paste(matched_lncs, collapse = ", ")
        })
        
        top_lncrna_labels <- sapply(names(pathway_hits), function(pw) {
          genes <- gene_sets[[pw]]
          matched_lncs <- intersect(toupper(genes), toupper(lnc_list))
          matched_lncs <- unique(matched_lncs)
          top5 <- head(matched_lncs, 5)
          paste(top5, collapse = ", ")
        })
        
        df_enrich <- data.frame(
          Pathway = names(pathway_hits),
          Matched_lncRNA_Count = pathway_hits,
          All_Matched_lncRNAs = all_matched_labels,
          Top_5_lncRNAs = top_lncrna_labels
        )
        
        if (selected_file != "Cancer_Hallmark_LncSEA.rds") {
          df_enrich <- head(df_enrich, 30)
        }
        
        enrichment_results(df_enrich)
        plot_ready(FALSE) # Ensure plot will generate only after explicit click
        showNotification("Pathway enrichment completed.", type = "message")
      }, error = function(e) {
        showNotification(paste("Enrichment error:", e$message), type = "error")
      })
    })
    
    observeEvent(input$plot_lnc_pathway, {
      req(enrichment_results())
      plot_ready(TRUE)
    })
    
    output$download_lnc_pathway <- downloadHandler(
      filename = function() {
        paste0(tools::file_path_sans_ext(input$selected_dataset), "_Pathway_Enrichment_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(enrichment_results())
        write.csv(enrichment_results(), file, row.names = FALSE)
      }
    )
    
    output$lnc_pathway_plot <- renderPlot({
      req(plot_ready())
      req(enrichment_results())
      
      selected_file <- input$selected_dataset
      plot_title <- "Top Pathways"
      plot_type <- if (selected_file %in% c("Cancer_Hallmark_LncSEA.rds", "Survival_LncSEA.rds")) "bar"
      else if (selected_file %in% c("Disease_Type_LncSEA.rds", "NPInter_lncRNA_interaction_clean.rds")) "lollipop"
      else "dot"
      
      df_plot <- data.frame(
        Pathway = enrichment_results()$Pathway,
        Count = enrichment_results()$Matched_lncRNA_Count,
        Top_lncRNAs = enrichment_results()$Top_5_lncRNAs
      )
      df_plot$Top_lncRNAs_wrapped <- stringr::str_wrap(df_plot$Top_lncRNAs, width = 30)
      
      clean_theme <- theme_bw() +
        theme(
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14)
        )
      
      if (plot_type == "bar") {
        ggplot(df_plot, aes(x = Count, y = reorder(Pathway, Count))) +
          geom_bar(stat = "identity", fill = "steelblue") +
          geom_text(aes(label = Top_lncRNAs_wrapped), hjust = -0.1, size = 3) +
          labs(
            title = paste(plot_title, "- Bar Plot"),
            x = "Matched lncRNA Count",
            y = "Pathway",
            caption = "Top 5 lncRNAs beside each bar"
          ) +
          clean_theme +
          xlim(0, max(df_plot$Count, na.rm = TRUE) * 1.4)
      } else if (plot_type == "lollipop") {
        ggplot(df_plot, aes(x = Count, y = reorder(Pathway, Count))) +
          geom_segment(aes(x = 0, xend = Count, yend = Pathway), color = "grey") +
          geom_point(color = "red", size = 4) +
          geom_text(aes(label = Top_lncRNAs_wrapped), hjust = -0.1, size = 3) +
          labs(
            title = paste(plot_title, "- Lollipop Plot"),
            x = "Matched lncRNA Count",
            y = "Pathway",
            caption = "Top 5 lncRNAs beside each point"
          ) +
          clean_theme +
          xlim(0, max(df_plot$Count, na.rm = TRUE) * 1.4)
      } else if (plot_type == "dot") {
        ggplot(df_plot, aes(x = Count, y = reorder(Pathway, Count))) +
          geom_point(color = "purple", size = 5, alpha = 0.7) +
          geom_text(aes(label = Top_lncRNAs_wrapped), hjust = -0.1, size = 3) +
          labs(
            title = paste(plot_title, "- Dot Plot"),
            x = "Matched lncRNA Count",
            y = "Pathway",
            caption = "Top 5 lncRNAs beside each dot"
          ) +
          clean_theme +
          xlim(0, max(df_plot$Count, na.rm = TRUE) * 1.4)
      }
    })
    
    output$download_lnc_pathway_plot <- downloadHandler(
      filename = function() {
        type <- if (input$selected_dataset == "Cancer_Hallmark_LncSEA.rds") "Bar_Plot"
        else if (input$selected_dataset == "Disease_Type_LncSEA.rds") "Lollipop_Plot"
        else if (input$selected_dataset == "Drug_LncSEA.rds") "Dot_Plot"
        else if (input$selected_dataset == "Survival_LncSEA.rds") "Bar_Plot"
        else if (input$selected_dataset == "lncRNADisease_clean.rds") "Dot_Plot"
        else if (input$selected_dataset == "NPInter_lncRNA_interaction_clean.rds") "Lollipop_Plot"
        else "Plot"
        paste0(tools::file_path_sans_ext(input$selected_dataset), "_", type, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        png(file, width = 12, height = 8, units = "in", res = 300)
        print(eval(parse(text = deparse(substitute(output$lnc_pathway_plot())))))
        dev.off()
      }
    )
  })
}
