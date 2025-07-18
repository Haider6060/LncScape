options(shiny.maxRequestSize = 1024 * 1024 ^ 2)

# ==== Libraries ====
library(shiny)
library(shinydashboard)
library(Seurat)
library(tools)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(DT)
library(tidyr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(pheatmap)
library(tibble)
library(slingshot)
library(gtools)       # Added for cluster mixedsort
library(GSEABase)
library(heatmaply)
library("shinybusy")
# ==== New Libraries for TF–lncRNA Dynamics Module ====
library(BSgenome.Hsapiens.UCSC.hg38)   # For promoter extraction
library(universalmotif)                # For motif import and handling
library(motifmatchr)                   # For motif scanning (optional, used in TF-lncRNA pairing)
# library(shinybusy)                   # Optional: For advanced loading spinners (if needed)

# ==== Global Settings ====
set.seed(1234)

# ==== Load GENCODE v44 GTF for lncRNA Discovery ====
gtf_path <- "E:/New_Tool_Split_App/Data/gencode.v44.annotation.gtf"

if (file.exists(gtf_path)) {
  
  message(paste0("✅ Loading GTF file: ", gtf_path))
  
  gtf_data <- rtracklayer::import(gtf_path)
  gtf_genes <- gtf_data[gtf_data$type == "gene"]
  
  gtf_df <- data.frame(
    gene_id = mcols(gtf_genes)$gene_id,
    gene_name = mcols(gtf_genes)$gene_name,
    gene_type = mcols(gtf_genes)$gene_type,
    chromosome = as.character(seqnames(gtf_genes)),
    start = start(gtf_genes),
    end = end(gtf_genes),
    strand = as.character(strand(gtf_genes)),
    stringsAsFactors = FALSE
  )
  
  lncRNA_list <- gtf_df %>%
    filter(gene_type %in% c(
      "lncRNA",
      "3prime_overlapping_ncRNA",
      "antisense",
      "non_coding",
      "sense_intronic",
      "sense_overlapping",
      "macro_lncRNA",
      "bidirectional_promoter_lncRNA"
    )) %>%
    distinct(gene_name, .keep_all = TRUE)
  
  message(paste0("✅ GTF loaded. lncRNAs identified: ", nrow(lncRNA_list)))
  
} else {
  stop(paste0("❌ GTF file not found at path: ", gtf_path, "\nPlease ensure the GTF file is present in your app directory."))
}
