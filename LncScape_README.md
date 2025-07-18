
# LncScape: lncRNA Discovery and Dynamics Analysis from Single-Cell RNA-Seq Data

## Overview

**LncScape** is a **Shiny application** designed for researchers and bioinformaticians working on **lncRNA dynamics, TF–lncRNA interactions, and pathway enrichment analysis from single-cell RNA-seq data**.

This tool allows users to analyze **both GEO datasets and 10x Genomics datasets** in a **generalized pipeline** without requiring advanced coding skills.

## Key Features

- ✅ Data Input & Clustering (UMAP-based)
- ✅ lncRNA Discovery (GENCODE v44 reference)
- ✅ Pseudotime Analysis (Trajectory Inference)
- ✅ LDS Method (lncRNA Dynamics Score - Variance + Pseudotime Correlation)
- ✅ TF–lncRNA Dynamics Prioritization (TLD Score)
- ✅ Pathway Enrichment & GSVA Analysis

## Input Data

LncScape supports:

- **Single-cell RNA-seq data from GEO or 10x Genomics**
- **Seurat v5 `.rds` files**

### How to Prepare Your Data

- Convert your dataset to **`.rds` format in R** using Seurat (v5 recommended).
- **Any dataset can be used** after this conversion.

## Demo Data

We provide **small demo datasets** for testing:

| **File** | **Description** |
|----------|----------------|
| `small_demo_seurat.rds` | Demo GEO dataset (200 cells × all genes) |
| `small_10x_demo_seurat.rds` | Demo 10x Genomics dataset (200 cells × all genes) |

## Usage Instructions (For Non-Coders)

- **Compatible Datasets:**  
You can use **any scRNA-seq dataset** from **GEO or 10x Genomics** by converting it to `.rds` format.

- **Important:**  
If you run locally (not via the app interface), **edit the file path in your code to point to your `.rds` dataset location**.

- **Analysis Workflow:**  
1️⃣ **Run Clustering First** (this step initializes downstream analysis)  
2️⃣ After clustering, proceed with **lncRNA discovery, LDS, TLD, pseudotime, and pathway enrichment**

## How to Run LncScape

### Option 1: Run in RStudio

1. Open `app.R` in RStudio  
2. Click **Run App**

### Option 2: Run in R Console

```r
shiny::runApp("your_app_folder_path")
```

Replace `"your_app_folder_path"` with the directory where you saved the app.

## Full Data Download (Pathways, GTF, TF Files)

Due to GitHub file size limitations, the following large files are provided via Google Drive:

- 6 Pathway `.rds` files  
- 6 GSVA `.rds` files  
- `gencode.v44.annotation.gtf`  
- `TF_target.rds`

### **Download Link:**

[Download Large Data Files from Google Drive](https://drive.google.com/drive/folders/1sIHNTeEyzXJc9gMayGmIOxeev1wwSOp6?usp=sharing)

After downloading, place these files into your local `Data/` folder in the app directory.

## License

MIT License

## Contact

For questions or support, contact:  
**haider@emails.bjut.edu.cn**

## Citation

If you use **LncScape** in your research, please cite this tool appropriately in your manuscript.
