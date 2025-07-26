# LncScape: lncRNA Discovery and Dynamics Analysis from Single-Cell RNA-Seq Data

## Overview

**LncScape** is a user-friendly **Shiny application** designed for researchers and bioinformaticians working on **lncRNA discovery**, **TF–lncRNA interaction analysis**, and **pathway enrichment** using **single-cell RNA-seq (scRNA-seq)** data.

The app enables automated analysis of **both GEO and 10x Genomics datasets** in a **generalized pipeline**, with **no coding required**.

---

## ✅ Key Features

- Data input via `.rds` files (GEO/10x supported)
- UMAP-based clustering
- lncRNA detection using **GENCODE v44**
- LDS score: identifies dynamic lncRNAs over pseudotime
- TLD score: transcription factor–lncRNA prioritization
- GSVA-based pathway enrichment
- High-resolution output plots and CSVs

---

## 📂 Input Requirements

- `.rds` file (Seurat v5 format)
- Can be generated from **GEO or 10x Genomics** count matrices

### ➤ Preparing Your Data

```r
# Example (from 10x or GEO matrix):
seurat_obj <- CreateSeuratObject(counts = your_matrix)
saveRDS(seurat_obj, file = "your_dataset.rds")
```

---

## 🧪 Datasets Tested

LncScape has been tested on **five real-world datasets**:

| Dataset Type | Source | Description |
|--------------|--------|-------------|
| Lung cancer (GEO) | GEO | Public scRNA-seq dataset |
| Breast cancer (GEO) | GEO | Validated on multiple samples |
| Pancreatic cancer (GEO) | GEO | Human CAR T-cell response |
| Brain tumor (10x) | 10x Genomics | Commercial single-cell dataset |
| Lung cancer (10x) | 10x Genomics | Public dataset, high complexity |

All datasets processed successfully with full outputs: **lncRNA files**, **LDS scores**, **TLD prioritization**, and **cluster plots**.

---

## 🎯 Quick Start

### Option 1: Run in RStudio

1. Open `app.R`
2. Click **Run App**

### Option 2: R Console

```r
shiny::runApp("your_app_folder_path")
```

---

## 📦 Demo Files Included

| File Name | Description |
|-----------|-------------|
| `small_demo_seurat.rds` | Lung cancer (GEO) |
| `small_demo_seurat_breast.rds` | Breast cancer |
| `small_demo_seurat_pancreas.rds` | Pancreatic cancer |
| `small_demo_seurat_brain.rds` | Brain tumor (10x) |
| `small_10x_demo_seurat_lung.rds` | Lung cancer (10x) |

> ⚠️ These are 200-cell subsamples to meet GitHub limits. Full datasets used for testing are available on request.

---

## 🔗 External Files

Some files are too large for GitHub and must be downloaded manually:

- 6 Pathway `.rds` files
- GSVA `.rds` files
- `gencode.v44.annotation.gtf`
- `TF_target.rds`

👉 [**Download from Google Drive**](https://drive.google.com/drive/folders/1sIHNTeEyzXJc9gMayGmIOxeev1wwSOp6?usp=sharing)

Place them into the app’s `Data/` folder after download.

---

## 📄 License

MIT License

---

## 📬 Contact

**Developer:** Haider  
📧 haider@emails.bjut.edu.cn

---

## 📚 Citation

If LncScape contributes to your research, please cite this tool in your publication.
