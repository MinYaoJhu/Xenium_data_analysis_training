# Xenium_data_analysis_training

Welcome to the **Xenium Data Analysis Training Repository**!  
This repository provides a structured 4-week program designed to introduce new researchers to **Xenium in situ** transcriptomics data analysis using **R** and **Xenium Explorer**.  

The goal is to build confidence in handling raw Xenium output, visualizing spatial gene expression, performing clustering, and identifying marker genes — all within an intuitive, reproducible R workflow.

---

## 📘 Overview

**Instructor:** Dr. Min-Yao Jhu  
**Format:** 4 weeks (1 hour guided session + ~3 days practice per week)  
**Tools:** R (Seurat v5, Xenium Ranger), Xenium Explorer  
**Level:** Beginner to intermediate (familiarity with single-cell analysis in R is helpful)

### 🧩 Learning Objectives

By the end of this training, you will be able to:
1. Understand the structure and content of Xenium output files  
2. Load and explore Xenium datasets in R  
3. Visualize transcript expression spatially using R and Xenium Explorer  
4. Perform PCA, clustering, and identify marker genes  
5. Interpret spatial patterns and document analyses reproducibly

---

## 🗂️ Repository Structure

```
xenium-training/
├── README.md ← overview and syllabus
├── week1_data_overview/ ← raw data loading & QC
├── week2_visualization/ ← visualization & Xenium Explorer
├── week3_clustering/ ← clustering & dimension reduction
├── week4_markers/ ← marker gene detection & interpretation
└── docs/ ← optional folder for GitHub Pages website
```

---

## 🧠 Weekly Outline

### **Week 1 — Data Overview & Loading**
- Understand Xenium output folder structure (`cells.csv.gz`, `transcripts.csv.gz`, `metadata.json`, etc.)  
- Explore metadata and perform basic quality checks  
- **Start with `analysis_summary.html`** to QC decoding and cell segmentation before diving into R.

### **Week 2 — Visualization & Exploration**
- Explore data interactively using **Xenium Explorer** (Basic)
- Explore data interactively using **Xenium Explorer** (Advanced)
- Export high-quality figures  

### **Week 3 — Clustering & Dimension Reduction**
- Load data using `LoadXenium()` in R  
- Visualize gene expression
- Normalize data with `SCTransform()`  
- Perform PCA and UMAP embedding  
- Identify clusters and visualize them spatially  

### **Week 4 — Marker Gene Detection & Interpretation**
- Use `FindAllMarkers()` to identify cluster-specific genes  
- Visualize marker genes with `FeaturePlot()` and `DotPlot()`  
- Interpret biological relevance of spatial clusters  

---

## 🧰 Installation & Setup

1. **Install R (≥4.3) and RStudio**
   - [Download R](https://cran.r-project.org/)  
   - [Download RStudio](https://posit.co/download/rstudio-desktop/)

2. **Install Required Packages**
   ```r
   install.packages("Seurat")
   install.packages("tidyverse")
   install.packages("patchwork")
   ```

---

## 📚 Recommended Resources

https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest

---

## GitHub website link: https://minyaojhu.github.io/Xenium_data_analysis_training/

---

© 2025–2026 Min-Yao Jhu
Donald Danforth Plant Science Center (course first developed at the Crop Science Centre, University of Cambridge)
For educational and research training use.
   