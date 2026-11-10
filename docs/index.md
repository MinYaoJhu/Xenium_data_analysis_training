---
title: Home
nav_order: 1
---

<img src="{{ site.baseurl }}/assets/logo.svg" alt="Course logo" width="140" align="right"/>

# 🌿 Xenium Data Analysis Training

Welcome to the **Xenium Data Analysis Training Course**!  
This 4-week program introduces the fundamentals of analyzing **10x Genomics Xenium in situ** transcriptomics data using **R** and **Xenium Explorer**.

The course is designed for new researchers and PhD students to gain hands-on experience with Xenium data — from loading raw output to visualizing gene expression, clustering cells, and identifying marker genes.

---

## 📘 Course Overview

**Instructor:** Dr. Min-Yao Jhu  
**Duration:** 4 weeks (1-hour guided session + ~3 days of independent practice per week)  
**Tools:** R (Seurat v5, Xenium Ranger) and Xenium Explorer  
**Prerequisites:** Basic R and single-cell analysis experience  

---

## 🎯 Learning Objectives

By the end of the course, you will be able to:

1. Understand the structure and components of Xenium output files  
2. Load Xenium data into R using Seurat functions  
3. Visualize spatial gene expression using R and Xenium Explorer  
4. Perform clustering and dimensionality reduction  
5. Detect and interpret marker genes for specific spatial domains  

---

## 🗓️ Weekly Syllabus

### **Week 1 – Data Overview & Loading**
Focus: Learn Xenium file structure 
➡️ [Go to Week 1 Lesson]({{ site.baseurl }}/week1.html)

References:

- ➡️ [Xenium Output Files at a Glance](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-at-a-glance)

- ➡️ [Understanding Xenium Outputs](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs)


---

### **Week 2 – Visualization & Exploration**
Focus: explore interactively in Xenium Explorer  
➡️ [Go to Week 2 Lesson]({{ site.baseurl }}/week2.html)

References:

- ➡️ [Xenium Explorer Input Files](https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/inputs)

- ➡️ [Navigating Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-getting-started-with-xenium-explorer)

---

### **Week 3 – Clustering & Dimension Reduction**
Focus: how to import data into R using `LoadXenium()`, Normalize data (`SCTransform()`), perform PCA/UMAP, and identify spatial clusters  
➡️ [Go to Week 3 Lesson]({{ site.baseurl }}/week3.html)

---

### **Week 4 – Marker Detection & Interpretation**
Focus: Identify cluster-specific marker genes with `FindAllMarkers()` and visualize them using `FeaturePlot()` and `DotPlot()`  
➡️ [Go to Week 4 Lesson]({{ site.baseurl }}/week4.html)


---

## 🧰 Setup

See **[Setup & data]({{ site.baseurl }}/setup.html)** for the course dataset, the tested software environment and how to render the notebooks.

## 📚 Recommended Resources

- [Xenium Onboard Analysis documentation](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest)
- [Analysis of Image-based Spatial Data in Seurat](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2) (Satija Lab), on which Weeks 3 and 4 are based


