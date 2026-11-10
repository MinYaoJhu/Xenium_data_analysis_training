---
title: "Week 4 · Marker genes"
nav_order: 5
---
# 🌱 Week 4 – Marker Gene Detection and Interpretation

**Learning Goals**
- Identify cluster markers
- Visualize marker expression
- Interpret biological meaning

## 🧭 Identify Markers
```r
markers <- FindAllMarkers(xen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
```

## 🎨 Visualize Markers
Well-known mouse brain markers in this panel: `Cux2` (upper cortical layers), `Bcl11b` (deep cortical layers), `Gad1` (inhibitory neurons) and `Gfap` (astrocytes).
```r
FeaturePlot(xen, features = c("Cux2", "Bcl11b", "Gad1", "Gfap"))           # on the UMAP
DotPlot(xen, features = c("Cux2", "Bcl11b", "Gad1", "Gfap")) + RotatedAxis()
ImageFeaturePlot(xen, features = c("Cux2", "Gad1"), max.cutoff = "q90")   # on the tissue
```

## 🧪 Practice
1. Save top 5 markers per cluster to CSV.  
2. Visualize two markers on the morphology image.  
3. Write a short interpretation (3–5 sentences).

---

## 📓 Full worked notebook

- [Open the rendered notebook]({{ site.baseurl }}/notebooks/week4_notebook.html) (code, output and figures)
- [R Markdown source](https://github.com/MinYaoJhu/Xenium_data_analysis_training/blob/main/week4_markers/Min-Yao_edited_seurat5_spatial_vignette-2.Rmd)


<sub>Source: [`week4_markers/lesson.md`](https://github.com/MinYaoJhu/Xenium_data_analysis_training/blob/main/week4_markers/lesson.md)</sub>
