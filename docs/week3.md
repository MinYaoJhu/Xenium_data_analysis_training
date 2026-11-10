---
title: "Week 3 · Clustering"
nav_order: 4
---
# 🔬 Week 3 – Clustering and Dimension Reduction

**Learning Goals**
- Normalize/scale data
- Run PCA/UMAP
- Identify clusters and visualize them spatially

## ⚙️ Preprocessing and Clustering
```r
xen <- SCTransform(xen)
xen <- RunPCA(xen)
xen <- RunUMAP(xen, dims = 1:20)
xen <- FindNeighbors(xen, dims = 1:20)
xen <- FindClusters(xen, resolution = 0.5)
```

## 📊 Visualize Clusters
```r
DimPlot(xen, reduction = "umap", label = TRUE)                       # clusters on the UMAP
ImageDimPlot(xen, group.by = "seurat_clusters", size = 1.5)          # clusters on the tissue
```

> `ImageDimPlot()` is for imaging-based data such as Xenium. `SpatialDimPlot()` is for spot-based Visium data and does not work with Xenium objects.

## 🧪 Practice
- Inspect PCs (`VizDimLoadings`, `DimPlot`).  
- Try resolutions 0.3 vs 1.0 and compare results.  
- Relate clusters to visible tissue regions.

---

## 📓 Full worked notebook

- [Open the rendered notebook]({{ site.baseurl }}/notebooks/week3_notebook.html) (code, output and figures)
- [R Markdown source](https://github.com/MinYaoJhu/Xenium_data_analysis_training/blob/main/week3_clustering/Min-Yao_edited_seurat5_spatial_vignette.Rmd)


<sub>Source: [`week3_clustering/lesson.md`](https://github.com/MinYaoJhu/Xenium_data_analysis_training/blob/main/week3_clustering/lesson.md)</sub>
