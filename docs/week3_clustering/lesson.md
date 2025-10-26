
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
SpatialDimPlot(xen, group.by = "seurat_clusters")
```

## 🧪 Practice
- Inspect PCs (`VizDimLoadings`, `DimPlot`).  
- Try resolutions 0.3 vs 1.0 and compare results.  
- Relate clusters to visible tissue regions.
