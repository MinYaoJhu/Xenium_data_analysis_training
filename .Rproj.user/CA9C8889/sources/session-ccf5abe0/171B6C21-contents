
# 🧭 Week 1 – Xenium Data Overview and Loading

**Learning Goals**
- Understand Xenium output files
- Load raw Xenium data into R using `LoadXenium()`
- Explore metadata and perform basic QC

## 🧬 Xenium Output Structure
- `cells.csv.gz` — per-cell counts & metadata
- `transcripts.csv.gz` — individual transcript coordinates
- `morphology.ome.tif` — morphology image
- `metadata.json` — experiment metadata

## 🧰 Load Data in R
```r
library(Seurat)
xen <- LoadXenium(path = "data/demo_xenium/")  # update this path to your dataset
xen
head(xen@meta.data)
```

## 🧪 Practice
1. Count cells/genes in the dataset.  
2. Plot total transcripts per cell with `VlnPlot(xen, "nCount_Spatial")`.  
3. Summarize file roles in a short paragraph.
