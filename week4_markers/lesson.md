
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
```r
FeaturePlot(xen, features = c("NIN"))
DotPlot(xen, features = c("ENOD11", "NF-YA1")) + RotatedAxis()
```

## 🧪 Practice
1. Save top 5 markers per cluster to CSV.  
2. Visualize two markers on the morphology image.  
3. Write a short interpretation (3–5 sentences).
