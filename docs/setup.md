---
title: "Setup & data"
nav_order: 6
---

# 🧰 Setup & data

## 1) Course dataset

The course uses the public 10x Genomics **Fresh Frozen Mouse Brain, Tiny subset** (Xenium Onboard Analysis 1.0.2), from the [Fresh Frozen Mouse Brain for Xenium Explorer Demo](https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard) dataset.

```bash
cd week3_clustering
curl -O https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
md5sum Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip   # a39fa6d0a751db1f206c915b6419e329
unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip -d Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs
```

The zip file is 3.5 GB (3,479,544,713 bytes; checked on 10 November 2026). The data folder is ignored by Git.

## 2) Software

**Option A: tested environment (recommended)**

```bash
micromamba create -f environment.yml     # or: conda env create -f environment.yml
micromamba activate xenium-training
```

This installs R 4.4.3 and Seurat 5.5.1, the versions used to render the course notebooks (full details in [`session_info.txt`](https://github.com/MinYaoJhu/Xenium_data_analysis_training/blob/main/session_info.txt)).

**Option B: your own R (≥ 4.3) and RStudio**

```r
install.packages(c("Seurat", "tidyverse", "patchwork", "future", "arrow", "hdf5r"))
```

[Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer/downloads) is needed for Week 2.

## 3) Render the notebooks

Run Week 3 first: it saves the processed Seurat object that Week 4 loads.

```bash
cd week3_clustering && Rscript -e 'rmarkdown::render("Min-Yao_edited_seurat5_spatial_vignette.Rmd")'
cd ../week4_markers   && Rscript -e 'rmarkdown::render("Min-Yao_edited_seurat5_spatial_vignette-2.Rmd")'
```

Each notebook takes about one minute.
