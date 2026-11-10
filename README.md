# Xenium Data Analysis Training <img src="assets/logo.svg" align="right" width="140" alt="Spatial training course logo"/>

**A 4-week course introducing new researchers to 10x Genomics Xenium *in situ* spatial transcriptomics analysis with R and Xenium Explorer**

[![Course website](https://img.shields.io/badge/Course%20website-minyaojhu.github.io-1A1A2E)](https://minyaojhu.github.io/Xenium_data_analysis_training/)
![R](https://img.shields.io/badge/R-4.4.3-276DC3?logo=r&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-5.5.1-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

👉 **Start here: [course website](https://minyaojhu.github.io/Xenium_data_analysis_training/)**, with one page per week, a setup guide and the fully rendered notebooks.

The course builds confidence in handling raw Xenium output, visualizing spatial gene expression, clustering cells and identifying marker genes, using a reproducible R workflow and a public 10x Genomics dataset.

## 📘 Overview

**Instructor:** Dr. Min-Yao Jhu
**Format:** 4 weeks (1-hour guided session + ~3 days of practice per week)
**Tools:** R (Seurat v5), Xenium Explorer
**Level:** beginner to intermediate (some single-cell analysis in R is helpful)
**Data:** public 10x Genomics Fresh Frozen Mouse Brain, Tiny subset (see [Setup & data](https://minyaojhu.github.io/Xenium_data_analysis_training/setup.html))

### 🧩 Learning objectives

By the end of the course, you will be able to:
1. Understand the structure and content of Xenium output files
2. Load and explore Xenium datasets in R
3. Visualize transcript expression spatially using R and Xenium Explorer
4. Perform PCA, clustering and identify marker genes
5. Interpret spatial patterns and document analyses reproducibly

## 🧠 Weekly outline

| Week | Topic | Materials |
|---|---|---|
| 1 | **Data overview & loading:** Xenium output files, `analysis_summary.html`, metadata and basic QC | [lesson](https://minyaojhu.github.io/Xenium_data_analysis_training/week1.html) |
| 2 | **Visualization & exploration** in Xenium Explorer (basic and advanced), exporting figures | [lesson](https://minyaojhu.github.io/Xenium_data_analysis_training/week2.html) |
| 3 | **Clustering & dimension reduction:** `LoadXenium()`, QC, `SCTransform()`, PCA, UMAP, clusters on the tissue with `ImageDimPlot()` | [lesson](https://minyaojhu.github.io/Xenium_data_analysis_training/week3.html) · [notebook](https://minyaojhu.github.io/Xenium_data_analysis_training/notebooks/week3_notebook.html) |
| 4 | **Marker genes & interpretation:** `FindAllMarkers()`, `FeaturePlot()`, `DotPlot()`, `ImageFeaturePlot()` | [lesson](https://minyaojhu.github.io/Xenium_data_analysis_training/week4.html) · [notebook](https://minyaojhu.github.io/Xenium_data_analysis_training/notebooks/week4_notebook.html) |

## 🗂️ Repository structure

```
Xenium_data_analysis_training/
├── week1_data_overview/lesson.md
├── week2_visualization/lesson.md
├── week3_clustering/         # lesson.md + R Markdown notebook (rendered .html/.md and figures)
├── week4_markers/            # lesson.md + R Markdown notebook + marker gene table
├── docs/                     # course website (GitHub Pages, Just the Docs theme)
├── scripts/build_site.py     # builds the website pages from each week's lesson.md
├── assets/                   # course logo and its generator
├── environment.yml           # tested conda environment
├── session_info.txt          # full R session used to render the notebooks
├── CITATION.cff
└── LICENSE
```

## 🧰 Setup

```bash
micromamba create -f environment.yml     # or: conda env create -f environment.yml
micromamba activate xenium-training
```

Download the course data and render the notebooks as described on the [Setup & data](https://minyaojhu.github.io/Xenium_data_analysis_training/setup.html) page. Week 3 must be run before Week 4.

## ✅ Tested versions

The Week 3 and Week 4 notebooks were re-rendered in November 2026 with **R 4.4.3** and **Seurat 5.5.1** (SeuratObject 5.4.0) and run without errors. If code from an older Seurat tutorial fails, check the version first: for example, `GetAssayData()` now takes `layer =` instead of `slot =`. Full details: [`session_info.txt`](session_info.txt).

## 🛠️ Updating the website

Edit the lessons in the week folders, then run:

```bash
python3 scripts/build_site.py
```

GitHub Pages rebuilds the website from `docs/` after each push.

## 📚 Credits and resources

- Weeks 3 and 4 are adapted from the Satija Lab vignette [Analysis of Image-based Spatial Data in Seurat](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2) (Seurat, MIT License).
- Course data: [Fresh Frozen Mouse Brain for Xenium Explorer Demo](https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard), 10x Genomics.
- [Xenium Onboard Analysis documentation](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest)

## Related repositories

- [MedicagoXeniumSpatialViewer](https://github.com/MinYaoJhu/MedicagoXeniumSpatialViewer): interactive Shiny app for exploring Xenium data
- [Medicago_Xenium_relabeling_resegmentation](https://github.com/MinYaoJhu/Medicago_Xenium_relabeling_resegmentation): Xenium Ranger relabeling and H&E-based cell resegmentation pipelines

## Citation

If you use these materials, please cite them with GitHub's **"Cite this repository"** button (from [`CITATION.cff`](CITATION.cff)).

## AI assistance

Updates to the code, website, figures and documentation in this repository were made with the assistance of Claude Opus 5.5 (Anthropic), used through Claude Code. The course materials were written by Min-Yao Jhu. All changes were directed, reviewed and approved by the author, who takes full responsibility for the content.

## Contact

**Min-Yao Jhu**, Senior Research Scientist, Donald Danforth Plant Science Center
Email: [MJhu@danforthcenter.org](mailto:MJhu@danforthcenter.org) · Website: [minyaojhu.com](https://www.minyaojhu.com)

© 2025–2026 Min-Yao Jhu · course first developed at the Crop Science Centre, University of Cambridge · MIT License
