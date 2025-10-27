
# 🧭 Week 1 – Xenium Data Overview and Loading

**Learning Goals**
- Understand Xenium output files
- Understand `analysis_summary.html`
- Explore metadata and perform basic QC

---

## 1) Xenium Output Files at a Glance

The **Xenium Onboard Analysis (XOA)** pipeline generates a rich set of decoded transcript data and morphology images for each selected region.  
These outputs provide the foundation for spatial analysis in R and visualization in **Xenium Explorer**.

---

### 🧬 Overview

Raw Xenium output bundles contain:
- **Decoded transcript coordinates and counts**
- **Morphology images** (DAPI / cell-segmentation channels)
- **Cell and nucleus segmentation polygons**
- **Derived matrices and QC summaries**

Low-level image sensor data are compressed into these standardized files so that essential spatial and quality information is preserved.

---

### Where your data live (folder layout)

Xenium run folders follow this pattern (one subfolder per selected region):

```
output/
└── <yyyymmdd><hhmmss><runName>/
└── output-<instrumentSN><slideID><regionName><yyyymmdd><hhmmss>/
├── experiment.xenium
├── analysis_summary.html
├── morphology.ome.tif
├── morphology_focus/
├── cells.csv.gz / cells.parquet
├── cells.zarr.zip
├── cell_boundaries.* / nucleus_boundaries.*
├── transcripts.parquet / transcripts.zarr.zip
├── cell_feature_matrix/ (MEX)
├── cell_feature_matrix.h5
├── cell_feature_matrix.zarr.zip
├── metrics_summary.csv
├── analysis/ (secondary analysis results)
├── analysis.zarr.zip
├── gene_panel.json
├── protein_panel.json (if protein assay used)
└── aux_outputs/
```

---

### 📁 Output Bundle Structure

Each Xenium experiment directory contains one folder per region of interest (ROI).  
Inside each ROI folder you’ll find the following groups of files:

| Category | File(s) | What it provides | Explorer readable? |
|---|---|---|---|
| **Manifest** | `experiment.xenium` | JSON with experiment metadata + relative paths to all files. | ✅ (Required) |
| **QC dashboard** | `analysis_summary.html` | Interactive run summary, decoding + segmentation + analysis tabs. | Openable via Explorer “Sample Information” |
| **Morphology images** | `morphology.ome.tif` | 3D DAPI Z-stack (pyramidal OME-TIFF, 16-bit). | ✅ |
|  | `morphology_focus/` | 2D autofocus projections for DAPI and, if used, additional stains (channel-named files in v4.0). | ✅ |
| **Cells (summary)** | `cells.csv.gz`, `cells.parquet` | Per-cell QC: counts, areas, centroids, segmentation method. | ✅ |
| **Masks & polygons** | `cells.zarr.zip` | Nucleus + cell **masks** (used for transcript→cell assignment). | ✅ |
|  | `cell_boundaries.*`, `nucleus_boundaries.*` | Polygon **approximations** of masks for fast visualization. | ✅ |
| **Transcripts** | `transcripts.parquet`, `transcripts.zarr.zip` | One row per decoded transcript (x,y,z, gene, qv, cell_id, etc.). | ✅ |
| **Matrices** | `cell_feature_matrix/` (MEX) | Sparse matrix (features × cells) + `features.tsv.gz`, `barcodes.tsv.gz`. | ✅ |
|  | `cell_feature_matrix.h5` | Same matrix in HDF5 (fast loading). | — |
|  | `cell_feature_matrix.zarr.zip` | Same matrix in Zarr; Explorer can read. | ✅ |
| **Metrics** | `metrics_summary.csv` | Key run metrics (decoding + segmentation). | — |
| **Secondary analysis** | `analysis/`, `analysis.zarr.zip` | PCA/UMAP, clustering (graph + K-means), differential expression. | `analysis.zarr.zip` is ✅ |
| **Panels** | `gene_panel.json`, `protein_panel.json` | Targets, codewords, probe sets; assay meta. | — |
| **Auxiliary** | `aux_outputs/` | FOV maps, per-cycle RNA images, overview scan, background images, QC masks. | Some importable (e.g., focus QC masks) |

---

### 🧩 Auxiliary Output Files (`aux_outputs/`)

Additional QC and image files supporting troubleshooting and alignment.

| **File / Folder** | **Description** | **Readable in Xenium Explorer?** |
|--------------------|-----------------|----------------------------------|
| `morphology_fov_locations.json` | Field-of-view positions in microns. | — |
| `overview_scan_fov_locations.json` | Field-of-view positions in pixels. | — |
| `per_cycle_channel_images/` | Down-sampled RNA images (TIFF) for each cycle × channel. | — |
| `overview_scan.png` | Whole-slide overview image. | — |
| `background_qc_images/` | Autofluorescence background images used for correction. | — |
| `morphology_focus_qc_masks/` | 8-bit saturation QC masks; 255 = saturated pixel, 0 = not. | ✅ (importable) |

---

## 2) Read the QC: `analysis_summary.html` (5 tabs to know)

Open the HTML (on the instrument, in a browser, or from Explorer → **Sample Information → Analysis**):

- **Summary** — high-level run info (panel, software versions), quick metrics, and overview images.  
- **Decoding** — Q-score distributions, detection rates, control probes/codewords diagnostics.  
- **Cell Segmentation** — cell/nucleus counts, size distributions, transcript partitioning.  
- **Analysis** — if run: normalization, PCA/UMAP, graph-based and K-means clustering, top features.  
- **Image QC** — galleries: downsampled **RNA cycle/channel** images, morphology stains, autofluorescence/background images.


> If secondary analysis wasn’t generated (older runs or disabled setting), this page may be minimal. You can compute clustering later in R.

---

## 3) Morphology images (orientation, viewing)

- All morphology images use the **same image coordinate system** (origin at top-left).  
- `morphology.ome.tif` is a **3D DAPI stack** (useful for segmentation QA or resegmentation).  
- `morphology_focus/` contains **2D projections** per channel (v4.0 files are channel/marker-named).  
- View in **Xenium Explorer**; multi-file focus images can also be viewed in **QuPath**, **Napari**, or **Fiji/ImageJ** (open any one file; metadata links the set).

---

## 4) Cells & transcripts (what’s inside the tables)

**`cells.*` columns (core subset):**
- `cell_id`, `x_centroid`, `y_centroid`, `cell_area`, `nucleus_area`, `nucleus_count`  
- `transcript_counts` (Q≥20), control counts (neg. probes/codewords, genomic if Prime)  
- `total_counts` (sum of above), `segmentation_method`

**`nucleus/cell_boundaries.*` (polygons):**
- `cell_id`, `vertex_x`, `vertex_y`, `label_id`  
> Polygons are simplified outlines for visualization; **masks** live in `cells.zarr.zip`.

**`transcripts.parquet` columns (core subset):**
- `transcript_id`, `feature_name` (gene/control), `cell_id` (if assigned), `x_location`, `y_location`, `z_location`, `qv`  
- `overlaps_nucleus` (0/1), `fov_name`, `nucleus_distance`, `codeword_index`, `codeword_category`, `is_gene`

**Cell-feature matrix (three formats):**
- **MEX** (`cell_feature_matrix/`), **HDF5** (`.h5`), **Zarr** (`.zarr.zip`)  
- Features include: **genes**, **negative controls**, **unassigned/deprecated codewords**, and (if used) **protein features** (v4.0).

---

## 5) Practice 

**Open** `analysis_summary.html` and write down 2–3 QC observations (e.g., median transcripts per cell, top markers, segmentation quality).

---

### Next (Week 2)

- Open and navigate in **Xenium Explorer**

---

## 🧠 Practical Tips

- **Start with `analysis_summary.html`** to QC decoding and cell segmentation before diving into R.  
- The **`cells.zarr.zip`** file contains segmentation masks — you generally won’t open it manually but it’s crucial for spatial assignments.  
- Keep the **folder hierarchy intact** when moving data; Xenium Explorer expects relative paths as generated by the pipeline.  
- Archive large morphology images (`.ome.tif`) on shared storage; they can be > 10 GB per ROI.


