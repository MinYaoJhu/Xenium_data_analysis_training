
# 🧭 Week 1 – Xenium Data Overview and Loading

**Learning Goals**
- Understand Xenium output files
- Understand `analysis_summary.html`
- Explore metadata and perform basic QC
- Open and navigate in **Xenium Explorer**

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

### 📁 Output Bundle Structure

Each Xenium experiment directory contains one folder per region of interest (ROI).  
Inside each ROI folder you’ll find the following groups of files:

| **Category** | **File(s)** | **Description** | **Readable in Xenium Explorer?** |
|---------------|-------------|-----------------|----------------------------------|
| **Experiment manifest** | `experiment.xenium` | Defines experiment metadata and links all outputs. | ✅ Required |
| **Interactive summary** | `analysis_summary.html` | HTML dashboard with key QC metrics, plots, and images. | — |
| **Morphology images** | `morphology.ome.tif` | 3-D DAPI image stack in OME-TIFF format. | ✅ |
|  | `morphology_focus/` | Multi-focus projections (2-D) for each stain; includes DAPI + additional segmentation channels. | ✅ (for multimodal assay) |
| **Cell summary** | `cells.csv.gz`, `cells.parquet` | Per-cell metrics: transcript counts, area, morphology stats. | ✅ |
| **Cell segmentation masks** | `cells.zarr.zip` | Zarr archive containing nucleus + cell masks used for transcript assignment. | ✅ |
| **Cell boundaries** | `cell_boundaries.csv.gz`, `cell_boundaries.parquet` | Polygon coordinates outlining each cell. | ✅ |
| **Nucleus boundaries** | `nucleus_boundaries.csv.gz`, `nucleus_boundaries.parquet` | Polygon coordinates outlining each nucleus. | ✅ |
| **Transcript data** | `transcripts.parquet`, `transcripts.zarr.zip` | Coordinates and gene identity for every detected transcript. | ✅ |
| **Cell-feature matrix** | `cell_feature_matrix/`, `cell_feature_matrix.h5`, `cell_feature_matrix.zarr.zip` | Sparse matrix of gene × cell counts in multiple formats (MTX / HDF5 / Zarr). | ✅ |
| **Metric summary** | `metrics_summary.csv` | Global summary of sequencing and decoding metrics. | — |
| **Secondary analysis** | `analysis/`, `analysis.zarr.zip` | Processed data and derived results (e.g., clustering). | — |
| **Panels** | `gene_panel.json`, `protein_panel.json` | Gene and (if applicable) protein panel definitions. | — |

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

## 2) Read the QC: `analysis_summary.html`

Open the HTML file in a web browser. Use it to:
- **Verify run metadata:** panel name, panel design ID, software version
- **Check detection metrics:** transcripts per cell, % cells with transcripts, top genes
- **Inspect segmentation quality:** cell counts, size distributions, boundary overlays (if present)
- **Spot issues quickly:** unusually low transcripts/cell, poor segmentation, or missing files

> If secondary analysis wasn’t generated (older runs or disabled setting), this page may be minimal. You can compute clustering later in R.

---

## 3) Open a Dataset in Xenium Explorer

**Start the app** → choose one of:
- **Drag & drop** the `experiment.xenium` file
- Click **Open New File** and select the `.xenium` file
- **Open file from path** (paste the full path to `experiment.xenium`)

> ⚠️ The **XOA output files must be in the same directory** as referenced in the manifest (unless you’ve edited paths inside `experiment.xenium`).

**Recent files**
- The Home page lists up to 100 recent datasets (v4.0+).  
- Closing a dataset auto-saves the current view as an **unnamed saved view** only for the Recent Files entry.  
- **Save your view explicitly** if you want a persistent named view (recommended).

---


## 4) The Explorer Interface (What You’ll Use Most)

When the dataset loads, you’ll see the DAPI (nuclei) image by default.

### A. Sample Information
- Click **Sample Information** to see **Run**, **Panel**, and **Analysis** metadata.
- You can open the QC summary (`analysis_summary.html`) from here in a separate window.

### B. Layers (Images / Cells / Transcripts)
Enable/disable these from the side panel:
- **Images**: Adjust brightness/contrast; select z-planes or focus projections; import & align images if needed.
- **Cells**: Show nuclei/cell boundaries (filled, outlined, or both).
- **Transcripts**: Toggle gene(s); points represent decoded transcripts.

> 💡 Combine layers: e.g. **Cells + Transcripts** to see per-cell transcript distributions; or **Images + Transcripts** to compare morphology context.

### C. Lasso Tools (ROI selection & export)
- **Rectangle (R)** or **Freehand (L)** to select a region.
- Selection pane shows counts for **transcripts** (by selected genes) and **cells** in the ROI.  
- Very large selections may hide transcript counts for performance—zoom in or make smaller ROIs.
- You can **export ROI** info (see “Export Regions of Interest” from the help menu).

### D. Go to Location
- Jump to coordinates with **Go to location** (format: `x,y` in µm).
- Paste multiple ROI coordinates using the **Lasso → Paste coordinates** feature.

---

## 5) Export Publication-Quality Images

Click the **download** button (top-right):

- **Quick Export Image**  
  Exports what you see in the viewport as a PNG (fast snapshots).

- **High Resolution Image Export** (v3.2+)  
  Exports a stitched high-res PNG matching the viewport aspect ratio.  
  - When zoomed in, Explorer fetches tiles at **max zoom** to maximize detail.  
  - When zoomed out, it renders at ~32× current level to preserve content.  
  - **Tip:** Increase *Transcript point size scale* before export so points remain visible.

> UI menus are not included in exports. Use **Settings** to show/hide scale axes or the picture-in-picture navigator.

---

## 6) Advanced (Optional)

### Edit the manifest to use custom analyses
- `experiment.xenium` is JSON. You can edit paths under `images` or `xenium_explorer_files` to point to **custom outputs** (e.g., updated boundaries or matrices).
- Paths can be **absolute or relative**. Avoid `~` for home expansion in relative paths.

### Large dataset feature limits (older versions)
- v3.0–3.1 had limits for extremely large matrices (>200M features).  
- v3.2+ restores most functions for XOA ≥ v3.1 datasets.

---

## 7) Practice (10–20 min)

1. **Open** a dataset by dragging `experiment.xenium` into Explorer.  
2. **Find** and open **Sample Information**. Note the panel name, panel design ID, and software versions.  
3. **Toggle layers**: show Images + Cells + Transcripts; add 1–2 gene layers.  
4. **Use lasso** (freehand and rectangle) to select an ROI and **inspect transcript counts**.  
5. **Export** a **Quick** image and a **High-Res** image (include scale axes).  
6. **Open** `analysis_summary.html` and write down 2–3 QC observations (e.g., median transcripts per cell, top markers, segmentation quality).

---

## 8) Troubleshooting

- **Nothing shows after opening `.xenium`:** check that all referenced files still exist and that folder paths in `experiment.xenium` are valid.
- **No transcripts visible:** ensure **Transcripts** layer is enabled and at least one gene is selected; zoom to a region with tissue.
- **Slow exports or missing points in high-res export:** increase the **Transcript point size scale**; be patient while tiles render.
- **Empty `/analysis` or missing clustering:** it may not have been run onboard; you can compute clustering later (Week 3) or re-run **Xenium Ranger** to populate `/analysis`.

---

### Next (Week 2)
We’ll start **visualizing genes** programmatically (R), reproduce key Explorer views, and standardize figure exports.

---

## 🧠 Practical Tips

- **Start with `analysis_summary.html`** to QC decoding and cell segmentation before diving into R.  
- The **`cells.zarr.zip`** file contains segmentation masks — you generally won’t open it manually but it’s crucial for spatial assignments.  
- Keep the **folder hierarchy intact** when moving data; Xenium Explorer expects relative paths as generated by the pipeline.  
- Archive large morphology images (`.ome.tif`) on shared storage; they can be > 10 GB per ROI.


