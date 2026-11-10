"""Build the GitHub Pages site in docs/ from the course materials.

The week folders are the single source of truth. This script copies each
lesson.md into docs/ with Jekyll front matter, copies the rendered Week 3
and Week 4 notebooks (self-contained HTML) and the logo, and leaves the
hand-written pages (index.md, setup.md, _config.yml) untouched.

Usage (from the repository root): python3 scripts/build_site.py
"""
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
REPO = "https://github.com/MinYaoJhu/Xenium_data_analysis_training"

WEEKS = [
    # folder, page, nav title, notebook (rendered HTML) or None
    ("week1_data_overview", "week1.md", "Week 1 · Data overview", None),
    ("week2_visualization", "week2.md", "Week 2 · Visualization", None),
    ("week3_clustering", "week3.md", "Week 3 · Clustering",
     "Min-Yao_edited_seurat5_spatial_vignette"),
    ("week4_markers", "week4.md", "Week 4 · Marker genes",
     "Min-Yao_edited_seurat5_spatial_vignette-2"),
]


def main():
    (DOCS / "notebooks").mkdir(parents=True, exist_ok=True)
    (DOCS / "assets").mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "assets" / "logo.svg", DOCS / "assets" / "logo.svg")

    for order, (folder, page, title, notebook) in enumerate(WEEKS, start=2):
        lesson = (ROOT / folder / "lesson.md").read_text().strip()
        parts = [f"---\ntitle: \"{title}\"\nnav_order: {order}\n---\n", lesson]
        if notebook:
            html = f"{page[:-3]}_notebook.html"
            shutil.copy(ROOT / folder / f"{notebook}.html", DOCS / "notebooks" / html)
            parts.append(
                "\n\n---\n\n## 📓 Full worked notebook\n\n"
                f"- [Open the rendered notebook]({{{{ site.baseurl }}}}/notebooks/{html}) "
                "(code, output and figures)\n"
                f"- [R Markdown source]({REPO}/blob/main/{folder}/{notebook}.Rmd)\n"
            )
        parts.append(f"\n\n<sub>Source: [`{folder}/lesson.md`]({REPO}/blob/main/{folder}/lesson.md)</sub>\n")
        (DOCS / page).write_text("".join(parts))
        print("wrote docs/" + page)


if __name__ == "__main__":
    main()
