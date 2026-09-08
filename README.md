# Gerasymchuk2026 — single-cell RNA-seq analysis

Analysis code for the single-cell RNA-seq panels of **Figure 4** in:

> Gerasymchuk *et al.* (2026), *Expression of genes involved in fatty acid
> metabolism in control and FXS astrocytes and their functional subtypes.*

Human iPSC-derived forebrain astrocytes from three lines — **PO2** (control),
**HEL69-5** (FXS) and **HEL100-1** (FXS with epilepsy) — profiled with 10x
Genomics 3' scRNA-seq and aggregated with `cellranger aggr`.

## Contents

| File | Purpose |
|---|---|
| `01_analysis.R` | Full pipeline: CellRanger matrix → QC → normalisation → doublet removal → Harmony → Leiden clustering → gene modules → functional enrichment. Writes the tables deposited on Zenodo. |
| `02_figures.R` | Regenerates the Figure 4 panels **from the deposited tables only**. Does not re-cluster. |
| `functions.R` | Helper functions used by both scripts. |
| `sessionInfo_validation.txt` | R and package versions `02_figures.R` was last verified against (R 4.4.3, Bioconductor 3.20). |

## Data

Neither script contains data. Both read from a `data/` directory that you
populate from the two public deposits:

| Source | Accession | What |
|---|---|---|
| ENA | [PRJEB110867](https://www.ebi.ac.uk/ena/browser/view/PRJEB110867) | Raw FASTQ reads (samples ERS29600475–ERS29600477) |
| Zenodo | [10.5281/zenodo.19322443](https://doi.org/10.5281/zenodo.19322443) | Count matrices, cell/gene metadata, cluster labels, TriMap coordinates, marker and enrichment tables |

```bash
mkdir -p data
# download and unpack the Zenodo record into data/
```

After unpacking, `data/` should contain `counts_filtered/`,
`counts_raw_cellranger/`, `cell_metadata.tsv.gz`, `gene_metadata.tsv.gz`,
`cluster_annotation.tsv`, `cell_line_annotation.tsv` and
`fig4A_gene_selection.tsv`. The Zenodo `README.md` documents every file.

## Running

```r
# Redraw the published Figure 4 panels (~5 min, ~8 GB RAM)
source("02_figures.R")

# Re-run the whole pipeline from the CellRanger matrix (~2 h, ~40 GB RAM)
source("01_analysis.R")
```

`02_figures.R` is the one to start with: it reproduces the published panels
directly and needs only the small deposited tables.

## Reproducibility

`01_analysis.R` is a faithful record of how the published clustering was
produced, but **it will not reproduce it exactly**. Doublet calling, PCA,
Harmony correction, the UMAP and TriMap embeddings, and Leiden clustering are
all stochastic. Seeds are set throughout, yet results still depend on the
BLAS/LAPACK build, the number of threads and the package versions, so a re-run
gives an equivalent — not identical — partition of the cells.

Two practical consequences:

- Step 7 of `01_analysis.R` maps curated biological names onto cluster IDs
  `k01`–`k10`. That mapping is only valid for the published run. If you re-run
  the pipeline, check the marker table before reusing it.
- The published cluster labels and TriMap coordinates are deposited on Zenodo
  (`cell_metadata.tsv.gz`). `02_figures.R` reads those, which is why it
  reproduces the figures and `01_analysis.R` does not.

## Dependencies

R ≥ 4.4, Bioconductor ≥ 3.20.

**CRAN** — `Matrix`, `ggplot2`, `matrixStats`, `stringr`, `doParallel`,
`aricode`, `harmony`, `RSQLite`, `igraph`, `uwot`

**Bioconductor** — `SingleCellExperiment`, `SummarizedExperiment`, `scuttle`,
`scater`, `scran`, `scDblFinder`, `ComplexHeatmap`, `BiocParallel`,
`JASPAR2024`, `GO.db`, `org.Hs.eg.db`, `reactome.db`

**GitHub** — [`oob`](https://github.com/DimitriMeistermann/oob) (TriMap/UMAP
wrappers, Leiden clustering, `heatmap.DM`, `proj2d`) and
[`GSDS`](https://github.com/DimitriMeistermann/GSDS) (module activation scores,
over-representation analysis).

```r
BiocManager::install(c("SingleCellExperiment", "scuttle", "scater", "scran",
                       "scDblFinder", "ComplexHeatmap", "JASPAR2024",
                       "GO.db", "org.Hs.eg.db", "reactome.db"))
remotes::install_github(c("DimitriMeistermann/oob",
                          "DimitriMeistermann/GSDS"))
```

Note that `GO.db`, `org.Hs.eg.db`, `reactome.db` and `GenomeInfoDbData` are
annotation *data* packages. Some conda builds ship them as empty stubs; install
them from Bioconductor if package loading fails with a missing-data error.

## Figure panels produced

| Panel | Output file | Content |
|---|---|---|
| 4A | `Fig4A_moduleExpressionHeatmap.pdf` | Module expression heatmap, 15 genes per module, 500 cells per astrocyte cluster |
| 4B | `Fig4B_trimapClusters.pdf` | TriMap embedding coloured by cluster |
| 4C | `Fig4C_clusterAbundance.pdf`, `Fig4C_oddsRatioHeatmap.pdf` | Cluster abundance per cell line; log2(OR) cluster × line association |
| 4E | `Fig4E_FADS_ELOVL.pdf` | *FADS1*, *FADS2*, *ELOVL2*, *ELOVL5* by supercluster |
| 4J | `Fig4J_ACOX3_ACSL.pdf` | *ACOX3*, *ACSL3*, *ACSL4* by supercluster |
| 4K | `Fig4K_FABP_SLC27A1.pdf` | *FABP7*, *FABP5*, *SLC27A1* by supercluster |

Panels 4D and 4F–4I, 4L–4O are not single-cell data and are not produced here.

## Licence

Code is released under the MIT licence (`LICENSE`). Data in the Zenodo record
are released under CC-BY-4.0.
