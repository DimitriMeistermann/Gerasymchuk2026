################################################################################
## 02_figures.R — regenerate the single-cell panels of Figure 4.
##
## Reads ONLY the Zenodo record (see README.md); no re-clustering, no raw reads.
## The cluster labels and TriMap coordinates it uses are the published ones, so
## the panels it draws match the paper.
##
## Input   data/counts_filtered/            (unpack counts_filtered.tar first)
##                                          matrix.mtx.gz, genes.tsv.gz,
##                                          barcodes.tsv.gz
##         data/cell_metadata.tsv.gz
##         data/gene_metadata.tsv.gz
##         data/cluster_annotation.tsv
##         data/cell_line_annotation.tsv
##         data/fig4A_gene_selection.tsv
##
## Output  figures/Fig4A_moduleExpressionHeatmap.pdf
##         figures/Fig4B_trimapClusters.pdf
##         figures/Fig4C_clusterAbundance.pdf
##         figures/Fig4C_oddsRatioHeatmap.pdf
##         figures/Fig4E_FADS_ELOVL.pdf
##         figures/Fig4J_ACOX3_ACSL.pdf
##         figures/Fig4K_FABP_SLC27A1.pdf
##
## Runtime ~5 min, peak memory ~8 GB.
################################################################################

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
  library(scuttle); library(ggplot2); library(ComplexHeatmap)
  library(oob)    # github.com/DimitriMeistermann/oob
  library(GSDS)   # github.com/DimitriMeistermann/GSDS
})

set.seed(666)
DATA <- "data"
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

source("functions.R")   # enrich2vector.V2

################################################################################
## 1. Load the deposited matrix and metadata
################################################################################
mdir <- file.path(DATA, "counts_filtered")
cts  <- Matrix::readMM(file.path(mdir, "matrix.mtx.gz")) |> as("dgCMatrix")
rownames(cts) <- readLines(gzfile(file.path(mdir, "genes.tsv.gz")))
colnames(cts) <- readLines(gzfile(file.path(mdir, "barcodes.tsv.gz")))

cellMeta <- read.delim(gzfile(file.path(DATA, "cell_metadata.tsv.gz")))
geneMeta <- read.delim(gzfile(file.path(DATA, "gene_metadata.tsv.gz")))
clustAnn <- read.delim(file.path(DATA, "cluster_annotation.tsv"))
lineAnn  <- read.delim(file.path(DATA, "cell_line_annotation.tsv"))
figGenes <- read.delim(file.path(DATA, "fig4A_gene_selection.tsv"))

stopifnot(identical(cellMeta$barcode, colnames(cts)),
          identical(geneMeta$gene_symbol, rownames(cts)))

## Palettes and factor orders exactly as published.
colorScale <- list(
  run    = setNames(lineAnn$colour,  lineAnn$cell_line),
  sample = setNames(lineAnn$colour,  lineAnn$genotype),
  anno   = setNames(clustAnn$colour, clustAnn$cluster_name)
)
lvlOrder <- clustAnn$cluster_name[order(clustAnn$plot_order)]

rownames(cellMeta)  <- cellMeta$barcode
cellMeta$anno       <- factor(cellMeta$cluster_name, levels = lvlOrder)
cellMeta$run        <- factor(cellMeta$cell_line,    levels = lineAnn$cell_line)
cellMeta$sample     <- factor(cellMeta$genotype,     levels = lineAnn$genotype)
cellMeta$superClust <- factor(cellMeta$supercluster,
                              levels = unique(sub("\\..*$", "", lvlOrder)))

sce <- SingleCellExperiment(list(counts = cts),
                            colData = DataFrame(cellMeta, row.names = rownames(cellMeta)))
sce <- scuttle::logNormCounts(sce, size.factors = cellMeta$size_factor)
rm(cts); gc()

## Subsets used below: the astrocyte subset (Fig 4E/J/K) and the 500-cells-per-
## cluster draw the published heatmap was rendered from (Fig 4A).
sceAstro     <- sce[, sce$in_astro_subset]
sceAstroDraw <- sce[, sce$in_fig4A_draw]
drawAstroAnnot <- setNames(as.character(sceAstroDraw$anno), colnames(sceAstroDraw))
trimapRes      <- as.matrix(cellMeta[, c("trimap_1", "trimap_2")])

################################################################################
## Figure 4A — module expression heatmap
################################################################################
## Rows: the deposited 15-genes-per-module selection plus the `other` group,
## in deposited order. Columns: 500 cells per astrocyte cluster, split by
## cluster, annotated by cell line.
genePerModSelVect <- factor(figGenes$row_group,
                            levels = unique(figGenes$row_group))
names(genePerModSelVect) <- figGenes$gene_symbol
stopifnot(all(names(genePerModSelVect) %in% rownames(sceAstroDraw)))

pdf(file.path(FIGS, "Fig4A_moduleExpressionHeatmap.pdf"), width = 16, height = 12)
heatmap.DM(
  assay(sceAstroDraw[names(genePerModSelVect), ], "logcounts"),
  column_split       = drawAstroAnnot,
  row_split          = genePerModSelVect,
  cluster_column_slice = FALSE,
  cluster_row_slices = FALSE,
  colData            = sceAstroDraw$run,
  colorAnnot         = colorScale["run"],
  row_title_rot      = 0,
  column_title_rot   = 45,
  show_column_names  = FALSE,
  use_raster         = TRUE,
  raster_by_magick   = FALSE,
  returnHeatmap      = FALSE
)
dev.off()

################################################################################
## Figure 4B — TriMap embedding coloured by cluster
################################################################################
drawProj <- sample.int(nrow(trimapRes))   # plot order only, avoids overplotting

pdf(file.path(FIGS, "Fig4B_trimapClusters.pdf"), width = 10, height = 9)
proj2d(trimapRes[drawProj, ], colorBy = sce$anno[drawProj],
       plotFactorsCentroids = TRUE, useScatterMore = FALSE,
       pointSize = 2, colorScale = colorScale$anno, legendTitle = "Cluster")
dev.off()

################################################################################
## Figure 4C — cluster abundance per cell line, and enrichment odds ratios
################################################################################
ggAbund <- ggplot(as.data.frame(colData(sce)), aes(x = run, fill = anno)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = colorScale$anno) +
  labs(x = "Cell line", y = "Number of cells", fill = "Cluster") +
  theme_bw()
ggsave(file.path(FIGS, "Fig4C_clusterAbundance.pdf"), plot = ggAbund,
       width = 5, height = 7)

## Odds ratio of each cluster in each cell line, relative to the whole dataset.
resEnrich <- enrich2vector.V2(sce$anno, sce$run,
                              varnames = c("clusters", "sample"))
ORmat <- matrixFromDimnames(levels(sce$anno), levels(sce$run))
for (ctype in rownames(ORmat))
  for (cond in colnames(ORmat))
    ORmat[ctype, cond] <- resEnrich[resEnrich$clusters == ctype &
                                    resEnrich$sample   == cond, "OR"]
ORmat[ORmat == 0] <- 1e-4     # keep log2 finite for absent combinations

pdf(file.path(FIGS, "Fig4C_oddsRatioHeatmap.pdf"), width = 10, height = 10)
heatmap.DM(log2(ORmat), preSet = "vanilla", showValues = TRUE,
           name = "log2(Odds-Ratio)",
           cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

################################################################################
## Figures 4E / 4J / 4K — expression by supercluster and cell line
################################################################################
## Violin (distribution) + inner boxplot (median and IQR), dodged by cell line
## within each supercluster. Astrocyte subset only.
vioboxPanel <- function(genes, file, width = 11, height = 5) {
  genes <- intersect(genes, rownames(sceAstro))
  stopifnot(length(genes) > 0)
  ggdat <- do.call(rbind, lapply(genes, function(g) data.frame(
    gene         = g,
    expr         = as.matrix(assay(sceAstro[g, ], "logcounts"))[1, ],
    superCluster = sceAstro$superClust,
    run          = sceAstro$run
  )))
  ggdat$gene <- factor(ggdat$gene, levels = genes)

  p <- ggplot(ggdat, aes(x = superCluster, y = expr, fill = run)) +
    geom_violin(position = position_dodge(width = 0.9, preserve = "total"),
                scale = "width") +
    geom_boxplot(aes(group = interaction(superCluster, run)),
                 fill = "black", width = 0.08, color = "white", linewidth = 0.3,
                 position = position_dodge(width = 0.9, preserve = "total"),
                 outlier.shape = NA) +
    facet_wrap(~ gene, ncol = 1, strip.position = "right", scales = "free_y") +
    scale_fill_manual(values = colorScale$run, name = "Cell line") +
    labs(x = "Super cluster", y = "logcount") +
    theme_bw()

  pdf(file.path(FIGS, file), width = width,
      height = height * max(1, length(genes) / 2))
  print(ggBorderedFactors(p))
  dev.off()
}

## 4E — fatty-acid desaturases and elongases
vioboxPanel(c("FADS1", "FADS2", "ELOVL2", "ELOVL5"), "Fig4E_FADS_ELOVL.pdf")
## 4J — peroxisomal / acyl-CoA synthetases
vioboxPanel(c("ACOX3", "ACSL3", "ACSL4"),            "Fig4J_ACOX3_ACSL.pdf")
## 4K — fatty-acid binding and transport
vioboxPanel(c("FABP7", "FABP5", "SLC27A1"),          "Fig4K_FABP_SLC27A1.pdf")

writeLines(capture.output(sessionInfo()), file.path(FIGS, "sessionInfo_02.txt"))
message("done — see ", FIGS)
