################################################################################
## 01_analysis.R — single-cell RNA-seq pipeline: CellRanger matrix -> clusters,
##                 gene modules and functional enrichment.
##
## hiPSC-derived astrocytes: PO2 (WT), HEL69-5 (FXS), HEL100-1 (FXS-epilepsy).
##
## Input   data/counts_raw_cellranger/    (Zenodo record — unpack
##                                        counts_raw_cellranger.tar; see README.md)
## Output  results/sce.rds                annotated SingleCellExperiment
##         results/cell_metadata.tsv.gz   per-cell table (as deposited)
##         results/gene_metadata.tsv.gz   per-gene table (as deposited)
##         results/cluster_markers.tsv.gz
##         results/enrichment/<module>.tsv
##
## ---------------------------------------------------------------------------
## REPRODUCIBILITY
##   Doublet detection, PCA, Harmony correction, the UMAP/TriMap embeddings and
##   Leiden clustering are stochastic. Seeds are set, but cluster IDs and
##   embedding coordinates still depend on the BLAS/LAPACK build, thread count
##   and package versions, so re-running this script gives an equivalent rather
##   than an identical clustering. Step 7 assigns curated names to clusters by
##   their `k01`..`k10` IDs; if you re-run the pipeline, verify the marker genes
##   before trusting that mapping.
##
##   The cluster labels and embedding coordinates behind the published figures
##   are deposited on Zenodo. To redraw the published panels, run 02_figures.R,
##   which reads those tables and does not re-cluster.
##
## Runtime ~2 h on 20 cores; peak memory ~40 GB.
################################################################################

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
  library(scuttle); library(scater); library(scran); library(scDblFinder)
  library(harmony); library(matrixStats); library(stringr)
  library(BiocParallel); library(doParallel)
  library(oob)    # github.com/DimitriMeistermann/oob
  library(GSDS)   # github.com/DimitriMeistermann/GSDS
})

set.seed(666)
BPPARAM <- MulticoreParam(workers = max(1, parallel::detectCores() - 2))
register(BPPARAM)
registerDoParallel(cores = BPPARAM$workers)

DATA    <- "data"
RESULTS <- "results"
dir.create(RESULTS, showWarnings = FALSE)
dir.create(file.path(RESULTS, "enrichment"), showWarnings = FALSE)

source("functions.R")   # fast_group_agg, renameModulePerTF, enrich2vector.V2, ...

################################################################################
## 1. Load the CellRanger aggregate matrix, collapse to gene symbols
################################################################################
mdir <- file.path(DATA, "counts_raw_cellranger")

counts_mat <- Matrix::readMM(file.path(mdir, "matrix.mtx.gz")) |> as("dgCMatrix")
cell_annot <- fastRead(file.path(mdir, "barcodes.tsv.gz"),
                       row.names = NULL, header = FALSE)
colnames(cell_annot) <- "agg_barcode"

## The CellRanger `aggr` suffix (-1/-2/-3) is the row order of aggregation.csv.
aggr <- read.csv(file.path(mdir, "aggregation.csv"))
cell_annot$run_no <- as.integer(sub(".*-", "", cell_annot$agg_barcode))
cell_annot$run    <- aggr$cell_line[cell_annot$run_no]
cell_annot$sample <- aggr$genotype [cell_annot$run_no]
rn(cell_annot)    <- cell_annot$agg_barcode

feature_annot <- fastRead(file.path(mdir, "features.tsv.gz"),
                          row.names = NULL, header = FALSE)[, c(1, 2)]
colnames(feature_annot) <- c("ensembl_id", "gene_symbol")

cn(counts_mat) <- rn(cell_annot)
## Several Ensembl IDs share a gene symbol; sum their counts so every downstream
## lookup is by symbol.
counts_mat <- fast_group_agg(x = counts_mat, by = feature_annot$gene_symbol) |>
  as("dgCMatrix")

sce <- SingleCellExperiment(list(counts = counts_mat), colData = cell_annot)
rm(cell_annot, counts_mat); gc()

################################################################################
## 2. Cell QC
################################################################################
per.cell <- scater::perCellQCMetrics(sce)
colData(sce) <- cbind(colData(sce), per.cell)

## Fixed floor on detected genes, plus scuttle's adaptive lower outlier cut on
## the same statistic. A fixed threshold alone would keep low-complexity cells
## in the deepest library; the adaptive cut alone would remove different
## fractions from each of the three libraries.
discard <- (per.cell$detected < 300) |
  isOutlier(per.cell$detected, type = "lower", log = TRUE)
sce <- sce[, !discard]

################################################################################
## 3. Normalisation and highly variable genes
################################################################################
sce <- scuttle::logNormCounts(sce)
dec <- scran::modelGeneVar(sce, assay.type = "logcounts")
hvg <- scran::getTopHVGs(dec, n = 3000)

################################################################################
## 4. Doublet removal
################################################################################
set.seed(666)
sce <- scDblFinder::scDblFinder(sce, samples = colData(sce)$sample)
sce <- sce[, sce$scDblFinder.class == "singlet"]; gc()

################################################################################
## 5. PCA, Harmony batch correction, embeddings
################################################################################
## Each cell line was sequenced as its own library, so line and batch are
## confounded by design. Harmony is run on `batch` (= library) to align shared
## cell states across libraries while keeping the genotype contrast readable.
sce$batch <- sce$run

sce <- scater::runPCA(sce, subset_row = hvg, ncomponents = 50)
reducedDim(sce, "reducedSpace") <- harmony::HarmonyMatrix(
  data_mat  = reducedDim(sce, "PCA"),      # cells x PCs
  meta_data = as.data.frame(colData(sce)),
  vars_use  = "batch",
  verbose   = TRUE
)

################################################################################
## 6. Clustering (Leiden on the UMAP k-NN graph)
################################################################################
## First pass identifies a low-quality cluster (k12: low complexity, no
## distinguishing markers); it is removed and the clustering redone on the
## remaining cells.
set.seed(666)
umapRes      <- oob::UMAP(reducedDim(sce, "reducedSpace"), n_neighbors = 200,
                          ret_nn = TRUE, transpose = FALSE)
sce$clusters <- leidenFromUMAP(umapRes, n_neighbors = 30,
                               resolution_parameter = 1)
sce <- sce[, sce$clusters != "k12"]

set.seed(666)
umapRes      <- oob::UMAP(reducedDim(sce, "reducedSpace"), n_neighbors = 300,
                          ret_nn = TRUE, transpose = FALSE)
sce$clusters <- leidenFromUMAP(umapRes, n_neighbors = 30,
                               resolution_parameter = 0.8)

set.seed(666)
trimapRes <- oob::TRIMAP(reducedDim(sce, "reducedSpace"), transpose = FALSE)
rownames(trimapRes) <- cn(sce)
reducedDim(sce, "TRIMAP") <- trimapRes

################################################################################
## 7. Markers and curated cluster names
################################################################################
markers <- getMarkers2(sce, sce$clusters, assay = "logcounts")
fastWrite(markers, file.path(RESULTS, "cluster_markers.tsv"))

## Names assigned by reading the marker table together with the per-module
## enrichment of step 9 — see the manuscript for the reasoning behind each.
## The component before the dot is the "supercluster" used in Fig 4E/J/K.
hashNameCluster <- c(
  "k01" = "MetaboHigh.highTranslation",
  "k02" = "Cycling",
  "k03" = "ImmuneDev.matureAstro",
  "k04" = "Morph.FXS_neurogen",
  "k05" = "MetaboHigh.neurogen",
  "k06" = "ImmuneDev.immature",
  "k07" = "Morph.TGFB",
  "k08" = "MetaboHigh.stressed",
  "k09" = "Morph.WT_neurogen",
  "k10" = "Mesenchymal"
)
colorScale <- list(
  run    = c("PO2" = "darkblue", "HEL69-5" = "brown", "HEL100-1" = "orange"),
  sample = c("WT"  = "darkblue", "FXS"     = "brown",
             "FXS-epilepsy" = "orange"),
  anno   = c("Cycling"                    = "#EE2F5F",
             "ImmuneDev.matureAstro"      = "#0A6B37",
             "ImmuneDev.immature"         = "#A2CD3A",
             "MetaboHigh.highTranslation" = "#E3AE27",
             "MetaboHigh.neurogen"        = "#FFFFA0",
             "MetaboHigh.stressed"        = "#A47844",
             "Morph.FXS_neurogen"         = "#67469C",
             "Morph.TGFB"                 = "#4AC8EE",
             "Morph.WT_neurogen"          = "#273A93",
             "Mesenchymal"                = "#B2B2B2")
)
lvlOrder <- names(colorScale$anno)

stopifnot(setequal(unique(sce$clusters), names(hashNameCluster)))
sce$anno         <- factor(hashNameCluster[sce$clusters], levels = lvlOrder)
sce$run          <- factor(sce$run,    levels = names(colorScale$run))
sce$sample       <- factor(sce$sample, levels = names(colorScale$sample))
sce$superCluster <- strsplitNth(as.character(sce$anno), split = ".",
                                fixed = TRUE, n = 1)

################################################################################
## 8. Gene modules
################################################################################
## Modules are Leiden communities in a gene-gene correlation k-NN graph, built
## on a 100-cell-per-cluster subsample (the full cell x cell correlation is not
## needed and the subsample keeps the graph tractable).
set.seed(666)
draw    <- oob::drawSamplePerGroup(cn(sce), sce$clusters, maxDrawSize = 100)
sceDraw <- sce[, names(draw)]

k <- 20
set.seed(666)
geneCluster <- UMAP(assay(sceDraw, "logcounts"), metric = "correlation",
                    transpose = FALSE, ret_nn = TRUE, n_neighbors = k) |>
  leidenFromUMAP(n_neighbors = k, resolution_parameter = 1.5)
geneCluster <- str_replace_all(geneCluster, "k", "M")
names(geneCluster) <- rn(sce)

## Modules with fewer than 10 genes are pooled into the unassigned group M0.
nGenePerModule <- table(geneCluster)
okSizeModule   <- names(nGenePerModule)[nGenePerModule >= 10]
genePerModule  <- sapply(okSizeModule, function(m) names(geneCluster)[geneCluster == m])

## Activation score = PC1 of the module's expression submatrix.
moduleActivationScoreDt <- GSDS::activeScorePCAlist(
  exprMatrix = as.matrix(assay(sceDraw, "logcounts")),
  geneList   = genePerModule[okSizeModule]
)
moduleActivationScore <- moduleActivationScoreDt$activScoreMat

geneAnnot         <- rowData(sce)
geneAnnot$Module  <- "M0"
geneAnnot[names(geneCluster), "Module"] <- geneCluster
geneAnnot[!geneAnnot$Module %in% okSizeModule, "Module"] <- "M0"
geneAnnot$Module  <- as.factor(geneAnnot$Module)

## Membership = correlation of the gene with its module's activation score.
## Contribution = the gene's loading on that PC1.
moduleMembership <- suppressWarnings(
  cor(t(as.matrix(assay(sceDraw, "logcounts"))), moduleActivationScore))
colnames(moduleMembership) <- colnames(moduleActivationScore)

geneAnnot$Membership   <- NA
geneAnnot$Contribution <- NA
for (gene in rn(geneAnnot)[geneAnnot$Module != "M0"]) {
  module <- as.character(geneAnnot[gene, "Module"])
  geneAnnot[gene, "Membership"]   <- moduleMembership[gene, module]
  geneAnnot[gene, "Contribution"] <- moduleActivationScoreDt$contributionList[[module]][gene]
}
rm(moduleActivationScoreDt); gc()

## Rename M1, M2, ... after the highest-contributing transcription factor in
## each module, giving the TCF12.Mod / JUND.Mod / ... labels used in Figure 4A.
geneAnnot <- renameModulePerTF(geneAnnot, dl_JAFAR_geneSyms_Human())
colnames(moduleActivationScore) <- levels(geneAnnot$Module)
names(genePerModule)            <- levels(geneAnnot$Module)

################################################################################
## 9. Astrocyte subset and a module recomputation restricted to it
################################################################################
## k10 (Mesenchymal, 179 cells) is a contaminating non-astrocyte population; it
## is excluded from the module analysis and from Figure 4A.
sceAstro     <- sce[, sce$clusters != "k10"]
set.seed(666)
drawAstro    <- oob::drawSamplePerGroup(cn(sceAstro), sceAstro$clusters,
                                        maxDrawSize = 500)
sceAstroDraw <- sceAstro[, names(drawAstro)]

## Recomputing modules on the astrocyte subset alone resolves one additional
## module (HES6.Mod) that the whole-dataset run splits across others. Only that
## module is carried over; the rest of the attribution stays as in step 8.
goodGenes <- rowSums2(assay(sceAstroDraw, "logcounts")) > 100
TFgynSyms <- dl_JAFAR_geneSyms_Human()
set.seed(666)
newGeneModules <- GSDS::computegenemodules(
  as.matrix(assay(sceAstroDraw, "logcounts")[goodGenes, ]),
  TF_list = TFgynSyms, resolution_parameter = 0.8
)
geneAnnot2 <- newGeneModules$gene_annot
HES6genes  <- rn(geneAnnot2)[geneAnnot2$Module_newName == "HES6.Mod"]
geneAnnot[HES6genes, "Module"]       <- "HES6.Mod"
geneAnnot[HES6genes, "Membership"]   <- geneAnnot2[HES6genes, "Membership"]
geneAnnot[HES6genes, "Contribution"] <- geneAnnot2[HES6genes, "Contribution"]

fastWrite(as.data.frame(geneAnnot), file.path(RESULTS, "module_attribution.tsv"))

################################################################################
## 10. Over-representation analysis per gene module
################################################################################
species.data     <- GSDS::getSpeciesData("Human")
geneSetsDatabase <- getDBterms(rn(sce), species = "Human",
                               returnGenesSymbol = TRUE,
                               database = c("kegg", "goBP", "goCC", "goMF", "reactom"))

for (module in levels(as.factor(geneAnnot$Module))) {
  isInModule        <- geneAnnot$Module == module
  names(isInModule) <- rn(geneAnnot)
  res <- enrich.ora(isInModule, speciesData = species.data,
                    db_terms = geneSetsDatabase, returnGenes = TRUE)
  exportEnrich(res, file.path(RESULTS, "enrichment", paste0(module, ".tsv")))
}

################################################################################
## 11. Export the tables deposited on Zenodo
################################################################################
saveRDS(sce, file.path(RESULTS, "sce.rds"))

cellMeta <- data.frame(
  barcode           = cn(sce),
  cell_line         = as.character(sce$run),
  genotype          = as.character(sce$sample),
  library_index     = sce$run_no,
  n_umi             = sce$sum,
  n_gene_detected   = sce$detected,
  size_factor       = sizeFactors(sce),
  scDblFinder_class = as.character(sce$scDblFinder.class),
  scDblFinder_score = sce$scDblFinder.score,
  cluster           = as.character(sce$clusters),
  cluster_name      = as.character(sce$anno),
  supercluster      = sce$superCluster,
  trimap_1          = trimapRes[cn(sce), 1],
  trimap_2          = trimapRes[cn(sce), 2],
  in_astro_subset   = cn(sce) %in% cn(sceAstro),
  in_fig4A_draw     = cn(sce) %in% cn(sceAstroDraw)
)
write.table(cellMeta, gzfile(file.path(RESULTS, "cell_metadata.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

ens <- tapply(feature_annot$ensembl_id, feature_annot$gene_symbol,
              function(x) paste(sort(unique(x)), collapse = ";"))
geneMeta <- data.frame(
  gene_symbol         = rn(sce),
  ensembl_id          = unname(ens[rn(sce)]),
  module              = as.character(geneAnnot[rn(sce), "Module"]),
  module_membership   = geneAnnot[rn(sce), "Membership"],
  module_contribution = geneAnnot[rn(sce), "Contribution"],
  is_hvg              = rn(sce) %in% hvg
)
write.table(geneMeta, gzfile(file.path(RESULTS, "gene_metadata.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(RESULTS, "sessionInfo_01.txt"))
message("done — see ", RESULTS)
