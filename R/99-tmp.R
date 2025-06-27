library(scPAS)
library(Matrix)
library(Seurat)

load(
  file = '~/Data/scPAS_tutorial_data-20250627T025840Z-1-001/scPAS_tutorial_data/data_survival.rda'
)
dim(sc_dataset)
#> [1] 21287  8853

head(bulk_dataset[, 1:4])
#>           TCGA-2Y-A9GS-01 TCGA-2Y-A9GT-01 TCGA-2Y-A9GU-01 TCGA-2Y-A9GV-01
#> ARHGEF10L         11.5199         12.2544         12.6434         11.9114
#> HIF3A              2.4425          4.6559          1.8083          3.7751
#> RNF17              3.7881          0.0000          0.0000          0.0000
#> RNF10             11.5109         11.7863         12.3187         11.7253
#> RNF11             10.6574         11.0133         10.5722         10.9586
#> RNF13             10.2623         10.4574          9.9284         10.5863
#>           TCGA-2Y-A9GW-01 TCGA-2Y-A9GX-01 TCGA-2Y-A9GY-01 TCGA-2Y-A9GZ-01
#> ARHGEF10L         11.8108         10.6616         10.6321         12.4218
#> HIF3A              2.8340          3.2456          7.1429          5.7972
#> RNF17              0.0000          0.0000          7.9222          0.7507
#> RNF10             11.8382         11.5872         11.5770         11.2699
#> RNF11             10.3071         10.7500          9.9313         10.7820
#> RNF13             10.1785         10.4880         10.0815         10.5819
#>           TCGA-2Y-A9H0-01 TCGA-2Y-A9H1-01
#> ARHGEF10L         10.8164         10.2715
#> HIF3A              3.7561          4.3687
#> RNF17              0.0000          0.0000
#> RNF10             12.7789         11.2257
#> RNF11             10.0653         10.1549
#> RNF13             10.6751          9.9593
head(phenotype)
#>                 time status
#> TCGA-2Y-A9GS-01  724      1
#> TCGA-2Y-A9GT-01 1624      1
#> TCGA-2Y-A9GU-01 1939      0
#> TCGA-2Y-A9GV-01 2532      1
#> TCGA-2Y-A9GW-01 1271      1
#> TCGA-2Y-A9GX-01 2442      0

keep_genes <- rownames(sc_dataset)[
  rowSums(sc_dataset@assays$RNA@counts > 1) >= 20
]
sc_dataset <- subset(sc_dataset, features = keep_genes)

dim(sc_dataset)
# > dim(sc_dataset)
# [1] 10775  8853

sc_dataset0 <- run_Seurat(sc_dataset, verbose = FALSE)
# # sc_dataset1 = SCPreProcess(sc_dataset, verbose = FALSE)

# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
# col_vector = unlist(mapply(
#   brewer.pal,
#   qual_col_pals$maxcolors,
#   rownames(qual_col_pals)
# ))
# set.seed(seed = 1234)
# cellType_col <- sample(col_vector, 7)
# names(cellType_col) <- unique(sc_dataset0$celltype)

library(randomcoloR)
cellType_col <- randomcoloR::distinctColorPalette(
  length(unique(sc_dataset0$celltype)),
  runTsne = TRUE
) %>%
  setNames(unique(sc_dataset0$celltype))

# * use my FetchUMAp
# UMAP_celltype <- DimPlot(
#   sc_dataset0,
#   reduction = "umap",
#   group.by = "celltype",
#   label = T,
#   cols = cellType_col
# )
UMAP_celltype = FetchUMAP(
  sc_dataset0,
  plot_color = cellType_col,
  group_by = "celltype",
  order = NULL,
  plot_show = T
)

# *bulk

bulk_dataset <- as.matrix(bulk_dataset)
class(bulk_dataset)
#> [1] "matrix" "array"
dim(bulk_dataset)
#> [1] 18807   370

# ?只保留那些0值数量少于总列数25%的行
bulk_dataset = bulk_dataset[
  apply(bulk_dataset, 1, function(x) sum(x == 0) < 0.25 * ncol(bulk_dataset)),
]
dim(bulk_dataset)
#> [1] 15874   370

# Check
all(rownames(phenotype) == colnames(bulk_dataset))

# The input provided above is survival data; therefore, scPAS needs to fit a Cox regression model (family = ‘cox’).
# imputation = F , without imputation.
# nfeature = 3000 indicates that the top 3000 highly variable genes are selected for model training.
# network_class = ‘SC’ indicates gene-gene similarity networks derived from single-cell data.
scPAS_result <- scPAS::scPAS(
  bulk_dataset = bulk_dataset,
  sc_dataset = sc_dataset0,
  assay = 'RNA',
  phenotype = phenotype,
  imputation = F,
  nfeature = 3000,
  alpha = 0.01,
  network_class = 'SC',
  family = 'cox'
)

# The scPAS provides both quantitative (scPAS_NRS) and qualitative (scPAS) prediction results:
library(ggplot2)
UMAP_scPAS <- DimPlot(
  sc_dataset,
  reduction = 'umap',
  group.by = 'scPAS',
  raster = FALSE,
  cols = c("0" = 'grey', "scPAS+" = 'indianred1', "scPAS-" = 'royalblue'),
  pt.size = 0.001,
  order = c("scPAS-", "scPAS+")
) +
  ggtitle(label = NULL) +
  labs(col = "scPAS")

UMAP_RS <- FeaturePlot(
  object = sc_dataset,
  raster = FALSE,
  features = 'scPAS_NRS',
  max.cutoff = 2,
  min.cutoff = -2
) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) +
  ggtitle(label = NULL) +
  labs(col = "Risk score")

UMAP_celltype | UMAP_RS | UMAP_scPAS
