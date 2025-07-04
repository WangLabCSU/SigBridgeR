library(Seurat)
library(preprocessCore)
library(scAB)


data("data_survival")
dim(sc_dataset)

head(bulk_dataset[, 1:10])
head(phenotype)

# sc_dataset <- run_seurat(sc_dataset, verbose = FALSE)
# UMAP_celltype <- DimPlot(sc_dataset, reduction = "umap", group.by = "celltype")
# UMAP_celltype

scAB_data <- create_scAB(sc_dataset, bulk_dataset, phenotype)

K <- select_K(scAB_data)
K

scAB_result <- scAB(Object = scAB_data, K = K)
sc_dataset <- findSubset(sc_dataset, scAB_Object = scAB_result, tred = 2)

UMAP_scAB <- DimPlot(
    sc_dataset,
    group.by = "scAB_select",
    cols = c("#80b1d3", "red"),
    pt.size = 0.001,
    order = c("scAB+ cells", "Other cells")
)
patchwork::wrap_plots(plots = list(UMAP_celltype, UMAP_scAB), ncol = 2)

UMAP_subset3 <- DimPlot(
    sc_dataset,
    group.by = "scAB_Subset3",
    cols = c("#80b1d3", "red"),
    pt.size = 0.001,
    order = c("scAB+ cells", "Other cells")
)
UMAP_subset5 <- DimPlot(
    sc_dataset,
    group.by = "scAB_Subset4",
    cols = c("#80b1d3", "red"),
    pt.size = 0.001,
    order = c("scAB+ cells", "Other cells")
)
UMAP_subset <- patchwork::wrap_plots(
    plots = list(UMAP_subset3, UMAP_subset5),
    ncol = 2
)
UMAP_subset

markers <- FindMarkers(
    sc_dataset,
    ident.1 = "scAB+ cells",
    group.by = 'scAB_select',
    logfc.threshold = 1
)
markers <- markers[which(markers$p_val_adj < 0.05), ]
head(markers)
