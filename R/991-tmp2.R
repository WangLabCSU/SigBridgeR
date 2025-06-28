# ---- binary ----

library(ScPP)
load(system.file("data/binary.RData", package = "ScPP"))

sc = sc_Preprocess(sc_count)
geneList = marker_Binary(bulk, binary, ref_group = "Normal")
res = ScPP(sc, geneList)
Matrix::head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
SeuratObject::Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey", "blue", "red"))


# ---- continuous ----

library(ScPP)
load(system.file("data/continuous.RData", package = "ScPP"))

sc = sc_Preprocess(sc_count)
geneList = marker_Continuous(bulk, continuous$TMB_non_silent)
res = ScPP(sc, geneList)
Matrix::head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
SeuratObject::Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey", "blue", "red"))

# ---- survival ----

library(ScPP)
load(system.file("data/survival.RData", package = "ScPP"))

sc = sc_Preprocess(sc_count)
geneList = marker_Survival(bulk, survival)
res = ScPP(sc, geneList)
Matrix::head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
SeuratObject::Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey", "blue", "red"))
