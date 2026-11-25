<a href="https://lifecycle.r-lib.org/articles/stages.html#deprecated"><img src="https://img.shields.io/badge/lifecycle-deprecated-orange.svg" alt="Lifecycle: deprecated" align="right"></a>

> **Deprecated.** This vignette is deprecated as of SigBridgeR 4.0.0 and
> will be removed in a future release. `SCPreProcess()` and the
> registry-based preprocessing pipeline are superseded

    library(SigBridgeR)

# Extend Seurat Methods

In general, any Seurat function can be registered for use within
`SCPreProcess`. However, due to the character-length limitation of the
pipeline code string, at most **52** functions can be registered-one for
each uppercase and lowercase letter (A-Z, a-z).

By default, the registry is defined as follows:

**Pipeline Code Table:**

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 30%" />
<col style="width: 62%" />
</colgroup>
<thead>
<tr>
<th>Code</th>
<th>Function</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>o</strong></td>
<td><code>CreateSeuratObject</code></td>
<td><strong>Required.</strong> Must be the first step.</td>
</tr>
<tr>
<td><strong>n</strong></td>
<td><code>NormalizeData</code></td>
<td>Standard normalization.</td>
</tr>
<tr>
<td><strong>s</strong></td>
<td><code>ScaleData</code></td>
<td>Scales data for PCA.</td>
</tr>
<tr>
<td><strong>v</strong></td>
<td><code>FindVariableFeatures</code></td>
<td>Selects highly variable genes.</td>
</tr>
<tr>
<td><strong>p</strong></td>
<td><code>RunPCA</code></td>
<td>Principal Component Analysis.</td>
</tr>
<tr>
<td><strong>e</strong></td>
<td><code>FindNeighbors</code></td>
<td>Computes SNN graph.</td>
</tr>
<tr>
<td><strong>c</strong></td>
<td><code>FindClusters</code></td>
<td>Louvain algorithm clustering.</td>
</tr>
<tr>
<td><strong>t</strong></td>
<td><code>RunTSNE</code></td>
<td>t-SNE reduction.</td>
</tr>
<tr>
<td><strong>u</strong></td>
<td><code>RunUMAP</code></td>
<td>UMAP reduction.</td>
</tr>
<tr>
<td><strong>r</strong></td>
<td><code>SCTransform</code></td>
<td><strong>SCT workflow.</strong> Replaces n, s, v.</td>
</tr>
</tbody>
</table>

For example, if we need to use `Seurat::LoadXenium` to load Xenium data
for screening, register this function first.

    Register(x = Seurat::LoadXenium)

To check:

    names(SCPreProcessStrategy)

    InterceptStrategy("SCPreProcessStrategy")

overwrite:

    Register(x = Seurat::LoadXenium, overwrite = TRUE)

## Session Info
