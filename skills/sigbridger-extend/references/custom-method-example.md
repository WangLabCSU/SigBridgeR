# Worked example: a custom screening method named `Foo`

A deliberately simple method — bulk differential genes drive a per-cell score — that
still satisfies every SigBridgeR contract.

## 1. The executor

```r
#' @title Foo: bulk-marker based cell scoring
#'
#' @param matched_bulk Genes x samples bulk expression matrix.
#' @param sc_data A preprocessed Seurat object.
#' @param phenotype Named vector, names matching `colnames(matched_bulk)`.
#' @param label_type Label written to `sc_data@meta.data`.
#' @param phenotype_class One of "binary", "survival", "continuous".
#' @param ... `verbose`, `assay`, `top_n`.
#'
#' @return A list with `scRNA_data` (Seurat) plus `foo_genes` and `foo_score`.
#' @export
Foo <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  ...
) {
  # -- contract: argument checks -------------------------------------------
  phenotype_class <- rlang::arg_match(phenotype_class)
  if (phenotype_class != "binary") {
    rlang::abort("Foo only supports binary phenotypes.")
  }

  dots <- rlang::list2(...)
  verbose <- if (is.null(dots$verbose)) TRUE else dots$verbose
  assay <- if (is.null(dots$assay)) "RNA" else dots$assay
  top_n <- if (is.null(dots$top_n)) 50L else as.integer(dots$top_n)
  label_type <- if (is.null(label_type)) "Foo" else label_type

  # -- contract: sample-id alignment ---------------------------------------
  if (is.null(names(phenotype))) {
    rlang::abort("`phenotype` must be a named vector for binary phenotypes.")
  }
  if (!setequal(names(phenotype), colnames(matched_bulk))) {
    rlang::abort("`names(phenotype)` must match `colnames(matched_bulk)`.")
  }
  grp <- as.character(phenotype[colnames(matched_bulk)])
  levels_grp <- unique(grp)
  if (length(levels_grp) != 2L) {
    rlang::abort("Foo needs exactly two phenotype groups.")
  }

  if (verbose) {
    cli::cli_alert_info("Foo: scoring cells with the top {top_n} bulk markers")
  }

  # -- algorithm: Welch t-test per gene on bulk ----------------------------
  pval <- apply(matched_bulk, 1L, function(x) {
    stats::t.test(x[grp == levels_grp[1L]], x[grp == levels_grp[2L]])$p.value
  })
  pval <- pval[!is.na(pval)]
  genes <- names(sort(pval))[seq_len(min(top_n, length(pval)))]

  # -- algorithm: per-cell mean expression of those markers ----------------
  expr <- Seurat::GetAssayData(sc_data, assay = assay, layer = "data")
  genes <- intersect(genes, rownames(expr))
  if (length(genes) == 0L) {
    rlang::abort("No bulk markers are present in the single-cell data.")
  }
  score <- Matrix::colMeans(expr[genes, , drop = FALSE])

  # -- contract: write labels into meta.data -------------------------------
  q <- stats::quantile(score, probs = c(0.2, 0.8))
  labels <- ifelse(
    score >= q[[2L]], "Positive",
    ifelse(score <= q[[1L]], "Negative", "Neutral")
  )
  sc_data[[label_type]] <- labels

  # -- contract: return a list whose `scRNA_data` is a Seurat object -------
  list(
    scRNA_data = sc_data,
    foo_genes = genes,
    foo_score = score
  )
}
```

## 2. Validate

```r
ValidateScreenFunc(Foo)   # TRUE, or an abort with a per-check report
```

## 3. Record the parameters

```r
config <- ScreenMethodConfig(
  method_name     = "Foo",
  param           = list(top_n = 50L, assay = "RNA"),
  method_version  = "0.1.0",
  phenotype_class = "binary",
  label_type      = "TP53"
)
```

## 4. Register

```r
Register(Foo = ScreenMethod(
  method_name     = "Foo",
  method_version  = "0.1.0",
  executor        = Foo,
  phenotype_class = "binary",
  mapper          = NULL
))

names(ScreenStrategy)                          # "Foo" is now listed
```

## 5. Use it exactly like a built-in method

```r
ref <- LoadRefData(data_type = "binary")

res <- Screen(
  matched_bulk    = ref[[2]],
  sc_data         = ref[[1]],
  phenotype       = ref[[3]],
  label_type      = "TP53",
  phenotype_class = "binary",
  screen_method   = "Foo",
  top_n           = 100L,        # forwarded to Foo() via ...
  verbose         = TRUE
)

print(res)                       # ScreenMethodResult
table(res@scRNA_data$TP53)       # Positive / Negative / Neutral
res$extra$foo_genes              # extras captured from the returned list
```

## 6. Combine and plot like any other method

```r
other <- Screen(ref[[2]], ref[[1]], ref[[3]],
                phenotype_class = "binary", screen_method = "scPP")
merged <- MergeResult(res, other)

ScreenFractionPlot(merged, group_by = "seurat_clusters",
                   screen_type = c("TP53", "scPP"))
```

Weighted voting (`weights =`) keys must match **both** a registered method name
and the corresponding metadata column, so it only applies when the label column
keeps the method name — i.e. when `label_type` is not overridden.

## Notes on the `mapper`

`mapper` receives the parameter list built by `Screen()` and must return a list.
Use it when your executor uses different argument names than the public API:

```r
mapper <- function(params) {
  params$top_n <- params$n_markers
  params$n_markers <- NULL
  params
}

Register(Foo = ScreenMethod(
  method_name = "Foo", method_version = "0.1.0", executor = Foo,
  phenotype_class = "binary", mapper = mapper
))
```
