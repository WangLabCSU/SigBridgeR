    library(SigBridgeR)
    library(Seurat)

# Overview

`BenchAcc()` compute the classification accuracy of screening methods
according to a known ground truth

It treats `"Positive"` as the positive class and every other supported
label (`"Negative"`, `"Neutral"`, `"Other"`) as the negative class. For
each screening method found in the object, it computes a full set of
binary classification metrics, so the methods can be ranked or compared
on exactly the same cells.

The function accepts three kinds of input:

- a `Seurat` object whose `meta.data` contains one column per screening
  method;
- a `ScreenMethodResult`, which is unwrapped to its `scRNA_data` slot;
- a `BenchAccConfig`, which runs every configured screening method first
  and then evaluates the merged labels in one pass.

# A small example

The example below uses a tiny object with six cells, a ground-truth
column and two simulated screening results. In a real analysis the
`Scissor` and `scPAS` columns would be produced by `Screen()` or
`MergeResult()`.

    counts <- Matrix::Matrix(
      matrix(
        seq_len(24L),
        nrow = 4L,
        dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:6))
      ),
      sparse = TRUE
    )

    sc <- CreateSeuratObject(counts = counts)

    sc$truth <- c(
      "Positive",
      "Negative",
      "Positive",
      "Other",
      "Positive",
      "Negative"
    )
    sc$Scissor <- c(
      "Positive",
      "Negative",
      "Positive",
      "Other",
      "Negative",
      "Negative"
    )
    sc$scPAS <- c(
      "Positive",
      "Positive",
      "Negative",
      "Other",
      "Negative",
      "Negative"
    )

    metrics <- BenchAcc(sc, bench_config = "truth", verbose = FALSE)

    round(do.call(rbind, metrics), 3)

Each element of the returned list is a one-row matrix whose row name is
the screening method. The same numbers are easier to compare after
binding them:

    metrics_dt <- data.table::rbindlist(
      lapply(metrics, as.data.frame),
      idcol = "method"
    )

    metrics_dt[, .(method, TPR, FPR, F1, Accuracy, Precision, Recall)]

# Specifying the ground truth

`bench_config` can be given in two equivalent forms:

    # * a length-one column name in sc@meta.data
    from_column <- BenchAcc(sc, bench_config = "truth", verbose = FALSE)

    # * a per-cell label vector, one label per cell
    from_vector <- BenchAcc(sc, bench_config = sc$truth, verbose = FALSE)

    identical(from_column, from_vector)

The labels must contain `"Positive"`, and the vector form must have one
element per cell. The function stops with an informative error when the
column is missing, when the labels have no positive class, or when the
vector length does not match the number of cells:

    BenchAcc(sc, bench_config = "missing_column", verbose = FALSE)

# Reading the metrics

With *P* positive and *N* negative ground-truth cells, and *T**P*,
*F**P*, *T**N*, *F**N* the usual confusion-matrix counts, `BenchAcc()`
reports:

<table>
<colgroup>
<col style="width: 27%" />
<col style="width: 41%" />
<col style="width: 31%" />
</colgroup>
<thead>
<tr>
<th>Metric</th>
<th>Definition</th>
<th>Meaning</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>TPR</code> / <code>Recall</code></td>
<td><span
class="math inline"><em>T</em><em>P</em>/(<em>T</em><em>P</em> + <em>F</em><em>N</em>)</span></td>
<td>Fraction of truly positive cells recovered</td>
</tr>
<tr>
<td><code>FPR</code></td>
<td><span
class="math inline"><em>F</em><em>P</em>/(<em>F</em><em>P</em> + <em>T</em><em>N</em>)</span></td>
<td>Fraction of negative cells incorrectly called positive</td>
</tr>
<tr>
<td><code>TNR</code></td>
<td><span
class="math inline"><em>T</em><em>N</em>/(<em>T</em><em>N</em> + <em>F</em><em>P</em>)</span></td>
<td>Fraction of negative cells correctly left out</td>
</tr>
<tr>
<td><code>FNR</code></td>
<td><span
class="math inline"><em>F</em><em>N</em>/(<em>F</em><em>N</em> + <em>T</em><em>P</em>)</span></td>
<td>Fraction of positive cells missed</td>
</tr>
<tr>
<td><code>Precision</code></td>
<td><span
class="math inline"><em>T</em><em>P</em>/(<em>T</em><em>P</em> + <em>F</em><em>P</em>)</span></td>
<td>Fraction of called cells that are truly positive</td>
</tr>
<tr>
<td><code>F1</code></td>
<td><span
class="math inline">2<em>T</em><em>P</em>/(2<em>T</em><em>P</em> + <em>F</em><em>P</em> + <em>F</em><em>N</em>)</span></td>
<td>Harmonic mean of precision and recall</td>
</tr>
<tr>
<td><code>Accuracy</code></td>
<td><span
class="math inline">(<em>T</em><em>P</em> + <em>T</em><em>N</em>)/(<em>P</em> + <em>N</em>)</span></td>
<td>Overall agreement with the ground truth</td>
</tr>
</tbody>
</table>

Because screening methods usually return only a small positive set,
`Accuracy` tends to be high and similar across methods. `TPR`/`Recall`,
`Precision` and `F1` are usually more informative when comparing
methods. A method that selects very few cells can have perfect precision
but low recall, while a permissive method can have high recall but many
false positives.

`BenchAcc()` only scores screening method columns. Other metadata
columns, including the ground-truth column itself, are ignored. Method
names are matched case-insensitively against the methods registered in
`ScreenStrategy`:

    sc$SCISSOR <- sc$Scissor
    sc$Scissor <- NULL

    names(BenchAcc(sc, bench_config = "truth", verbose = FALSE))

# Benchmarking the screening step

If the screening has not been run yet, wrap the inputs in a
`BenchAccConfig`. The configuration stores the bulk data, the
single-cell object, the phenotype and the methods to run. `BenchAcc()`
then screens the data with each method and evaluates the accumulated
labels:

    config <- BenchAccConfig(
      matched_bulk = bulk,
      sc_data = sc,
      phenotype = phenotype,
      label_type = "{method_name}",
      phenotype_class = "binary",
      screen_method = c("Scissor", "scPAS"),
      method_param = MethodParamSet(
        Scissor_param = list(alpha = 0.05),
        scPAS_param = list(alpha = 0.01)
      ),
      bench_config = "truth"
    )

    acc <- BenchAcc(config, verbose = FALSE)

The `bench_config` of a `BenchAccConfig` is validated when the object is
created, so an invalid ground truth fails early instead of after a long
screening run. `method_param` is the place for method-specific
arguments; see `?MethodParamSet` for the full list of fields.

# Practical notes

- **Use the same cells for every method.** `BenchAcc()` compares the
  labels cell by cell. Merge the screening results first (for example
  with `MergeResult()`) so that every method shares the same cell names
  and order.
- **Use an independent ground truth.** A label derived from the same
  bulk phenotype that drove the screening is not an independent
  validation; use marker-based annotation, an external cell-type
  reference or held-out experimental labels when possible.
- **Report more than accuracy.** `TPR`, `Precision` and `F1` describe
  the trade-off between missing positive cells and including negative
  ones, which is usually the relevant question for phenotype-associated
  cell screening.
- **Only the four supported labels are accepted.** Predicted labels must
  be one of `"Positive"`, `"Negative"`, `"Neutral"` or `"Other"`; typos
  and `NA` abort with the offending column and row.

# See also

- `BenchGSEA()` for comparing the functional signal of the selected
  cells.
- `BenchGeneScore()` for comparing pathway activity at the single-cell
  level.
- `BenchMemory()` for comparing run time and memory.
- `BenchAccConfig()`, `MethodParamSet()` and `?Screen` for the
  configuration and screening interfaces.
