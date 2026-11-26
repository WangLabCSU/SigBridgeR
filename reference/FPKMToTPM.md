# Convert bulk RNA-seq FPKM to TPM

Converts a gene-level FPKM matrix to TPM (Transcripts Per Million).
Because the per-gene length factor is constant within a sample, the
FPKM-to-TPM transformation requires no gene lengths: each column is
scaled so that its values sum to \\10^6\\.

For a given sample, every FPKM value is multiplied by \\10^6 /
\sum(\text{FPKM})\\. Samples whose FPKM column sums to zero produce a
zero column in the result and trigger a warning.

Missing values (`NA`) are either replaced by zero (`na_as_zero = TRUE`,
the default) or trigger an error (`na_as_zero = FALSE`).

## Usage

``` r
FPKMToTPM(data, na_as_zero = TRUE, verbose = TRUE, ...)
```

## Arguments

- data:

  A two-dimensional matrix-like object of FPKM values (a base matrix, a
  data frame, or an S4 `Matrix`). Rows are genes, columns are samples.
  Data frames are coerced and must contain only numeric columns. Values
  must be non-negative and finite.

- na_as_zero:

  A logical flag. If `TRUE` (default), missing values are replaced by
  zero before normalization. If `FALSE`, an error is raised when `data`
  contains missing values.

- verbose:

  A logical flag. If `TRUE` (default), diagnostic messages
  (missing-value report and completion summary) are shown.

- ...:

  Additional arguments passed to
  [`CheckNA()`](https://wanglabcsu.github.io/sigbridger/reference/CheckNA.md)
  when `data` contains missing values and `verbose` is `TRUE` (for
  example `max_print`).

## Value

A numeric TPM matrix with the same dimensions and `dimnames` as `data`.
Each column sums to \\10^6\\ unless the corresponding FPKM column sums
to zero.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create an FPKM matrix (3 genes x 2 samples)
fpkm <- matrix(
  c(10, 20, 30, 15, 25, 35),
  nrow = 3, ncol = 2,
  dimnames = list(c("GENE1", "GENE2", "GENE3"), c("Sample1", "Sample2"))
)

# Convert FPKM to TPM
tpm <- FPKMToTPM(fpkm)
colSums(tpm)  # Each column should sum to 1e6

# Error if missing values are present and na_as_zero = FALSE
fpkm_na <- fpkm
fpkm_na[1, 1] <- NA
tpm2 <- FPKMToTPM(fpkm_na, na_as_zero = FALSE)
} # }
```
