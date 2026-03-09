# Preprocess Phenotype Data

Aligns and validates phenotype data with bulk RNA-seq expression matrix.
Supports three phenotype types (binary, continuous, survival) with
automatic sample matching, type validation, and optional conditional
mapping via
[`PhenoMap`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md).

## Usage

``` r
PhenoPreProcess(
  bulk,
  phenotype,
  phenotype_class = c("binary", "continuous", "survival"),
  ...,
  select = NULL,
  verbose = getFuncOption("verbose")
)
```

## Arguments

- bulk:

  A two-dimensional matrix or data frame of bulk RNA-seq expression data
  with genes as rows and samples as columns. Must have at least 2
  samples.

- phenotype:

  Phenotype data, either:

  - Named numeric vector (for binary/continuous phenotypes)

  - Data frame/matrix with row names matching `colnames(bulk)` (for
    survival or multi-column phenotypes)

- phenotype_class:

  Character. Type of phenotype:

  `"binary"`

  :   Two-class categorical outcome (e.g., case/control). Must have
      exactly 2 unique values.

  `"continuous"`

  :   Continuous measurement (e.g., age, expression). Must have \>2
      unique values.

  `"survival"`

  :   Survival data with time and status columns. Automatically guesses
      columns named "time" and "status/censor" if `select` is `NULL`.

  Partial matching is supported.

- ...:

  Conditional mapping rules passed to
  [`PhenoMap`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md)
  for transforming phenotype values before validation. Format:
  `condition ~ value`. Example:
  `status == "death" ~ 1, status == "alive" ~ 0`.

- select:

  Character. Column name(s) to select from `phenotype` when it is
  two-dimensional. Required for survival data with multiple columns, or
  to specify which column to use for binary/continuous phenotypes.

- verbose:

  Logical. Whether to print diagnostic messages. Default: inherits from
  `getOption("SigBridgeRUtils.verbose")`.

## Value

Preprocessed phenotype data:

- For binary/continuous: Named numeric vector with sample names as names

- For survival: Data frame with two columns (time, status) and sample
  names as row names

Only samples present in both `bulk` and `phenotype` are retained.

## Validation Rules

- **Binary**: Must have exactly 2 unique values (e.g., 0/1, TRUE/FALSE)

- **Continuous**: Must have \>2 unique values (to avoid confusion with
  binary)

- **Survival**: Must have exactly 2 unique values in status column
  (typically 0 = censored, 1 = event)

- **Sample matching**: Common samples identified via
  `intersect(colnames(bulk), names/rownames(phenotype))`

## Automatic Features

- **Survival column guessing**: When `phenotype_class = "survival"` and
  `select = NULL`, automatically detects columns containing "time" and
  "status"/"censor" in their names (case-insensitive).

- **Type conversion**: Logical values are converted to numeric
  (TRUE-\>1, FALSE-\>0)

- **Sample alignment**: Returns only samples present in both `bulk` and
  `phenotype`

## See also

[`PhenoMap()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md)

Other input_preprocess:
[`BulkPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/BulkPreProcess.md),
[`PhenoMap()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Binary phenotype
bulk_data <- matrix(rpois(100, 10), nrow = 20, ncol = 5)
colnames(bulk_data) <- paste0("Sample", 1:5)

pheno_binary <- c(Sample1 = 1, Sample2 = 0, Sample3 = 1, Sample4 = 0, Sample5 = 1)

result <- PhenoPreProcess(
  bulk = bulk_data,
  phenotype = pheno_binary,
  phenotype_class = "binary"
)

# Example 2: Continuous phenotype with discretization
pheno_age <- c(Sample1 = 25, Sample2 = 35, Sample3 = 45, Sample4 = 55, Sample5 = 65)

result <- PhenoPreProcess(
  bulk = bulk_data,
  phenotype = pheno_age,
  phenotype_class = "continuous",
  .x < 40 ~ "Young", .x >= 40 ~ "Old"
)

# Example 3: Survival data with automatic column detection
pheno_surv <- data.frame(
  time = c(12, 24, 18, 36, 30),
  status = c(1, 0, 1, 1, 0),
  row.names = paste0("Sample", 1:5)
)

result <- PhenoPreProcess(
  bulk = bulk_data,
  phenotype = pheno_surv,
  phenotype_class = "survival"
)

# Example 4: Survival data with explicit column names
pheno_surv_custom <- data.frame(
  follow_up = c(12, 24, 18, 36, 30),
  event = c(1, 0, 1, 1, 0),
  row.names = paste0("Sample", 1:5)
)

result <- PhenoPreProcess(
  bulk = bulk_data,
  phenotype = pheno_surv_custom,
  phenotype_class = "survival",
  select = c("follow_up", "event")
)
} # }
```
