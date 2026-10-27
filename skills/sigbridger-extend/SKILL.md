---
name: sigbridger-extend
description: Add custom algorithms to SigBridgeR without forking it — scaffold a screening function with TemplateScreenFunc(), validate it with ValidateScreenFunc(), record parameters as a ScreenMethodConfig, wrap it in ScreenMethod and register it via Register() so Screen() can dispatch to it; also register annotation methods (AnnotationMethod into SCAnnotateStrategy) and Seurat preprocessing steps (functions into SCPreProcessStrategy). Use when the user wants to plug a new screening/annotation method into SigBridgeR or asks how Register()/ScreenStrategy work.
---

# Extending SigBridgeR

## Purpose

SigBridgeR dispatches through three S7-backed registries, all writable at runtime:

| Registry | Holds | Dispatched by |
| --- | --- | --- |
| `ScreenStrategy` | `ScreenMethod` objects | `Screen()` |
| `SCAnnotateStrategy` | `AnnotationMethod` objects | `SCAnnotate()` |
| `SCPreProcessStrategy` | plain functions | `SCIntegrate(pipeline=)` / deprecated `SCPreProcess()` |

`Register()` is the single entry point; it chooses the registry from the class of
the value it receives.

## Adding a screening method

### 1. Scaffold

```r
TemplateScreenFunc(filename = "R/foo.R", func_name = "Foo", documentation = TRUE)
```

`filename = "current"` appends to the active Positron/RStudio script (requires
`rstudioapi`). `append = FALSE` overwrites after confirmation.

### 2. Implement the fixed contract

```r
Foo <- function(matched_bulk, sc_data, phenotype,
                label_type = NULL,
                phenotype_class = c("binary", "survival", "continuous"),
                ...) {
  phenotype_class <- rlang::arg_match(phenotype_class)
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% TRUE

  # ... algorithm ...

  list(
    scRNA_data = sc_data,          # REQUIRED: Seurat with labels in meta.data
    foo_scores = scores            # optional extras -> ScreenMethodResult$extra
  )
}
```

Non-negotiable rules:

- The first five argument names must not change — `ScreenMethod`'s validator
  checks the executor's formals.
- `phenotype` is a named vector for `binary` / `continuous`, and a
  `data.frame(time, status)` for `survival`; names/row names must equal
  `colnames(matched_bulk)`.
- The returned list **must** contain `scRNA_data`, a `Seurat` object.
- Prefix new metadata columns with the method name to avoid collisions.
- Take `label_type` and store it so downstream merge/plot code can find the
  label column.

### 3. Validate

```r
ValidateScreenFunc(Foo)
```

Runs five checks (inputs, verbose handling, syntax, return value, directory
naming) and either returns `TRUE` or aborts with a per-check report.

### 4. Record parameters

```r
config <- ScreenMethodConfig(
  method_name     = "Foo",                  # after registration: a ScreenStrategy key
  param           = list(cutoff = 1L),
  method_version  = "0.1.0",
  phenotype_class = c("binary", "continuous", "survival"),
  label_type      = "TP53"
)
```

`param` must be a non-empty named list. Store the config in `@misc` so the run is
reproducible (see `sigbridger-cache-config`).

### 5. Wrap and register

```r
Register(Foo = ScreenMethod(
  method_name     = "Foo",
  method_version  = "0.1.0",
  executor        = Foo,
  phenotype_class = c("binary", "continuous", "survival"),
  mapper          = NULL            # optional list -> list parameter transform
))
```

Then `Screen(..., screen_method = "Foo")` works immediately;
`names(ScreenStrategy)` lists it. Re-register with `overwrite = TRUE`; otherwise
duplicate keys abort. `mapper` is applied to the parameter list before the
executor and must return a list — use it to rename/derive method arguments.

## Adding an annotation method

```r
Register(MyAnno = AnnotationMethod(
  method_name = "MyAnno",
  executor    = MyAnnoFn,     # function(sc, ...) -> Seurat
  pkg_name    = "myannopkg",  # version resolved from R, or from Python if py_env set
  py_env      = NULL
))

SCAnnotate(sc, method = "MyAnno", ...)
```

`AnnotationMethod` extends `Method` with `py_env`: when a conda environment name
or Python path is given, `method_version` is resolved with `py_pkg_version()`
instead of `r_pkg_version()`. Built-in methods (`SingleR`, `CellTypist`,
`mLLMCelltype`) are registered the same way — inspect them with
`names(SCAnnotateStrategy)`.

## Adding a Seurat preprocessing step

```r
Register(h = Seurat::RunHarmony)      # key "h" becomes usable in a pipeline string
```

The value is a plain function; `RegisterImpl` validates it as a Seurat-level
function and stores it in `SCPreProcessStrategy`, whose built-in keys are:

`o` CreateSeuratObject · `n` NormalizeData · `v` FindVariableFeatures ·
`s` ScaleData · `p` RunPCA · `c` FindClusters · `e` FindNeighbors · `t` RunTSNE ·
`u` RunUMAP · `r` SCTransform · `i` IntegrateLayers

Pipelines are strings, e.g. `"nsvpiectu"` for `SCIntegrate()`.

## Deprecated entry points

`RegisterScreenMethod()`, `RegisterAnnoMethod()` and `RegisterSeuratMethod()` are
deprecated since 4.0.0 and emit lifecycle warnings. Migrate to `Register()`:
the name of each argument becomes the registry key.

## Verification checklist

```r
names(ScreenStrategy)          # method visible
ValidateScreenFunc(Foo)        # contract satisfied
Screen(bulk, sc, pheno, phenotype_class = "binary", screen_method = "Foo")
```

## Reference

- `?SigBridgeR::TemplateScreenFunc`, `?SigBridgeR::ValidateScreenFunc`,
  `?SigBridgeR::ScreenMethod`, `?SigBridgeR::AnnotationMethod`,
  `?SigBridgeR::Register`, `?SigBridgeR::ScreenStrategy`,
  `?SigBridgeR::SCPreProcessStrategy`
- Vignettes `Extending_ScreenMethod`, `Extending_Annotation`, `Extending_Seurat`
- Worked end-to-end example: `references/custom-method-example.md`
