# ScreenMethod: Screening Method Class for Phenotype-Associated Cell Detection

`ScreenMethod` is the S7 class for phenotype-associated cell screening
methods. It extends
[Method](https://wanglabcsu.github.io/sigbridger/reference/Method.md)
with phenotype classification and mapper function support, and enforces
a standardized function signature for all screening executors.

## Usage

``` r
ScreenMethod(
  method_name,
  method_version,
  executor,
  phenotype_class,
  mapper = NULL,
  sigbridger_version = get_pkg_version()
)
```

## Arguments

- method_name:

  `character`. (Inherited from
  [Method](https://wanglabcsu.github.io/sigbridger/reference/Method.md))
  A string identifying the screening method.

- method_version:

  `character`. (Inherited from
  [Method](https://wanglabcsu.github.io/sigbridger/reference/Method.md))
  A string indicating the version of the underlying software package.

- executor:

  `function`. (Inherited from
  [Method](https://wanglabcsu.github.io/sigbridger/reference/Method.md))
  The screening function. Must accept `matched_bulk`, `sc_data`,
  `phenotype`, `label_type`, and `phenotype_class` as arguments.

- phenotype_class:

  `character`. A character vector specifying the phenotype
  classification type(s) the method supports. Valid values: `"binary"`,
  `"continuous"`, `"survival"`. Default: all three.

- mapper:

  `function` or `NULL`. An optional function used to transform or map
  parameters before screening. When `NULL`, no mapper is applied.

- sigbridger_version:

  Package Version of SigBridgeR

## Value

A `ScreenMethod` S7 object.

## Details

This class inherits from
[Method](https://wanglabcsu.github.io/sigbridger/reference/Method.md)
and adds screening-specific properties. It validates that the executor
function accepts the required arguments: `matched_bulk`, `sc_data`,
`phenotype`, `label_type`, and `phenotype_class`. This ensures a
consistent interface across all screening algorithms (e.g., Scissor,
scPAS, scPP, TiRank).

## See also

Other SigBridgeR-Classes:
[`Method()`](https://wanglabcsu.github.io/sigbridger/reference/Method.md),
[`SigBridgeRBase()`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md)
