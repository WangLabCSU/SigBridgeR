# Registry of Cell Type Annotation Methods

An environment storing methods for annotating cell types.

## Usage

``` r
SCAnnotateStrategy
```

## Format

An object of class `environment` of length 3.

## Details

Storing structure - named list

- `key`: method name

- `value`: list

  - `method_name`: method name

  - `executor`: function implementation of the method

## See also

Other Single_Cell_Annotation_Method:
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md),
[`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)
