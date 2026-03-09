# Registry of Phenotype-Associated Cell Screening Methods

An environment storing methods for screening phenotype-associated cells.

## Usage

``` r
ScreenStrategy
```

## Format

An object of class `environment` of length 9.

## Details

Storing structure - named list

- `key`: method name

- `value`: list

  - `method_name`: method name

  - `executor`: function implementation of the method

  - `phenotypes`: The phenotype types supported by this method

  - `mapper`: A function that transforms the parameters passed to the
    `Screen` function before forwarding them to the executor; both input
    and output must be of type `list`.

## See also

Other Add_Screen_method:
[`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md),
[`TemplateScreenFunc()`](https://wanglabcsu.github.io/sigbridger/reference/TemplateScreenFunc.md),
[`ValidateScreenFunc()`](https://wanglabcsu.github.io/sigbridger/reference/ValidateScreenFunc.md)
