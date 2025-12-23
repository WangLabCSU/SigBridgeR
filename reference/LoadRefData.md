# Download & Load Reference Data

This function checks if the data already exists in a local cache. If
not, it downloads the file from the remote repository using multiple
sources with fallback.

## Usage

``` r
LoadRefData(
  data_type = c("survival", "binary", "continuous"),
  path = NULL,
  timeout = SigBridgeRUtils::getFuncOption("timeout"),
  ...
)
```

## Arguments

- data_type:

  The type of data to download. Must be one of "survival", "binary", or
  "continuous".

- path:

  Optional path to save the downloaded file, default: NULL, saving in
  package.

- timeout:

  Integer. Connection timeout in seconds (default: 60)

- ...:

  Additional arguments. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

## Value

The requested datasets, stored in a list.

## Examples

``` r
if (FALSE) { # \dontrun{
ref_data <- LoadRefData(data_type = c("survival"))
} # }
```
