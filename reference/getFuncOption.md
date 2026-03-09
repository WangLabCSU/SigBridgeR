# Configuration Functions for SigBridgeR Package

These functions provide a centralized configuration system for the
SigBridgeR package, allowing users to set and retrieve package-specific
options with automatic naming conventions.

setFunOption sets one or more configuration options for the SigBridgeR
package. Options are automatically prefixed with "SigBridgeR." if not
already present, ensuring proper namespace isolation.

getFuncOption retrieves configuration options for the SigBridgeR
package. The function automatically handles the "SigBridgeR." prefix,
allowing users to reference options with or without the explicit prefix.

checkFuncOption validates configuration options for the SigBridgeR
package. This function ensures that all options meet type and value
requirements before they are set in the global options.

## Usage

``` r
getFuncOption(option = NULL, default = NULL)
```

## Arguments

- option:

  Character string specifying the option name to check

- default:

  The default value to return if the option is not set. Defaults to
  NULL.

## Value

Invisibly returns NULL. The function is called for its side effects of
setting package options.

The value of the specified option if it exists, otherwise the provided
default value.

No return value. The function throws an error if the value doesn't meet
the required specifications for the given option.

## Details

This function performs the following validations:

- verbose:

  Must be single logical values (TRUE/FALSE)

- timeout, seed:

  Must be single integer values

The function is automatically called by setFuncOption to ensure all
configuration options are valid before they are set.

## Available Opions

- `verbose`: A logical value indicating whether to print verbose
  messages, defaults to `TRUE`

- `timeout`: An integer specifying the timeout in seconds for parallel
  processing, defaults to \`180L“

- `seed`: An integer specifying the random seed for reproducible
  results, defaults to \`123L“

## Examples

``` r
# Set options with automatic prefixing
setFuncOption(verbose = TRUE)

# Retrieve options with automatic prefixing
getFuncOption("verbose")
#> [1] TRUE
```
