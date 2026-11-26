# ScreenMethodCache: Screening Method Cache Container

`ScreenMethodCache` bundles cached screening data with its configuration
and metadata. It serves as the primary container for reading from and
writing to the cache system.

## Usage

``` r
ScreenMethodCache(
  cache_path,
  cache_config_path,
  cache_data,
  screen_method_config,
  sigbridger_version = get_pkg_version()
)
```

## Arguments

- cache_path:

  `character`. Path to the cache directory where data files are stored.

- cache_config_path:

  `character`. Path to the `cache_config.json` metadata file. Must end
  with `"cache_config.json"`.

- cache_data:

  `list`. A named list of cached data objects. Can be `NULL` or empty
  (for load-only use).

- screen_method_config:

  A `ScreenMethodConfig` S7 object containing the screening method
  configuration.

- sigbridger_version:

  Package Version of SigBridgeR

## Value

A `ScreenMethodCache` S7 object.

## Details

This class is defined using the S7 object-oriented system and inherits
from
[SigBridgeRBase](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md).
It wraps a `ScreenMethodConfig` together with cached data objects, file
paths, and a read-only SigBridgeR version stamp.

### Validator

The class validator checks that:

- `cache_path` is an existing directory.

- The parent directory of `cache_config_path` exists.

- `cache_config_path` ends with `"cache_config.json"`.

## See also

Other S7-Classes:
[`AnnotationMethod()`](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md),
[`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
[`BenchConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md),
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md),
[`ScreenMethodConfig()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md),
[`ScreenMethodResult()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)

## Additional properties

- `@sigbridger_version`:

  `character` (read-only). The SigBridgeR package version at the time
  the cache was created. Automatically set from
  [`get_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/get_pkg_version.md)
  and cannot be modified by the user.
