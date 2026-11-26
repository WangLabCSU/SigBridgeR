# Stage Benchmark Inputs in a Temporary Directory

Serialises the three shared benchmark inputs into a fresh temporary
directory with `qs2`, so that each measurement child process can load
its own copies instead of sharing mutable objects with the parent.

## Usage

``` r
bench_tmp_dir(sc_data, bulk_data, pheno_data)
```

## Arguments

- sc_data:

  A `Seurat` object holding the cells to screen.

- bulk_data:

  The matched bulk data passed to
  [`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md).

- pheno_data:

  The phenotype data passed to
  [`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md).

## Value

The path of the created directory, containing `sc_data.qs2`,
`matched_bulk.qdata` and `phenotype.qdata` and the number of CPU cores.
The caller is responsible for removing it, as
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
does via [`on.exit()`](https://rdrr.io/r/base/on.exit.html).

## Details

Serialisation uses half of the detected CPU cores (at least one) as the
`qs2` thread count.
