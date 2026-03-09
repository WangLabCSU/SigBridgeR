# Is 1 Much More Than 0

A function to check whether Y is highly skewed. The expected ideal
distribution is p1 : p0 = cutoff : (1 - cutoff). If the deviation
exceeds n_sd times the standard deviation, it is considered skewed,
i.e., the input phenotype raw data is skewed, which means the
reliability test may be unreliable or error-prone.

## Usage

``` r
is_skewed_dynamic(x, target = 0, expected_p = 0.8, n_sd = 4L)
```
