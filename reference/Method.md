# Method: Abstract Base Class for Analysis Methods

`Method` is the abstract base class for all analysis methods (screening
and annotation) in the SigBridgeR package. It defines the common
interface that every method must implement, including a method name,
version, and executor function.

## Usage

``` r
Method(
  method_name,
  executor,
  method_version,
  sigbridger_version = get_pkg_version()
)
```

## Arguments

- method_name:

  `character`. A string identifying the method (e.g., `"Scissor"`,
  `"scPAS"`).

- executor:

  `function`. The function that performs the actual computation. This is
  the core workhorse of each method.

- method_version:

  `character`. A string indicating the version of the underlying
  software package.

- sigbridger_version:

  Package Version of SigBridgeR

## Value

A `Method` S7 object (abstract; cannot be instantiated directly).

## Details

This class is defined using the S7 object-oriented system and inherits
from
[SigBridgeRBase](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md).
It is marked as abstract, meaning it cannot be instantiated directly.
Concrete subclasses such as
[ScreenMethod](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethod.md)
and `AnnotationMethod` provide specific method categories.

## See also

Other SigBridgeR-Classes:
[`ScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethod.md),
[`SigBridgeRBase()`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md)
