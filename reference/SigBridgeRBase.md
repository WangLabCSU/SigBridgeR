# SigBridgeRBase: Abstract Base Class for SigBridgeR

`SigBridgeRBase` is the abstract base class for all S7 classes in the
SigBridgeR package. It provides a unified foundation for S3 method
dispatch and ensures consistent behavior across all subclasses.

## Usage

``` r
SigBridgeRBase(sigbridger_version = "4.0.0")
```

## Arguments

- sigbridger_version:

  Package Version of SigBridgeR

## Details

This class is defined using the S7 object-oriented system and is marked
as abstract, meaning it cannot be instantiated directly. All other
classes in the SigBridgeR package should inherit from `SigBridgeRBase`.

The primary purpose of this base class is to serve as a dispatch target
for base R S3 generic functions, enabling polymorphic behavior across
all SigBridgeR S7 objects through S3 method dispatch.

## See also

Other SigBridgeR-Classes:
[`Method()`](https://wanglabcsu.github.io/sigbridger/reference/Method.md),
[`ScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethod.md)
