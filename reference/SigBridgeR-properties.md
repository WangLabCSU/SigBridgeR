# SigBridgeR Base Class Properties

These are S7 property definitions used within the `SigBridgeRBase` class
and its subclasses. They define the types, validators, and default
values for class attributes. These properties are exported for community
developers to extend and subclass `SigBridgeRBase`.

## Usage

``` r
property_sigbridger_verison

property_verison

property_chr

property_phenotype_class

property_fn

property_list

property_mapper_fn

property_py_env

property_data_list

property_seurat

property_phenotype

property_matrix
```

## Format

An S7 property of class `class_character`.

An S7 property of class `class_character`.

An S7 property of class `class_character`.

An S7 property of class `class_character`.

An S7 property of class `class_function`.

An S7 property of class `class_list`.

An S7 property of class `class_any`.

An S7 property of class `class_any`.

An S7 property of class `class_list`.

An S7 property of class `class_any`.

An S7 property of class `class_any`.

An S7 property of class `class_any`.

## Details

A read-only character property that returns the current version of the
`SigBridgeR` package. Attempting to set this property will raise an
error.

The version string is validated to contain at least one dot (`.`). The
default value is obtained from
[`get_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/get_pkg_version.md).

A read-only character property that returns the current version of the
`SigBridgeR` package. Attempting to set this property will raise an
error.

The version string is validated to contain at least one dot (`.`). The
default value is obtained from
[`get_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/get_pkg_version.md).

A general-purpose character property with no additional validation. Can
be used for any string-typed attribute in subclasses.

A character property that specifies the phenotype classification
type(s). Valid values are `"binary"`, `"continuous"`, and `"survival"`.
Multiple values can be provided as a character vector.

Default: `c("binary", "continuous", "survival")`.

A property that accepts a function value. Used to store callback or
transformation functions in class attributes.

A named list property. The value must be a named list (not `NULL` or
`NA`) with unique names. Duplicate names are not allowed.

A property that accepts either a function or `NULL`. Used for storing
mapper/transformation functions that may be optionally specified.

If a value is provided, it must be a function; otherwise, `NULL` is
accepted to indicate no mapper is configured.

A property representing a Python environment configuration. Accepts
either a character string (e.g., a conda environment name or path) or
`NULL` (when no Python environment is specified).

A named list property for storing multiple data objects. The list must
be named (if not `NULL`) and names must be unique.

Unlike `property_list`, this property allows `NULL`, `NA`, and empty
lists ([`list()`](https://rdrr.io/r/base/list.html)) as valid values,
making it suitable for optional data storage.

A property that requires a
[Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object. The value must inherit from the `"Seurat"` class.

This property is used to store single-cell data in Seurat format for
downstream analysis.

A property that accepts a phenotype data

A property that accepts a matrix.

## See also

[SeuratObject::Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
