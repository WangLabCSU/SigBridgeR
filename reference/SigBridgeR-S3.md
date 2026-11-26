# Base R Methods for SigBridgeR Objects

S3 methods that make every `SigBridgeRBase` object – and therefore every
SigBridgeR object, such as `ScreenMethod`, `ScreenMethodResult`, and
`InterceptStrategy` – behave like a standard R object. They support the
base generics [`print()`](https://rdrr.io/r/base/print.html),
[`length()`](https://rdrr.io/r/base/length.html),
[`names()`](https://rdrr.io/r/base/names.html), `names<-()`,
[`format()`](https://rdrr.io/r/base/format.html),
[`as.list()`](https://rdrr.io/r/base/list.html), and
[`as.character()`](https://rdrr.io/r/base/character.html).

Pretty-print a `SigBridgeR` object using the `cli` package.

The class name is displayed in bold blue. Slots are displayed as an
unordered list, with slot names in regular blue and each value preceded
by its type.

Returns the number of properties (slots) stored in a `SigBridgeRBase`
object.

Returns the names of all properties (slots) stored in a `SigBridgeRBase`
object.

Returns the names of all properties (slots) stored in a `SigBridgeRBase`
object.

Format a `SigBridgeRBase` object as a compact character vector.

Converts a `SigBridgeRBase` object into a named list of its properties.

Coercion of `SigBridgeRBase` to character is intentionally unsupported.
Use [`format()`](https://rdrr.io/r/base/format.html) to obtain a
character representation instead.

## Usage

``` r
# S3 method for class '`SigBridgeR::SigBridgeRBase`'
print(x, ...)

# S3 method for class '`SigBridgeR::SigBridgeRBase`'
summary(object, ...)

# S3 method for class '`SigBridgeR::SigBridgeRBase`'
length(x)

# S3 method for class '`SigBridgeR::SigBridgeRBase`'
names(x)

# S3 method for class '`SigBridgeR::SigBridgeRBase`'
names(x) <- value

# S3 method for class '`SigBridgeR::SigBridgeRBase`'
format(x, ...)

# S3 method for class '`SigBridgeR::SigBridgeRBase`'
as.list(x, ...)

# S3 method for class '`SigBridgeR::SigBridgeRBase`'
as.character(x, ...)
```

## Arguments

- x:

  A `SigBridgeRBase` object. Basically all SigBridgeR objects.

- ...:

  Additional arguments passed to the method. Currently unused.

- object:

  A `SigBridgeRBase` object. Basically all SigBridgeR objects.

- value:

  A character vector of property names.

## Value

Invisibly returns `x`.

Invisibly returns `object`.

An integer: the number of properties.

A character vector of property names.

A character vector of property names.

A character vector.

A named list of properties.

Throws an error; never returns.

## Details

All SigBridgeR objects inherit from the abstract base class
[SigBridgeRBase](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md)
and store their data in S7 properties. The methods documented on this
page expose those properties through familiar base R interfaces:

- [`print()`](https://rdrr.io/r/base/print.html) displays a
  human-readable summary of the object. The class name is shown in bold
  blue and each property is listed with its name, value type, and value.

- [`length()`](https://rdrr.io/r/base/length.html) returns the number of
  properties stored in the object.

- [`names()`](https://rdrr.io/r/base/names.html) returns the names of
  the object's properties.

- `names<-()` is restricted: properties are validated by their S7
  definitions, so modifying them in place is not allowed. Reconstruct a
  new object instead.

- [`format()`](https://rdrr.io/r/base/format.html) converts the object
  to a compact character vector, starting with a header of the form
  `<Class: SigBridgeR <version>>`.

- [`as.list()`](https://rdrr.io/r/base/list.html) converts the object to
  a named list of its properties.

[`as.character()`](https://rdrr.io/r/base/character.html) is
intentionally unsupported: coercing a `SigBridgeRBase` object to a
character vector is not meaningful and always throws an error. Use
[`format()`](https://rdrr.io/r/base/format.html) instead.

The class name is displayed in bold blue. Slots are displayed as an
unordered list, with slot names in regular blue and each value preceded
by its type.

The first line is displayed as `<Class: name version>` when the object
has `method_name` and `method_version` properties, for example
`<ScreenMethod: SigBridgeR 4.0.0>`. Each following line contains one
property name and its value on the same line. Property types are not
shown.

## See also

- [SigBridgeRBase](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md)
  for the abstract base class on which these methods dispatch.

- [`S7::props()`](https://rconsortium.github.io/S7/reference/props.html)
  for the underlying property accessor.
