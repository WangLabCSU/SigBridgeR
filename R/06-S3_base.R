# nocov start

#' Base R Methods for SigBridgeR Objects
#'
#' @description
#' S3 methods that make every `SigBridgeRBase` object -- and therefore every
#' SigBridgeR object, such as `ScreenMethod`, `ScreenMethodResult`, and
#' `InterceptStrategy` -- behave like a standard R object. They support the
#' base generics `print()`, `length()`, `names()`, `names<-()`, `format()`,
#' `as.list()`, and `as.character()`.
#'
#' @details
#' All SigBridgeR objects inherit from the abstract base class
#' [SigBridgeRBase] and store their data in S7 properties. The methods
#' documented on this page expose those properties through familiar base R
#' interfaces:
#'
#' - `print()` displays a human-readable summary of the object. The class
#'   name is shown in bold blue and each property is listed with its name,
#'   value type, and value.
#' - `length()` returns the number of properties stored in the object.
#' - `names()` returns the names of the object's properties.
#' - `names<-()` is restricted: properties are validated by their S7
#'   definitions, so modifying them in place is not allowed. Reconstruct a
#'   new object instead.
#' - `format()` converts the object to a compact character vector, starting
#'   with a header of the form `<Class: SigBridgeR <version>>`.
#' - `as.list()` converts the object to a named list of its properties.
#'
#' `as.character()` is intentionally unsupported: coercing a
#' `SigBridgeRBase` object to a character vector is not meaningful and
#' always throws an error. Use [format()] instead.
#'
#' @seealso
#' - [SigBridgeRBase] for the abstract base class on which these methods
#'   dispatch.
#' - [S7::props()] for the underlying property accessor.
#'
#' @name SigBridgeR-S3
NULL

#' Print a SigBridgeR object
#'
#' Pretty-print a `SigBridgeR` object using the `cli` package.
#'
#' The class name is displayed in bold blue. Slots are displayed as an
#' unordered list, with slot names in regular blue and each value preceded by
#' its type.
#'
#' @param x A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#' @param ... Additional arguments passed to the print method. Currently unused.
#'
#' @return Invisibly returns `x`.
#'
#' @rdname SigBridgeR-S3
#' @rawNamespace S3method(print,"SigBridgeR::SigBridgeRBase")
#' @export
`print.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  cls_x <- gsub(".*:", "", class(x)[1L])
  properties <- S7::props(x)
  nms_properties <- names(properties)

  # src/list_to_chr.cpp
  print_sigbridger_base_cpp(cls_x, nms_properties, properties)

  invisible(x)
}

#' Summary of a SigBridgeR object
#'
#' The class name is displayed in bold blue. Slots are displayed as an
#' unordered list, with slot names in regular blue and each value preceded by
#' its type.
#'
#' @param object A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#' @param ... Additional arguments passed to the print method. Currently unused.
#'
#' @return Invisibly returns `object`.
#'
#' @rdname SigBridgeR-S3
#' @rawNamespace S3method(summary,"SigBridgeR::SigBridgeRBase")
#' @export
`summary.SigBridgeR::SigBridgeRBase` <- function(object, ...) {
  `print.SigBridgeR::SigBridgeRBase`(x = object, ...)
}


#' Number of Properties in a SigBridgeR Object
#'
#' @description
#' Returns the number of properties (slots) stored in a `SigBridgeRBase` object.
#'
#' @param x A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#'
#' @return An integer: the number of properties.
#'
#' @rdname SigBridgeR-S3
#' @export
#' @rawNamespace S3method(length,"SigBridgeR::SigBridgeRBase")
`length.SigBridgeR::SigBridgeRBase` <- function(x) {
  length(S7::props(x))
}

#' Names of Properties in a SigBridgeR Object
#'
#' @description
#' Returns the names of all properties (slots) stored in a `SigBridgeRBase`
#' object.
#'
#' @param x A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#'
#' @return A character vector of property names.
#'
#' @rdname SigBridgeR-S3
#' @export
#' @rawNamespace S3method(names,"SigBridgeR::SigBridgeRBase")
`names.SigBridgeR::SigBridgeRBase` <- function(x) {
  S7::prop_names(x)
}

#' Names of Properties in a SigBridgeR Object
#'
#' @description
#' Returns the names of all properties (slots) stored in a `SigBridgeRBase`
#' object.
#'
#' @param x A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#' @param value A character vector of property names.
#'
#' @return A character vector of property names.
#'
#' @rdname SigBridgeR-S3
#' @export
#' @rawNamespace S3method(`names<-`,"SigBridgeR::SigBridgeRBase")
`names<-.SigBridgeR::SigBridgeRBase` <- function(x, value) {
  Abort(
    "Modification is restricted",
    "Please reconstruct a new object",
    type = "[RESTRICTED]"
  )
}

#' Format a SigBridgeR base object
#'
#' Format a `SigBridgeRBase` object as a compact character vector.
#'
#' The first line is displayed as `<Class: name version>` when the object has
#' `method_name` and `method_version` properties, for example
#' `<ScreenMethod: SigBridgeR 4.0.0>`. Each following line contains one property name
#' and its value on the same line. Property types are not shown.
#'
#' @param x A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#' @param ... Additional arguments passed to the format method. Currently unused.
#'
#' @return A character vector.
#'
#' @rdname SigBridgeR-S3
#' @export
#' @rawNamespace S3method(format,"SigBridgeR::SigBridgeRBase")
`format.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  cls_x <- sub("^.*::", "", class(x)[1L])

  properties <- S7::props(x)
  nms_properties <- names(properties)

  value_to_chr <- function(value) {
    if (is.function(value)) {
      toString(rlang::fn_fmls_names(value))
    } else if (is.null(value)) {
      "NULL"
    } else if (length(value) == 0L) {
      ""
    } else if (is.atomic(value)) {
      toString(as.character(value))
    } else {
      toString(class(value))
    }
  }

  header <- paste0("<", cls_x, ":", " SigBridgeR ", get_pkg_version(), ">")

  property_lines <- vapply(
    nms_properties,
    function(nm) {
      paste0("- ", nm, ": ", value_to_chr(properties[[nm]]))
    },
    character(1L)
  )

  paste(c(header, unname(property_lines)), sep = "\n", collapse = "\n")
}

#' Coerce a SigBridgeR Object to a List
#'
#' @description
#' Converts a `SigBridgeRBase` object into a named list of its properties.
#'
#' @param x A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#' @param ... Additional arguments passed to the method. Currently unused.
#'
#' @return A named list of properties.
#'
#' @rdname SigBridgeR-S3
#' @export
#' @rawNamespace S3method(as.list,"SigBridgeR::SigBridgeRBase")
`as.list.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  S7::props(x)
}

#' Coerce a SigBridgeR Object to a Character Vector
#'
#' @description
#' Coercion of `SigBridgeRBase` to character is intentionally unsupported.
#' Use [format()] to obtain a character representation instead.
#'
#' @param x A `SigBridgeRBase` object. Basically all SigBridgeR objects.
#' @param ... Additional arguments passed to the method. Currently unused.
#'
#' @return Throws an error; never returns.
#'
#' @rdname SigBridgeR-S3
#' @export
#' @rawNamespace S3method(as.character,"SigBridgeR::SigBridgeRBase")
`as.character.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  Abort(
    "cannot coerce type {.cls SigBridgeRBase} to vector of type {.cls character}",
    "Maybe you want to use {.fun format}?",
    type = "[TYPE ERROR]"
  )
}

# nocov end
