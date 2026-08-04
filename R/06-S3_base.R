#' Print a ScreenMethod object
#'
#' Pretty-print a `ScreenMethod` object using the `cli` package.
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
#' @family SigBridgeR-S3
#' @rawNamespace S3method(print,"SigBridgeR::SigBridgeRBase")
#' @export
`print.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  cls_x <- gsub(".*:", "", class(x)[1])
  properties <- props(x)
  nms_properties <- names(properties)

  cli::cli_div(
    theme = list(
      span.red = list(color = "red"),
      span.blue = list(color = "blue"),
      span.orange = list(color = "orange"),
      span.purple = list(color = "purple"),
      span.green = list(color = "green"),
      span.magenta = list(color = "magenta"),
      span.cyan = list(color = "cyan"),
      span.bold = list("font-weight" = "bold"),
      span.typeof = list(color = "grey")
    )
  )

  cli::cli_text("{.cls {cls_x}}")

  for (i in seq_along(properties)) {
    name <- nms_properties[i]
    key <- properties[[i]]

    cls_key <- class(key)[1] # class(NULL) -> NULL

    key <- if (is_function(key)) {
      fn_fmls_names(key)
    } else if (is.list(key)) {
      as.character.list(key)
    } else if (typeof(key) %chin% c("S3", "S4", "S7", "R6")) {
      class(key)
    } else {
      key
    }

    cli::cli_li(
      "{.blue {name}} {.typeof [{cls_key}]}: {key}"
    )
  }

  cli::cli_end()
  cli::cli_end()

  invisible(x)
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
#' @family SigBridgeR-S3
#' @export
#' @rawNamespace S3method(length,"SigBridgeR::SigBridgeRBase")
`length.SigBridgeR::SigBridgeRBase` <- function(x) {
  length(props(x))
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
#' @family SigBridgeR-S3
#' @export
#' @rawNamespace S3method(names,"SigBridgeR::SigBridgeRBase")
`names.SigBridgeR::SigBridgeRBase` <- function(x) {
  prop_names(x)
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
#' @family SigBridgeR-S3
#' @export
#' @rawNamespace S3method(`names<-`,"SigBridgeR::SigBridgeRBase")
`names<-.SigBridgeR::SigBridgeRBase` <- function(x, value) {
  Abort("Modification is restricted", "Please reconstruct a new object")
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
#' @family SigBridgeR-S3
#' @export
#' @rawNamespace S3method(format,"SigBridgeR::SigBridgeRBase")
`format.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  cls_x <- sub("^.*::", "", class(x)[1])

  properties <- props(x)
  nms_properties <- names(properties)

  value_to_chr <- function(value) {
    if (is.function(value)) {
      toString(fn_fmls_names(value))
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
#' @family SigBridgeR-S3
#' @export
#' @rawNamespace S3method(as.list,"SigBridgeR::SigBridgeRBase")
`as.list.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  props(x)
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
#' @family SigBridgeR-S3
#' @export
#' @rawNamespace S3method(as.character,"SigBridgeR::SigBridgeRBase")
`as.character.SigBridgeR::SigBridgeRBase` <- function(x, ...) {
  Abort(
    "cannot coerce type {.cls SigBridgeRBase} to vector of type {.cls character}",
    "Maybe you want to use {.fun format}?"
  )
}

#' Convert a List to a Compact Character Representation
#'
#' Converts an R `list` into a single character string using the format
#' `"key: value"`. Named elements use their names as keys, while unnamed
#' elements use their 1-based positions. Nested lists are converted
#' recursively and enclosed in parentheses.
#'
#' This method is registered as `as.character.list`, so it is called
#' automatically when `as.character()` is applied to a list.
#'
#' @param x A list to be converted.
#' @param ... Additional arguments passed to or from other methods.
#'   Currently ignored.
#'
#' @return A character vector of length one.
#'
#' @examples
#' x <- list(key1 = "val1", key2 = "val2")
#' as.character(x)
#' # "key1: val1, key2: val2"
#'
#' x <- list(key1 = list(key2 = list(key3 = "val3")))
#' as.character(x)
#' # "key1: (key2: (key3: val3))"
#'
#' x <- list(
#'   a = 1,
#'   b = list(
#'     c = 2,
#'     d = list(e = "hello")
#'   ),
#'   f = c(1, 2, 3)
#' )
#' as.character(x)
#' # "a: 1, b: (c: 2, d: (e: hello)), f: 1, 2, 3"
#'
#' x <- list("foo", bar = NA, baz = NULL)
#' as.character(x)
#' # "1: foo, bar: NA, baz: NULL"
#'
#' @export
as.character.list <- function(x, ...) {
  list_to_character_recursive(x)
}
