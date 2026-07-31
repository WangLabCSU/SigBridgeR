#' @title SigBridgeR Base Class Properties
#'
#' @description
#' These are S7 property definitions used within the `SigBridgeRBase` class
#' and its subclasses. They define the types, validators, and default values
#' for class attributes. These properties are exported for community developers
#' to extend and subclass `SigBridgeRBase`.
#'
#' @name SigBridgeR-properties
#' @rdname SigBridgeR-properties
#' @docType data
NULL

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_character`.
#' @details
#' A read-only character property that returns the current version of the
#' `SigBridgeR` package. Attempting to set this property will raise an error.
#'
#' The version string is validated to contain at least one dot (`.`).
#' The default value is obtained from `get_pkg_version()`.
property_sigbridger_verison <- new_property(
  class = class_character,
  validator = \(value) {
    if (!grepl(pattern = "\\.", value)) {
      "Package version does not contain a dot."
    }
  },
  default = get_pkg_version()
)

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_character`.
#' @details
#' A general-purpose character property with no additional validation.
#' Can be used for any string-typed attribute in subclasses.
property_chr <- new_property(class = class_character)

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_character`.
#' @details
#' A character property that specifies the phenotype classification type(s).
#' Valid values are `"binary"`, `"continuous"`, and `"survival"`.
#' Multiple values can be provided as a character vector.
#'
#' Default: `c("binary", "continuous", "survival")`.
#'
#' @examples
#' # Valid usage
#' property_phenotype_class@validator("binary")
#' property_phenotype_class@validator(c("binary", "survival"))
property_phenotype_class <- new_property(
  class = class_character,
  validator = \(value) {
    valid_val <- c("binary", "continuous", "survival")

    if (!all(value %chin% valid_val)) {
      cli::cli_fmt(cli::cli_text(
        "`phenotype_class` must be one or more of {.val {valid_val}}"
      ))
    }
  },
  default = c("binary", "continuous", "survival")
)

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_function`.
#' @details
#' A property that accepts a function value. Used to store callback or
#' transformation functions in class attributes.
property_fn <- new_property(class = class_function)


#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_list`.
#' @details
#' A named list property. The value must be a named list (not `NULL` or
#' `NA`) with unique names. Duplicate names are not allowed.
property_list <- new_property(class = class_list, validator = \(value) {
  # NULL, NA & list() are invalid
  if (!is_named(value)) {
    "`list` must be named or not NULL."
  }

  nms <- names(value)
  if (unique(nms) != nms) {
    "`list` names must be unique."
  }
})

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_any`.
#' @details
#' A property that accepts either a function or `NULL`. Used for storing
#' mapper/transformation functions that may be optionally specified.
#'
#' If a value is provided, it must be a function; otherwise, `NULL` is
#' accepted to indicate no mapper is configured.
property_mapper_fn <- new_property(class = class_any, validator = \(value) {
  if (!is.function(value) && !is.null(value)) {
    cli::cli_fmt(cli::cli_text(
      "mapper function must be a {.cls function} or {.val NULL}."
    ))
  }
})

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_any`.
#' @details
#' A property representing a Python environment configuration.
#' Accepts either a character string (e.g., a conda environment name or path)
#' or `NULL` (when no Python environment is specified).
property_py_env <- new_property(
  class = class_any,
  validator = \(value) {
    if (!is.null(value) && !is.character(value)) {
      cli::cli_fmt(cli::cli_text(
        "`value` must be a {.cls character} or {.val NULL}."
      ))
    }
  }
)

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_list`.
#' @details
#' A named list property for storing multiple data objects.
#' The list must be named (if not `NULL`) and names must be unique.
#'
#' Unlike `property_list`, this property allows `NULL`, `NA`, and empty
#' lists (`list()`) as valid values, making it suitable for optional
#' data storage.
property_data_list <- new_property(class = class_list, validator = \(value) {
  # NULL, NA & list() are valid
  if (!is_named2(value)) {
    "`data_list` must be named."
  }

  nms <- names(value) # OK with NULL
  if (unique(nms) != nms) {
    "`data_list` names must be unique."
  }
})

#' @rdname SigBridgeR-properties
#' @format An S7 property of class `class_any`.
#' @details
#' A property that requires a [Seurat][SeuratObject::Seurat-class] object.
#' The value must inherit from the `"Seurat"` class.
#'
#' This property is used to store single-cell data in Seurat format
#' for downstream analysis.
#'
#' @seealso [SeuratObject::Seurat-class]
property_seurat <- new_property(
  class = class_any,
  validator = \(value) {
    if (!inherits(value, "Seurat")) {
      cli::cli_fmt(cli::cli_text("value must be a {.cls Seurat} object"))
    }
  }
)

# --------------------------------------------------------------------------------

#' @title SigBridgeRBase: Abstract Base Class for SigBridgeR
#'
#' @description
#' `SigBridgeRBase` is the abstract base class for all S7 classes in the
#' SigBridgeR package. It provides a unified foundation for S3 method
#' dispatch and ensures consistent behavior across all subclasses.
#'
#' @details
#' This class is defined using the S7 object-oriented system and is marked
#' as abstract, meaning it cannot be instantiated directly. All other
#' classes in the SigBridgeR package should inherit from `SigBridgeRBase`.
#'
#' The primary purpose of this base class is to serve as a dispatch target
#' for base R S3 generic functions, enabling polymorphic behavior across
#' all SigBridgeR S7 objects through S3 method dispatch.
#' @family SigBridgeR-Classes
#' @export
SigBridgeRBase <- new_class(
  name = "SigBridgeRBase",
  properties = list(sigbridger_version = property_sigbridger_verison),
  abstract = TRUE
)
