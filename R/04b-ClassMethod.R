# --------------------------------------------------------------------------------

#' @title Method: Abstract Base Class for Analysis Methods
#'
#' @description
#' `Method` is the abstract base class for all analysis methods (screening and
#' annotation) in the SigBridgeR package. It defines the common interface that
#' every method must implement, including a method name, version, and executor
#' function.
#'
#' @details
#' This class is defined using the S7 object-oriented system and inherits from
#' [SigBridgeRBase]. It is marked as abstract, meaning it cannot be instantiated
#' directly. Concrete subclasses such as [ScreenMethod] and `AnnotationMethod`
#' provide specific method categories.
#'
#' @param method_name `character`. A string identifying the method (e.g.,
#'   `"Scissor"`, `"scPAS"`).
#' @param executor `function`. The function that performs the actual
#'   computation. This is the core workhorse of each method.
#' @param method_version `character`. A string indicating the version of the
#'   underlying software package.
#'
#' @returns A `Method` S7 object (abstract; cannot be instantiated directly).
#'
#' @family SigBridgeR-Classes
#' @export
Method <- new_class(
  name = "Method",
  parent = SigBridgeRBase,
  properties = list(
    method_name = property_chr,
    method_version = property_chr,
    executor = property_fn
  ),
  abstract = TRUE,
  constructor = \(method_name, executor, method_version) {
    new_object(
      S7_object(),
      method_name = method_name,
      method_version = method_version,
      executor = func
    )
  }
)

# --------------------------------------------------------------------------------

#' @title ScreenMethod: Screening Method Class for Phenotype-Associated Cell Detection
#'
#' @description
#' `ScreenMethod` is the S7 class for phenotype-associated cell screening
#' methods. It extends [Method] with phenotype classification and mapper
#' function support, and enforces a standardized function signature for all
#' screening executors.
#'
#' @details
#' This class inherits from [Method] and adds screening-specific properties.
#' It validates that the executor function accepts the required arguments:
#' `matched_bulk`, `sc_data`, `phenotype`, `label_type`, and
#' `phenotype_class`. This ensures a consistent interface across all screening
#' algorithms (e.g., Scissor, scPAS, scPP, TiRank).
#'
#' @param method_name `character`. (Inherited from [Method]) A string
#'   identifying the screening method.
#' @param executor `function`. (Inherited from [Method]) The screening
#'   function. Must accept `matched_bulk`, `sc_data`, `phenotype`,
#'   `label_type`, and `phenotype_class` as arguments.
#' @param method_version `character`. (Inherited from [Method]) A string
#'   indicating the version of the underlying software package.
#' @param phenotype_class `character`. A character vector specifying the
#'   phenotype classification type(s) the method supports. Valid values:
#'   `"binary"`, `"continuous"`, `"survival"`. Default: all three.
#' @param mapper `function` or `NULL`. An optional function used to transform
#'   or map parameters before screening. When `NULL`, no mapper is applied.
#'
#' @returns A `ScreenMethod` S7 object.
#'
#' @family SigBridgeR-Classes
#' @export
ScreenMethod <- new_class(
  name = "ScreenMethod",
  parent = Method,
  properties = list(
    phenotype_class = property_phenotype_class,
    mapper = property_mapper_fn
  ),
  validator = \(self) {
    param_names <- fn_fmls_names(self@executor)
    expected_args <- c(
      "matched_bulk",
      "sc_data",
      "phenotype",
      "label_type",
      "phenotype_class"
    )
    if (!all(expected_args %chin% param_names)) {
      missing_args <- setdiff(expected_args, param_names)

      cli::cli_fmt(cli::cli_text(
        "Missing arguments in screening function: {.arg {missing_args}}"
      ))
    }
  }
)
