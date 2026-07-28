#' @title Register a Custom Screening Method for Phenotype-Driven Analysis
#' @family Add_Screen_method
#' @family Registering
#' @description
#' `r lifecycle::badge('deprecated')` Registers one or more user-defined screening functions into a shared registry
#' (i.e., \code{ScreenStrategy}), enabling dynamic dispatch based on phenotype type
#' (binary, survival, or continuous). Each method is stored with metadata including
#' supported phenotypes and an optional parameter remapping function.
#'
#' This function supports community-driven extension of analysis pipelines—ideal for
#' integrating novel single cell phenotypic screening methods while
#' maintaining compatibility with downstream validation and execution frameworks.
#'
#' @param ... Named arguments where each value is a function implementing a screening method.
#'   The name becomes the method identifier (e.g., \code{Scissor = DoScissor} registers
#'   \code{DoScissor} under the key \code{"Scissor"}). Unnamed arguments are auto-named
#'   using their expression (e.g., \code{Scissor} → name \code{"Scissor"}).
#' @param supported_phenotypes Character vector specifying which phenotype types the method supports.
#'   Must be a subset of \code{c("binary", "survival", "continuous")}. Default: all three.
#' @param parameter_mapper A function that transforms the input parameter list before
#'   passing it to the executor. Useful for changing parameters from interface function.
#'   Receives a named list \code{params} and must return a modified list.
#'   Default: \code{NULL}.
#' @param registry An environment used as a method registry (e.g., \code{ScreenStrategy}).
#'   New methods are stored as \code{registry[["MethodName"]]} = metadata list.
#'   Default: \code{ScreenStrategy}.
#' @param overwrite Logical. If \code{FALSE} (default), throws an error when attempting
#'   to replace an existing method. Set to \code{TRUE} to allow updates.
#' @param verbose Logical. Whether to print a success message upon registration.
#'   Default: inherits from package option.
#'
#' @section Registry Entry Structure:
#' Each registered method is stored as a list with four elements:
#' \describe{
#'   \item{\code{method_name}}{Character: the method key (e.g., "Scissor")}
#'   \item{\code{executor}}{Function: the actual screening implementation}
#'   \item{\code{phenotypes}}{Character vector: supported phenotype classes}
#'   \item{\code{mapper}}{Function: parameter transformation hook}
#' }
#'
#' @return Invisibly returns \code{TRUE} on successful registration.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Register a custom screen method for survival data
#' DoMyScreen <- function(...) {
#'   stop("Not implemented")
#' }
#'
#' RegisterScreenMethod(
#'   # Key = function
#'   MyCox = DoMyScreen,
#'   # If `continuous` phenotype is not supported
#'   supported_phenotypes = c("survival", "binary"),
#'   parameter_mapper = function(params) {
#'     # If I prefer to use `quiet` instead of `verbose`
#'     params$quiet <- !params$verbose
#'     params
#'   }
#' )
#' }
RegisterScreenMethod <- function(
  # ? Register a new screen method using named functions:
  # ? Usage:
  # (Key = function())
  # Scissor = DoScissor(...)
  ...,
  supported_phenotypes = c("binary", "survival", "continuous"),
  # ? Use to change the name of parameters when passing to specific method
  # ? Usage:
  # parameter_mapper = function(params) {
  #   # original -> renamed
  #   params$scissor_family <- switch(
  #     params$phenotype_class,
  #     "binary" = "binomial",
  #     "survival" = "cox",
  #     "continuous" = "gaussian"
  #   )
  #   params
  # }
  parameter_mapper = NULL,
  # ? an environment to store the screen methods
  registry = ScreenStrategy,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  lifecycle::deprecate_warn("4.0.0", "RegisterScreenMethod()", "Register()")
  # * check the input
  if (!is.null(parameter_mapper)) {
    chk::chk_function(parameter_mapper)
  }
  if (!all(supported_phenotypes %chin% c("binary", "survival", "continuous"))) {
    cli::cli_abort(c(
      "x" = "unsupported phenotype class when registering screen method",
      ">" = "Current supported phenotypes are: binary, survival, continuous",
      ">" = "Provided phenotypes are: {supported_phenotypes}"
    ))
  }
  chk::chk_logical(verbose)
  chk::chk_logical(overwrite)
  chk::chk_environment(registry)

  # * detect where are functions
  dots <- rlang::list2(...)
  is_fun <- vapply(
    X = dots,
    FUN = \(x) rlang::is_function(x),
    FUN.VALUE = logical(1)
  )
  method_names <- get_names_4_ids(..., .quoses = rlang::enquos(...))[is_fun]

  for (i in seq_len(sum(is_fun))) {
    method_name <- method_names[i]
    executor <- dots[is_fun][[i]]

    # Check for existing method
    exists_already <- method_name %chin% names(registry)
    if (exists_already && !overwrite) {
      cli::cli_abort(c(
        "x" = "Method already exists: {.val {method_name}}",
        "i" = "Registered methods: {.val {names(registry)}}",
        "i" = "Use {.code overwrite = TRUE} to force replacement"
      ))
    }

    registry[[method_name]] <- rlang::list2(
      method_name = method_name,
      executor = executor,
      phenotypes = supported_phenotypes,
      mapper = parameter_mapper
    )

    if (verbose) {
      if (exists_already) {
        cli::cli_alert_warning("Updated {.field {method_name}}")
      } else {
        cli::cli_alert_success("Registered {.field {method_name}}")
      }
    }
  }

  invisible(TRUE)
}
