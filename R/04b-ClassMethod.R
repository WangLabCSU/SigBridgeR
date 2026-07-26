#' @export
Method <- new_class(
  name = "Method",
  parent = SigBridgeRBase,
  properties = list(
    method_name = property_chr,
    method_version = property_chr,
    executor = property_fn,
  ),
  abstract = TRUE,
  constructor = \(method_name, executor, algorithm_version) {
    new_object(
      S7_object(),
      method_name = method_name,
      method_version = algorithm_version,
      executor = func
    )
  }
)

#' @export
ScreenMethod <- new_class(
  name = "ScreenMethod",
  parent = Method,
  properties = list(
    phenotype_class = property_phenotype_class,
    mapper = property_fn
  ),
  validator = \(self) {
    param_names <- fn_fmls_names(self$executor)
    expected_args <- c(
      "matched_bulk",
      "sc_data",
      "phenotype",
      "label_type",
      "phenotype_class"
    )
    if (!all(expected_args %in% param_names)) {
      missing_args <- setdiff(expected_args, param_names)

      cli::cli_fmt(cli::cli_alert_danger(
        "Missing arguments in function: {.arg {missing_args}}"
      ))
    }
  }
)
