Algorithm <- new_class(
  name = "Algorithm",
  parent = SigBridgeRBase,
  properties = list(
    method_name = class_character,
    executor = class_function,
    algorithm_version = class_character
  ),
  abstract = TRUE,
  constructor = \(method_name, executor, algorithm_version) {
    new_object(
      S7_object(),
      method_name = method_name,
      executor = func,
      algorithm_version = algorithm_version
    )
  }
)

ScreenMethod <- new_class(
  name = "ScreenMethod",
  parent = Algorithm,
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
      return(
        cli::cli_fmt(cli::cli_alert_danger(
          "Missing arguments in function: {.arg {missing_args}}"
        ))
      )
    }
  }
)

AnnotationMethod <- new_class(
  name = "AnnotationMethod",
  parent = Algorithm
)
