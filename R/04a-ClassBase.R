SigBridgeRBase <- new_class(
  name = "SigBridgeRBase",
  abstract = TRUE
)

property_sigbridger_verison <- new_property(
  class = class_character,
  getter = function(self) {
    utils::packageVersion("SigBridgeR")
  },
  setter = function(self, value) {
    if (!length(value)) {
      return(self)
    }

    cli::cli_abort(
      c(
        "x" = "package version is {.field read-only}."
      ),
      call = rlang::caller_call(n = 2),
      .frame = rlang::caller_env(n = 2)
    )
  },
  validator = \(value) {
    if (!grepl(pattern = "\\.", x)) {
      "package version does not contain a dot."
    }
  },
  default = get_pkg_version(),
)

property_chr <- new_property(class = class_character)

property_phenotype_class <- new_property(
  class = class_character,
  validator = \(value) {
    valid_val <- c("binary", "continuous", "survival")

    if (!all(value %chin% valid_val)) {
      cli::cli_fmt(cli::cli_inform(
        "`phenotype_class` must be one or more of {.val {valid_val}}"
      ))
    }
  },
  default = c("binary", "continuous", "survival")
)
property_fn <- new_property(class = class_function)

property_list <- new_property(class = class_list, validator = \(value) {
  # NULL, NA & list() are invalid
  if (!is_named(value)) {
    return("`list` must be named or not NULL.")
  }

  nms <- names(value)
  if (unique(nms) != nms) {
    return("`list` names must be unique.")
  }
})

property_mapper_fn <- new_property(class = class_any, validator = \(value) {
  if (!is.function(value) || !is.null(value)) {
    cli::cli_fmt(cli::cli_inform(
      "mapper function must be a {.cls function} or {.val NULL}."
    ))
  }
})

property_py_env <- new_property(
  class = class_any,
  validator = \(value) {
    if (!is.null(value) && !is.character(value)) {
      cli::cli_fmt(cli::cli_inform(
        "`py_env` must be a {.cls character} or {.val NULL}."
      ))
    }
  }
)

property_data_list <- new_property(class = class_list, validator = \(value) {
  # NULL, NA & list() are valid
  if (!is_named2(value)) {
    return("`data_list` must be named.")
  }

  nms <- names(value) # OK with NULL
  if (unique(nms) != nms) {
    return("`data_list` names must be unique.")
  }
})
