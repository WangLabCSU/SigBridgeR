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
      cli::cli_fmt(cli::cli_alert_danger(
        "package version does not contain a dot."
      ))
    }
  },
  default = quote(utils::packageVersion("SigBridgeR")),
)

property_chr <- new_property(class = class_character)

property_phenotype_class <- new_property(
  class = class_character,
  validator = \(value) {
    valid_val <- c("binary", "continuous", "survival")

    if (!all(value %chin% valid_val)) {
      cli::cli_fmt(cli::cli_alert_danger(
        "`phenotype_class` must be one or more of {.val {valid_val}}"
      ))
    }
  },
  default = c("binary", "continuous", "survival")
)
property_fn <- new_property(class = class_function)

property_list <- new_property(class = class_list, validator = \(value) {
  if (!is_named(value)) {
    cli::cli_fmt(cli::cli_alert_danger(
      "`list` must be named."
    ))
  }
})
