# 同时用作screen result和cache config
ScreenMethodConfig <- new_class(
  name = "ScreenMethodConfig",
  parent = SigBridgeRBase,
  properties = list(
    method_name = property_chr,
    method_version = property_chr,
    screen_method = property_chr,
    phenotype_class = property_chr,
    label_type = property_chr,
    params = property_list
  ),
  validator = \(self) {
    if (length(self@params) == 0L) {
      return("@params is empty")
    }

    valid_phenotype_class <- c("binary", "continuous", "survival")
    if (!self@phenotype_class %chin% valid_phenotype_class) {
      return(cli::cli_fmt(cli::cli_inform(
        "@phenotype_class must be one of {.val {valid_phenotype_class}}, but got {.val {self@phenotype_class}}"
      )))
    }

    valid_method_name <- names(ScreenStrategy)
    if (!self@method_name %chin% valid_method_name) {
      return(cli::cli_fmt(cli::cli_inform(
        "@method_name must be one of {.val {valid_method_name}}, but got {.val {self@method_name}}"
      )))
    }
  }
)
