# 同时用作screen result和cache config
ScreenMethodConfig <- new_class(
  name = "ScreenMethodConfig",
  parent = SigBridgeRBase,
  properties = list(
    method_name = property_chr,
    method_version = property_chr,
    phenotype_class = property_chr,
    label_type = property_chr,
    param = property_data_list # NULL is OK; stores the parameters
  ),
  validator = \(self) {
    if (length(self@params) == 0L) {
      Abort("@params is empty")
    }

    valid_phenotype_class <- c("binary", "continuous", "survival")
    if (!self@phenotype_class %chin% valid_phenotype_class) {
      Abort(
        "@phenotype_class must be one of {.val {valid_phenotype_class}}, but got {.val {self@phenotype_class}}"
      )
    }

    valid_method_name <- names(ScreenStrategy)
    if (!self@method_name %chin% valid_method_name) {
      Abort(
        "@method_name must be one of {.val {valid_method_name}}, but got {.val {self@method_name}}"
      )
    }
  },
  constructor = \(
    method_name,
    param,
    method_version = NULL,
    phenotype_class = c("binary", "continuous", "survival"),
    label_type = character(0L)
  ) {
    new_object(
      S7_object(),
      method_name = method_name,
      method_version = method_version %||% r_pkg_version(method_name),
      phenotype_class = phenotype_class,
      label_type = label_type,
      param = param
    )
  }
)


ScreenMethodCache <- new_class(
  name = "ScreenMethodCache",
  parent = SigBridgeRBase,
  properties = list(
    cache_path = property_chr,
    cache_config_path = property_chr,
    cache_data = property_data_list, # NULL is OK; stores the data
    sigbridger_version = property_sigbridger_verison, # chr
    screen_method_config = new_property(class = ScreenMethodConfig)
  ),
  validator = \(self) {
    if (!dir.exists(self@cache_path)) {
      return(cli::cli_fmt(cli::cli_text(
        "@cache_path ({.file {self@cache_path}}) does not exist"
      )))
    }
    EnsureParentDir(self@cache_config_path)

    if (!endsWith(path, "cache_config.json")) {
      cli::cli_fmt(cli::cli_text(
        "@cache_config_path ({.file {self@cache_config_path}}) must end with {.file cache_config.json}"
      ))
    }
  },
  constructor = \(
    cache_path,
    cache_config_path,
    cache_data,
    screen_method_config
  ) {
    new_object(
      S7_object(),
      cache_path = cache_path,
      cache_config_path = cache_config_path,
      cache_data = cache_data, # NULL is OK; stores the data
      sigbridger_version = get_pkg_version(), # chr
      screen_method_config = screen_method_config
    )
  }
)
