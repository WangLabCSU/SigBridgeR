# --------------------------------------------------------------------------------

#' @title ScreenMethodConfig: Screening Method Configuration
#'
#' @description
#' `ScreenMethodConfig` stores the configuration for a phenotype-associated
#' cell screening method. It is used both to define screening parameters and
#' as a cache result descriptor.
#'
#' @details
#' This class is defined using the S7 object-oriented system and inherits from
#' [SigBridgeRBase]. It validates that `method_name` is a registered screening
#' method in [ScreenStrategy] and that `phenotype_class` is one of the
#' supported types.
#'
#' ## Validator
#'
#' The class validator enforces:
#'
#' * `param` must be a non-empty named list.
#' * `phenotype_class` must be one or more of `"binary"`,
#'   `"continuous"`, `"survival"`.
#' * `method_name` must match a key in [ScreenStrategy].
#'
#' @param sigbridger_version Package Version of SigBridgeR
#' @param method_name `character`. The name of the screening method (e.g.,
#'   `"Scissor"`, `"scPAS"`). Must match a key in [ScreenStrategy].
#' @param param `list`. A named list of parameters passed to the screening
#'   function. Must be non-empty.
#' @param method_version `character` or `NULL`. The version of the
#'   underlying software package. When `NULL`, automatically resolved via
#'   `r_pkg_version(method_name)`.
#' @param phenotype_class `character`. The phenotype classification type(s).
#'   One or more of `"binary"`, `"continuous"`, `"survival"`. Default: all
#'   three.
#' @param label_type `character`. The label type for the screening method.
#'   Default: `character(0L)`.
#'
#' @returns A `ScreenMethodConfig` S7 object.
#'
#' @family S7-Classes
#' @export
ScreenMethodConfig <- new_class(
  name = "ScreenMethodConfig",
  parent = SigBridgeRBase,
  properties = list(
    method_name = property_chr,
    method_version = property_verison,
    phenotype_class = property_chr,
    label_type = property_chr,
    param = property_data_list # NULL is OK; stores the parameters
  ),
  validator = \(self) {
    if (length(self@param) == 0L) {
      return("@param is empty")
    }

    valid_phenotype_class <- c("binary", "continuous", "survival")
    if (!any(self@phenotype_class %chin% valid_phenotype_class)) {
      return(cli::cli_fmt(cli::cli_text(
        "@phenotype_class must be one of {.val {valid_phenotype_class}},\
         but got {.val {self@phenotype_class}}"
      )))
    }

    valid_method_name <- names(ScreenStrategy)
    if (!any(self@method_name %chin% valid_method_name)) {
      return(cli::cli_fmt(cli::cli_text(
        "@method_name must be one of {.val {valid_method_name}},\
         but got {.val {self@method_name}}"
      )))
    }
    invisible()
  },
  constructor = \(
    method_name,
    param,
    method_version = NULL,
    phenotype_class = c("binary", "continuous", "survival"),
    label_type = character(0L),
    sigbridger_version = get_pkg_version()
  ) {
    new_object(
      S7_object(),
      method_name = method_name,
      method_version = method_version %||% r_pkg_version(method_name),
      phenotype_class = phenotype_class,
      label_type = label_type,
      param = param,
      sigbridger_version = sigbridger_version
    )
  }
)

# --------------------------------------------------------------------------------

#' @title ScreenMethodCache: Screening Method Cache Container
#'
#' @description
#' `ScreenMethodCache` bundles cached screening data with its configuration
#' and metadata. It serves as the primary container for reading from and
#' writing to the cache system.
#'
#' @details
#' This class is defined using the S7 object-oriented system and inherits from
#' [SigBridgeRBase]. It wraps a `ScreenMethodConfig` together with cached data
#' objects, file paths, and a read-only SigBridgeR version stamp.
#'
#' ## Validator
#'
#' The class validator checks that:
#'
#' * `cache_path` is an existing directory.
#' * The parent directory of `cache_config_path` exists.
#' * `cache_config_path` ends with `"cache_config.json"`.
#'
#' @param sigbridger_version Package Version of SigBridgeR
#' @param cache_path `character`. Path to the cache directory where data
#'   files are stored.
#' @param cache_config_path `character`. Path to the `cache_config.json`
#'   metadata file. Must end with `"cache_config.json"`.
#' @param cache_data `list`. A named list of cached data objects. Can be
#'   `NULL` or empty (for load-only use).
#' @param screen_method_config A `ScreenMethodConfig` S7 object containing
#'   the screening method configuration.
#'
#' @prop sigbridger_version `character` (read-only). The SigBridgeR package
#'   version at the time the cache was created. Automatically set from
#'   `get_pkg_version()` and cannot be modified by the user.
#'
#' @returns A `ScreenMethodCache` S7 object.
#'
#' @family S7-Classes
#' @export
ScreenMethodCache <- new_class(
  name = "ScreenMethodCache",
  parent = SigBridgeRBase,
  properties = list(
    cache_path = property_chr,
    cache_config_path = property_chr,
    cache_data = property_data_list, # NULL is OK; stores the data
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
      return(cli::cli_fmt(cli::cli_text(
        "@cache_config_path ({.file {self@cache_config_path}}) must end with {.file cache_config.json}"
      )))
    }
    invisible()
  },
  constructor = \(
    cache_path,
    cache_config_path,
    cache_data,
    screen_method_config,
    sigbridger_version = get_pkg_version()
  ) {
    new_object(
      S7_object(),
      cache_path = cache_path,
      cache_config_path = cache_config_path,
      cache_data = cache_data, # NULL is OK; stores the data
      screen_method_config = screen_method_config,
      sigbridger_version = sigbridger_version
    )
  }
)
