# --------------------------------------------------------------------------------

#' @title BenchAccConfig: Accuracy Benchmark Configuration
#'
#' @description
#' `BenchAccConfig` stores the configuration of an accuracy benchmark run
#' performed by [BenchAcc()]. It extends [BenchConfig] with the ground-truth
#' cell labels collected in `bench_config`; the inputs shared by every
#' screening method (`matched_bulk`, `sc_data`, `phenotype`, `label_type`,
#' `phenotype_class`, `screen_method`) and the method specific arguments
#' collected in `method_param` are inherited from [BenchConfig].
#'
#' @details
#' `bench_config` supplies the reference ("ground-truth") labels that the
#' screening results are scored against. It is a `character` vector interpreted
#' in one of two ways:
#'
#' * Length one: the name of a column in `sc_data@meta.data` holding the
#'   per-cell labels. The column must exist and contain a `"Positive"` level.
#' * Length greater than one: a per-cell label vector, one label per cell. It
#'   must contain a `"Positive"` level and have the same length as the number
#'   of cells in `sc_data`.
#'
#' In both cases a `"Positive"` label is required, since [BenchAcc()] computes
#' binary classification metrics (TPR, FPR, TNR, FNR, F1, Accuracy) treating
#' `"Positive"` as the positive class. An empty `bench_config` is rejected by
#' the validator.
#'
#' @inheritParams BenchConfig
#' @param bench_config `character`. Ground-truth labels used to score the
#'   screening results: either a length-one column name in `sc_data@meta.data`,
#'   or a per-cell label vector of length `ncol(sc_data)`. Must contain a
#'   `"Positive"` label. Default: `character()`, which is rejected by the
#'   validator.
#' @param ... Must be empty. Method specific arguments belong in `method_param`.
#'
#' @returns A `BenchAccConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchAccConfig <- S7::new_class(
  name = "BenchAccConfig",
  parent = BenchConfig,
  properties = list(
    bench_config = property_chr
  ),
  validator = \(self) {
    if (length(self@bench_config) == 0L) {
      return(cli::format_inline(
        "{.val @bench_config} is not specified"
      ))
    }

    if (length(self@bench_config) == 1L) {
      exists_bench_label <- self@bench_config %chin% colnames(self@sc_data[[]])
      if (!exists_bench_label) {
        return(cli::format_inline(
          "{.val {self@bench_config}} is not a valid column name in sc_data@meta.data"
        ))
      }
      bench_label <- self@sc_data[[]][[self@bench_config]] |> table() |> names()
      if (!"Positive" %chin% bench_label) {
        return(cli::format_inline(
          "{.val {self@bench_config}} in sc_data@meta.data doesn't have a `Positive` label"
        ))
      }
    }

    bench_label <- self@bench_config |> table() |> names()
    if (!"Positive" %chin% bench_label) {
      return(cli::format_inline(
        "{.val @bench_config} doesn't have a `Positive` label"
      ))
    }
    if (length(self@bench_config) != ncol(self@sc_data)) {
      n_cell <- ncol(self@sc_data)
      n_bench_config <- length(self@bench_config)
      return(cli::format_inline(
        "{.val bench_config} doesn't have the same length as cells. \
        Number of cells: {.val {n_cell}}, length of bench_config: {.val {n_bench_config)}}"
      ))
    }
  },
  constructor = \(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "{method_name}",
    phenotype_class = c("binary", "survival", "continuous"),
    screen_method = names(ScreenStrategy),
    method_param = MethodParamSet(),
    bench_config = character(),
    ...
  ) {
    rlang::check_dots_empty()
    S7::new_object(
      BenchConfig(
        matched_bulk = matched_bulk,
        sc_data = sc_data,
        phenotype = phenotype,
        label_type = label_type,
        phenotype_class = phenotype_class,
        screen_method = screen_method,
        method_param = method_param,
        sigbridger_version = get_pkg_version()
      ),
      bench_config = bench_config
    )
  }
)
