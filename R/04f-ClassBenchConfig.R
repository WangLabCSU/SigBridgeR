#' @returns A `BenchConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchConfig <- new_class(
  name = "BenchConfig",
  parent = SigBridgeRBase,
  properties = list(
    param = property_data_list # NULL is OK; stores the parameters
  ),
  constructor = \(
    sigbridger_version = get_pkg_version()
  ) {
    new_object(
      S7_object(),
      sigbridger_version = sigbridger_version
    )
  }
)


#' @title BenchAccConfig: Accuracy Benchmark Configuration
#'
#' @description
#' `BenchAccConfig` stores the configuration of an accuracy benchmark run
#' performed by [BenchAcc()]. It extends [BenchConfig] with the benchmarked
#' object, and forwards every additional argument to [BenchAcc()] as `param`.
#'
#' @param obj The object to be benchmarked, either a
#'   [Seurat][SeuratObject::Seurat-class] object or a [ScreenMethodResult]
#'   object.
#' @param ... Additional named arguments forwarded to [BenchAcc()] (e.g., the
#'   ground-truth labels or the metric to compute). Stored in `param`.
#' @param sigbridger_version Package Version of SigBridgeR
#'
#' @returns A `BenchAccConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchAccConfig <- new_class(
  name = "BenchAccConfig",
  parent = BenchConfig,
  properties = list(
    obj = property_bench_obj
  ),
  constructor = \(obj, ..., sigbridger_version = get_pkg_version()) {
    new_object(
      S7_object(),
      obj = obj,
      param = rlang::list2(...),
      sigbridger_version = sigbridger_version
    )
  }
)

# --------------------------------------------------------------------------------

#' @title BenchGSEAConfig: GSEA Benchmark Configuration
#'
#' @description
#' `BenchGSEAConfig` stores the configuration of a GSEA benchmark run
#' performed by [BenchGSEA()]. It extends [BenchConfig] with the benchmarked
#' object, and forwards every additional argument to [BenchGSEA()] as `param`.
#'
#' @inheritParams BenchAccConfig
#'
#' @returns A `BenchGSEAConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchGSEAConfig <- new_class(
  name = "BenchGSEAConfig",
  parent = BenchConfig,
  properties = list(
    obj = property_bench_obj
  ),
  constructor = \(obj, ..., sigbridger_version = get_pkg_version()) {
    new_object(
      S7_object(),
      obj = obj,
      param = rlang::list2(...),
      sigbridger_version = sigbridger_version
    )
  }
)


# --------------------------------------------------------------------------------

#' @title BenchSSGSEAConfig: ssGSEA Benchmark Configuration
#'
#' @description
#' `BenchSSGSEAConfig` stores the configuration of an ssGSEA benchmark run
#' performed by [BenchSSGSEA()]. It extends [BenchConfig] with the benchmarked
#' object, and forwards every additional argument to [BenchSSGSEA()] as
#' `param`.
#'
#' @inheritParams BenchAccConfig
#'
#' @returns A `BenchSSGSEAConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchSSGSEAConfig <- new_class(
  name = "BenchSSGSEAConfig",
  parent = BenchConfig,
  properties = list(
    obj = property_bench_obj
  ),
  constructor = \(obj, ..., sigbridger_version = get_pkg_version()) {
    new_object(
      S7_object(),
      obj = obj,
      param = rlang::list2(...),
      sigbridger_version = sigbridger_version
    )
  }
)
