# nocov start

#' @title Registry of Phenotype-Associated Cell Screening Methods
#' @family Add_Screen_method
#' @description
#' An environment storing methods for screening phenotype-associated cells.
#'
#' @details
#' Storing structure - named list
#'  - `key`: method name
#'  - `value`: list
#'     - `method_name`: method name
#'     - `executor`: function implementation of the method
#'     - `phenotype_class`: The phenotype types supported by this method
#'     - `mapper`: A function that transforms the parameters passed to the `Screen` function before forwarding them to the executor; both input and output must be of type `list`.
#'
#' @export
ScreenStrategy <- new_environment(
  list(
    Scissor = ScreenMethod(
      method_name = "Scissor",
      method_version = r_pkg_version("Scissor"),
      executor = DoScissor,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    scPAS = ScreenMethod(
      method_name = "scPAS",
      method_version = r_pkg_version("scPAS"),
      executor = DoscPAS,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    scAB = ScreenMethod(
      method_name = "scAB",
      method_version = r_pkg_version("scAB"),
      executor = DoscAB,
      phenotype_class = c("binary", "survival"),
      mapper = NULL
    ),
    scPP = ScreenMethod(
      method_name = "scPP",
      method_version = r_pkg_version("ScPP"),
      executor = DoscPP,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    DEGAS = ScreenMethod(
      method_name = "DEGAS",
      method_version = r_pkg_version("DEGAS"),
      executor = DoDEGAS,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    DEGASv2 = ScreenMethod(
      method_name = "DEGASv2",
      method_version = r_pkg_version("DEGASv2"),
      executor = DoDEGASv2,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    LP_SGL = ScreenMethod(
      method_name = "LP_SGL",
      method_version = r_pkg_version("LPSGL"),
      executor = DoLP_SGL,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    PIPET = ScreenMethod(
      method_name = "PIPET",
      method_version = r_pkg_version("PIPET"),
      executor = DoPIPET,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    SIDISH = ScreenMethod(
      method_name = "SIDISH",
      method_version = r_pkg_version("rSIDISH"),
      executor = DoSIDISH,
      phenotype_class = "survival",
      mapper = NULL
    ),
    SCIPAC = ScreenMethod(
      method_name = "SCIPAC",
      method_version = r_pkg_version("SCIPAC"),
      executor = DoSCIPAC,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    TiRank = ScreenMethod(
      method_name = "TiRank",
      method_version = r_pkg_version("rTiRank"),
      executor = DoTiRank,
      phenotype_class = c("binary", "survival", "continuous"),
      mapper = NULL
    )
  )
)
