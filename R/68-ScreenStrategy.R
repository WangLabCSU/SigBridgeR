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
#'     - `phenotypes`: The phenotype types supported by this method
#'     - `mapper`: A function that transforms the parameters passed to the `Screen` function before forwarding them to the executor; both input and output must be of type `list`.
#'
#' @export
ScreenStrategy <- rlang::new_environment(
  list(
    Scissor = list(
      method_name = "Scissor",
      executor = DoScissor,
      phenotypes = c("binary", "survival", "continuous"),
      mapper = function(params) {
        params$family <- switch(
          params$phenotype_class,
          "binary" = "binomial",
          "survival" = "cox",
          "continuous" = "gaussian"
        )
        params
      }
    ),
    scPAS = list(
      method_name = "scPAS",
      executor = DoscPAS,
      phenotypes = c("binary", "survival", "continuous"),
      mapper = function(params) {
        params$family <- switch(
          params$phenotype_class,
          "binary" = "binomial",
          "survival" = "cox",
          "continuous" = "gaussian"
        )
        params
      }
    ),
    scAB = list(
      method_name = "scAB",
      executor = DoscAB,
      phenotypes = c("binary", "survival"),
      mapper = NULL
    ),
    scPP = list(
      method_name = "scPP",
      executor = DoscPP,
      phenotypes = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    DEGAS = list(
      method_name = "DEGAS",
      executor = DoDEGAS,
      phenotypes = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    LP_SGL = list(
      method_name = "LP_SGL",
      executor = DoLP_SGL,
      phenotypes = c("binary", "survival", "continuous"),
      mapper = function(params) {
        params$family <- switch(
          params$phenotype_class,
          "binary" = "logit",
          "survival" = "cox",
          "continuous" = "linear"
        )
        params
      }
    ),
    PIPET = list(
      method_name = "PIPET",
      executor = DoPIPET,
      phenotypes = c("binary", "survival", "continuous"),
      mapper = NULL
    ),
    SIDISH = list(
        method_name = "SIDISH",
        executor = DoSIDISH,
        phenotypes = "survival",
        mapper = NULL
    )
  )
)
