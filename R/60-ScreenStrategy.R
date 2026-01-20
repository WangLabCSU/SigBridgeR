#' @title Registry of Phenotype-Associated Cell Screening Methods
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
ScreenStrategy <- rlang::new_environment()

#' @keywords internal
GetExistingStrategy <- function(
  # ? an environment
  registry = ScreenStrategy,
  ...
) {
  chk::chk_environment(registry)
  names(ScreenStrategy)
}
