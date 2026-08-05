#' @title Map Phenotype Values Using Conditional Rules
#'
#' @description
#' Efficiently transforms phenotype values based on user-defined conditional rules.
#' Similar to \code{dplyr::case_when} but optimized for performance on large datasets.
#'
#' @param data A vector, data frame, or data.table containing the phenotype data to transform.
#' @param ... Named or unnamed formulas of the form \code{condition ~ value}.
#'   Conditions are evaluated sequentially; the first matching condition determines the output value.
#'   Example: \code{col > 10 ~ "High", col <= 10 ~ "Low"}
#' @param .default Value to use when no conditions match. Default: \code{NA}.
#'
#' @return A data.table with the transformed column. If input is a vector, returns a vector.
#' @export
#' @family input_preprocess
#'
#' @examples
#' \dontrun{
#' # Example 1: Discretize a continuous phenotype
#' scores <- rnorm(100, mean = 50, sd = 10)
#' scores <- PhenoMap(scores, scores > 60 ~ "High", scores > 40 ~ "Medium", .default = "Low")
#'
#' # Example 2: Transform a column in a data frame
#' df <- data.frame(
#'   age = c(25, 35, 45, 55, 65),
#'   gender = c("M", "F", "M", "F", "M")
#' )
#' df <- PhenoMap(df, age < 30 ~ "Young", age < 50 ~ "Middle", .default = "Senior")
#'
#' # Example 3: Multiple conditions with complex logic
#' df <- data.frame(value = c(5, 15, 25, 35, 45))
#' df <- PhenoMap(
#'   df,
#'   value < 10 ~ "Very Low",
#'   value < 20 ~ "Low",
#'   value < 30 ~ "Medium",
#'   value < 40 ~ "High",
#'   .default = "Very High"
#' )
#' }
#'
#' @seealso \code{\link[data.table]{fcase}}, \code{\link[dplyr]{case_when}}
PhenoMap <- function(data, ..., .default = NA) {
  check_installed("dplyr")
  rules <- list2(...)

  if (length(rules) == 0) {
    Abort(
      "Condition is empty",
      tips = "Format e.g.: {.code col > 10 ~ 1, col <= 10 ~ 0}"
    )
  }

  if (!all(vapply(X = rules, FUN = is.call, FUN.VALUE = logical(1)))) {
    Abort(
      "Not all conditions are formula",
      tips = "Use e.g.: {.code col > 10 ~ 1, col <= 10 ~ 0}"
    )
  }

  conditions <- lapply(rules, `[[`, 2)
  values <- lapply(rules, `[[`, 3)

  col <- all.vars(conditions[[1]])[1]

  if (!is_2d(data)) {
    original_names <- names(data)
    res <- dplyr::case_when(
      ...,
      .default = .default,
      .unmatched = "default",
      .ptype = NULL,
      .size = NULL
    )
    names(res) <- original_names
    return(res)
  }

  dt <- data.table::as.data.table(data)

  col <- all.vars(conditions[[1L]])[1L]

  if (is.na(col) || !col %in% names(dt)) {
    Abort(
      "Cannot determine target column from the first condition",
      tips = "Use e.g.: {.code mpg > 15 ~ 1, mpg <= 15 ~ 0}"
    )
  }

  args <- vector("list", length(rules) * 2L)

  for (i in seq_along(rules)) {
    env <- environment(rules[[i]])

    args[[2L * i - 1L]] <- eval(
      conditions[[i]],
      envir = dt,
      enclos = env
    )

    args[[2L * i]] <- eval(
      values[[i]],
      envir = dt,
      enclos = env
    )
  }

  args$default <- .default

  res <- do.call(data.table::fcase, args)

  data.table::set(dt, j = col, value = res)

  dt
}

is_2d <- function(data) {
  if (length(dim(data)) == 2L) TRUE else FALSE
}
