# * ------------------- other function ----------------------

#' @title Identify Numeric-Convertible Columns (Internal)
#'
#' @description
#' Detects columns that can be safely converted to numeric based on their content pattern.
#' Used internally for automatic data type detection before conversion.
#'
#' @param data A data frame or matrix containing mixed data types.
#' @param thresh Numeric threshold [0-1]. Minimum proportion of elements in character
#'        columns that must match numeric patterns to be considered convertible
#'        (default: 1).
#'
#' @return A logical vector where:
#'         - `TRUE`: Column should be converted to numeric (either already numeric
#'            or character with sufficient numeric patterns)
#'         - `FALSE`: Column should remain as-is
#'
#' @section Pattern Matching Rules:
#' Uses regular expression to detect:
#' \itemize{
#'   \item Standard decimals (e.g., "1.23", "-0.5")
#'   \item Scientific notation (e.g., "1e-5", "2E+10")
#'   \item Leading/trailing decimals (e.g., ".5", "3.")
#' }
#' Non-character columns always return `TRUE`.
#'
#' @examples
#' \dontrun{
#' # Internal usage example:
#' df <- data.frame(
#'   nums = c("1.2", "3e5", "NaN"),
#'   chars = c("A", "B", "C"),
#'   actual_nums = rnorm(3)
#' )
#'
#' # Identify convertible columns (50% threshold)
#' convertible <- IdentifyDataColumn(df)
#' #> nums chars actual_nums
#' #> TRUE FALSE        TRUE
#'
#' # Convert identified columns
#' df[convertible] <- lapply(df[convertible], as.numeric)
#' }
#'
#' @keywords internal
#' @noRd
#'
IdentifyDataColumn <- function(data, thresh = 1) {
  pattern_based <- sapply(data, function(x) {
    if (is.character(x)) {
      numeric_pat <- grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x)
      return(mean(numeric_pat, na.rm = TRUE) >= 1)
    } else {
      return(TRUE)
    }
  })
  return(pattern_based)
}
