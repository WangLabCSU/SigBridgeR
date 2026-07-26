#' @title Preprocess Phenotype Data
#'
#' @description
#' Aligns and validates phenotype data with bulk RNA-seq expression matrix.
#' Supports three phenotype types (binary, continuous, survival) with automatic
#' sample matching, type validation, and optional conditional mapping via
#' \code{\link{PhenoMap}}.
#'
#' @param bulk A two-dimensional matrix or data frame of bulk RNA-seq expression
#'   data with genes as rows and samples as columns. Must have at least 2 samples.
#' @param phenotype Phenotype data, either:
#'   \itemize{
#'     \item Named numeric vector (for binary/continuous phenotypes)
#'     \item Data frame/matrix with row names matching \code{colnames(bulk)}
#'           (for survival or multi-column phenotypes)
#'   }
#' @param phenotype_class Character. Type of phenotype:
#'   \describe{
#'     \item{\code{"binary"}}{Two-class categorical outcome (e.g., case/control).
#'           Must have exactly 2 unique values.}
#'     \item{\code{"continuous"}}{Continuous measurement (e.g., age, expression).
#'           Must have >2 unique values.}
#'     \item{\code{"survival"}}{Survival data with time and status columns.
#'           Automatically guesses columns named "time" and "status/censor" if
#'           \code{select} is \code{NULL}.}
#'   }
#'   Partial matching is supported.
#' @param ... Conditional mapping rules passed to \code{\link{PhenoMap}} for
#'   transforming phenotype values before validation. Format: \code{condition ~ value}.
#'   Example: \code{status == "death" ~ 1, status == "alive" ~ 0}.
#' @param select Character. Column name(s) to select from \code{phenotype} when
#'   it is two-dimensional. Required for survival data with multiple columns, or
#'   to specify which column to use for binary/continuous phenotypes.
#' @param verbose Logical. Whether to print diagnostic messages. Default: inherits
#'   from \code{getOption("SigBridgeRUtils.verbose")}.
#'
#' @return Preprocessed phenotype data:
#'   \itemize{
#'     \item For binary/continuous: Named numeric vector with sample names as names
#'     \item For survival: Data frame with two columns (time, status) and sample
#'           names as row names
#'   }
#'   Only samples present in both \code{bulk} and \code{phenotype} are retained.
#'
#' @section Validation Rules:
#' \itemize{
#'   \item \strong{Binary}: Must have exactly 2 unique values (e.g., 0/1, TRUE/FALSE)
#'   \item \strong{Continuous}: Must have >2 unique values (to avoid confusion with binary)
#'   \item \strong{Survival}: Must have exactly 2 unique values in status column
#'         (typically 0 = censored, 1 = event)
#'   \item \strong{Sample matching}: Common samples identified via \code{intersect(colnames(bulk), names/rownames(phenotype))}
#' }
#'
#' @section Automatic Features:
#' \itemize{
#'   \item \strong{Survival column guessing}: When \code{phenotype_class = "survival"}
#'         and \code{select = NULL}, automatically detects columns containing "time"
#'         and "status"/"censor" in their names (case-insensitive).
#'   \item \strong{Type conversion}: Logical values are converted to numeric (TRUE->1, FALSE->0)
#'   \item \strong{Sample alignment}: Returns only samples present in both \code{bulk} and \code{phenotype}
#' }
#'
#' @family input_preprocess
#' @export
#' @seealso [PhenoMap()]
#'
#' @examples
#' \dontrun{
#' # Example 1: Binary phenotype
#' bulk_data <- matrix(rpois(100, 10), nrow = 20, ncol = 5)
#' colnames(bulk_data) <- paste0("Sample", 1:5)
#'
#' pheno_binary <- c(Sample1 = 1, Sample2 = 0, Sample3 = 1, Sample4 = 0, Sample5 = 1)
#'
#' result <- PhenoPreProcess(
#'   bulk = bulk_data,
#'   phenotype = pheno_binary,
#'   phenotype_class = "binary"
#' )
#'
#' # Example 2: Continuous phenotype with discretization
#' pheno_age <- c(Sample1 = 25, Sample2 = 35, Sample3 = 45, Sample4 = 55, Sample5 = 65)
#'
#' result <- PhenoPreProcess(
#'   bulk = bulk_data,
#'   phenotype = pheno_age,
#'   phenotype_class = "continuous",
#'   .x < 40 ~ "Young", .x >= 40 ~ "Old"
#' )
#'
#' # Example 3: Survival data with automatic column detection
#' pheno_surv <- data.frame(
#'   time = c(12, 24, 18, 36, 30),
#'   status = c(1, 0, 1, 1, 0),
#'   row.names = paste0("Sample", 1:5)
#' )
#'
#' result <- PhenoPreProcess(
#'   bulk = bulk_data,
#'   phenotype = pheno_surv,
#'   phenotype_class = "survival"
#' )
#'
#' # Example 4: Survival data with explicit column names
#' pheno_surv_custom <- data.frame(
#'   follow_up = c(12, 24, 18, 36, 30),
#'   event = c(1, 0, 1, 1, 0),
#'   row.names = paste0("Sample", 1:5)
#' )
#'
#' result <- PhenoPreProcess(
#'   bulk = bulk_data,
#'   phenotype = pheno_surv_custom,
#'   phenotype_class = "survival",
#'   select = c("follow_up", "event")
#' )
#' }
PhenoPreProcess <- function(
  bulk,
  phenotype,
  phenotype_class = c("binary", "continuous", "survival"),
  ...,
  select = NULL,
  verbose = getFuncOption("verbose")
) {
  lifecycle::deprecate_warn(
    "4.0.0",
    "PhenoPreProcess()",
    "PhenoMap()",
    details = "PhenoPreProcess() will be deprecated because of ambiguous functionality."
  )
  check_installed("tibble")
  if (!is.null(select)) {
    chk::chk_character(select)
  }

  is_phenotype_2d <- is_2d(phenotype)

  dots <- list2(...)
  if (length(dots) > 0) {
    phenotype <- PhenoMap(phenotype, ...)

    if (is_phenotype_2d) {
      conditions <- lapply(dots, `[[`, 2)
      select <- all.vars(conditions[[1]])[1]
    }
  }

  # Validate phenotype_class
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("binary", "continuous", "survival"),
    NULL
  )

  # Validate bulk is a two-dimensional matrix with genes as rows and samples as columns
  if (!is_2d(bulk)) {
    cli::cli_abort(c(
      "x" = "`bulk` must be a 2-dimensional matrix",
      ">" = "Current type is {.cls {class(bulk)}}"
    ))
  }
  n_samples <- ncol(bulk)
  if (n_samples < 2L) {
    cli::cli_abort(c(
      "x" = "`bulk` must have at least 2 samples (columns)"
    ))
  }

  bulk <- as.matrix(bulk)
  sample_names <- colnames(bulk)

  if (verbose) {
    b_dimension <- dim(bulk)

    cli::cli_alert_info(
      "Bulk dimension: {.pkg ({b_dimension[1]},{b_dimension[2]})}"
    )

    if (is_phenotype_2d) {
      p_dimension <- dim(phenotype)

      cli::cli_alert_info(
        "Phenotype dimension: {.pkg ({p_dimension[1]},{p_dimension[2]})}"
      )
    } else {
      p_length <- length(phenotype)
      cli::cli_alert_info(
        "Phenotype length: {.pkg {p_length}}"
      )
    }
  }

  if (is_phenotype_2d) {
    common_samples <- if (n_samples != nrow(phenotype)) {
      if (verbose) {
        cli::cli_alert_info("Matching `bulk` and `phenotype`")
      }
      intersect(sample_names, rownames(phenotype))
    } else {
      sample_names
    }
    if (length(common_samples) == 0L) {
      cli::cli_abort(c(
        "x" = "No common sample names between `bulk` and `phenotype`",
        ">" = "Check the colnames of `bulk` and rownames of `phenotype`"
      ))
    }
  } else {
    # vector
    common_samples <- if (n_samples != length(phenotype)) {
      if (verbose) {
        cli::cli_alert_info("Matching `bulk` and `phenotype`")
      }
      intersect(sample_names, names(phenotype))
    } else {
      sample_names
    }
    if (length(common_samples) == 0L) {
      cli::cli_abort(c(
        "x" = "No common sample names between `bulk` and `phenotype`",
        ">" = "Check the colnames of `bulk` and names of `phenotype`"
      ))
    }
  }

  # Determine if phenotype is two-dimensional
  result <- handle_phenotype_select(
    is_phenotype_2d = is_phenotype_2d,
    is_select_null = is.null(select),
    phenotype_class = phenotype_class,
    phenotype = phenotype,
    sample_names = common_samples,
    select = select
  )

  # Convert phenotype to proper format based on phenotype_class

  if (verbose) {
    cli::cli_alert_success("Phenotype preprocessing completed")
  }

  result
}

handle_phenotype_select <- function(
  is_phenotype_2d,
  is_select_null,
  phenotype_class = c("binary", "continuous", "survival"),
  phenotype,
  sample_names,
  select = NULL,
  ...
) {
  # 分派执行
  if (is_phenotype_2d) {
    if (is.null(select)) {
      handle_case_1(
        phenotype_class = phenotype_class,
        phenotype = phenotype,
        sample_names = sample_names,
        select = select
      )
    } else {
      handle_case_2(
        phenotype_class = phenotype_class,
        phenotype = phenotype,
        sample_names = sample_names,
        select = select
      )
    }
  } else {
    if (is.null(select)) {
      handle_case_3(
        phenotype_class = phenotype_class,
        phenotype = phenotype,
        sample_names = sample_names,
        select = select
      )
    } else {
      handle_case_4(
        phenotype_class = phenotype_class,
        phenotype = phenotype,
        sample_names = sample_names,
        select = select
      )
    }
  }
}

#' @keywords internal
col_count <- function(data, col) {
  length(table(data[[col]]))
}

#' @keywords internal
col_class <- function(data, col) {
  unique(vapply(X = data[[col]], FUN = base::class, FUN.VALUE = character(1L)))
}

#' @keywords internal
handle_case_1 <- function(
  phenotype_class = c("binary", "continuous", "survival"),
  phenotype,
  sample_names,
  ...
) {
  # TRUE, NULL: 处理逻辑
  if (phenotype_class == "survival") {
    cli::cli_alert_info("Guess survival time and status from data")
    possible_cols <- tolower(colnames(phenotype))
    time_col <- grepv("time", possible_cols)
    status_col <- grepv("status|censor|event", possible_cols)

    if (length(time_col) > 1 || length(status_col) > 1) {
      cli::cli_abort(c(
        "x" = "Unable to guess time and status columns, multiple columns found",
        ">" = "Try specify them"
      ))
    }

    select <- c(time_col, status_col)

    if (length(unique(unlist(phenotype[status_col]))) != 2L) {
      cli::cli_abort(c(
        "x" = "Status column must have exactly 2 unique values",
        ">" = "Current guessed columns: {.val {time_col}} {.val {status_col}}"
      ))
    }

    return(phenotype[sample_names, select]) # df
  }

  # "binary", "continuous"
  if (ncol(phenotype) == 1) {
    col_val_count <- col_count(phenotype, 1L)
    if (phenotype_class == "binary" && col_val_count > 2L) {
      cli::cli_abort(c(
        "x" = "Binary phenotype must have exactly 2 unique values",
        ">" = "Currenly has {col_val_count} unique values"
      ))
    }

    if (phenotype_class == "continuous" && col_val_count == 2L) {
      cli::cli_abort(c(
        "x" = "Continuous phenotype must have more than 2 unique values",
        ">" = "Currenly has exactly 2 unique values, \
          if data is correct, specify {.code phenotype_class = \"binary\"}"
      ))
    }

    if (col_val_count == 1L) {
      cli::cli_abort(c("x" = "Only one unique value found"))
    }

    col_val_class <- col_class(phenotype, 1)
    if (col_val_class == "character") {
      cli::cli_abort(c("x" = "Must be numeric"))
    } else if (col_val_class == "logical") {
      return(stats::setNames(
        as.numeric(unlist(phenotype[sample_names, 1])),
        sample_names
      ))
    }

    # vec
    return(stats::setNames(
      unlist(phenotype[sample_names, 1]),
      sample_names
    ))
  } else {
    cli::cli_abort(c(
      "x" = "Multiple columns found but `select` is NULL",
      ">" = "Use `select = <colname>` to clearly specified"
    ))
  }
}

#' @keywords internal
handle_case_2 <- function(
  phenotype_class = c("binary", "continuous", "survival"),
  phenotype,
  sample_names,
  select,
  ...
) {
  # TRUE, 非NULL: 处理逻辑
  if (phenotype_class %chin% c("binary", "continuous")) {
    cli::cli_warn(
      "`select` is only used when `phenotype` is two-dimensional. Ignoring `select`."
    )
    return(handle_case_1(
      phenotype_class = phenotype_class,
      phenotype = phenotype,
      sample_names = sample_names,
    ))
  }

  # * survival
  chk::chk_length(select, 2)
  col_exists <- select %chin% colnames(phenotype)
  if (!all(col_exists)) {
    cli::cli_abort(c(
      "x" = "Column {.val {colnames(phenotype)[!col_exists]}} not found in `phenotype`"
    ))
  }

  return(phenotype[sample_names, select])
}

#' @keywords internal
handle_case_3 <- function(
  phenotype_class = c("binary", "continuous", "survival"),
  phenotype,
  sample_names,
  ...
) {
  # FALSE, NULL: 处理逻辑
  if (phenotype_class == "survival") {
    cli::cli_abort(c(
      "x" = "Invalid type of `phenotype`",
      ">" = "Current type: {.cls {class(phenotype)}}, expect: {.cls {c('surv', 'data.frame')}}"
    ))
  }

  val_class <- unique(vapply(phenotype, class, character(1)))
  if (length(val_class) > 1 || val_class == "character") {
    cli::cli_abort(c(
      "x" = "Invalid type of `phenotype`",
      ">" = "Current type: {.cls {val_class}}, expect: {.cls {c('numeric', 'integer')}}"
    ))
  }

  # continuous & binary
  return(phenotype[sample_names])
}

#' @keywords internal
handle_case_4 <- function(
  phenotype_class = c("binary", "continuous", "survival"),
  phenotype,
  sample_names,
  select,
  ...
) {
  # FALSE, 非NULL: 处理逻辑
  if (phenotype_class == "survival") {
    cli::cli_abort(c(
      "x" = "Invalid type of `phenotype`",
      ">" = "Current type: {.cls {class(phenotype)}}"
    ))
  }

  if (phenotype_class %chin% c("binary", "continuous")) {
    cli::cli_warn(
      "`select` is only used when `phenotype` is two-dimensional. Ignoring `select`."
    )
    return(handle_case_3(
      phenotype_class = phenotype_class,
      phenotype = phenotype,
      sample_names = sample_names,
    ))
  }
}
