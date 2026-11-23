# ---
# repo: WangLabCSU/SigBridgeR
# file: standalone-args_to_roxy_doc.R
# last-updated: 2026-11-22
# license: https://unlicense.org
# ---
#
#
# nocov start

#' Extract a function's roxygen2 argument documentation (`@param`) as Markdown
#'
#' Reads the `\arguments` section of a function's help page (Rd) and returns a
#' character vector whose elements are Markdown list items of the form
#' "- `arg1`: desc1" and "- `arg2`: desc2".
#' Supports both packages under development that have been documented with
#' `devtools::document()` (read from `man/*.Rd` in the source directory) and
#' packages installed on the system.
#'
#' @param fn Function name (a string) or a function object, e.g. `"lm"`,
#'   `stats::lm`, `my_fun`.
#' @param package Package name. When `NULL`, it is inferred automatically:
#'   first from the function object, then from the package under development,
#'   and finally from the help database of loaded/installed packages.
#' @param path Root directory of the package under development (containing
#'   `DESCRIPTION` and `man/`); defaults to the current directory.
#'
#' @return A character vector with one Markdown unordered-list item per
#'   element; returns `character(0)` when no documentation is found.
#' @noRd
#'
#' @examples
#' args_to_roxy_doc("lm", package = "stats")
#' args_to_roxy_doc(stats::sd)
args_to_roxy_doc <- function(fn, package = NULL, path = ".") {
  # ---- Resolve the function name and package name ----
  if (is.function(fn)) {
    fn_name <- sub("^.*::+", "", paste(deparse(substitute(fn)), collapse = ""))
    if (is.null(package)) {
      ns <- environment(fn)
      if (!is.null(ns)) {
        top <- topenv(ns)
        if (isNamespace(top)) package <- getNamespaceName(top)[[1L]]
      }
    }
  } else if (is.character(fn) && length(fn) == 1L) {
    fn_name <- sub("^.*::+", "", fn)
    if (is.null(package) && grepl("::", fn, fixed = TRUE)) {
      package <- sub("::+.*$", "", fn)
    }
  } else {
    stop("`fn` must be a function or a single string.", call. = FALSE)
  }

  # ---- Candidate Rd source 1: package under development (man/*.Rd) ----
  rd_list <- NULL
  desc_file <- file.path(path, "DESCRIPTION")
  man_dir <- file.path(path, "man")
  if (file.exists(desc_file) && dir.exists(man_dir)) {
    dev_pkg <- read.dcf(desc_file, fields = "Package")[1L, 1L]
    if (is.null(package) || identical(package, dev_pkg)) {
      rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
      rd_list <- lapply(rd_files, function(f) {
        tryCatch(
          tools::parse_Rd(f, encoding = "UTF-8"),
          error = function(e) NULL
        )
      })
      rd_list <- Filter(Negate(is.null), rd_list)
      rd <- find_rd_by_alias(rd_list, fn_name)
      if (!is.null(rd)) return(rd_arguments_to_md(rd))
    }
  }

  # ---- Candidate Rd source 2: installed package ----
  if (is.null(package)) {
    # Look up the owning package via help()
    hp <- suppressWarnings(utils::help(fn_name))
    if (length(hp) > 0L) {
      package <- basename(dirname(dirname(as.character(hp)[1L])))
    }
  }
  if (is.null(package)) {
    warning("Cannot determine the package of `", fn_name, "`.", call. = FALSE)
    return(character(0L))
  }
  # After load_all() on a development package, its install-path info points
  # back to the source directory
  if (
    !requireNamespace(package, quietly = TRUE) &&
      !nzchar(system.file(package = package))
  ) {
    warning("Package `", package, "` is not installed.", call. = FALSE)
    return(character(0L))
  }
  db <- get_rd_db(package)
  if (is.null(db)) {
    warning("Cannot read Rd documentation of package `", package, "`.", call. = FALSE)
    return(character(0L))
  }
  rd <- find_rd_by_alias(db, fn_name)
  if (is.null(rd)) {
    warning("No Rd documentation found for `", fn_name, "` in package `", package, "`.", call. = FALSE)
    return(character(0L))
  }
  rd_arguments_to_md(rd)
}

# Find an Rd object by its \alias or \name
find_rd_by_alias <- function(rd_list, topic) {
  for (rd in rd_list) {
    tags <- vapply(rd, function(x) attr(x, "Rd_tag") %||% "", character(1L))
    keys <- vapply(
      rd[tags %in% c("\\alias", "\\name")],
      function(x) trimws(paste(unlist(x), collapse = "")),
      character(1L)
    )
    if (topic %in% keys) return(rd)
  }
  NULL
}

# Extract \arguments from an Rd object and return a Markdown vector
rd_arguments_to_md <- function(rd) {
  tags <- vapply(rd, function(x) attr(x, "Rd_tag") %||% "", character(1L))
  args_sec <- rd[tags == "\\arguments"]
  if (length(args_sec) == 0L) {
    return(character(0L))
  }

  out <- character(0L)
  for (item in args_sec[[1L]]) {
    if (!identical(attr(item, "Rd_tag"), "\\item") || length(item) < 2L) {
      next
    }
    arg_names <- rd_to_md(item[[1L]])
    desc <- rd_to_md(item[[2L]])
    desc <- trimws(gsub("\\s+", " ", desc))
    # A single \item may document several arguments: "x, y"
    arg_names <- trimws(strsplit(arg_names, ",", fixed = TRUE)[[1L]])
    arg_names <- arg_names[nzchar(arg_names)]
    out <- c(out, sprintf("- `%s`: %s", arg_names, desc))
  }
  out
}

# Recursively convert an Rd fragment to Markdown text
rd_to_md <- function(x) {
  tag <- attr(x, "Rd_tag") %||% ""
  if (tag %in% c("\\dots", "\\ldots")) {
    return("...")
  }
  if (!is.list(x)) {
    return(paste(as.character(x), collapse = ""))
  }
  inner <- paste(vapply(x, rd_to_md, character(1L)), collapse = "")
  switch(
    tag,
    "\\code" = ,
    "\\verb" = ,
    "\\samp" = ,
    "\\env" = ,
    "\\file" = ,
    "\\pkg" = ,
    "\\command" = ,
    "\\option" = paste0("`", inner, "`"),
    "\\emph" = ,
    "\\em" = paste0("*", inner, "*"),
    "\\strong" = ,
    "\\bold" = paste0("**", inner, "**"),
    "\\href" = {
      if (length(x) >= 2L) {
        sprintf("[%s](%s)", rd_to_md(x[[2L]]), rd_to_md(x[[1L]]))
      } else {
        inner
      }
    },
    "\\item" = paste0(" ", inner),
    inner
  )
}

get_rd_db <- function(package, path = ".") {
  loaded <- tryCatch(find.package(package), error = function(e) NULL)
  for (d in unique(c(loaded, if (file.exists(file.path(path, "DESCRIPTION"))) normalizePath(path)))) {
    if (!is.null(d) && dir.exists(file.path(d, "man"))) {
      db <- tryCatch(tools::Rd_db(dir = d), error = function(e) NULL)
      if (!is.null(db)) return(db)
    }
  }

  inst <- tryCatch(
    find.package(package, lib.loc = .libPaths()),
    error = function(e) NULL
  )
  if (!is.null(inst) && file.exists(file.path(inst, "help", paste0(package, ".rdx")))) {
    db <- tryCatch(
      tools::Rd_db(package, lib.loc = dirname(inst)),
      error = function(e) NULL
    )
    if (!is.null(db)) return(db)
  }
  NULL
}

# nocov end
