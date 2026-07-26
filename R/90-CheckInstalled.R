#' @keywords internal
CheckInstalled <- function(
  pkg = character(),
  where = c("cran", "bioc", "github"),
  reason = "this function",
  abort = TRUE
) {
  chk::chk_chr(pkg)
  where <- if (grepl("/", pkg)) {
    "github"
  } else {
    SigBridgeRUtils::MatchArg(where, c("cran", "bioc", "github"), NULL)
  }
  pkg_name <- gsub(".*/", "", pkg)
  if (is_installed(pkg_name)) {
    return(invisible(TRUE))
  }

  choice <- utils::menu(
    c("Yes", "Cancel"),
    title = cli::cli_fmt(cli::cli_inform(
      "{cli::symbol$menu} The following package are required for {reason}: {.pkg {pkg_name}}, Do you want to install?"
    ))
  )

  if (choice == 2L) {
    if (abort) {
      Abort("Installation aborted by user, processing canceled")
    }
    return(invisible(FALSE))
  }

  # cran and bioc
  if (where == "cran") {
    utils::install.packages(pkg)
    return(invisible(TRUE))
  }

  if (where == "bioc") {
    check_installed("BiocManager")
    BiocManager::install(pkg)
    return(invisible(TRUE))
  }

  # github
  check_installed("pak")
  pak::pkg_install(pkg)
  invisible(TRUE)
}
