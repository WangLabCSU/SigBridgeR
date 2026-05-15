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
  if (rlang::is_installed(pkg_name)) {
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
      cli::cli_abort(c(
        "x" = "Installation aborted by user, processing canceled"
      ))
    }
    return(invisible(FALSE))
  }

  # cran and bioc
  if (where == "cran") {
    utils::install.packages(pkg)
    return(invisible(TRUE))
  }

  if (where == "bioc") {
    rlang::check_installed("BiocManager")
    BiocManager::install(pkg)
    return(invisible(TRUE))
  }

  # github
  if (rlang::is_installed("pak")) {
    pak::pkg_install(pkg)
    return(invisible(TRUE))
  }
  if (rlang::is_installed("remotes")) {
    remotes::install_github(pkg)
    return(invisible(TRUE))
  }
  cli::cli_abort(c(
    "x" = "Package {.pkg pak}, {.pkg BiocManager} or {.pkg remotes} is required for installation!"
  ))
}
