#' @keywords internal
CheckInstalled <- function(
  pkg = character(),
  where = c("cran", "bioc", "github")
) {
  chk::chk_chr(pkg)
  where <- SigBridgeRUtils::MatchArg(where, c("cran", "bioc", "github"))
  pkg_name <- gsub(".*/", "", pkg)
  if (rlang::is_installed(pkg_name)) {
    return(invisible(TRUE))
  }

  choice <- utils::menu(
    c("Yes", "Cancel"),
    title = cli::cli_fmt(cli::cli_alert_info(
      "The following package{?s} are required for this function: {.pkg {pkg_name}}, Do you want to install?"
    ))
  )

  if (choice == 2L) {
    return(invisible(FALSE))
  }

  # cran and bioc
  if (where == "cran") {
    install.packages(pkg)
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
