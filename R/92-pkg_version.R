r_pkg_version <- function(pkg_name) {
  tryCatch(utils::packageVersion(pkg_name), error = NULL)
}

py_pkg_version <- function(pkg_name, conda_env) {
  cmd <- sprintf(
    "conda run -n %s python -c \"import importlib.metadata as m; print(m.version('%s'))\"",
    shQuote(conda_env),
    pkg_name
  )

  out <- suppressWarnings(
    system(cmd, intern = TRUE, ignore.stderr = TRUE)
  )

  if (!is.null(attr(out, "status"))) {
    return(NULL)
  }
  out[1]
}
