r_pkg_version <- function(pkg_name) {
  tryCatch(as.character(utils::packageVersion(pkg_name)), error = NULL)
}

get_pkg_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("SigBridgeR")),
    error = function(e) {
      # devtools::load_all() 会将工作目录切换至包根目录
      read.dcf("DESCRIPTION")[, "Version"]
    }
  )
}

py_pkg_version <- function(pkg_name, conda_env) {
  # 1. Prefer conda metadata: conda list -n ENV PKG
  conda_out <- tryCatch(
    system2(
      command = "conda",
      args = c(
        "list",
        "-n",
        shQuote(conda_env),
        shQuote(pkg_name)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character(0L)
  )

  conda_out <- conda_out[!grepl("^#", trimws(conda_out))]

  if (length(conda_out) == 0L) {
    return(NULL)
  }

  fields <- strsplit(trimws(conda_out), "\\s+") # list

  hit <- fields[vapply(
    X = fields,
    FUN = function(x) {
      length(x) >= 2L && identical(x[[1L]], pkg_name)
    },
    FUN.VALUE = logical(1L)
  )]

  if (length(hit) > 0L) {
    return(hit[[1L]][[2L]])
  }
  NULL
}
