#' @inherit zeallot::`%<-%`
#' @export
`%<-%` <- function(lhs, rhs) {
  lhs_expr <- substitute(lhs)
  env <- parent.frame()

  plan <- compile_lhs(lhs_expr)

  fast_assign_plan(plan, rhs, env)

  invisible(rhs)
}

compile_lhs <- function(expr, path = integer()) {
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    if (nm %in% c("_", ".")) {
      return(list(list(type = "ignore", path = path)))
    }

    return(list(list(type = "assign", name = nm, path = path)))
  }

  if (is.call(expr) && identical(expr[[1]], as.name("c"))) {
    args <- as.list(expr[-1])
    out <- list()

    for (i in seq_along(args)) {
      out <- c(out, compile_lhs(args[[i]], c(path, i)))
    }

    return(out)
  }

  stop("Invalid destructuring target")
}
