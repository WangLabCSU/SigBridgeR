# rexpand.R -- build-time R source expander (R-src/ -> R/), rule based.
# Base R only; `names_from = "eval"` additionally needs pkgload.
#
#   source("dev/rexpand.R")
#   rexpand("R-src", out = "R", min_len = 2)        # write R/
#   rexpand("R-src")                                # dry run, print edits
#   rexpand("R-src", names_from = "eval")           # names via pkgload
#   devtools::document()
#
# Built-in rule `rule_extend_vec()`:
#   1. collect the names of the environment `ScreenStrategy`
#   2. every `c("a", "b", ...)` made only of string literals, all of which are
#      names of that environment (case-sensitive), with length >= min_len,
#      is rewritten in place to the original elements followed by the missing
#      names. Already complete vectors are left alone (idempotent).
#   3. put `# rexpand: keep` on the line of a `c(` to exclude it.
#
# Extend: a rule is list(prepare = function(files), edits = function(ctx)).
# `files` is a named list of source text; `ctx` holds path, text, pd (parse
# data). `edits` returns a data.frame(start, end, new) of 1-based character
# offsets into `text`. `min_len`, `names_from` and `pkg` only configure the
# default rule; pass `rules =` to use your own.

# ---------------------------------------------------------------- engine ----
rexpand <- function(
  src = "R-src",
  out = NULL,
  rules = NULL,
  min_len = 2L,
  names_from = c("static", "eval"),
  pkg = "."
) {
  names_from <- match.arg(names_from)
  if (is.null(rules)) {
    rules <- list(rule_extend_vec(min_len = min_len, names_from = names_from, pkg = pkg))
  }
  paths <- list.files(src, "\\.[Rr]$", full.names = TRUE)
  if (!length(paths)) {
    stop("no .R files found under `", src, "`", call. = FALSE)
  }
  files <- stats::setNames(
    lapply(paths, function(p) paste(readLines(p, encoding = "UTF-8", warn = FALSE), collapse = "\n")),
    paths
  )
  for (r in rules) {
    r$prepare(files)
  }

  n_edit <- 0L
  if (!is.null(out)) {
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
  }
  for (p in paths) {
    text <- files[[p]]
    pd <- utils::getParseData(parse(text = text, keep.source = TRUE))
    if (is.null(pd)) {
      next
    }
    ctx <- list(path = p, text = text, pd = pd, lines = strsplit(text, "\n", fixed = TRUE)[[1]])
    ed <- do.call(rbind, lapply(rules, function(r) r$edits(ctx)))
    if (!is.null(ed) && nrow(ed)) {
      ed <- ed[order(-ed$start), ]
      if (is.null(out)) {
        for (i in rev(seq_len(nrow(ed)))) {
          cat(sprintf("%s\n  - %s\n  + %s\n", basename(p), substr(text, ed$start[i], ed$end[i]), ed$new[i]))
        }
      }
      for (i in seq_len(nrow(ed))) {
        text <- paste0(substr(text, 1L, ed$start[i] - 1L), ed$new[i], substr(text, ed$end[i] + 1L, nchar(text)))
      }
      n_edit <- n_edit + nrow(ed)
    }
    if (!is.null(out)) {
      writeLines(
        text,
        file.path(out, basename(p)),
        useBytes = TRUE
      )
    }
  }
  message(n_edit, " replacement(s)")
  invisible(n_edit)
}

# ------------------------------------------------------------- utilities ----
# (line, col) from getParseData -> 1-based character offset in the whole text.
# getParseData columns expand tabs like the R lexer, so map them back.
pos_to_offset <- function(ctx, line, col) {
  starts <- cumsum(c(0L, nchar(ctx$lines) + 1L))
  ln <- ctx$lines[line]
  if (is.na(ln)) {
    return(NA_integer_)
  }
  if (grepl("\t", ln, fixed = TRUE)) {
    ch <- strsplit(ln, "")[[1]]
    cur <- 1L
    for (k in seq_along(ch)) {
      if (cur >= col) {
        col <- k
        break
      }
      cur <- if (ch[k] == "\t") {
        ((cur - 1L) %/% 8L + 1L) * 8L + 1L
      } else {
        cur + 1L
      }
    }
  }
  starts[line] + col
}

# Static names of `var`: var$x <- , var[["x"]] <- , assign("x", , envir = var),
# and constructor calls such as var <- new_environment(list(a = , b = )).
static_env_names <- function(pd, var) {
  if (is.null(pd) || !nrow(pd)) {
    return(character())
  }
  out <- character()
  term <- pd[pd$terminal, ]
  term <- term[order(term$line1, term$col1), ]
  unq <- function(s) gsub("^[\"'`]|[\"'`]$", "", s)
  is_assign <- function(i) {
    # token after row i is `<-` / `=`
    i + 1L <= nrow(term) && term$token[i + 1L] %in% c("LEFT_ASSIGN", "EQ_ASSIGN")
  }
  kids_of <- function(id) {
    k <- pd[pd$parent == id, , drop = FALSE]
    k[order(k$line1, k$col1), , drop = FALSE]
  }
  # name of the function called by the call node `id`
  call_name <- function(id) {
    k <- kids_of(id)
    fn <- k[k$token == "expr", , drop = FALSE]
    if (!nrow(fn)) {
      return("")
    }
    fk <- kids_of(fn$id[1L])
    sym <- fk[fk$token == "SYMBOL_FUNCTION_CALL", , drop = FALSE]
    if (!nrow(sym)) {
      return("")
    }
    sym$text[nrow(sym)]
  }
  call_args <- function(id) {
    k <- kids_of(id)
    k$id[k$token == "expr"]
  }
  # names written as `name = value` directly in a call
  call_arg_names <- function(id) {
    k <- kids_of(id)
    nms <- character()
    for (j in seq_len(nrow(k))) {
      if (k$token[j] %in% c("SYMBOL_SUB", "STR_CONST") && j < nrow(k) && k$token[j + 1L] == "EQ_SUB") {
        nms <- c(nms, unq(k$text[j]))
      }
    }
    unique(nms)
  }
  # `new_environment(list(a = , b = ))` wraps the names one level down; only
  # descend into list-like constructors, not into nested method constructors.
  ctor_names <- function(id, depth = 0L) {
    fun <- call_name(id)
    nms <- character()
    if (fun %in% c("list", "c")) {
      nms <- call_arg_names(id)
    } else if (fun %in% c("new_environment", "new.env", "env")) {
      nms <- setdiff(call_arg_names(id), c("data", "parent", "hash", "size"))
      if (depth < 2L) {
        for (a in call_args(id)) {
          if (call_name(a) %in% c("list", "c")) {
            nms <- c(nms, ctor_names(a, depth + 1L))
          }
        }
      }
    }
    unique(nms)
  }
  for (i in seq_len(nrow(term))) {
    if (term$token[i] != "SYMBOL" || term$text[i] != var) {
      next
    }
    t1 <- if (i + 1L <= nrow(term)) term$text[i + 1L] else NA_character_
    t2 <- if (i + 2L <= nrow(term)) term$text[i + 2L] else NA_character_
    if (
      identical(t1, "$") &&
        i + 2L <= nrow(term) &&
        term$token[i + 2L] %in% c("SYMBOL", "STR_CONST") &&
        is_assign(i + 2L)
    ) {
      out <- c(out, unq(t2))
    } else if (
      identical(t1, "[[") &&
        i + 4L <= nrow(term) &&
        term$token[i + 2L] == "STR_CONST" &&
        is_assign(i + 4L)
    ) {
      out <- c(out, unq(t2))
    }
  }
  # var <- list(...) / var <- new_environment(list(...))
  for (i in seq_len(nrow(term))) {
    if (term$token[i] != "SYMBOL" || term$text[i] != var || !is_assign(i)) {
      next
    }
    sym_expr <- pd$parent[pd$id == term$id[i]]
    if (length(sym_expr) != 1L || is.na(sym_expr)) {
      next
    }
    asg <- pd$parent[pd$id == sym_expr]
    if (length(asg) != 1L || is.na(asg)) {
      next
    }
    k <- kids_of(asg)
    op <- which(k$token %in% c("LEFT_ASSIGN", "EQ_ASSIGN"))
    if (!length(op)) {
      next
    }
    rhs <- k[
      k$token == "expr" &
        (k$line1 > k$line2[op[1L]] |
          (k$line1 == k$line2[op[1L]] & k$col1 > k$col2[op[1L]])),
      ,
      drop = FALSE
    ]
    if (!nrow(rhs)) {
      next
    }
    out <- c(out, ctor_names(rhs$id[1L]))
  }
  # assign("x", value, envir = var)
  for (i in which(term$token == "SYMBOL_FUNCTION_CALL" & term$text == "assign")) {
    if (i + 2L > nrow(term) || term$token[i + 2L] != "STR_CONST") {
      next
    }
    j <- i + 3L
    depth <- 1L
    hit <- FALSE
    while (j <= nrow(term) && depth > 0L) {
      if (term$text[j] == "(") {
        depth <- depth + 1L
      }
      if (term$text[j] == ")") {
        depth <- depth - 1L
      }
      if (term$token[j] == "SYMBOL" && term$text[j] == var) {
        hit <- TRUE
      }
      j <- j + 1L
    }
    if (hit) out <- c(out, unq(term$text[i + 2L]))
  }
  unique(out)
}

wrap_vec <- function(vals, indent, width = 80L) {
  q <- paste0("\"", gsub("([\\\\\"])", "\\\\\\1", vals), "\"")
  one <- paste0("c(", paste(q, collapse = ", "), ")")
  if (nchar(one) + indent <= width) {
    return(one)
  }
  pad <- strrep(" ", indent + 2L)
  lines <- character()
  cur <- ""
  for (i in seq_along(q)) {
    piece <- paste0(q[i], if (i < length(q)) "," else "")
    if (nzchar(cur) && nchar(pad) + nchar(cur) + 1L + nchar(piece) > width) {
      lines <- c(lines, cur)
      cur <- piece
    } else {
      cur <- if (nzchar(cur)) paste(cur, piece) else piece
    }
  }
  paste0("c(\n", paste0(pad, c(lines, cur), collapse = "\n"), "\n", strrep(" ", indent), ")")
}

# ------------------------------------------------------------------ rule ----
rule_extend_vec <- function(var = "ScreenStrategy", min_len = 2L, names_from = c("static", "eval"), pkg = ".") {
  names_from <- match.arg(names_from)
  nm <- character()
  list(
    prepare = function(files) {
      nm <<- if (names_from == "eval") {
        pkgload::load_all(pkg, quiet = TRUE)
        ls(get(var, envir = asNamespace(read.dcf(file.path(pkg, "DESCRIPTION"), "Package")[1])), all.names = TRUE)
      } else {
        unique(unlist(lapply(files, function(t) {
          static_env_names(utils::getParseData(parse(text = t, keep.source = TRUE)), var)
        })))
      }
      if (!length(nm)) {
        stop("no names found for `", var, "`; try names_from = 'eval'", call. = FALSE)
      }
      message("[", var, "] ", length(nm), " names: ", paste(nm, collapse = ", "))
    },
    edits = function(ctx) {
      pd <- ctx$pd
      calls <- pd$parent[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "c"]
      rows <- list()
      for (e in pd$parent[pd$id %in% calls]) {
        # e: the call's expr node
        kids <- pd[pd$parent == e, ]
        kids <- kids[order(kids$line1, kids$col1), ]
        if (kids$token[2] != "'('" || kids$token[nrow(kids)] != "')'") {
          next
        }
        args <- kids[-c(1:2, nrow(kids)), ]
        isarg <- args$token == "expr"
        if (!all(args$token[!isarg] == "','")) {
          next
        }
        vals <- character()
        ok <- sum(isarg) > 0
        for (a in args$id[isarg]) {
          k <- pd[pd$parent == a, ]
          if (nrow(k) != 1L || k$token != "STR_CONST" || grepl("\\\\", k$text) || grepl("^[rR]", k$text)) {
            ok <- FALSE
            break
          }
          vals <- c(vals, substr(k$text, 2L, nchar(k$text) - 1L))
        }
        if (!ok || length(vals) < min_len || !all(vals %in% nm)) {
          next
        }
        add <- setdiff(nm, vals)
        if (!length(add)) {
          next
        } # already complete
        row <- pd[pd$id == e, ]
        if (grepl("rexpand: keep", ctx$lines[row$line1], fixed = TRUE)) {
          next
        }
        s <- pos_to_offset(ctx, row$line1, row$col1)
        t <- pos_to_offset(ctx, row$line2, row$col2)
        chunk <- substr(ctx$text, s, t) # safety check on offsets
        if (!startsWith(chunk, "c(") || !endsWith(chunk, ")")) {
          next
        }
        indent <- nchar(sub("^(\\s*).*", "\\1", ctx$lines[row$line1]))
        rows[[length(rows) + 1L]] <- data.frame(
          start = s,
          end = t,
          new = wrap_vec(c(unique(vals), add), indent)
        )
      }
      if (length(rows)) do.call(rbind, rows) else NULL
    }
  )
}
