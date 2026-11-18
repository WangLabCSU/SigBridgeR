# rexpand.R -- build-time R source expander (R-src/ -> R/), rule based. Base R only.
#
#   source("dev/rexpand.R")
#   rexpand("R-src", out = "R", min_len = 2)        # write R/
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
# offsets into `text`.

# ---------------------------------------------------------------- engine ----
rexpand <- function(src = "R-src", out = NULL, rules = list(rule_extend_vec()), min_len = NULL) {
  if (!is.null(min_len)) {
    rules <- list(rule_extend_vec(min_len = min_len))
  }
  paths <- list.files(src, "\\.[Rr]$", full.names = TRUE)
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
        c("# Generated from R-src/ by rexpand.R, do not edit", text),
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
  if (grepl("\t", ln, fixed = TRUE)) {
    ch <- strsplit(ln, "")[[1]]
    cur <- 0L
    for (k in seq_along(ch)) {
      cur <- cur + 1L
      if (ch[k] == "\t") {
        cur <- bitwAnd(cur + 6L, bitwNot(7L)) + 1L
      }
      if (cur >= col) {
        col <- k
        break
      }
    }
  }
  starts[line] + col
}

# Static names of `var`: var$x <- , var[["x"]] <- , assign("x", , envir = var)
static_env_names <- function(pd, var) {
  out <- character()
  term <- pd[pd$terminal, ]
  term <- term[order(term$line1, term$col1), ]
  unq <- function(s) gsub("^[\"'`]|[\"'`]$", "", s)
  is_assign <- function(i) {
    # token after row i is `<-` / `=`
    i + 1L <= nrow(term) && term$token[i + 1L] %in% c("LEFT_ASSIGN", "EQ_ASSIGN")
  }
  for (i in seq_len(nrow(term))) {
    if (term$token[i] != "SYMBOL" || term$text[i] != var) {
      next
    }
    t1 <- term$text[i + 1L]
    t2 <- term$text[i + 2L]
    if (identical(t1, "$") && term$token[i + 2L] %in% c("SYMBOL", "STR_CONST") && is_assign(i + 2L)) {
      out <- c(out, unq(t2))
    } else if (identical(t1, "[[") && term$token[i + 2L] == "STR_CONST" && is_assign(i + 4L)) {
      out <- c(out, unq(t2))
    }
  }
  # assign("x", value, envir = var)
  for (i in which(term$token == "SYMBOL_FUNCTION_CALL" & term$text == "assign")) {
    if (term$token[i + 2L] != "STR_CONST") {
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
        stop("no names found for `", var, "`; try names_from = 'eval'")
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
