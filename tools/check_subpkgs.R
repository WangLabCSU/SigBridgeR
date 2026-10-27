# tools/check_subpkgs.R
# Install-time, NON-FATAL check for optional sub-packages.
# Called from `configure` / `configure.win`; the working directory is the package root.

owner <- "Exceret"

# Single source of truth: sub-packages listed in DESCRIPTION `Remotes` as "Exceret/<pkg>".
# (Repo name == package name.)  Fallback to a hard-coded vector if the field is absent.
get_subpkgs <- function() {
  fallback <- character() # e.g. c("pkgA", "pkgB")
  if (!file.exists("DESCRIPTION")) {
    return(fallback)
  }

  remotes <- tryCatch(
    read.dcf("DESCRIPTION", fields = "Remotes")[1, 1],
    error = function(e) NA_character_
  )
  if (is.na(remotes)) {
    return(fallback)
  }

  x <- trimws(strsplit(remotes, "[,\n]")[[1]])
  x <- sub("^github::", "", x) # allow "github::owner/repo"
  x <- sub("@.*$", "", x) # drop @ref
  x <- sub("#.*$", "", x) # drop #pr
  x <- x[startsWith(x, paste0(owner, "/"))]
  unique(sub(paste0("^", owner, "/"), "", x))
}

is_installed <- function(pkg) length(find.package(pkg, quiet = TRUE)) > 0L

pkgs <- get_subpkgs()
if (length(pkgs) == 0L) {
  quit(save = "no", status = 0L)
}

missing <- pkgs[!vapply(pkgs, is_installed, logical(1L))]


if (length(missing) > 0L) {
  repos <- sprintf("\"%s/%s\"", owner, missing)
  message(
    "\n* Optional sub-packages not found: ",
    toString(missing),
    "\n",
    "  They are not required to install this package, but the related features\n",
    "  will be unavailable until they are installed. To install them:\n\n",
    "    pak::pak(c(",
    toString(repos),
    "))\n",
    "    # or: remotes::install_github(",
    repos[1L],
    ")\n"
  )
} else {
  message(
    "* All optional sub-packages are installed: ",
    paste(pkgs, collapse = ", ")
  )
}

quit(save = "no", status = 0L) # never fail the installation
