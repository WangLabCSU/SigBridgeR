# Install every CRAN / Bioconductor dependency declared in DESCRIPTION.
#
# The R-CMD-check workflow used to hard-code this list, so every DESCRIPTION
# update also required editing the workflow. This script derives the list
# instead: GitHub-only packages are dropped (they cannot be resolved from
# CRAN/Bioconductor) and everything else is installed with pak.
#
# Two exclusion layers keep the list self-maintaining:
#   1. `github_only` below -- packages known to live on GitHub only. Extend it
#      when you add such a dependency; the comments document where it came from.
#   2. Repository lookup -- whatever is absent from CRAN/Bioconductor is
#      skipped automatically, so an unlisted GitHub dependency never breaks the
#      workflow.
#
# Usage:
#   Rscript tools/install_cran_bioc_deps.R            # install missing deps
#   Rscript tools/install_cran_bioc_deps.R --dry-run  # report only, no install
#   Rscript tools/install_cran_bioc_deps.R --list     # print package names only
#   Rscript tools/install_cran_bioc_deps.R --no-strict
#       # exit with status 0 even if some package fails to install
#
# Environment variables:
#   DESCRIPTION_PATH  path of DESCRIPTION (default: "DESCRIPTION")
#   CRAN_REPO         CRAN mirror (default: getOption("repos")[["CRAN"]], then
#                     https://cloud.r-project.org)
#   BIOC_REPO         Bioconductor repository URL (default: the Posit Package
#                     Manager mirror, falling back to bioconductor.org)
#   EXTRA_REPOS       comma-separated extra repositories, e.g. r-universe

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(x) x %in% args

list_only <- has_flag("--list")
dry_run <- list_only || has_flag("--dry-run")
strict <- !has_flag("--no-strict")

desc_path <- Sys.getenv("DESCRIPTION_PATH", unset = "DESCRIPTION")
if (!file.exists(desc_path)) {
  stop("DESCRIPTION not found: ", desc_path, call. = FALSE)
}

## ---- declared dependencies ------------------------------------------------

read_deps <- function(path) {
  dcf <- read.dcf(path)
  fields <- intersect(
    c("Depends", "Imports", "LinkingTo", "Suggests"),
    colnames(dcf)
  )
  raw <- unlist(lapply(fields, function(f) strsplit(dcf[1L, f], ",")[[1L]]))
  raw <- trimws(sub("\\(.*", "", raw)) # drop version constraints
  unique(raw[nzchar(raw) & raw != "R"])
}

deps <- read_deps(desc_path)

# Packages shipped with R itself are never installable from a repository.
base_pkgs <- rownames(utils::installed.packages(priority = "base"))
deps <- setdiff(deps, base_pkgs)

## ---- exclusion layer 1: GitHub-only packages ------------------------------

# Packages that live on GitHub only, i.e. are invisible to CRAN and to the
# Bioconductor release/devel channels. Add every new GitHub dependency here; the
# workflow must keep working when these cannot be installed.
github_only <- c(
  "DEGAS",
  "DEGASv2",
  "LPSGL",
  "PIPET",
  "rSIDISH",
  "rTiRank",
  "scAB",
  "SCIPAC",
  "Scissor",
  "scPAS",
  "ScPP"
)

skipped_github <- intersect(deps, github_only)
candidates <- setdiff(deps, github_only)

## ---- repositories ---------------------------------------------------------

# CRAN: honour whatever the environment already configured (setup-r-dependencies
# points this at a binary mirror), falling back to the canonical mirror.
cran_repo <- Sys.getenv("CRAN_REPO", unset = "")
if (!nzchar(cran_repo)) {
  cran_repo <- getOption("repos")[["CRAN"]]
}
if (
  is.null(cran_repo) || !nzchar(cran_repo) || identical(cran_repo, "@CRAN@")
) {
  cran_repo <- "https://cloud.r-project.org"
}

# Bioconductor: "release" for a released R, "devel" while R is unstable.
bioc_channel <- "release"
if (grepl("unstable|devel", R.version$status)) {
  bioc_channel <- "devel"
}

# Posit Package Manager mirrors Bioconductor, with binary packages on every
# platform, at /bioconductor/<snapshot>/packages/<bioc-version>/<sub-repo>; the
# release version paired with the running R therefore has to be spelled out.
bioc_versions <- c(
  "4.7" = "3.23",
  "4.6" = "3.22",
  "4.5" = "3.21",
  "4.4" = "3.20",
  "4.3" = "3.19",
  "4.2" = "3.18"
)
r_series <- paste(
  R.version$major,
  strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L],
  sep = "."
)
p3m_version <- unname(bioc_versions[r_series])

bioc_repos <- character()
if (bioc_channel == "release" && !is.na(p3m_version)) {
  p3m_root <- paste0(
    "https://packagemanager.posit.co/bioconductor/latest/packages/",
    p3m_version
  )
  bioc_repos <- c(
    BioCsoft = paste0(p3m_root, "/bioc"),
    BioCann = paste0(p3m_root, "/data/annotation"),
    BioCexp = paste0(p3m_root, "/data/experiment")
  )
}

# bioconductor.org stays behind the mirror: pak installs from the first
# repository that provides a package, so the mirror wins while the upstream copy
# fills in the packages it has not synced yet.
upstream <- file.path("https://bioconductor.org/packages", bioc_channel)
bioc_repos <- c(
  bioc_repos,
  BioCsoftUpstream = paste0(upstream, "/bioc"),
  BioCannUpstream = paste0(upstream, "/data/annotation"),
  BioCexpUpstream = paste0(upstream, "/data/experiment")
)

bioc_override <- Sys.getenv("BIOC_REPO", unset = "")
if (nzchar(bioc_override)) {
  bioc_repos <- c(BioC = bioc_override, bioc_repos)
}

## ---- exclusion layer 2: keep only packages a repository can provide --------

repos <- c(CRAN = cran_repo, bioc_repos)

extra_repos <- Sys.getenv("EXTRA_REPOS", unset = "")
if (nzchar(extra_repos)) {
  extra <- trimws(strsplit(extra_repos, ",")[[1L]])
  repos <- c(repos, extra[nzchar(extra)])
}

if (is.na(p3m_version) || bioc_channel != "release") {
  message("* Bioconductor source: bioconductor.org (", bioc_channel, ")")
} else {
  message(
    "* Bioconductor source: P3M (",
    p3m_version,
    ") + bioconductor.org"
  )
}

# Names advertised by the repository set; a mirror that answers 404 is skipped
# by utils::available.packages, leaving the remaining repositories in place.
idx <- suppressWarnings(tryCatch(
  utils::available.packages(repos = unname(repos), type = "source"),
  error = function(e) NULL
))
if (!is.null(idx) && nrow(idx) > 0L) {
  available <- unique(rownames(idx))
} else {
  warning(
    "Could not query any repository, skipping the availability filter.",
    call. = FALSE
  )
  available <- character()
}

if (length(available) > 0L) {
  not_in_repo <- setdiff(candidates, available)
  candidates <- intersect(candidates, available)
} else {
  not_in_repo <- character()
}

## ---- report ---------------------------------------------------------------

message(
  "Declared dependencies in ",
  desc_path,
  ": ",
  length(deps),
  " (",
  length(candidates),
  " from CRAN/Bioconductor)"
)
if (length(skipped_github) > 0L) {
  message("* GitHub-only, skipped: ", toString(skipped_github))
}
if (length(not_in_repo) > 0L) {
  message("* Not on CRAN/Bioconductor, skipped: ", toString(not_in_repo))
}
published <- intersect(github_only, available)
if (length(published) > 0L) {
  message(
    "* Now on CRAN/Bioconductor, drop from `github_only`: ",
    toString(published)
  )
}

if (list_only) {
  cat(candidates, sep = "\n")
  quit(save = "no", status = 0L)
}

## ---- install --------------------------------------------------------------

is_installed <- function(pkg) length(find.package(pkg, quiet = TRUE)) > 0L

missing <- candidates[!vapply(candidates, is_installed, logical(1L))]

if (length(missing) == 0L) {
  message("* All CRAN/Bioconductor dependencies are already installed.")
  quit(save = "no", status = 0L)
}

message("* To install (", length(missing), "): ", toString(missing))

if (dry_run) {
  quit(save = "no", status = 0L)
}

if (!requireNamespace("pak", quietly = TRUE)) {
  utils::install.packages("pak", repos = unname(repos))
}
# Keep the names: pak needs them to recognise a binary CRAN mirror configured by
# setup-r (`Rprofile.site`) instead of treating it as a generic source repo.
options(repos = repos)

ok <- tryCatch(
  {
    # pak resolves the whole set at once and reports every package it could not
    # install, so a single call is enough.
    pak::pkg_install(missing, upgrade = FALSE)
    TRUE
  },
  error = function(e) {
    message("* Installation failed: ", conditionMessage(e))
    FALSE
  }
)

if (!ok) {
  if (strict) {
    quit(save = "no", status = 1L)
  }
  message("* Continuing anyway (--no-strict).")
} else {
  message("* Done.")
}

quit(save = "no", status = 0L)
