# ---- qs2-backed workspace image helpers -------------------------------------
#
# Drop-in replacements for base::save.image() / base::load() that use the
# `qs2` package instead of base serialization. Intended for interactive use:
#
#   source("dev/save_image2.R")
#   save.image2()            # writes .image.qs2 with every object in .GlobalEnv
#   load.image2()            # restores them into .GlobalEnv
#
# Threads default to half of the detected logical CPUs (at least 1).

#' @title Default qs2 Thread Count
#'
#' @description
#' Half of the logical CPU count reported by [parallel::detectCores()], with a
#' lower bound of one. Falls back to `1L` when the core count cannot be
#' detected (`NA`).
#'
#' @return `integer(1)`.
#' @noRd
image2_default_nthreads <- function() {
  cores <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(e) NA_integer_,
    warning = function(w) NA_integer_
  )

  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    return(1L)
  }

  max(1L, as.integer(cores %/% 2L))
}

#' @title Validate a qs2 Thread Count
#'
#' @param nthreads A single positive number.
#'
#' @return `integer(1)`.
#' @noRd
image2_check_nthreads <- function(nthreads) {
  if (
    !is.numeric(nthreads) ||
      length(nthreads) != 1L ||
      is.na(nthreads) ||
      nthreads < 1
  ) {
    stop("`nthreads` must be a single positive number.", call. = FALSE)
  }

  as.integer(nthreads)
}

#' @title Validate a Workspace Image Path
#'
#' @param file A single non-empty character path.
#' @param must_exist Logical. Whether the file has to exist already.
#'
#' @return The path with `~` expanded, as `character(1)`.
#' @noRd
image2_check_file <- function(file, must_exist = FALSE) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a single non-empty path.", call. = FALSE)
  }

  file <- path.expand(file)

  if (isTRUE(must_exist) && !file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  file
}

#' @title Write a qs2 Workspace Image
#'
#' @description
#' Writes `objects` to `file` with [qs2::qs_save()]. When `safe = TRUE` the
#' data is first written to a temporary file in the same directory and then
#' moved into place, so an interrupted save cannot leave a truncated image.
#'
#' @param objects A named list of objects to save.
#' @param file `character(1)`. Destination path.
#' @param nthreads `integer(1)`. Threads passed to [qs2::qs_save()].
#' @param compress_level `integer(1)`. Compression level passed to
#'   [qs2::qs_save()].
#' @param shuffle `logical(1)`. Shuffle bytes before compression.
#' @param safe `logical(1)`. Write atomically via a temporary file.
#'
#' @return Invisible `file`.
#' @noRd
image2_qs_save <- function(
  objects,
  file,
  nthreads,
  compress_level = 3L,
  shuffle = TRUE,
  safe = TRUE
) {
  dir <- dirname(file)
  if (!dir.exists(dir)) {
    stop(
      "Directory does not exist: ",
      dir,
      "\nCreate it first or pass a path in an existing directory.",
      call. = FALSE
    )
  }

  if (!isTRUE(safe)) {
    qs2::qs_save(
      objects,
      file,
      compress_level = compress_level,
      shuffle = shuffle,
      nthreads = nthreads
    )
    return(invisible(file))
  }

  tmp <- tempfile(pattern = paste0(".", basename(file), "-"), tmpdir = dir)
  on.exit(unlink(tmp), add = TRUE)

  qs2::qs_save(
    objects,
    tmp,
    compress_level = compress_level,
    shuffle = shuffle,
    nthreads = nthreads
  )

  if (!file.rename(tmp, file)) {
    # file.rename() cannot overwrite on some platforms (e.g. Windows)
    if (!file.copy(tmp, file, overwrite = TRUE)) {
      stop("Failed to move temporary image into place: ", file, call. = FALSE)
    }
  }

  invisible(file)
}

#' @title Save the Current Environment as a qs2 Image
#'
#' @description
#' qs2-backed counterpart of [base::save.image()]. Collects every object in
#' `envir` (by default the calling environment, i.e. `.GlobalEnv` when called
#' from the top level) and writes them to a single file with [qs2::qs_save()].
#' The file is read back with [load.image2()].
#'
#' @param file `character(1)`. Output path. Default: `".image.qs2"`.
#' @param envir Environment whose objects are saved. Default: the calling
#'   environment.
#' @param nthreads `integer(1)`. Threads used by `qs2`. Default: half of the
#'   detected logical CPUs (at least one).
#' @param compress_level `integer(1)`. Compression level, passed to
#'   [qs2::qs_save()]. Default: `3L`.
#' @param shuffle `logical(1)`. Shuffle bytes before compression, passed to
#'   [qs2::qs_save()]. Default: `TRUE`.
#' @param all.names `logical(1)`. Include objects whose names start with a dot,
#'   as in [base::ls()]. Default: `TRUE`.
#' @param exclude `character`. Object names to skip. Default: `NULL`.
#' @param verbose `logical(1)`. Report the number of saved objects.
#'   Default: `TRUE`.
#' @param safe `logical(1)`. Write atomically via a temporary file.
#'   Default: `TRUE`.
#'
#' @return Invisibly, the path of the written file.
#'
#' @seealso [load.image2()]
#' @export
#'
#' @examples
#' \dontrun{
#' source("dev/save_image2.R")
#'
#' x <- 1:10
#' y <- data.frame(a = 1:3)
#' save.image2()                 # -> .image.qs2
#' save.image2("analysis.qs2")   # custom path
#' }
save.image2 <- function(
  file = ".image.qs2",
  envir = parent.frame(),
  nthreads = image2_default_nthreads(),
  compress_level = 3L,
  shuffle = TRUE,
  all.names = TRUE,
  exclude = NULL,
  verbose = TRUE,
  safe = TRUE
) {
  file <- image2_check_file(file)
  nthreads <- image2_check_nthreads(nthreads)

  if (!is.environment(envir)) {
    stop("`envir` must be an environment.", call. = FALSE)
  }

  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("Package `qs2` is required to save images.", call. = FALSE)
  }

  obj_names <- ls(envir = envir, all.names = all.names)

  if (length(exclude)) {
    obj_names <- setdiff(obj_names, exclude)
  }

  objects <- if (length(obj_names)) {
    mget(obj_names, envir = envir, inherits = FALSE)
  } else {
    list()
  }

  image2_qs_save(
    objects,
    file,
    nthreads = nthreads,
    compress_level = compress_level,
    shuffle = shuffle,
    safe = safe
  )

  if (isTRUE(verbose)) {
    message(sprintf(
      "Saved %d object%s to %s (qs2, %d thread%s)",
      length(obj_names),
      if (length(obj_names) == 1L) "" else "s",
      file,
      nthreads,
      if (nthreads == 1L) "" else "s"
    ))
  }

  invisible(file)
}

#' @title Load a qs2 Workspace Image
#'
#' @description
#' Companion of [save.image2()]. Reads a file written by [save.image2()] with
#' [qs2::qs_read()] and assigns the stored objects into `envir`, overwriting
#' existing bindings with the same names (like [base::load()]).
#'
#' @param file `character(1)`. Image written by [save.image2()].
#'   Default: `".image.qs2"`.
#' @param envir Environment to restore into. Default: the calling environment.
#' @param nthreads `integer(1)`. Threads used by `qs2`. Default: half of the
#'   detected logical CPUs (at least one).
#' @param validate_checksum `logical(1)`. Verify the stored checksum, passed to
#'   [qs2::qs_read()]. Default: `TRUE`.
#' @param verbose `logical(1)`. Report the number of loaded objects.
#'   Default: `TRUE`.
#'
#' @return Invisibly, a character vector with the names of the restored
#'   objects.
#'
#' @seealso [save.image2()]
#' @export
#'
#' @examples
#' \dontrun{
#' source("dev/save_image2.R")
#' load.image2()                 # restores .image.qs2 into .GlobalEnv
#' load.image2("analysis.qs2")   # custom path
#' }
load.image2 <- function(
  file = ".image.qs2",
  envir = parent.frame(),
  nthreads = image2_default_nthreads(),
  validate_checksum = TRUE,
  verbose = TRUE
) {
  file <- image2_check_file(file, must_exist = TRUE)
  nthreads <- image2_check_nthreads(nthreads)

  if (!is.environment(envir)) {
    stop("`envir` must be an environment.", call. = FALSE)
  }

  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("Package `qs2` is required to load images.", call. = FALSE)
  }

  objects <- qs2::qs_read(
    file,
    validate_checksum = validate_checksum,
    nthreads = nthreads
  )

  if (!is.list(objects) || (length(objects) && is.null(names(objects)))) {
    stop(
      "`",
      file,
      "` is not a workspace image written by `save.image2()`; ",
      "read it with `qs2::qs_read()` instead.",
      call. = FALSE
    )
  }

  obj_names <- names(objects)
  if (is.null(obj_names)) {
    obj_names <- character()
  }

  if (length(objects)) {
    list2env(objects, envir = envir)
  }

  if (isTRUE(verbose)) {
    message(sprintf(
      "Loaded %d object%s from %s (qs2, %d thread%s)",
      length(obj_names),
      if (length(obj_names) == 1L) "" else "s",
      file,
      nthreads,
      if (nthreads == 1L) "" else "s"
    ))
  }

  invisible(obj_names)
}
