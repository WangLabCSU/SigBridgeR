# Cache system (R/43-49) --------------------------------------------------

# Shared constructors ------------------------------------------------------

new_cache_config <- function(
  method_name = "Scissor",
  param = list(alpha = 0.05),
  phenotype_class = "survival",
  label_type = "survival",
  method_version = "1.0.0"
) {
  ScreenMethodConfig(
    method_name = method_name,
    param = param,
    method_version = method_version,
    phenotype_class = phenotype_class,
    label_type = label_type
  )
}

new_cache_root <- function(method_name = "Scissor") {
  file.path(withr::local_tempdir(), paste0(method_name, "_res"))
}

# 43-WriteCacheMeta.R ------------------------------------------------------

describe("WriteCacheMeta - ScreenMethodConfig", {
  it("writes a JSON file and returns the metadata invisibly", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config()

    expect_invisible(WriteCacheMeta(config, file = json, verbose = FALSE))
    expect_true(file.exists(json))
  })

  it("writes the expected top-level sections", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config()

    WriteCacheMeta(config, file = json, verbose = FALSE)
    meta <- jsonlite::read_json(json)

    expect_named(meta, c("general", "config", "description"))
  })

  it("records general system information", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config()

    meta <- WriteCacheMeta(config, file = json, verbose = FALSE)

    expect_true(all(
      c(
        "file",
        "os",
        "time",
        "r_version",
        "arch",
        "platform",
        "SigBridgeR_version"
      ) %in%
        names(meta$general)
    ))
    expect_identical(meta$general$r_version, R.version.string)
    expect_match(
      meta$general$time,
      "^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}$"
    )
  })

  it("round-trips the config properties", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config(param = list(alpha = 0.05, nfold = 5L))

    WriteCacheMeta(config, file = json, verbose = FALSE)
    cached <- jsonlite::read_json(json)$config

    expect_identical(cached$method_name, "Scissor")
    expect_identical(cached$method_version, "1.0.0")
    expect_identical(cached$phenotype_class, "survival")
    expect_identical(cached$label_type, "survival")
    expect_equal(cached$param$alpha, 0.05)
    expect_equal(cached$param$nfold, 5L)
  })

  it("stores additional_description", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config()

    WriteCacheMeta(
      config,
      file = json,
      additional_description = "Scissor screening",
      verbose = FALSE
    )

    expect_identical(
      jsonlite::read_json(json)$description,
      "Scissor screening"
    )
  })

  it("creates missing parent directories", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "a", "b", "cache_config.json")
    config <- new_cache_config()

    WriteCacheMeta(config, file = json, verbose = FALSE)

    expect_true(file.exists(json))
  })

  it("rejects a non-character file argument", {
    config <- new_cache_config()

    expect_error(
      WriteCacheMeta(config, file = 1L, verbose = FALSE),
      class = "chk_error"
    )
  })

  it("rejects an unsupported object class", {
    expect_error(
      WriteCacheMeta(list(a = 1L), file = tempfile()),
      class = "rlang_error"
    )
  })
})

describe("WriteCacheMeta - helpers", {
  it("get_time() returns a formatted timestamp", {
    expect_match(get_time(), "^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}$")
  })

  it("get_r_info() exposes the R build information", {
    info <- get_r_info()

    expect_named(info, c("r_version", "arch", "platform"))
    expect_identical(info$r_version, R.version.string)
  })

  it("get_os_info() returns a single non-empty string", {
    os <- get_os_info()

    expect_type(os, "character")
    expect_length(os, 1L)
    expect_true(nzchar(os))
  })
})

# 44-CheckCache.R ----------------------------------------------------------

describe("CheckCache - ScreenMethodConfig", {
  it("returns TRUE when the cache config matches", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config()
    WriteCacheMeta(config, file = json, verbose = FALSE)

    expect_invisible(CheckCache(config, path = dir))
    expect_true(CheckCache(config, path = dir))
  })

  it("accepts the path to cache_config.json directly", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config()
    WriteCacheMeta(config, file = json, verbose = FALSE)

    expect_true(CheckCache(config, path = json))
  })

  it("aborts when a parameter differs", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    WriteCacheMeta(new_cache_config(), file = json, verbose = FALSE)

    expect_error(
      CheckCache(new_cache_config(param = list(alpha = 0.01)), path = dir),
      "not consistent",
      class = "rlang_error"
    )
  })

  it("aborts when the method name differs", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    WriteCacheMeta(new_cache_config(), file = json, verbose = FALSE)

    expect_error(
      CheckCache(new_cache_config(method_name = "scPAS"), path = dir),
      class = "rlang_error"
    )
  })

  it("aborts on an unsupported path", {
    config <- new_cache_config()

    expect_error(
      CheckCache(config, path = tempfile()),
      class = "rlang_error"
    )
  })

  it("rejects an unsupported object class", {
    expect_error(
      CheckCache("not a config", path = tempdir()),
      class = "rlang_error"
    )
  })
})

# 45-ChooseCache.R ---------------------------------------------------------

describe("ChooseCache", {
  it("aborts when the directory name is not a cache root", {
    dir <- withr::local_tempdir()

    expect_error(ChooseCache(dir), class = "rlang_error")
  })

  it("aborts when the cache root does not exist", {
    root <- file.path(withr::local_tempdir(), "Scissor_res")

    expect_error(ChooseCache(root), class = "rlang_error")
  })

  it("aborts when the cache root is empty", {
    root <- new_cache_root()
    dir.create(root, recursive = TRUE)

    expect_error(ChooseCache(root), class = "rlang_error")
  })

  it("returns the only cache directory without prompting", {
    root <- new_cache_root()
    cache <- file.path(root, "survival_20250101")
    dir.create(cache, recursive = TRUE)

    expect_identical(ChooseCache(root), cache)
  })

  it("asks the user to pick when several caches exist", {
    root <- new_cache_root()
    first <- file.path(root, "survival_20250101")
    second <- file.path(root, "survival_20250102")
    dir.create(first, recursive = TRUE)
    dir.create(second, recursive = TRUE)
    local_mocked_bindings(menu = function(...) 2L, .package = "utils")

    expect_identical(ChooseCache(root), second)
  })

  it("aborts on an out-of-range menu choice", {
    root <- new_cache_root()
    dir.create(file.path(root, "survival_20250101"), recursive = TRUE)
    dir.create(file.path(root, "survival_20250102"), recursive = TRUE)
    local_mocked_bindings(menu = function(...) 99L, .package = "utils")

    expect_error(ChooseCache(root), class = "rlang_error")
  })
})

# 46-WriteCache.R ----------------------------------------------------------

describe("WriteCache - ScreenMethodConfig", {
  it("writes the cache metadata to the requested file", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    config <- new_cache_config()

    expect_no_error(WriteCache(config, file = json, verbose = FALSE))
    expect_true(file.exists(json))
  })

  it("rejects an unsupported object class", {
    expect_error(WriteCache(list(a = 1L)), class = "rlang_error")
  })
})

describe("detect_format", {
  it("uses csv for small 2D objects", {
    expect_identical(detect_format(data.frame(a = 1:3), 1000L, 20L), "csv")
  })

  it("does not use csv when the object exceeds max_rows/max_cols", {
    df <- data.frame(a = seq_len(50L))

    expect_false(detect_format(df, 10L, 20L) == "csv")
    expect_false(detect_format(mtcars, 1000L, 5L) == "csv")
  })

  it("falls back to rdata for non-qdata objects when qs2 is absent", {
    skip_if(is_installed("qs2"))
    expect_identical(detect_format(mean, 1000L, 20L), "rdata")
  })
})

describe("WriteSingleCacheData", {
  it("writes and reads back a csv file", {
    file <- withr::local_tempfile(fileext = ".csv")
    df <- data.frame(a = 1:3, b = letters[1:3])

    WriteSingleCacheData(df, file = file, format = "csv")

    expect_equal(LoadCache(file), data.table::as.data.table(df))
  })

  it("writes and reads back an rds file", {
    file <- withr::local_tempfile(fileext = ".rds")
    obj <- list(a = 1:3)

    WriteSingleCacheData(obj, file = file, format = "rdata")

    expect_identical(LoadCache(file), obj)
  })

  it("writes and reads back a qs2 file", {
    skip_if_not_installed("qs2")
    file <- withr::local_tempfile(fileext = ".qs2")
    obj <- list(a = 1:3)

    WriteSingleCacheData(obj, file = file, format = "qs2")

    expect_identical(LoadCache(file), obj)
  })

  it("creates missing parent directories", {
    dir <- withr::local_tempdir()
    file <- file.path(dir, "nested", "x.rds")

    WriteSingleCacheData(list(a = 1L), file = file, format = "rdata")

    expect_true(file.exists(file))
  })
})

describe("WriteCache - ScreenMethodCache", {
  it("writes each cached object with a readable file extension", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    cache <- ScreenMethodCache(
      cache_path = dir,
      cache_config_path = json,
      cache_data = list(result = data.frame(a = 1:3, b = letters[1:3])),
      screen_method_config = new_cache_config()
    )

    files <- WriteCache(cache, verbose = FALSE)

    expect_length(files, 1L)
    expect_true(file.exists(file.path(dir, "result.csv")))
    expect_true(file.exists(json))
    expect_equal(LoadCache(file.path(dir, "result.csv"))$a, 1:3)
  })

  it("honours an explicit format", {
    dir <- withr::local_tempdir()
    json <- file.path(dir, "cache_config.json")
    cache <- ScreenMethodCache(
      cache_path = dir,
      cache_config_path = json,
      cache_data = list(result = data.frame(a = 1:3)),
      screen_method_config = new_cache_config()
    )

    files <- WriteCache(cache, format = "rdata", verbose = FALSE)

    expect_match(files, "\\.rds$")
    expect_true(all(file.exists(files)))
  })
})

# 47-LoadCache.R -----------------------------------------------------------

describe("LoadCache", {
  it("loads a csv file", {
    file <- withr::local_tempfile(fileext = ".csv")
    data.table::fwrite(data.frame(x = 1:3), file)

    expect_equal(LoadCache(file)$x, 1:3)
  })

  it("loads an rds file", {
    file <- withr::local_tempfile(fileext = ".rds")
    saveRDS(list(x = 1:3), file)

    expect_identical(LoadCache(file), list(x = 1:3))
  })

  it("loads an RData file", {
    file <- withr::local_tempfile(fileext = ".RData")
    saveRDS(list(x = 1:3), file)

    expect_identical(LoadCache(file), list(x = 1:3))
  })

  it("loads a qs2 file", {
    skip_if_not_installed("qs2")
    file <- withr::local_tempfile(fileext = ".qs2")
    qs2::qs_save(list(y = 1:2), file)

    expect_identical(LoadCache(file), list(y = 1:2))
  })

  it("loads a qdata file", {
    skip_if_not_installed("qs2")
    file <- withr::local_tempfile(fileext = ".qdata")
    qs2::qs_save(list(z = 1:2), file)

    expect_identical(LoadCache(file), list(z = 1:2))
  })

  it("aborts on an unsupported extension", {
    file <- withr::local_tempfile(fileext = ".txt")
    writeLines("x", file)

    expect_error(
      LoadCache(file),
      "Unsupported file extension",
      class = "rlang_error"
    )
  })

  it("aborts when the file does not exist", {
    expect_error(LoadCache(tempfile(fileext = ".csv")), class = "chk_error")
  })

  it("aborts on a non-string file argument", {
    expect_error(LoadCache(1L), class = "chk_error")
  })
})

# 48-CacheSetHere.R --------------------------------------------------------

describe("CacheSetHere - save mode", {
  it("creates a cache under a parent layer path", {
    parent <- withr::local_tempdir()
    config <- new_cache_config()

    cache_dir <- CacheSetHere(
      parent,
      config,
      mode = "save",
      timestamp = "20250101"
    )

    expect_identical(
      cache_dir,
      normalizePath(
        file.path(parent, "Scissor_res", "survival_20250101"),
        mustWork = FALSE
      )
    )
    expect_true(dir.exists(cache_dir))
  })

  it("creates a cache under a root layer path", {
    root <- new_cache_root()
    dir.create(root, recursive = TRUE)
    config <- new_cache_config()

    cache_dir <- CacheSetHere(
      root,
      config,
      mode = "save",
      timestamp = "20250101"
    )

    expect_identical(basename(cache_dir), "survival_20250101")
    expect_identical(dirname(cache_dir), normalizePath(root, mustWork = FALSE))
  })

  it("aborts when the target path is already a cache", {
    dir <- withr::local_tempdir()
    WriteCacheMeta(
      new_cache_config(),
      file = file.path(dir, "cache_config.json"),
      verbose = FALSE
    )

    expect_error(
      CacheSetHere(dir, new_cache_config(), mode = "save"),
      "Recursive caching",
      class = "rlang_error"
    )
  })

  it("aborts when the cache directory already exists", {
    parent <- withr::local_tempdir()
    config <- new_cache_config()
    CacheSetHere(parent, config, mode = "save", timestamp = "20250101")

    expect_error(
      CacheSetHere(parent, config, mode = "save", timestamp = "20250101"),
      "already exists",
      class = "rlang_error"
    )
  })
})

describe("CacheSetHere - load mode", {
  it("returns the cache path unchanged for a cache layer", {
    dir <- withr::local_tempdir()
    WriteCacheMeta(
      new_cache_config(),
      file = file.path(dir, "cache_config.json"),
      verbose = FALSE
    )

    expect_identical(
      CacheSetHere(dir, new_cache_config(), mode = "load"),
      normalizePath(dir, mustWork = FALSE)
    )
  })

  it("resolves a root layer path holding a single cache", {
    root <- new_cache_root()
    cache <- file.path(root, "survival_20250101")
    dir.create(cache, recursive = TRUE)

    expect_identical(
      CacheSetHere(root, new_cache_config(), mode = "load"),
      cache
    )
  })

  it("resolves a parent layer path holding a cache root", {
    parent <- withr::local_tempdir()
    cache <- file.path(parent, "Scissor_res", "survival_20250101")
    dir.create(cache, recursive = TRUE)

    expect_identical(
      CacheSetHere(parent, new_cache_config(), mode = "load"),
      cache
    )
  })

  it("aborts when no cache root exists under the parent", {
    parent <- withr::local_tempdir()

    expect_error(
      CacheSetHere(parent, new_cache_config(), mode = "load"),
      class = "rlang_error"
    )
  })
})

describe("CacheSetHere - argument validation", {
  it("aborts on an invalid mode", {
    expect_error(
      CacheSetHere(withr::local_tempdir(), new_cache_config(), mode = "auto"),
      class = "rlang_error"
    )
  })

  it("aborts on a non-string path", {
    expect_error(CacheSetHere(1L, new_cache_config()), class = "chk_error")
  })

  it("aborts on unused dots", {
    expect_error(
      CacheSetHere(withr::local_tempdir(), new_cache_config(), foo = 1),
      class = "rlang_error"
    )
  })
})

describe("generate_cache_name", {
  it("joins the phenotype class and the timestamp", {
    expect_identical(
      generate_cache_name("survival", "20250101"),
      "survival_20250101"
    )
  })

  it("falls back to the current time", {
    expect_match(
      generate_cache_name("binary"),
      "^binary_\\d{14}$"
    )
  })
})

# 49-CacheSysCall.R -------------------------------------------------------

describe("CacheSysCall - argument validation", {
  it("aborts when a ScreenMethodConfig is used in save mode", {
    expect_error(
      CacheSysCall("save", cache = new_cache_config()),
      class = "rlang_error"
    )
  })

  it("aborts when cache is neither config nor cache", {
    expect_error(
      CacheSysCall("save", cache = list(a = 1L)),
      class = "rlang_error"
    )
  })

  it("aborts on an invalid mode", {
    expect_error(
      CacheSysCall("bogus", cache = new_cache_config()),
      class = "rlang_error"
    )
  })

  it("aborts when path is missing in load mode", {
    expect_error(
      CacheSysCall("load", cache = new_cache_config()),
      "must be supplied",
      class = "rlang_error"
    )
  })

  it("aborts when path is missing in load mode with a cache object", {
    dir <- withr::local_tempdir()
    cache <- ScreenMethodCache(
      cache_path = dir,
      cache_config_path = file.path(dir, "cache_config.json"),
      cache_data = list(),
      screen_method_config = new_cache_config()
    )

    expect_error(
      CacheSysCall("load", cache = cache),
      "must be supplied",
      class = "rlang_error"
    )
  })
})

describe("CacheSysCall - round trip", {
  it("falls back to cache@cache_path in save mode", {
    parent <- withr::local_tempdir()
    cache <- ScreenMethodCache(
      cache_path = parent,
      cache_config_path = file.path(parent, "cache_config.json"),
      cache_data = list(result = data.frame(a = 1:3)),
      screen_method_config = new_cache_config()
    )

    CacheSysCall("save", cache = cache, verbose = FALSE, timestamp = "20250101")

    expect_true(
      file.exists(file.path(
        parent,
        "Scissor_res",
        "survival_20250101",
        "result.csv"
      ))
    )
  })

  it("writes into the cache directory resolved by CacheSetHere()", {
    parent <- withr::local_tempdir()
    cache <- ScreenMethodCache(
      cache_path = parent,
      cache_config_path = file.path(parent, "cache_config.json"),
      cache_data = list(result = data.frame(a = 1:3, b = letters[1:3])),
      screen_method_config = new_cache_config()
    )

    expect_no_error(CacheSysCall(
      "save",
      path = parent,
      cache = cache,
      verbose = FALSE,
      timestamp = "20250101"
    ))

    cache_dir <- file.path(parent, "Scissor_res", "survival_20250101")
    expect_true(file.exists(file.path(cache_dir, "cache_config.json")))
    expect_true(file.exists(file.path(cache_dir, "result.csv")))
  })

  it("loads a previously saved cache as a named list", {
    parent <- withr::local_tempdir()
    config <- new_cache_config()
    cache <- ScreenMethodCache(
      cache_path = parent,
      cache_config_path = file.path(parent, "cache_config.json"),
      cache_data = list(result = data.frame(a = 1:3, b = letters[1:3])),
      screen_method_config = config
    )
    CacheSysCall(
      "save",
      path = parent,
      cache = cache,
      verbose = FALSE,
      timestamp = "20250101"
    )

    loaded <- CacheSysCall(
      "load",
      path = file.path(parent, "Scissor_res", "survival_20250101"),
      cache = config,
      verbose = FALSE
    )

    expect_type(loaded, "list")
    expect_named(loaded, "result")
    expect_equal(loaded$result$a, 1:3)
  })

  it("aborts when the cache directory holds no data", {
    dir <- withr::local_tempdir()
    WriteCacheMeta(
      new_cache_config(),
      file = file.path(dir, "cache_config.json"),
      verbose = FALSE
    )

    expect_error(
      CacheSysCall(
        "load",
        path = dir,
        cache = new_cache_config(),
        verbose = FALSE
      ),
      "No cached data found",
      class = "rlang_error"
    )
  })

  it("aborts when loading an inconsistent cache", {
    dir <- withr::local_tempdir()
    WriteCacheMeta(
      new_cache_config(),
      file = file.path(dir, "cache_config.json"),
      verbose = FALSE
    )

    expect_error(
      CacheSysCall(
        "load",
        path = dir,
        cache = new_cache_config(param = list(alpha = 0.01)),
        verbose = FALSE
      ),
      class = "rlang_error"
    )
  })
})
