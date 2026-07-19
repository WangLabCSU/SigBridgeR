utils::globalVariables(name = c("class_seurat"))

SigBridgeRBase <- new_class(
  name = "SigBridgeRBase",
  abstract = TRUE
)

# --------------------------------------------------------------------------
# ! Inputs in Screen()

class_seurat <- new_property(
  class = class_any,
  getter = NULL,
  setter = NULL,
  validator = \(value) {
    cls_value <- class(value)
    if (cls_value != "Seurat") {
      cli::cli_fmt(cli::cli_inform(
        "epxected a {.cls Seurat} object, got {.cls cls_value}"
      ))
    }
  },
  default = quote(cli::cli_abort(c(
    "x" = "a {.cls Seurat} object is required in {.field @data}"
  )))
)

class_InMemoryAnnData <- new_property(class = class_any, validator = \(seld) {
  cls_value <- class(value)
  if (!"InMemoryAnnData" %in% cls_value) {
    cli::cli_fmt(cli::cli_inform(
      "epxected a {.cls InMemoryAnnData} object, got {.cls cls_value}"
    ))
  }
})

class_AnnDataR6 <- new_property(class = class_any, validator = \(seld) {
  cls_value <- class(value)
  if (!"AnnDataR6" %in% cls_value) {
    cli::cli_fmt(cli::cli_inform(
      "epxected a {.cls AnnDataR6} object, got {.cls cls_value}"
    ))
  }
})

class_matrix <- new_property(
  class = class_any,
  validator = \(value) {
    cls_value <- class(value)
    if (!"matrix" %in% cls_value) {
      cli::cli_fmt(cli::cli_inform(
        "epxected a {.cls matrix} object, got {.cls {cls_value}}"
      ))
    }
  }
)

class_datalike2d <- new_property(
  class = class_any,
  validator = \(value) {
    if (!is_2d(value)) {
      cli::cli_fmt(cli::cli_inform(
        "epxected a 2d data like object, got {.cls {cls_value}}"
      ))
    }
  }
)

# --------------------------------------------------------------------------
# ! Base Info

class_sigbridger_verison <- new_property(
  class = class_character,
  getter = function(self) {
    utils::packageVersion("SigBridgeR")
  },
  setter = function(self, value) {
    if (!length(value)) {
      return(self)
    }

    cli::cli_abort(
      c(
        "x" = "package version is {.field read-only}."
      ),
      call = rlang::caller_call(n = 2),
      .frame = rlang::caller_env(n = 2)
    )
  },
  default = quote(utils::packageVersion("SigBridgeR")),
)
