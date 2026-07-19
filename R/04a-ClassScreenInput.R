# ! Preprocessed Result of `PreProcess()``

utils::globalVariables(
  names = c("ScreenInput", "SCInput", "BulkInput", "PhenoInput")
)

ScreenInput <- new_class(
  name = "ScreenInput",
  parent = SigBridgeRBase,
  properties = list(
    version = class_sigbridger_verison,
    preprocessed = class_logical
  ),
  abstract = TRUE
)

# --------------------------------------------------------------------------------

SCInput <- new_class(
  name = "SCInput",
  parent = ScreenInput,
  properties = list(data = class_seurat),
  abstract = TRUE
)

# ------------------------------------------------------------------------------

BulkInput <- new_class(
  name = "BulkInput",
  parent = ScreenInput,
  properties = list(data = class_matrix),
  constructor = \(x) {
    new_object(
      .parent = S7_object(),
      data = x
    )
  },
  validator = function(self) {
    if (is.null(rownames(self@data))) {
      "@data must have rownames (gene names)"
    }
    if (is.null(rownames(self@data))) {
      "@data must have colnames (sample names)"
    }
  }
)

# --------------------------------------------------------------------------

PhenoInput <- new_generic(name = "PhenoInput", dispatch_args = "data")

method(PhenoInput, class_data.frame) <- function(data) {
  new_object(PhenotypeInput.survival, data = data)
}
method(PhenoInput, class_vector) <- function(data) {
  new_object(PhenotypeInput.binary, data = data)
}

PhenotypeInput <- new_class(
  name = "PhenotypeInput",
  parent = ScreenInput,
  properties = list(data = class_any),
  abstract = TRUE
)


PhenotypeInput.survival <- new_class(
  name = "PhenotypeInput.survival",
  parent = PhenotypeInput,
  properties = list(data = class_data.frame),
  constructor = \(x) {
    new_object(
      .parent = S7_object(),
      data = x
    )
  },
  validator = \(self) {
    if (is.null(rownames(self@data))) {
      "@data must have rownames (sample names)"
    }
    if (is.null(colnames(self@data))) {
      "@data must have colnames (phenotype names)"
    }
    if (anyNA(self@data) || any(is.null(self@data))) {
      cli::cli_fmt(cli::cli_inform("@data contains {.val NA} or {.val NULL}"))
    }
  }
)

PhenotypeInput.binary <- new_class(
  name = "PhenotypeInput.binary",
  parent = PhenotypeInput,
  properties = list(data = class_vector),
  abstract = FALSE,
  constructor = \(x) {
    new_object(
      .parent = S7_object(),
      data = x
    )
  },
  validator = \(self) {
    if (is.null(names(self@data))) {
      "@data must have names (sample names)"
    }
    if (length(table(self@data)) != 2) {
      "data must have two unique values"
    }
    if (anyNA(self@data)) {
      cli::cli_fmt(cli::cli_inform("@data contains {.val NA}"))
    }
  }
)
