ScreenMethodResult <- new_class(
  name = "ScreenMethodResult",
  parent = SigBridgeRBase,
  properties = list(
    scRNA_data = property_seurat
  ),
  validator = \(self) {
    if (!prop_exists(self, "scRNA_data")) {
      return(cli::cli_fmt(cli::cli_text(
        "{.val scRNA_data} is required"
      )))
    }

    if (!inherits(self@scRNA_data, "Seurat")) {
      return(cli::cli_fmt(cli::cli_text(
        "{.val scRNA_data} must be a {.cls Seurat} object"
      )))
    }
  },
  constructor = \(...) {
    dots <- list(...)

    obj <- new_object(
      S7_object(),
      scRNA_data = scRNA_data
    )
    props(obj) <- dots
    obj
  }
)
