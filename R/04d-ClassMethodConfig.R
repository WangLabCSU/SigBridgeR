# 同时用作screen result和cache config
ScreenMethodConfig <- new_class(
  name = "ScreenMethodConfig",
  parent = SigBridgeRBase,
  properties = list(
    method_name = property_chr,
    method_version = property_chr,
    config = property_list
  ),
  abstract = TRUE,
  validator = \(self) {
    if (length(self@config) == 0L) {
      return("@config is empty")
    }
  }
)
