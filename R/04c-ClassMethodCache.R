MethodCache <- new_class(
  name = "MethodCache",
  parent = SigBridgeRBase,
  properties = list(
    cache_path = class_character,
    cache_config_path = class_character,
    cache_config = class_list,
    data = ScreenMethodConfig
  ),
  validator = \(self) {
    if (!is_named(self@data)) {
      return("@data must have names (sample names)")
    }
    if (!is_named(self@cache_config)) {
      return("@cache_config must have names (sample names)")
    }
  }
)
