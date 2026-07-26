ScreenMethodCache <- new_class(
  name = "ScreenMethodCache",
  parent = SigBridgeRBase,
  properties = list(
    cache_path = property_chr,
    cache_config_path = property_chr,
    cache_config = property_list,
    screen_method_config = new_property(class = ScreenMethodConfig)
  )
)
