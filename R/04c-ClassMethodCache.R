MethodCache <- new_class(
  name = "MethodCache",
  parent = SigBridgeRBase,
  properties = list(
    cache_path = class_character,
    cache_config_path = class_character,
    cache_config = class_list,
    data = class_list
  )
)
