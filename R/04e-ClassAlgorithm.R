Algorithm <- new_class(
  name = "Algorithm",
  parent = SigBridgeRBase,
  properties = list(
    method_name = class_character,
    executor = class_function,
    algorithm_version = class_character
  ),
  abstract = TRUE,
  constructor = \(method_name, executor, algorithm_version) {
    new_object(
      S7_object(),
      method_name = method_name,
      executor = func,
      algorithm_version = algorithm_version
    )
  }
)

ScreenMethod <- new_class(
  name = "ScreenMethod",
  parent = Algorithm
)
