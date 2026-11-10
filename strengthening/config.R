# Fixed, independently reviewed exploratory specification. No model selection.
strengthening_config <- list(
  expected_R = "4.5.1", packages = c(survey = "4.5", foreign = "0.8.90", digest = "0.6.37", jsonlite = "2.0.0"),
  output_relative = "../outputs/strengthening/results", minimum_probability = 0.05,
  maximum_factor = 20, expected_n = c(parent = 39156L, eligible = 2424L, complete = 2320L, context = 2175L))
