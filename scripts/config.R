# Edit paths here. Relative paths resolve from reproducibility/scripts.
config <- list(
  data_root = "../data",
  output_root = "../outputs/original",
  secondary_output_root = "../outputs/secondary",
  r_library = "../.library",
  expected_R = "4.5.1",
  packages = c(survey = "4.5", jsonlite = "2.0.0", digest = "0.6.37", foreign = "0.8.90"),
  stage_elapsed_limit_seconds = 1800
)
