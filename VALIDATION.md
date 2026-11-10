# Reproduction validation

Validation date: 9 November 2026.

## Completed checks

- Native Windows with R 4.5.1 and the exact versions in `scripts/dependency_versions.csv`.
- All 32 raw XPORT input hashes matched the fixed public-data manifest. A separate downloader smoke test retrieved one missing SMQ file from its official CDC URL, verified its hash, and verified the other 31 pre-existing files. This was not a new download of all 32 files.
- The first complete run started with no generated outputs or copied model checkpoints. It reconstructed all three original models, all three earlier exploratory models, and every fixed additional sensitivity analysis from verified raw inputs.
- All 48 expected reference files matched at numeric tolerance `1e-9`: 26 original-model/descriptive files including three full coefficient covariance matrices, seven secondary-analysis files, and 15 additional-sensitivity files including response diagnostics.
- All 121 response-model/outcome-model jackknife replicates completed. Both original-model coefficient/covariance reconstruction checks passed. The seven synthetic scientific tests passed.
- A specification-text encoding correction changed no scientific code. The affected sensitivity stage was rebuilt in an empty output directory and again matched all fixed references. The original and secondary stages retained their valid completed results.
- Reference-verifier fixtures confirmed rejection of a missing expected output, a changed numerical value, and a renamed column. A valid fixture passed. CSV dimensions, column names/order and covariance term labels are checked before numeric comparison.
- The final `--verify` command checked source/output identities and all independent aggregate references. No expected reference value was generated from the validation run.
- The eight original scientific/helper modules and the additional sensitivity analysis/test modules retain the manuscript implementation. Distribution changes concern orchestration, installation, downloading, relative paths and validation against safe fixed references.

## Dependency installation scope

The default installer verified the already installed pinned versions. An isolated installation excluded the pre-existing user library and successfully downloaded the pinned Rcpp 1.1.0 source. It stopped at compilation because the test machine lacked `make` and Rtools45. It did not substitute another version. A complete isolated dependency installation was therefore not demonstrated on this machine. Install Rtools45 before building archived compiled packages on Windows; the README lists this prerequisite. Exact current binaries are used when CRAN provides the pinned version.

## Limits

The validated platform is Windows. macOS and Linux execution has not been validated. PNG byte equality depends on graphics devices and installed fonts; the plot's point estimates and interval endpoints are checked numerically. This release reproduces analysis values and the figure, not manually edited Word prose or pagination. Successful reproduction does not establish causal validity, remove selection/reporting error, or verify the untestable response assumptions.

Raw inputs, generated outputs, participant-level diagnostics and RDS checkpoints remain local and are excluded from version control. The repository contains the code, setup instructions, public input manifests and safe aggregate comparison values needed to repeat these checks.
