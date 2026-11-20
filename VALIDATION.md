# Reproduction validation

## Version 1.1.0: adjusted RQ1 energy partition (20 November 2026)

- Adds a fourth stage, `rq1_adjusted/`: five survey-weighted models of the paired energy components (ethanol, beverage-energy remainder, linked additions, other foods and beverages, total) on the recall-order sign and six weekday-indicator differences in the primary 3–10-day domain (n = 2,320). The intercept is the drinking-day difference at equal recall-order weighting with both recalls on the same weekday. Specified after all earlier results were known and before these models were fitted; post hoc sensitivity analysis.
- Native Windows, R 4.5.1, pinned dependency versions. The stage's three synthetic tests passed.
- The stage rebuilt the paired data from the 28 verified raw XPORT files. Before fitting, unadjusted means reproduced Table 2 within 4e-13 kcal (gate tolerance 1e-6).
- The adjusted components summed to the adjusted total (gap about 1e-12 kcal). Results agreed within 5e-13 with the first computation and matched the two new fixed references at tolerance 1e-9.
- The earlier three stages and their code are unchanged. `Rscript RUN_ALL.R --verify` matched all 50 fixed reference files (48 earlier and 2 new), and the reference-verifier fixtures passed.
- The new stage was run in the existing validated checkout, whose earlier-stage outputs came from the v1.0.0 fresh run below. A complete four-stage rerun in a new empty checkout was not repeated for v1.1.0.

## Version 1.0.0

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
