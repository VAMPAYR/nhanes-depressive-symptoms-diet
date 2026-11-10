# Depressive symptoms and dietary energy on drinking and comparison days

This repository reproduces the NHANES 2011–2018 analyses for **Depressive symptoms and dietary energy on days with and without alcoholic beverages in NHANES 2011 to 2018**. It rebuilds the data from official public XPORT files, fits all manuscript models, produces the main and supplementary numerical results and Figure 1, and compares the results with fixed aggregate reference values.

The primary analysis includes 2,320 adults who reported an alcoholic beverage on exactly one of two reliable dietary recalls 3–10 days apart. The outcome removes complete alcoholic-beverage items and their linked additions from total energy. The PHQ-9 association is observational and imprecise; the analysis does not establish causation or equivalence.

## Requirements

- **R 4.5.1**. Install this exact version from [CRAN](https://cran.r-project.org/). The package checks the version rather than silently accepting another runtime.
- The exact package versions in [`scripts/dependency_versions.csv`](scripts/dependency_versions.csv), including indirect dependencies. The installer checks all listed versions.
- **Windows was the validated platform.** The runner locates Rscript from the active R installation. macOS and Linux paths are supported by the entry point, but results on those platforms have not been validated.
- For source installation of compiled packages on Windows, install the [Rtools45 toolchain](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html). Other operating systems need their corresponding compiler toolchain. RStudio is optional.
- Internet access for package installation and 32 official NHANES files. The inputs require approximately 690 MiB; allow several GB for the R library and local outputs. Figure rendering uses Cairo and Times New Roman; font/device differences can change PNG bytes without changing the plotted values.
- Run commands in a terminal that can find `Rscript`, or use the Rscript executable from your own R installation. No author-specific drive, user folder, or data location is required.

## Quick start

Clone or download this repository. Open a terminal in its root directory, the directory containing `RUN_ALL.R`.

```sh
git clone https://github.com/VAMPAYR/nhanes-depressive-symptoms-diet.git
cd nhanes-depressive-symptoms-diet
Rscript scripts/install_dependencies.R
Rscript scripts/download_data.R --download
Rscript RUN_ALL.R --preflight
Rscript RUN_ALL.R --run
Rscript RUN_ALL.R --verify
```

The default paths are relative to the checkout: `data/`, `.library/`, and `outputs/`. These directories are ignored by Git. Downloaded public inputs and generated participant-level diagnostic files stay local.

For RStudio, open an RStudio project in this repository or set its working directory to the repository root. Complete installation and downloading using the Terminal commands above, then run:

```r
source("RUN_ALL.R")
```

The sourced entry point uses the same three-stage `--run` workflow. No separate manual ordering or path editing is required.

## What each command does

| Command | Operation |
| --- | --- |
| `Rscript scripts/install_dependencies.R` | Checks pinned versions and installs missing versions into the checkout's `.library/`; no participant data are read. |
| `Rscript scripts/install_dependencies.R --isolated` | Excludes the existing user library and installs/checks the pinned dependencies using this checkout's library and R's bundled packages. Archived compiled versions can require Rtools45. |
| `Rscript scripts/download_data.R --download` | Downloads missing official files and checks every SHA-256 checksum. Existing correct files are reused; an existing changed file causes an error. |
| `Rscript scripts/download_data.R --verify` | Checks the 32 local inputs without downloading. HTML/error responses cannot pass the expected hashes. |
| `Rscript RUN_ALL.R --preflight` | Checks R/packages, code/reference identities and input hashes; runs seven synthetic sensitivity tests. It fits no observed participant-outcome model. |
| `Rscript RUN_ALL.R --run` | Rebuilds all three analysis stages from raw files when outputs are empty. If a stage is already complete, verifies its recorded identity and outputs rather than silently refitting. |
| `Rscript RUN_ALL.R --verify` | Verifies completed local output hashes and compares 48 fixed aggregate reference files, including full coefficient covariance matrices for the original three models. It is not a fresh analysis. |
| `Rscript scripts/test_validation.R` | Tests that the reference verifier rejects missing files, changed numerical values and renamed columns, using synthetic temporary fixtures only. |

To conduct another **fresh** reproduction, use a new checkout with empty `outputs/`. Reuse the verified `data/` files by copying them into the same relative layout if desired. Do not delete another person's active run or edit expected values to make a failed comparison pass. The original and secondary stages can resume valid interrupted checkpoints. An incomplete sensitivity-stage output directory is preserved for diagnosis; start that stage in a new empty checkout after resolving its cause.

## Data

[`data_manifest.csv`](data_manifest.csv) lists all 32 input paths, official CDC download URLs, codebook URLs, byte counts and SHA-256 hashes. Files are organized by cycle (`2011-2012`, `2013-2014`, `2015-2016`, `2017-2018`). Each cycle needs:

| File family | Purpose |
| --- | --- |
| `DEMO` | Age, sex, survey design, income, education and race/ethnicity. |
| `DPQ` | Nine depression-screening items. |
| `DR1TOT`, `DR2TOT` | Daily totals, recall reliability, dietary weights and timing fields. |
| `DR1IFF`, `DR2IFF` | Individual-food energy, ethanol, food codes and combination fields. |
| `DRXFCD` | Cycle-specific food-code descriptions. |
| `SMQ` | Smoking categories for the fixed additional-adjustment analysis. |

The repository does not redistribute raw data, participant records, local model RDS objects, or copyrighted literature. NHANES is publicly accessible through [NCHS](https://wwwn.cdc.gov/nchs/nhanes/). The manifests retain the exact releases used for these estimates; a changed upstream file causes a checksum failure requiring investigation.

## Analysis stages and scripts

| File | Responsibility |
| --- | --- |
| `RUN_ALL.R` | Runs the original, secondary and additional sensitivity stages in order, then validates all fixed references. |
| `scripts/install_dependencies.R` | Installs/checks pinned R dependencies. |
| `scripts/download_data.R` | Retrieves/verifies the 32 public inputs. |
| `scripts/config.R` | Relative data, output and optional library locations; fixed R and direct dependency versions. |
| `scripts/run.R` | Original primary, all-interval and quadratic models; participant descriptions, energy accounting, diagnostics and figure. |
| `scripts/run_secondary.R` | PHQ-by-order, PHQ-by-sex and appetite-deleted models; two Holm-adjusted interaction tests and recall-pattern descriptions. |
| `scripts/code/master_builder.R` | Merges cycle-specific participant, symptom and dietary data. |
| `scripts/code/paired_data.R` | Classifies beverages and additions, constructs paired outcomes/domains, and calculates descriptive energy components. |
| `scripts/code/core_helpers.R` | Builds the full survey parent and computes standardized estimates and Taylor covariance. |
| `scripts/code/paired_engine.R` | Fits original models and derives profiles, contrasts, residual/leverage and survey-influence diagnostics. |
| `scripts/code/secondary_engine.R` | Defines and fits the three exploratory models and recall-pattern tables. |
| `scripts/code/manuscript_descriptives.R` | Produces manuscript participant and selection tables. |
| `scripts/code/reporting.R` | Coordinates descriptive exports and draws Figure 1. |
| `scripts/code/checkpoint_helpers.R` | Stores run identities, owned locks, resumable checkpoints and output receipts. |
| `scripts/public_validation.R` | Checks pinned dependencies and distribution identities; compares fixed aggregate reference values. |
| `scripts/test_validation.R` | Tests valid and deliberately invalid reference-verifier fixtures without touching scientific results. |
| `strengthening/run.R` | Rebuilds inputs and runs the fixed influence, missingness, context-adjustment and response-weighting checks. |
| `strengthening/analysis.R` | Scientific functions for those checks, including response and outcome refitting in 121 JKn replicates. |
| `strengthening/tests.R` | Seven synthetic tests of deletion identities, eligibility, recoding, response replication and same-sample comparisons. |
| `strengthening/config.R` | Fixed eligibility counts and operational response-weighting thresholds. |
| `scripts/SECONDARY_SPEC.md` | Specification of the three earlier exploratory models. |
| `strengthening/PLAN.md` | Scientific transcription of the later fixed exploratory specification and interpretation boundaries. |

`scripts/active_code_manifest.csv`, `scripts/secondary_code_manifest.csv`, and `distribution_manifest.csv` bind the supplied implementation. `reference_manifest.csv` identifies all 48 expected reference files and their generated counterparts. These checks intentionally fail after a scientific code or reference change; investigate and document such changes before creating a distinct version. The reproduction commands never rewrite the expected references.

The repository's `.gitattributes` disables automatic line-ending conversion because the SHA-256 manifests bind exact file bytes. Keep that file when cloning or copying the package, including on Windows with `core.autocrlf=true`.

## Mapping manuscript results to outputs

All paths below are relative to `outputs/` unless marked as a static repository file. Full-precision files are retained; the paper rounds values for presentation.

| Manuscript display | Generated evidence |
| --- | --- |
| Table 1 | `original/descriptives/table1_primary.csv` (and `table1_all_intervals.csv` for the sensitivity cohort). |
| Table 2 | `original/weighted_energy_decomposition.csv`, restricted to `primary_days3to10`, `All PHQ scores`; `total_energy_identity.csv` checks the partition. |
| Table 3 | `original/PHQ_contrasts.csv` (`PHQ10_minus0`) and `model_diagnostics.csv` for cohort sizes and degrees of freedom. |
| Figure 1 | `original/figures/paired_phq_profiles.png`, drawn from `original/PHQ_profiles.csv`. |
| Appendix Table A1 | Static `scripts/table_A1_crosswalk.csv`; implemented definitions are in `paired_data.R` and `paired_engine.R`. No new model fit is required. |
| Appendix Table B1 | `original/model_diagnostics.csv`, `leverage_summary.csv`, and `influence_summary.csv` for `contrast:PHQ10_minus0`. |
| Appendix Table B2 | `original/residual_summaries.csv`, primary model with `grouping=PHQ_group`. |
| Supplementary Table S1 | `secondary/secondary_contrasts.csv`; `interaction_tests.csv` contains original and Holm-adjusted p values; `contrast_covariance.csv` retains covariance. |
| Supplementary Table S2 | `secondary/recall_patterns_by_sex_PHQ.csv`; `subgroup_support.csv` documents model support. |
| Supplementary Table S3 | `strengthening/results/C_context_comparison.csv`, `D_response_comparison.csv`, and `original_Taylor_comparisons.csv`. |
| Supplementary Table S4 | `strengthening/results/A_summary.csv`; all individual and stratum deletion results are generated locally as `A_cases.csv` and `A_strata.csv`. |
| Supplementary Table S5 | `strengthening/results/B_description.csv`, supported by `B_flow.csv` and `B_item_patterns.csv`. |

Narrative sample flow comes from `original/flow.csv`, `analysis_domain_flow.csv`, and `descriptives/selection_age_sex.csv`. Interval composition comes from `timing_distribution.csv` and `timing_mean.csv`. The primary trend test is in `primary_test.csv`. Response-probability, weight and replicate checks are in `strengthening/results/D_response_diagnostics.json`, `D_weight_diagnostics.csv`, and `D_replicate_diagnostics.csv`.

The package reproduces the statistical analyses, full-precision tables and plot. It does **not** regenerate manually edited Word prose, red text, bibliography formatting, or final document pagination. The public reference PNG and newly rendered PNG can differ in label pixels because of font/device rendering; all plotted estimates and interval endpoints are independently checked in the numeric references.

## Scientific scope and provenance

The release was tested on Windows with R 4.5.1 and the pinned package versions. A fresh three-stage analysis rebuilt all outputs from verified raw inputs and matched all 48 fixed reference files. The seven sensitivity tests and four reference-verifier fixture checks passed. The default dependency check passed against installed pinned versions. An isolated installation successfully retrieved the pinned Rcpp source and stopped because the test machine lacked Rtools45; a complete isolated dependency installation has therefore not been demonstrated on that machine. Install the compiler prerequisite before using archived compiled packages. No alternate dependency version was accepted.

The first three models estimate a common PHQ slope or a quadratic PHQ relationship. The next three were specified after the original results were known. The later influence, incomplete-response, additional-adjustment and response-weighted checks were specified after all six earlier results were known and before their own estimates were inspected. All additions remain exploratory. Neither a later code release nor a successful rerun changes that timing.

The scientific analysis modules and synthetic tests retain the original manuscript implementation. Distribution changes concern generic paths, dependency/data setup, reference checking and removal of machine-bound history. Freshly generated original-model checkpoints support the later coefficient/covariance reconstruction checks. Independently saved aggregate references additionally verify the original fits, all later comparisons and response diagnostics. The fixed sensitivity rules still stop inference on failed convergence, rank, replication or response-probability checks; no alternative model is selected automatically.

The main limitation remains the between-person symptom comparison despite the within-person dietary outcome. Order and mode coincide in the analyzed cycles. Response weighting relies on untestable conditional-response assumptions and does not address the selection of people reporting alcohol on exactly one recall. These analyses do not establish an appetite mechanism, long-term energy compensation, weight change, or a causal dietary recommendation.

## Licensing

Original analysis code is provided under the MIT license in `LICENSE`. NHANES data and external R packages retain their own terms; the code license does not relicense them. Cite NHANES and the relevant statistical packages when using their data or methods. Until a manuscript citation is available, identify the repository URL and exact release tag or commit used for reproduction.
