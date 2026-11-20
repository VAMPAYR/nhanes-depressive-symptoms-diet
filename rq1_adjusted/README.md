# RQ1 adjusted energy partition (fourth stage)

This stage estimates the drinking-day energy difference and its beverage and nonbeverage components after adjustment for recall order and weekday. It answers RQ1 on the same adjustment basis that the manuscript requires for drinking-day contrasts and that Breslow et al. (2013) used. The unadjusted Table 2 and the six existing models are unchanged.

## Status

Post hoc sensitivity analysis. It was specified on 20 November 2026, after the original results were known and before any adjusted model was fitted, following an independent review that found the RQ1 partition unadjusted. `SPEC.md` is the frozen specification; the runner stops if its SHA-256 changes. The result is reported whatever it shows. No model selection, imputation or refitting occurs.

## Model

For each component (ethanol energy, calculated beverage-energy remainder, linked additions, other foods and beverages, total energy), a survey-weighted linear model regresses the drinking-day-minus-comparison-day difference on the recall-order sign (+1 drinking on recall 2, −1 on recall 1) and the six weekday-indicator differences (Sunday reference). The intercept is the difference at equal recall-order weighting with both recalls on the same weekday. The domain is the primary 3–10-day sample (n = 2,320) in the primary model's design: two-day weight / 4, released strata and nested PSUs, Taylor linearization. Intervals use t with the fitted residual degrees of freedom (55).

## Checks performed by the runner

1. R and package versions match the existing stages; the distribution specification hash and `distribution_manifest.csv` match.
2. The 28 raw XPORT files match `scripts/raw_source_manifest.csv`.
3. Original and secondary outputs match their fixed aggregate references before and after the run.
4. Synthetic tests: exact recovery of a known intercept and additivity of component intercepts.
5. Gate: unadjusted means rebuilt from raw sources reproduce Table 2 (`scripts/reference_outputs/paired/weighted_energy_decomposition.csv`) within 1e-6 kcal before any adjusted model is read.
6. The four adjusted components sum to the adjusted total.
7. Results agree within 1e-8 with the fixed aggregate reference `reference/rq1_adjusted/rq1_adjusted_results.csv`, and the results and summary match `reference_manifest.csv` at tolerance 1e-9.

## Commands

```text
Rscript --vanilla rq1_adjusted/run.R --preflight
Rscript --vanilla rq1_adjusted/run.R --run      # fresh run needs an empty output folder
Rscript --vanilla rq1_adjusted/run.R --verify
Rscript --vanilla rq1_adjusted/run.R --run --output=C:/absolute/empty/path
```

`RUN_ALL.R` calls this stage after the strengthening stage (release v1.1.0).

## Outputs (`outputs/rq1_adjusted/results/`)

| File | Contents |
|---|---|
| `rq1_adjusted_results.csv` | Adjusted intercept, SE, 95% CI, recall-order coefficient and interval, residual df and n for each component |
| `rq1_adjusted_summary.json` | Beverage sum and share, identity gap, degrees of freedom |
| `gate_table2_reproduction.csv` | Table 2 reproduction |
| `reference_agreement.csv` | Agreement with the first computation |
| `synthetic_tests.csv`, `input_manifest.csv`, `code_manifest.csv`, `output_manifest.csv`, `sessionInfo.txt`, `STARTED.json`, `COMPLETED.json` | Run records |

## Manuscript use

In the manuscript, Section 3.5 describes the model. Section 4.2 reports the adjusted total (339.7 kcal; 95% CI 253.2, 426.2), the components (ethanol 235.2, remainder 93.8, additions 4.9 kcal; 98.3% of the total) and the other-intake difference (5.8 kcal; 95% CI −71.0, 82.6). Sections 5.1 to 5.3 interpret them.
