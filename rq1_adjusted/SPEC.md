# RQ1 order- and weekday-adjusted energy partition: fixed specification

Distribution transcription of the specification fixed on 20 November 2026, before any adjusted model was fitted. The original results were known when this sensitivity analysis was specified, after an independent review found the RQ1 energy partition unadjusted. The original local specification SHA-256 was `6444e4e208826647d6301e6bb626ee720409fb47dd0fc2534e820b6f0f5e11ff`. This transcription has a different hash because a machine-specific data location is replaced by the public raw-file rebuild; the statistical content is unchanged. The result is reported whatever it shows.

## Question

How much higher is reported energy on the drinking day, and how much of that difference lies in alcoholic-beverage components, after adjustment for recall order and weekday?

## Data and domain

- Data: the paired master rebuilt from the 28 verified raw XPORT files listed in `scripts/raw_source_manifest.csv` (39,156 rows). The first computation used a hashed cached copy of the same master; the distribution rebuilds it from raw files.
- Domain: the primary 3–10-day domain from `pe_prepare()` in `scripts/code/paired_engine.R` (expected n = 2,320).
- Design: `parent_design(data, weight = "WTDR2D")` from `scripts/code/core_helpers.R` (two-day weight / 4, strata, nested PSUs), then `subset()` to the domain, exactly as in the primary model.

## Gate (must pass before the adjusted models are read)

An intercept-only survey-weighted model of `pd_delta_total_kcal` in the domain must reproduce Table 2: 396.3 kcal (95% CI 315.7, 476.9), to 0.1 kcal, using the design's degrees of freedom. The other four components must reproduce 244.6, 98.1, 5.8 and 47.8 kcal. If the gate fails, stop and report; fit nothing else.

## Models

For each component c in {beverage ethanol, calculated beverage-energy remainder, linked additions, other foods and beverages, total}:

`svyglm(pd_delta_c ~ pd_order_sign + pd_weekday2 + ... + pd_weekday7, family = gaussian(identity))`

- `pd_order_sign` is +1 when drinking is on recall 2 and −1 when on recall 1; the weekday terms are drinking-day-minus-comparison-day indicator differences (Sunday reference), as in the primary model.
- The intercept is the adjusted drinking-day difference at equal recall-order weighting (order sign 0) and for the same weekday on both days (weekday differences 0). This matches the within-person adjustment for recall order and day of week used by Breslow et al. (2013).
- No PHQ, sex, age or cycle terms: these do not vary within person, and RQ1 concerns the drinking-day difference itself.
- Inference: Taylor linearization; 95% CIs from t with the fitted model's residual degrees of freedom (`fit$df.residual`), the policy of the existing models.

## Reported quantities

1. Adjusted intercept and 95% CI for each component and the total.
2. Identity check: the four component intercepts sum to the total intercept (same design, linear models; tolerance 1e-6 kcal).
3. Beverage share: (ethanol + remainder + additions intercepts) / total intercept, point estimate only.
4. The `pd_order_sign` coefficient for other intake and total, for context.

## Use in the manuscript

The unadjusted Table 2 stays as it is. The adjusted total, other-intake difference and beverage share are added in text (Method 3.5, Results 4.2, Discussion 5.1 and 5.2), labeled as a sensitivity analysis added after review. No other model is fitted, and no specification is changed after seeing the results.
