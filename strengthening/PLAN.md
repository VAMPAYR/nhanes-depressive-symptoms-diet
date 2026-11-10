# Fixed exploratory sensitivity specification

Distribution transcription of the scientific specification fixed on 8 November 2026.
The six earlier model results were known when these checks were planned; the later
outcome estimates had not yet been inspected. This is exploratory work, not a
prospective public preregistration. The original local specification SHA-256 was
`d3802fb595b6dde92da7184ebbc7057ecfe3696e86c0c44134e2b9148e124e24`.
This distribution transcription has a different hash because internal manuscript
editing instructions and machine-specific workflow records are omitted. The
statistical functions are unchanged. The repository validates newly generated
outputs against immutable aggregate reference values instead of historical
machine-bound checkpoint identities.

## Fixed additions

### A. Exact deletion influence

Question: how much do the primary and quadratic PHQ 10-versus-0 point contrasts depend on one person or one survey stratum?

Reconstruct both existing fits with the exact current formulas, domain and weights. Require agreement with saved estimates and covariance before proceeding. Use all 2,320 leave-one-person-out weighted least-squares deletion identities for each fit. Compute every leave-one-stratum-out point contrast on the same full-rank coefficient space. Deleting a stratum avoids interpreting a lone-PSU variance adjustment as evidence about influence. These are point-estimate diagnostics, not exclusion rules or alternative inferential samples.

Report full-sample contrast, deletion estimate range, maximum absolute shift, maximum shift divided by original standard error, and rank failures. Preserve all deletion results in machine-readable output. Validate the exact identities against direct weighted least-squares refits for the first case, maximum leverage and maximum absolute contrast shift. Do not delete all observations above a leverage threshold; the flagged set carries substantial weight. Do not select exclusions for p values. Do not attach a confidence interval to the range of deletion estimates.

### B. PHQ-independent eligibility and response assessment

Reconstruct all original inclusion criteria except PHQ completeness. Verify that the complete-PHQ subset exactly recovers the original 2,320 participant keys. Distinguish entirely absent PHQ records, present records with all invalid items, partial completion and special response codes where the raw source permits this distinction. Report the 2,424/2,320/104 flow and item-completion pattern, with original dietary weights.

Describe complete versus incomplete responders by age, sex, cycle, drinking-day recall order, paired other-intake difference, mean daily energy and drinking-day beverage ethanol. This is a selection description, with no significance screen and no claim to reconstruct missing symptom scores. Include weighted summaries and unweighted denominators. It does not address participants reporting alcohol on both or neither day.

### C. One context-adjustment sensitivity

Primary model only. Add family income-to-poverty ratio continuously on its observed 0–5 top-coded scale; education in four categories (less than high school, high school/GED, some college/associate degree, college or higher); five-category RIDRETH1 as a social/context covariate; smoking as never/former/current using SMQ020/SMQ040 and documented skips. Refusal, unknown and missing values remain missing.

These contemporaneous variables describe an extended adjusted association. Their inclusion is not proof of causal confounding control. Race/ethnicity is not interpreted as a biological cause. BMI and medications can be related to symptoms, treatment or reporting and lack a sufficiently specified role for this fixed block. Dose adjustment changes the comparison and is not added to the primary model.

Education represents attained socioeconomic circumstances; income-to-poverty represents current material resources; race/ethnicity represents measured social/context differences; smoking represents a related health behavior and intake/reporting context. Each could relate to symptoms and the paired outcome, but their temporal roles are not established here and some may be consequences or intermediates. The purpose is to test sensitivity to this explicitly conditional association, not to claim an identified confounder set. Do not cross this block with response weighting or add interactions.

Fit the original baseline on the exact 2,175 complete-context participants, then the extended model on those same people, preserving the full survey parent for domain inference. Compare full-sample primary, restricted-sample baseline and extended estimates with 95% intervals, sample sizes, higher-symptom counts, residual degrees of freedom and model ranks. This separates sample restriction from adjustment. Report all three regardless of direction or significance. Stop rather than collapse categories after results if the prespecified matrix is rank-deficient.

### D. Conditional response-weighted sensitivity

Estimate the PHQ association for all 2,424 otherwise eligible pairs under conditional response exchangeability, using the 2,320 observed symptom scores. This is a broader response target, not the primary complete-case estimand.

Fixed survey-weighted logistic response model: complete PHQ indicator on centered age/10, its square, sex, survey cycle, drinking-day recall-order sign and the observed paired other-intake outcome per 1,000 kcal. The outcome is observed in the eligible nonresponders and is included because selection on a missing predictor may relate to the observed outcome. These are all observed predictors. The identifying assumption is that PHQ completeness is independent of the unobserved PHQ score conditional on these observed predictors within the 2,424-person target, with positive response probabilities and a correctly specified response model. This assumption is untestable from these data. Do not choose response terms using association results, significance or AIC. Audit source-record presence and documented administration exclusions before fitting; do not infer a person's reason for nonresponse from an all-missing row.

Multiply original dietary weights by inverse fitted response probability for complete responders. A common stabilization constant does not change the regression estimate. Retain the original outcome model and PHQ contrast. Refit BOTH response and outcome models within a survey JKn replication design built from the full positive-weight parent, with zero contribution outside eligibility. Use `as.svrepdesign(type='JKn', mse=TRUE)` and the package's documented scale and replicate scales, centered on the full-sample estimate. This accounts for estimated response weights. Use a t reference with the replicate design's computed degrees of freedom (`survey::degf`), explicitly reporting the value rather than silently inheriting the original regression residual df. Report a dietary-weight-only baseline computed using this identical JKn scheme and reference alongside the response-weighted result. Check baseline point agreement with the original Taylor result, without requiring identical variance estimators. Preserve the original Taylor interval as the primary result.

Report response rates and probability/factor diagnostics across all 2,424 eligible people, plus complete-responder weight concentration/effective sample size, rank/convergence and replicate failures. Check the response model and positivity BEFORE fitting/exporting the adjusted outcome estimate. No ad hoc truncation, alternate propensity model or outcome model follows an unfavorable result. If probabilities below 0.05 or factors above 20 occur, end this planned response sensitivity with diagnostics and no inferential correction claim; preserve any already computed intermediate results internally. These thresholds are operational extrapolation flags, not proof of sufficient overlap. Any failed replicate or nonconvergence stops inferential reporting until a defined computational cause is resolved and reviewed. A method change requires an amended plan before another result is examined.

No multiple imputation is planned: 100 of 104 eligible nonresponders have no valid symptom items, and adding a symptom imputation model would require additional assumptions beyond this bounded response sensitivity. Response weighting cannot resolve nonignorable missingness; report this explicitly.

## Outcomes and interpretation rules

- If estimates remain similar, report the amount of sensitivity actually tested; do not claim universal robustness, equivalence or no bias.
- If deletion shifts are substantial, report dependence and narrow the interpretation. Do not discard influential valid observations.
- If context adjustment changes the estimate, distinguish the sample-restriction change from the adjustment change; do not select one as the preferred causal answer.
- If response adjustment changes the estimate or increases uncertainty, retain the primary result and report the sensitivity with its assumptions.
- If a check fails its design assumptions, report the failure and retain the corresponding limitation. Failure is an informative result, not permission for unplanned substitutions.
- Broad intervals remain reportable. No new 100-kcal margin or retrospective equivalence claim will be introduced.

