# Source-only paired-recall construction. No participant read occurs on source.
# Official local sources: item_codebooks/DR1IFF_J.txt:83-85,134-177,240-269,
# 353-375,721-727; dietary_sources/codebooks/DR2TOT_J.txt:93-96,186-215;
# item_codebooks/DRXFCD_G.txt:132-152; sources/breslow_2013_fulltext.txt:59.
pp_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
pp_cycles <- function() c("2011-2012" = "G", "2013-2014" = "H", "2015-2016" = "I", "2017-2018" = "J")
pp_positive <- function(x) is.finite(x) & x > 0
pp_key <- function(release, seqn) paste(as.character(release), sprintf("%.0f", seqn), sep = ":")
pp_id_valid <- function(x) is.numeric(x) && all(is.finite(x) & x > 0 & x == floor(x))
pp_require <- function(d, fields, label) pp_assert(is.data.frame(d) && all(fields %in% names(d)),
  paste(label, "lacks required fields:", paste(setdiff(fields, names(d)), collapse = ", ")))
pp_equal <- function(a, b, tolerance = 1e-8) length(a) == length(b) &&
  all((is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= tolerance * pmax(1, abs(a), abs(b))))
pp_tot_fields <- function(day) c("SEQN", "WTDR2D", paste0("DR", day, c("DRSTZ", "DBIH", "DAY", "TKCAL", "TALCO")))
pp_iff_fields <- function(day) c("SEQN", paste0("DR", day, c("ILINE", "DRSTZ", "IFDCD", "IKCAL", "IALCO", "CCMNM", "CCMTX")))
pp_components <- function() c(beverage_ethanol_kcal = "Ethanol energy in coded beverages",
  beverage_nonethanol_kcal = "Calculated coded-beverage energy remainder", linked_addition_kcal = "Energy in linked beverage additions",
  other_kcal = "Energy from other reported foods and beverages", total_kcal = "Total reported energy")

pp_read_xpt <- function(path, fields) {
  pp_assert(requireNamespace("foreign", quietly = TRUE), "Installed native R package foreign is required; do not substitute haven.")
  pp_assert(file.exists(path), paste("Required XPORT source is absent:", path))
  con <- file(path, open = "rb")
  header <- tryCatch(readBin(con, what = "raw", n = 80L), finally = close(con))
  pp_assert(length(header) == 80L && startsWith(rawToChar(header), "HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!"),
    paste("Invalid XPORT library header:", path))
  raw <- foreign::read.xport(path)
  pp_require(raw, fields, basename(path))
  # Reduce immediately; do not retain the full individual-food file.
  raw[fields]
}

pp_initialize <- function(master, canonical = FALSE) {
  fields <- c("PERSON_KEY", "SEQN", "RELEASE_ID", "RIDAGEYR", "RIAGENDR", "WTDR2D", "SDMVSTRA", "SDMVPSU",
    "DR1DRSTZ", "DR2DRSTZ", "DR1TKCAL", "DR2TKCAL", "PHQ9_TOTAL", sprintf("DPQ%03d", seq(10, 90, 10)))
  pp_require(master, fields, "Canonical master")
  pp_assert(!canonical || nrow(master) == 39156L, "The canonical master must retain all 39,156 people.")
  pp_assert(!any(startsWith(names(master), "pd_")), "The paired-day builder cannot overwrite existing pd_ fields.")
  pp_assert(pp_id_valid(master$SEQN) && !anyDuplicated(master$SEQN), "Master SEQN is invalid, duplicated, or reused across releases.")
  pp_assert(all(as.character(master$RELEASE_ID) %in% names(pp_cycles())), "Master contains an unsupported release.")
  pp_assert(!anyNA(master$PERSON_KEY) && !anyDuplicated(master$PERSON_KEY), "Master person keys are invalid.")
  master
}

pp_validate_descriptions <- function(descriptions) {
  pp_require(descriptions, c("DRXFDCD", "DRXFCSD", "DRXFCLD"), "Food descriptions")
  pp_assert(pp_id_valid(descriptions$DRXFDCD) && !anyDuplicated(descriptions$DRXFDCD), "Food-description codes are invalid or duplicated.")
  pp_assert(is.character(descriptions$DRXFCSD) && is.character(descriptions$DRXFCLD), "Food descriptions must contain text labels.")
  invisible(TRUE)
}

# Pure reducer for one release/day. The synthetic tests use this same path.
pp_reduce_day <- function(master, tot, iff, descriptions, cycle, day) {
  pp_assert(cycle %in% names(pp_cycles()) && day %in% 1:2, "Unsupported release/day.")
  tf <- pp_tot_fields(day); inf <- pp_iff_fields(day)
  pp_require(tot, tf, "Total nutrient source"); pp_require(iff, inf, "Individual-food source")
  pp_validate_descriptions(descriptions)
  pp_assert(all(vapply(tot[tf], is.numeric, logical(1))) && all(vapply(iff[inf], is.numeric, logical(1))), "Recall source fields must be numeric.")
  pp_assert(pp_id_valid(tot$SEQN) && !anyDuplicated(tot$SEQN), "Total nutrient SEQN is invalid or duplicated.")
  names_tot <- c("SEQN", "WEIGHT", "STATUS", "DBIH", "WEEKDAY", "KCAL", "ALCO")
  names_iff <- c("SEQN", "LINE", "STATUS", "FOODCODE", "KCAL", "ALCO", "COMBONUM", "COMBOTYPE")
  tot <- tot[tf]; names(tot) <- names_tot
  iff <- iff[inf]; names(iff) <- names_iff
  pp_assert(pp_id_valid(iff$SEQN) && pp_id_valid(iff$LINE) && !anyDuplicated(iff[c("SEQN", "LINE")]), "Individual-food release/person/day/line keys are invalid or duplicated.")
  source_to_master <- match(tot$SEQN, master$SEQN)
  present <- !is.na(source_to_master)
  pp_assert(all(as.character(master$RELEASE_ID[source_to_master[present]]) == cycle), "A total nutrient source belongs to another release.")
  rows <- which(as.character(master$RELEASE_ID) == cycle)
  master_keys <- pp_key(master$RELEASE_ID[rows], master$SEQN[rows])
  idx <- match(master_keys, pp_key(cycle, tot$SEQN))
  numeric_idx <- match(master$SEQN[rows], tot$SEQN)
  pp_assert(identical(idx, numeric_idx), "Normalized release/SEQN keys disagree with numeric matching.")
  matched <- !is.na(idx)
  controls <- c(WTDR2D = "WEIGHT", setNames("STATUS", paste0("DR", day, "DRSTZ")), setNames("KCAL", paste0("DR", day, "TKCAL")))
  for (field in names(controls)) pp_assert(pp_equal(master[[field]][rows[matched]], tot[[controls[[field]]]][idx[matched]]),
    paste("Existing master/source control disagrees:", cycle, day, field))
  item_tot <- match(iff$SEQN, tot$SEQN)
  pp_assert(!anyNA(item_tot), "An individual-food person is absent from its release/day total source.")
  reliable_items <- iff$STATUS %in% 1
  pp_assert(all(tot$STATUS[item_tot[reliable_items]] %in% 1), "Reliable item status disagrees with total recall status.")
  all_item_n <- nrow(iff)
  iff <- iff[reliable_items, , drop = FALSE]
  pp_assert(all(is.finite(iff$KCAL) & iff$KCAL >= 0 & is.finite(iff$ALCO) & iff$ALCO >= 0), "Reliable item energy or ethanol is unavailable or negative.")
  pp_assert(pp_id_valid(iff$FOODCODE), "Reliable item food codes are invalid.")
  pp_assert(all(is.finite(iff$COMBONUM) & iff$COMBONUM >= 0 & iff$COMBONUM == floor(iff$COMBONUM)) &&
    all(is.finite(iff$COMBOTYPE) & iff$COMBOTYPE >= 0 & iff$COMBOTYPE == floor(iff$COMBOTYPE)), "Reliable combination fields are invalid.")
  di <- match(iff$FOODCODE, descriptions$DRXFDCD)
  pp_assert(!anyNA(di), "A reliable item lacks a matching-release food description.")
  category <- iff$FOODCODE >= 93100000 & iff$FOODCODE <= 93599999 & iff$FOODCODE != 93401300
  coded <- category & iff$ALCO > 0
  combo <- paste(sprintf("%.0f", iff$SEQN), sprintf("%.0f", iff$COMBONUM), sep = ":")
  positive_combo <- coded & iff$COMBONUM > 0 & iff$COMBOTYPE == 1
  linked <- !coded & iff$COMBONUM > 0 & iff$COMBOTYPE == 1 & combo %in% unique(combo[positive_combo])
  other <- !(coded | linked)
  pp_assert(!any(coded & linked) && all(as.integer(coded) + as.integer(linked) + as.integer(other) == 1L), "Beverage/addition exclusions overlap or fail to partition items.")
  values <- cbind(item_kcal = iff$KCAL, item_alcohol_g = iff$ALCO,
    beverage_kcal = ifelse(coded, iff$KCAL, 0), beverage_alcohol_g = ifelse(coded, iff$ALCO, 0),
    linked_addition_kcal = ifelse(linked, iff$KCAL, 0), other_kcal = ifelse(other, iff$KCAL, 0),
    beverage_item_n = as.integer(coded), linked_addition_item_n = as.integer(linked),
    zero_ethanol_category_item_n = as.integer(category & iff$ALCO == 0), item_n = rep(1, nrow(iff)))
  if (nrow(iff)) {
    sums <- rowsum(values, group = sprintf("%.0f", iff$SEQN), reorder = FALSE)
    gi <- match(sprintf("%.0f", tot$SEQN), rownames(sums))
    aggregated <- sums[gi, , drop = FALSE]
  } else { gi <- rep(NA_integer_, nrow(tot)); aggregated <- matrix(NA_real_, nrow(tot), ncol(values), dimnames = list(NULL, colnames(values))) }
  reliable <- tot$STATUS %in% 1
  valid_tot <- is.finite(tot$KCAL) & tot$KCAL >= 0 & is.finite(tot$ALCO) & tot$ALCO >= 0
  pp_assert(all(valid_tot[reliable]), "Reliable total energy or ethanol is unavailable or negative.")
  fasting <- reliable & is.na(gi) & tot$KCAL == 0 & tot$ALCO == 0
  aggregated[fasting, ] <- 0
  rownames(aggregated) <- NULL
  kd <- aggregated[, "item_kcal"] - tot$KCAL
  ad <- aggregated[, "item_alcohol_g"] - tot$ALCO
  pp_assert(all(is.finite(kd[reliable]) & abs(kd[reliable]) <= 1e-6) &&
    all(is.finite(ad[reliable]) & abs(ad[reliable]) <= 1e-6), "Reliable item kcal/ethanol sums do not reconcile to total nutrients.")
  component <- data.frame(total_kcal = tot$KCAL, total_alcohol_g = tot$ALCO,
    beverage_kcal = aggregated[, "beverage_kcal"], beverage_alcohol_g = aggregated[, "beverage_alcohol_g"],
    beverage_ethanol_kcal = 6.93 * aggregated[, "beverage_alcohol_g"],
    beverage_nonethanol_kcal = aggregated[, "beverage_kcal"] - 6.93 * aggregated[, "beverage_alcohol_g"],
    linked_addition_kcal = aggregated[, "linked_addition_kcal"], other_kcal = aggregated[, "other_kcal"],
    beverage_item_n = aggregated[, "beverage_item_n"], linked_addition_item_n = aggregated[, "linked_addition_item_n"],
    weekday = tot$WEEKDAY, dbih = tot$DBIH, status = tot$STATUS, valid = reliable & valid_tot)
  component$beverage_positive <- ifelse(component$valid, component$beverage_alcohol_g > 0, NA)
  identity <- rowSums(component[c("beverage_ethanol_kcal", "beverage_nonethanol_kcal", "linked_addition_kcal", "other_kcal")]) - component$total_kcal
  pp_assert(all(is.finite(identity[reliable]) & abs(identity[reliable]) <= 1e-6), "Daily energy decomposition does not sum to total energy.")
  attached <- component[idx, , drop = FALSE]
  attached$raw_matched <- matched
  attached$valid[!matched] <- FALSE
  names(attached) <- paste0("pd_day", day, "_", names(attached))
  source_unmatched <- sum(!(tot$SEQN %in% master$SEQN[rows]))
  negative_remainders <- component$beverage_nonethanol_kcal[reliable & component$beverage_nonethanol_kcal < 0]
  audit <- data.frame(cycle = cycle, day = day, master_n = length(rows), source_n = nrow(tot),
    matched_n = sum(matched), unmatched_master_n = sum(!matched), unmatched_source_n = source_unmatched,
    all_item_n = all_item_n, reliable_item_n = nrow(iff), reliable_person_n = sum(reliable),
    reliable_zero_without_items_n = sum(fasting), coded_beverage_item_n = sum(coded), linked_addition_item_n = sum(linked),
    zero_ethanol_category_item_n = sum(category & iff$ALCO == 0),
    negative_coded_nonethanol_energy_n = sum(component$beverage_nonethanol_kcal[reliable] < 0),
    minimum_negative_coded_beverage_remainder_kcal = if (length(negative_remainders)) min(negative_remainders) else 0,
    maximum_negative_coded_beverage_remainder_magnitude_kcal = if (length(negative_remainders)) max(abs(negative_remainders)) else 0,
    max_item_total_kcal_difference = if (any(reliable)) max(abs(kd[reliable])) else 0,
    max_item_total_alcohol_difference = if (any(reliable)) max(abs(ad[reliable])) else 0,
    max_decomposition_difference = if (any(reliable)) max(abs(identity[reliable])) else 0)
  pp_assert(audit$source_n == audit$matched_n + audit$unmatched_source_n &&
    audit$master_n == audit$matched_n + audit$unmatched_master_n, "Cycle/day join counts fail to reconcile.")
  codes <- sort(unique(iff$FOODCODE[coded | linked]))
  mapped <- match(codes, descriptions$DRXFDCD)
  list(rows = rows, data = attached, audit = audit,
    description_audit = data.frame(cycle = rep(cycle, length(codes)), day = rep(day, length(codes)), food_code = codes,
      short_description = descriptions$DRXFCSD[mapped], long_description = descriptions$DRXFCLD[mapped]))
}

pp_finish_data <- function(data, audits, description_audits = list()) {
  phq <- as.matrix(data[sprintf("DPQ%03d", seq(10, 90, 10))])
  strict <- rowSums(matrix(phq %in% 0:3, nrow = nrow(data))) == 9L
  total <- rep(NA_real_, nrow(data)); total[strict] <- rowSums(phq[strict, , drop = FALSE])
  pp_assert(pp_equal(data$PHQ9_TOTAL[strict], total[strict]), "Strict nine-item PHQ total disagrees with the master.")
  data$pd_phq_total <- total
  data$pd_phq_band <- factor(ifelse(!strict, NA_character_, ifelse(total < 5, "0-4", ifelse(total < 10, "5-9", ">=10"))), levels = c("0-4", "5-9", ">=10"))
  design_valid <- is.finite(data$SDMVSTRA) & is.finite(data$SDMVPSU)
  two_valid <- data$pd_day1_valid & data$pd_day2_valid
  data$pd_base_eligible <- is.finite(data$RIDAGEYR) & data$RIDAGEYR >= 20 & pp_positive(data$WTDR2D) & design_valid & strict & two_valid
  data$pd_discordant <- two_valid & xor(data$pd_day1_beverage_positive %in% TRUE, data$pd_day2_beverage_positive %in% TRUE)
  data$pd_pair_pattern <- rep("unavailable", nrow(data))
  for (a in c(FALSE, TRUE)) for (b in c(FALSE, TRUE)) {
    mask <- two_valid & data$pd_day1_beverage_positive %in% a & data$pd_day2_beverage_positive %in% b
    data$pd_pair_pattern[mask] <- if (a && b) "both_positive" else if (a) "day1_positive_day2_zero" else if (b) "day1_zero_day2_positive" else "both_zero"
  }
  data$pd_gap_days <- data$pd_day2_dbih - data$pd_day1_dbih
  data$pd_order_sign <- rep(NA_real_, nrow(data))
  data$pd_order_sign[data$pd_discordant] <- ifelse(data$pd_day2_beverage_positive[data$pd_discordant], 1, -1)
  weekday_valid <- data$pd_day1_weekday %in% 1:7 & data$pd_day2_weekday %in% 1:7
  for (day in 2:7) {
    contrast <- data$pd_order_sign * (as.numeric(data$pd_day2_weekday %in% day) - as.numeric(data$pd_day1_weekday %in% day))
    contrast[!weekday_valid] <- NA_real_
    data[[paste0("pd_weekday", day)]] <- contrast
  }
  for (component in names(pp_components())) data[[paste0("pd_delta_", component)]] <-
    data$pd_order_sign * (data[[paste0("pd_day2_", component)]] - data[[paste0("pd_day1_", component)]])
  delta_identity <- rowSums(data[paste0("pd_delta_", names(pp_components())[1:4])]) - data$pd_delta_total_kcal
  pp_assert(all(abs(delta_identity[data$pd_discordant]) <= 1e-6), "Oriented component differences do not sum to the total difference.")
  timing <- do.call(rbind, lapply(names(pp_cycles()), function(cycle) {
    selected <- as.character(data$RELEASE_ID) == cycle & two_valid & weekday_valid & is.finite(data$pd_gap_days)
    mismatch <- selected & ((data$pd_day2_weekday - data$pd_day1_weekday) %% 7 != data$pd_gap_days %% 7)
    data.frame(cycle = cycle, comparable_pair_n = sum(selected), weekday_gap_mismatch_n = sum(mismatch),
      nonpositive_observed_gap_n = sum(selected & data$pd_gap_days <= 0))
  }))
  list(data = data, audit = list(cycle_day = do.call(rbind, audits), food_descriptions = do.call(rbind, description_audits),
    timing_consistency = timing, canonical_master_overwritten = FALSE,
    beverage_rule = "Code 93100000-93599999 except 93401300 AND item ethanol >0; day reference is exactly zero coded-beverage ethanol",
    linked_addition_rule = "Disjoint noncoded items sharing release/person/day and positive CCMNM with a coded positive beverage; both CCMTX=1",
    energy_definition = "Other reported food/beverage energy excludes the union of coded positive-ethanol beverage items and their identified linked additions; unlinked mixers remain unresolved",
    coded_beverage_energy_remainder = "Signed tabulated coded-beverage item energy minus 6.93 kcal/g times coded-beverage ethanol grams; negative values are retained without attributing their cause",
    negative_remainder_audit = "Minimum negative value and maximum negative magnitude are both zero when no negative remainder occurs; beverage_nonethanol_kcal remains the stable machine field",
    energy_multiplier_kcal_per_g = 6.93, zero_energy_policy = "Reliable zero-energy days are retained, including zero totals without item rows",
    interval_definition = "DR2DBIH minus DR1DBIH; DRDINT is never used as elapsed time; weekday congruence is diagnostic only",
    order_mode = "Day1 MEC versus day2 telephone; recall order and mode cannot be separated",
    weight_parent = "All positive WTDR2D carriers with valid strata/PSUs; person-level WTDR2D/4; all original rows retained"))
}

pp_build_from_sources <- function(master, totals, items, descriptions) {
  data <- pp_initialize(master)
  pp_assert(all(vapply(list(totals, items, descriptions), function(collection) is.list(collection) &&
    length(collection) == 4L && !anyDuplicated(names(collection)) &&
    setequal(names(collection), names(pp_cycles())), logical(1))), "Synthetic/source collections must contain exactly four unique cycles.")
  pp_assert(all(vapply(c(totals, items), function(days) is.list(days) && length(days) == 2L, logical(1))),
    "Each recall collection must contain exactly two days.")
  audits <- mapping <- list()
  for (cycle in names(pp_cycles())) for (day in 1:2) {
    reduced <- pp_reduce_day(master, totals[[cycle]][[day]], items[[cycle]][[day]], descriptions[[cycle]], cycle, day)
    for (field in names(reduced$data)) {
      if (!field %in% names(data)) data[[field]] <- if (is.logical(reduced$data[[field]])) rep(FALSE, nrow(data)) else rep(NA_real_, nrow(data))
      data[[field]][reduced$rows] <- reduced$data[[field]]
    }
    audits[[length(audits) + 1L]] <- reduced$audit; mapping[[length(mapping) + 1L]] <- reduced$description_audit
  }
  pp_finish_data(data, audits, mapping)
}

pp_build_data <- function(master, root) {
  pp_assert(requireNamespace("digest", quietly = TRUE), "Installed R package digest is required for source hashes.")
  data <- pp_initialize(master, canonical = TRUE)
  audits <- mapping <- manifest <- list()
  for (cycle in names(pp_cycles())) {
    suffix <- pp_cycles()[[cycle]]
    description_path <- file.path(root, cycle, paste0("DRXFCD_", suffix, ".xpt"))
    descriptions <- pp_read_xpt(description_path, c("DRXFDCD", "DRXFCSD", "DRXFCLD"))
    pp_validate_descriptions(descriptions)
    manifest[[length(manifest) + 1L]] <- data.frame(cycle = cycle, day = NA_integer_, kind = "descriptions",
      path = normalizePath(description_path, winslash = "/"), sha256 = digest::digest(file = description_path, algo = "sha256", serialize = FALSE))
    for (day in 1:2) {
      paths <- file.path(root, cycle, paste0("DR", day, c("TOT_", "IFF_"), suffix, ".xpt"))
      tot <- pp_read_xpt(paths[1], pp_tot_fields(day))
      iff <- pp_read_xpt(paths[2], pp_iff_fields(day))
      reduced <- pp_reduce_day(master, tot, iff, descriptions, cycle, day)
      pp_assert(reduced$audit$unmatched_source_n == 0L, "A total nutrient source is outside the canonical master roster.")
      rm(tot, iff)
      for (field in names(reduced$data)) {
        if (!field %in% names(data)) data[[field]] <- if (is.logical(reduced$data[[field]])) rep(FALSE, nrow(data)) else rep(NA_real_, nrow(data))
        data[[field]][reduced$rows] <- reduced$data[[field]]
      }
      audits[[length(audits) + 1L]] <- reduced$audit; mapping[[length(mapping) + 1L]] <- reduced$description_audit
      for (i in 1:2) manifest[[length(manifest) + 1L]] <- data.frame(cycle = cycle, day = day,
        kind = c("totals", "items")[i], path = normalizePath(paths[i], winslash = "/"),
        sha256 = digest::digest(file = paths[i], algo = "sha256", serialize = FALSE))
    }
    rm(descriptions)
  }
  result <- pp_finish_data(data, audits, mapping)
  result$audit$source_files <- do.call(rbind, manifest)
  result$audit$reader <- paste0("foreign::read.xport ", utils::packageVersion("foreign"), "; each file parsed once; no haven or SAS/FPED reads")
  result
}

pp_descriptive_domains <- function(data) {
  pp_require(data, c("pd_base_eligible", "pd_discordant", "pd_gap_days", "pd_delta_other_kcal", "pd_order_sign",
    paste0("pd_weekday", 2:7), "RIDAGEYR", "RIAGENDR", "RELEASE_ID", "PHQ9_TOTAL"), "Paired data")
  valid <- data$pd_base_eligible & data$pd_discordant & is.finite(data$pd_delta_other_kcal) &
    data$pd_order_sign %in% c(-1, 1) & is.finite(data$RIDAGEYR) & data$RIAGENDR %in% 1:2 &
    as.character(data$RELEASE_ID) %in% names(pp_cycles()) & is.finite(data$PHQ9_TOTAL) &
    rowSums(is.finite(as.matrix(data[paste0("pd_weekday", 2:7)]))) == 6L
  list(all_intervals = valid, primary_days3to10 = valid & is.finite(data$pd_gap_days) &
    data$pd_gap_days >= 3 & data$pd_gap_days <= 10)
}

pp_support <- function(design, mask) {
  z <- design$variables[mask, c(".weight", ".stratum", ".psu"), drop = FALSE]
  psus <- nrow(unique(z[c(".stratum", ".psu")]))
  strata <- length(unique(z$.stratum))
  parent <- attr(design, "pp_parent_support")
  if (is.null(parent)) parent <- list(n = nrow(design$variables),
    psus = nrow(unique(design$variables[c(".stratum", ".psu")])),
    strata = length(unique(design$variables$.stratum)), df = survey::degf(design))
  data.frame(n = nrow(z), weight_denominator = sum(z$.weight), psus = psus, strata = strata,
    domain_df = psus - strata, parent_n = parent$n, parent_psus = parent$psus,
    parent_strata = parent$strata, parent_df = parent$df)
}
pp_cell_status <- function(s) if (s$n == 0L) "empty" else if (s$n < 5L) "n_below_5" else
  if (s$psus < 2L) "fewer_than_2_psus" else if (s$domain_df <= 0) "nonpositive_domain_df" else "reported"

pp_joint_mean <- function(design, mask, variables) {
  pp_assert(is.logical(mask) && length(mask) == nrow(design$variables) && !anyNA(mask), "Invalid descriptive domain.")
  x <- as.matrix(design$variables[variables])
  pp_assert(all(is.finite(x[mask, , drop = FALSE])), "Unavailable values entered an energy summary.")
  support <- pp_support(design, mask)
  status <- pp_cell_status(support)
  estimate <- setNames(rep(NA_real_, length(variables)), variables)
  covariance <- matrix(NA_real_, length(variables), length(variables), dimnames = list(variables, variables))
  if (support$n) estimate <- colSums(x[mask, , drop = FALSE] * design$variables$.weight[mask]) / support$weight_denominator
  if (status == "reported") {
    fit <- survey::svymean(stats::reformulate(variables), design[mask, ], na.rm = FALSE)
    estimate <- stats::coef(fit); covariance <- stats::vcov(fit)
    if (!all(is.finite(covariance))) { covariance[,] <- NA_real_; status <- "variance_unavailable" }
  }
  se <- sqrt(pmax(0, diag(covariance)))
  delta <- if (status == "reported") qt(.975, support$domain_df) * se else rep(NA_real_, length(variables))
  list(summary = cbind(data.frame(variable = variables, estimate = as.numeric(estimate), se = as.numeric(se),
    ci_lower = as.numeric(estimate - delta), ci_upper = as.numeric(estimate + delta), ci_status = status),
    support[rep(1L, length(variables)), , drop = FALSE]), covariance = covariance)
}

pp_proportion <- function(design, domain, membership) {
  pp_assert(length(domain) == length(membership) && !anyNA(domain) && !anyNA(membership), "Incomplete proportion domain.")
  working <- design; working$variables$.pp_member <- as.numeric(membership)
  result <- pp_joint_mean(working, domain, ".pp_member")$summary
  category <- pp_support(design, domain & membership)
  complement <- pp_support(design, domain & !membership)
  result$category_n <- category$n; result$category_psus <- category$psus; result$category_strata <- category$strata
  if (result$ci_status == "reported" && (category$n < 5L || category$psus < 2L || complement$n < 5L || complement$psus < 2L)) {
    result$se <- result$ci_lower <- result$ci_upper <- NA_real_
    result$ci_status <- "sparse_category_or_complement"
  }
  result$ci_lower <- pmax(0, result$ci_lower); result$ci_upper <- pmin(1, result$ci_upper)
  result$variable <- NULL
  result
}

pp_descriptives <- function(data, core = NULL) {
  pp_assert(requireNamespace("survey", quietly = TRUE), "Installed R package survey is required.")
  domains <- pp_descriptive_domains(data)
  fixed <- c("WTDR2D", "SDMVSTRA", "SDMVPSU", "RELEASE_ID", "RIDAGEYR", "RIAGENDR", "PHQ9_TOTAL",
    "pd_phq_band", "pd_phq_total", "pd_base_eligible", "pd_discordant", "pd_day1_valid", "pd_day2_valid", "pd_order_sign", "pd_gap_days")
  columns <- unique(c(fixed, paste0("pd_weekday", 2:7),
    unlist(lapply(c("day1", "day2", "delta"), function(prefix) paste0("pd_", prefix, "_", names(pp_components()))))))
  pp_require(data, columns, "Paired descriptive data")
  z <- data[columns]; z$.pp_row <- seq_len(nrow(z))
  z$.pp_all <- domains$all_intervals; z$.pp_primary <- domains$primary_days3to10
  if (!is.null(core)) {
    pp_assert((is.list(core) || is.environment(core)) && is.function(core$parent_design), "Invalid frozen survey core.")
    design <- core$parent_design(z, weight = "WTDR2D")
  } else {
    keep <- pp_positive(z$WTDR2D) & is.finite(z$SDMVSTRA) & is.finite(z$SDMVPSU)
    parent <- z[keep, , drop = FALSE]; parent$.weight <- parent$WTDR2D / 4
    parent$.stratum <- parent$SDMVSTRA; parent$.psu <- parent$SDMVPSU
    pp_assert(nrow(parent) > 0L, "No positive two-day survey carrier remains.")
    design <- survey::svydesign(ids = ~.psu, strata = ~.stratum, weights = ~.weight, data = parent, nest = TRUE)
  }
  expected <- pp_positive(z$WTDR2D) & is.finite(z$SDMVSTRA) & is.finite(z$SDMVPSU)
  pp_assert(identical(design$variables$.pp_row, z$.pp_row[expected]) &&
    pp_equal(design$variables$.weight, z$WTDR2D[expected] / 4), "Descriptives changed the full positive-WTDR2D parent.")
  support <- pp_support(design, rep(TRUE, nrow(design$variables)))
  pp_assert(support$domain_df == support$parent_df, "Expected the frozen single-stage nested Taylor design.")
  attr(design, "pp_parent_support") <- list(n = support$n, psus = support$psus, strata = support$strata, df = support$parent_df)
  v <- design$variables
  windows <- list(primary_days3to10 = v$.pp_primary, all_intervals = v$.pp_all)
  adult <- is.finite(v$RIDAGEYR) & v$RIDAGEYR >= 20
  flow_masks <- list(full_two_day_weight_parent = rep(TRUE, nrow(v)), adults = adult,
    adults_strict_PHQ = adult & is.finite(v$pd_phq_total),
    adults_strict_PHQ_two_reliable_valid_days = v$pd_base_eligible,
    eligible_beverage_discordant = v$pd_base_eligible & v$pd_discordant,
    model_valid_all_intervals = v$.pp_all, model_valid_days3to10 = v$.pp_primary)
  flow <- lapply(names(flow_masks), function(name) {
    mask <- flow_masks[[name]]; s <- pp_support(design, mask)
    status <- pp_cell_status(s); se <- lower <- upper <- NA_real_
    if (status == "reported") {
      working <- design; working$variables$.pp_member <- as.numeric(mask)
      fit <- survey::svytotal(~.pp_member, working)
      se <- as.numeric(survey::SE(fit)); delta <- qt(.975, s$parent_df) * se
      lower <- max(0, s$weight_denominator - delta); upper <- s$weight_denominator + delta
    }
    cbind(data.frame(domain = name, population_se = se, population_ci_lower = lower,
      population_ci_upper = upper, population_ci_status = status), s)
  })
  order <- summary <- covariance <- identities <- list()
  for (window in names(windows)) for (band in c("All PHQ scores", "0-4", "5-9", ">=10")) {
    mask <- windows[[window]] & if (band == "All PHQ scores") TRUE else v$pd_phq_band %in% band
    for (sign in c(-1, 1)) order[[length(order) + 1L]] <- cbind(
      data.frame(window = window, PHQ_band = band, beverage_positive_day = if (sign == -1) "Day 1 (MEC)" else "Day 2 (telephone)"),
      pp_proportion(design, mask, v$pd_order_sign %in% sign))
    for (state in c("beverage_positive", "beverage_zero", "oriented_difference")) {
      working <- design
      for (component in names(pp_components())) {
        a <- v[[paste0("pd_day1_", component)]]; b <- v[[paste0("pd_day2_", component)]]
        value <- if (state == "oriented_difference") v[[paste0("pd_delta_", component)]] else if
          (state == "beverage_positive") ifelse(v$pd_order_sign == 1, b, a) else ifelse(v$pd_order_sign == 1, a, b)
        working$variables[[paste0(".pp_", component)]] <- value
      }
      variables <- paste0(".pp_", names(pp_components()))
      joint <- pp_joint_mean(working, mask, variables)
      id <- data.frame(window = window, PHQ_band = band, state = state)
      table <- joint$summary; table$component <- unname(pp_components())
      summary[[length(summary) + 1L]] <- cbind(id, table)
      cov <- expand.grid(component_row = unname(pp_components()), component_column = unname(pp_components()), stringsAsFactors = FALSE)
      cov$covariance <- as.vector(joint$covariance)
      covariance[[length(covariance) + 1L]] <- cbind(id, cov)
      mean_difference <- sum(table$estimate[1:4]) - table$estimate[5]
      variance_difference <- sum(joint$covariance[1:4, 1:4, drop = FALSE]) - joint$covariance[5, 5]
      if (table$n[1] > 0L) pp_assert(is.finite(mean_difference) && abs(mean_difference) < 1e-6, "Weighted component means fail the total-energy identity.")
      if (table$ci_status[1] == "reported") pp_assert(is.finite(variance_difference) &&
        abs(variance_difference) < 1e-6 * max(1, abs(joint$covariance[5, 5])), "Joint component covariance fails the total-energy identity.")
      identities[[length(identities) + 1L]] <- cbind(id, data.frame(mean_sum_minus_total = mean_difference,
        component_covariance_sum_minus_total_variance = variance_difference, ci_status = table$ci_status[1]))
    }
  }
  gap <- v$pd_gap_days
  gap_groups <- list(days3to10 = is.finite(gap) & gap >= 3 & gap <= 10,
    above10 = is.finite(gap) & gap > 10, below3 = is.finite(gap) & gap < 3, unavailable = !is.finite(gap))
  timing <- lapply(names(gap_groups), function(name) cbind(data.frame(interval = name), pp_proportion(design, v$.pp_all, gap_groups[[name]])))
  gap_mean <- pp_joint_mean(design, v$.pp_all & is.finite(gap), "pd_gap_days")$summary
  list(flow = do.call(rbind, flow), paired_order_PHQ_support = do.call(rbind, order),
    timing_distribution = do.call(rbind, timing), timing_mean = gap_mean,
    weighted_energy_decomposition = do.call(rbind, summary), energy_joint_covariance = do.call(rbind, covariance),
    total_energy_identity = do.call(rbind, identities), metadata = list(
      weights = "Person-level WTDR2D/4; full positive-weight design parent retained before domain restriction",
      primary_domain = "Eligible beverage-discordant people with model-valid age, sex, cycle, weekdays and observed intake gap 3-10 days",
      broad_domain = "The same model-valid beverage-discordant people at all intervals, including unavailable gaps",
      units = "Energy kcal; one participant per paired record; ethanol energy 6.93 kcal/g",
      difference = "Beverage-positive day minus exactly-zero coded-beverage ethanol day",
      remainder = "Other reported food/beverage energy excludes coded positive-ethanol beverages and identified linked additions; unlinked mixers remain unresolved",
      calculated_beverage_remainder = "Signed difference between tabulated coded-beverage item energy and 6.93 kcal/g ethanol conversion; negative values are retained and their cause is not assigned",
      uncertainty = "Joint Taylor covariance for all energy components; t 95% intervals with domain PSU-minus-stratum df",
      cell_rule = "Point estimates retained; suppress intervals for n<5, PSUs<2 or domain df<=0; proportions require category/complement n>=5 and PSUs>=2",
      PHQ_bands = "0-4, 5-9 and >=10 describe support only; no model cutoffs or tests",
      order_mode = "Day1 MEC and day2 telephone are inseparable; these tables do not estimate separate mode effects"))
}
