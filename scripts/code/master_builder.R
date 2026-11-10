# Scoped upstream master construction. Reads four public components per release.
# The 22 retained values and row order reproduce the canonical master projection.
mb_fields <- function() c("PERSON_KEY", "SEQN", "RELEASE_ID", "RIDAGEYR", "RIAGENDR", "WTDR2D",
  "SDMVSTRA", "SDMVPSU", "DR1DRSTZ", "DR2DRSTZ", "DR1TKCAL", "DR2TKCAL", "PHQ9_TOTAL",
  sprintf("DPQ%03d", seq(10, 90, 10)))

mb_build <- function(data_root) {
  phq <- sprintf("DPQ%03d", seq(10, 90, 10))
  fields <- list(DEMO = c("SEQN", "RIDAGEYR", "RIAGENDR", "SDMVSTRA", "SDMVPSU"),
    DPQ = c("SEQN", phq), DR1TOT = c("SEQN", "DR1DRSTZ", "DR1TKCAL"),
    DR2TOT = c("SEQN", "WTDR2D", "DR2DRSTZ", "DR2TKCAL"))
  frames <- audit <- list()
  for (cycle in names(pp_cycles())) {
    suffix <- pp_cycles()[[cycle]]
    frame <- NULL
    for (component in names(fields)) {
      path <- file.path(data_root, cycle, paste0(component, "_", suffix, ".xpt"))
      source <- pp_read_xpt(path, fields[[component]])
      pp_assert(pp_id_valid(source$SEQN) && !anyDuplicated(source$SEQN), paste("Invalid master source IDs:", path))
      if (component == "DEMO") {
        frame <- source
        frame$RELEASE_ID <- cycle
      } else {
        # match is an order-preserving one-to-one left join after ID validation.
        index <- match(frame$SEQN, source$SEQN)
        additions <- setdiff(names(source), "SEQN")
        pp_assert(!any(additions %in% names(frame)), "Master source fields collide.")
        for (field in additions) frame[[field]] <- source[[field]][index]
        audit[[length(audit) + 1L]] <- data.frame(cycle = cycle, component = component,
          roster_n = nrow(frame), source_n = nrow(source), matched_n = sum(!is.na(index)),
          unmatched_n = sum(is.na(index)), source_outside_roster_n = sum(!source$SEQN %in% frame$SEQN),
          cardinality = "one_to_one_left_join", stringsAsFactors = FALSE)
      }
    }
    values <- as.matrix(frame[phq])
    valid <- is.finite(values) & values %in% 0:3
    strict <- rowSums(valid) == 9L
    frame$PHQ9_TOTAL <- rowSums(values)
    frame$PHQ9_TOTAL[!strict] <- NA_real_
    frame$PERSON_KEY <- pp_key(cycle, frame$SEQN)
    frames[[cycle]] <- frame[mb_fields()]
  }
  master <- do.call(rbind, frames)
  rownames(master) <- NULL
  pp_initialize(master, canonical = TRUE)
  list(data = master, join_audit = do.call(rbind, audit),
    metadata = list(rows = nrow(master), columns = names(master),
      raw_reader = paste0("foreign::read.xport ", utils::packageVersion("foreign")),
      roster = "DEMO source order within G/H/I/J releases; all 39156 rows retained",
      PHQ = "Sum of nine DPQ items only when all nine are exactly 0, 1, 2 or 3"))
}
