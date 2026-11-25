# Source this file to DEFINE run_table2(); no analysis runs on source().
# inputs is the list returned by export_inputs(), or its saved RDS path.

run_table2 <- function(inputs, output_dir = "table2_output", install_packages = FALSE) {
  packages <- c("geepack", "emmeans", "dplyr", "flextable", "officer")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && isTRUE(install_packages)) {
    # Use a writable personal library when the system library is protected.
    personal_lib <- path.expand(Sys.getenv("R_LIBS_USER"))
    if (!nzchar(personal_lib)) personal_lib <- file.path(path.expand("~"), "R", "table2-library")
    personal_lib <- strsplit(personal_lib, .Platform$path.sep, fixed = TRUE)[[1]][1]
    dir.create(personal_lib, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(personal_lib, .libPaths()))
    install.packages(missing, lib = personal_lib, repos = "https://cloud.r-project.org")
  }
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". Install them or call run_table2(..., install_packages = TRUE).")
  suppressPackageStartupMessages(library(dplyr))

  bundle <- if (is.character(inputs) && length(inputs) == 1L) readRDS(inputs) else inputs
  needed <- c("fig4_df", "fit_extub_gee", "table2_final", "p_interaction_gee")
  if (!is.list(bundle) || !all(needed %in% names(bundle))) {
    stop("Input must be the list saved by export_inputs.R; do not supply the raw dataset alone.")
  }
  dat <- as.data.frame(bundle$fig4_df)
  original_fit <- bundle$fit_extub_gee
  original_table <- as.data.frame(bundle$table2_final)
  interaction_p <- bundle$p_interaction_gee
  if (!inherits(original_fit, "geeglm")) stop("fit_extub_gee must be a fitted geeglm object.")
  if (length(interaction_p) != 1L || !is.numeric(interaction_p) ||
      !is.finite(interaction_p) || interaction_p < 0 || interaction_p > 1) {
    stop("p_interaction_gee must be a single numeric probability, not a rounded label.")
  }
  vars <- c("phenotype", "treat", "extub14_num", "GESTAGE", "male", "randcent", "rand_id")
  if (!all(vars %in% names(dat))) stop("fig4_df lacks required columns: ",
                                      paste(setdiff(vars, names(dat)), collapse = ", "))
  if (anyNA(dat[vars])) stop("fig4_df contains missing model inputs; export the actual complete-case cohort.")
  if (!all(dat$extub14_num %in% c(0, 1))) stop("extub14_num must be coded 0/1.")
  if (!identical(levels(dat$treat), c("Placebo", "Hydrocortisone"))) {
    stop("treat must be a factor with levels Placebo, Hydrocortisone in that order.")
  }
  expected_pheno <- paste("Phenotype", 1:4)
  if (!setequal(as.character(unique(dat$phenotype)), expected_pheno)) {
    stop("Expected exactly Phenotype 1 through Phenotype 4 in fig4_df.")
  }
  table_cols <- c("phenotype", "n_over_N", "pct_placebo", "pct_hc",
                  "adj_pct_placebo", "adj_pct_hc", "RR_CI", "ARD_CI", "OR_CI", "p_fmt")
  if (!all(table_cols %in% names(original_table)) || nrow(original_table) != 4L ||
      !setequal(as.character(original_table$phenotype), expected_pheno)) {
    stop("table2_final must include the four phenotype rows and adjusted probability columns. Rerun the workflow with the updated export_inputs.R to refresh the saved inputs.")
  }
  if (nobs(original_fit) != nrow(dat)) stop("Saved GEE and fig4_df have different sample sizes.")
  mf <- model.frame(original_fit)
  check_vars <- c("phenotype", "treat", "GESTAGE", "male", "randcent", "extub14_num")
  if (!all(check_vars %in% names(mf)) ||
      !all(vapply(check_vars, function(v) {
        isTRUE(all.equal(unname(mf[[v]]), unname(dat[[v]]), check.attributes = FALSE))
      }, logical(1)))) {
    stop("Saved GEE does not match fig4_df row order or model variables. Export matching objects.")
  }
  if (original_fit$family$family != "binomial" || original_fit$family$link != "logit" ||
      !any(grepl("phenotype:treat|treat:phenotype", attr(terms(original_fit), "term.labels")))) {
    stop("Expected the original binomial-logit GEE containing phenotype-by-treatment interaction.")
  }

  # Keep the saved model, table rows, and interaction P intact.
  dat <- dat[order(dat$rand_id), , drop = FALSE]
  fit_all <- geepack::geeglm(
    extub14_num ~ treat + GESTAGE + male + randcent,
    id = rand_id, data = dat,
    family = binomial(link = "logit"), corstr = "exchangeable",
    std.err = "san.se"
  )
  if (!is.null(fit_all$geese$error) && any(fit_all$geese$error != 0)) {
    stop("Pooled GEE did not converge; inspect the data/model before reporting results.")
  }
  if (any(!is.finite(coef(fit_all))) || any(!is.finite(vcov(fit_all)))) {
    stop("Pooled GEE returned non-finite coefficients/covariance; inspect the model.")
  }

  emm <- emmeans::emmeans(fit_all, ~ treat, vcov. = vcov(fit_all))
  # NEW: probabilities from the same overall grid used for its RR and RD.
  overall_probabilities <- as.data.frame(summary(emm, type = "response"))
  if (!all(c("treat", "prob") %in% names(overall_probabilities)) ||
      nrow(overall_probabilities) != 2L ||
      !setequal(as.character(overall_probabilities$treat), c("Placebo", "Hydrocortisone")) ||
      anyDuplicated(overall_probabilities$treat) ||
      any(!is.finite(overall_probabilities$prob)) ||
      any(overall_probabilities$prob <= 0 | overall_probabilities$prob >= 1)) {
    stop("Expected two estimable overall treatment probabilities.")
  }
  overall_pbo_prob <- overall_probabilities$prob[
    as.character(overall_probabilities$treat) == "Placebo"]
  overall_hc_prob <- overall_probabilities$prob[
    as.character(overall_probabilities$treat) == "Hydrocortisone"]
  contrast_def <- list("HC vs Placebo" = c(-1, 1))
  summarize_contrast <- function(grid, response = FALSE) {
    as.data.frame(summary(
      emmeans::contrast(grid, method = contrast_def, adjust = "none"),
      type = if (response) "response" else "link",
      infer = c(TRUE, TRUE), level = 0.95, adjust = "none"
    ))
  }
  or <- summarize_contrast(emm, TRUE)
  rr <- summarize_contrast(emmeans::regrid(emm, transform = "log"), TRUE)
  rd <- summarize_contrast(emmeans::regrid(emm, transform = "response"))
  pick <- function(x, choices) {
    hit <- intersect(choices, names(x))
    if (!length(hit)) stop("Unexpected emmeans columns: ", paste(names(x), collapse = ", "))
    x[[hit[1]]]
  }
  ci_text <- function(x, est_names, mult = 1, digits = 2, suffix = "") {
    values <- mult * c(pick(x, est_names), pick(x, c("asymp.LCL", "lower.CL")),
                      pick(x, c("asymp.UCL", "upper.CL")))
    if (length(values) != 3L || any(!is.finite(values))) stop("Non-estimable pooled treatment contrast.")
    txt <- formatC(values, digits = digits, format = "f")
    paste0(txt[1], suffix, " (", txt[2], " to ", txt[3], ")")
  }
  counts <- dat %>% group_by(treat) %>%
    summarise(events = sum(extub14_num), total = n(), .groups = "drop")
  pbo <- filter(counts, treat == "Placebo")
  hc <- filter(counts, treat == "Hydrocortisone")
  overall <- data.frame(
    phenotype = "All infants",
    n_over_N = sprintf("%d/%d vs %d/%d", pbo$events, pbo$total, hc$events, hc$total),
    pct_placebo = sprintf("%.0f%%", 100 * pbo$events / pbo$total),
    pct_hc = sprintf("%.0f%%", 100 * hc$events / hc$total),
    adj_pct_placebo = sprintf("%.2f%%", 100 * overall_pbo_prob),
    adj_pct_hc = sprintf("%.2f%%", 100 * overall_hc_prob),
    RR_CI = ci_text(rr, c("ratio", "response")),
    ARD_CI = ci_text(rd, "estimate", 100, 1),
    OR_CI = ci_text(or, "odds.ratio"),
    p_fmt = ifelse(or$p.value < 0.001, "<0.001", sprintf("%.2f", or$p.value))
  )
  original_table$phenotype <- as.character(original_table$phenotype)
  final <- bind_rows(overall, original_table[table_cols])
  stopifnot(nrow(final) == 5L)
  # Check saved phenotype values without depending on row-name attributes.
  stopifnot(all(vapply(table_cols, function(v) {
    identical(as.character(final[[v]][-1]), as.character(original_table[[v]]))
  }, logical(1))))

  output <- output_dir
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output)) stop("Could not create output directory: ", output)
  interaction_label <- sprintf("Global treatment-by-phenotype interaction P = %.2f", interaction_p)
  footnote <- paste0(
    "Counts and observed percentages are unadjusted. Adjusted percentages are model-estimated ",
    "extubation probabilities from binomial-logit generalized estimating equations adjusted for ",
    "gestational age, sex, and study center, with exchangeable correlation and robust standard ",
    "errors clustered by randomization unit. The emmeans reference grid evaluates gestational ",
    "age at its analytic-cohort mean and averages equally over sex and study-center levels on ",
    "the logit scale before back-transformation, rather than using individual-level empirical ",
    "standardization. Adjusted risk ratios and risk differences are calculated by dividing and ",
    "subtracting, respectively, the adjusted treatment-specific probabilities, not the observed ",
    "percentages. Effects are calculated before rounding; risk differences are expressed in ",
    "percentage points. All infants represents the same phenotype-assigned analytic cohort ",
    "(n = ", nrow(dat), "), not the full parent-trial population, and uses a separate adjusted ",
    "model without phenotype terms. Phenotype-specific estimates and the global interaction ",
    "P-value come from the adjusted model containing the treatment-by-phenotype interaction. ",
    "Row-specific P-values test treatment odds ratios; the global interaction test is on the ",
    "log-odds scale. HC, hydrocortisone; RR, risk ratio; OR, odds ratio; CI, confidence interval."
  )
  ft <- flextable::flextable(final)
  ft <- flextable::set_header_labels(ft,
    phenotype = "Phenotype", n_over_N = "n/N (Placebo vs HC)",
    pct_placebo = "Observed\nPlacebo %", pct_hc = "Observed\nHC %",
    adj_pct_placebo = "Adjusted\nPlacebo %", adj_pct_hc = "Adjusted\nHC %",
    RR_CI = "Adjusted RR\n(95% CI)",
    ARD_CI = "Adjusted risk difference\nPercentage points (95% CI)",
    OR_CI = "Adjusted OR\n(95% CI)", p_fmt = "P-value")
  ft <- flextable::add_header_lines(ft, values = interaction_label)
  ft <- flextable::add_header_lines(ft,
    values = "Table 2. Hydrocortisone vs placebo by baseline phenotype.")
  ft <- flextable::bold(ft, i = 1, part = "body")
  ft <- flextable::fontsize(ft, size = 8, part = "all")
  ft <- flextable::add_footer_lines(ft, values = footnote)
  ft <- flextable::width(ft, width = c(0.78, 1.25, 0.65, 0.60, 0.65, 0.60, 1.16, 1.67, 1.16, 0.55))
  ft <- flextable::set_table_properties(ft, layout = "fixed")
  section <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = 0.5, bottom = 0.5, left = 0.5, right = 0.5))
  flextable::save_as_docx(ft, path = file.path(output, "Table_2_with_All_Infants.docx"), pr_section = section)
  write.csv(final, file.path(output, "Table_2_with_All_Infants.csv"), row.names = FALSE, na = "")
  # Export numeric adjusted probabilities without display rounding, including CIs.
  if (is.null(bundle$adjusted_probabilities)) {
    phenotype_probabilities <- as.data.frame(summary(
      emmeans::emmeans(original_fit, ~ treat | phenotype), type = "response"))
  } else {
    phenotype_probabilities <- as.data.frame(bundle$adjusted_probabilities)
  }
  overall_probabilities$phenotype <- "All infants"
  phenotype_probabilities$phenotype <- as.character(phenotype_probabilities$phenotype)
  all_probabilities <- dplyr::bind_rows(overall_probabilities, phenotype_probabilities)
  write.csv(all_probabilities, file.path(output, "Table_2_adjusted_probabilities.csv"),
            row.names = FALSE, na = "")
  saveRDS(list(table = final, overall_model = fit_all, overall_OR = or,
               overall_RR = rr, overall_RD = rd, counts = counts,
               adjusted_probabilities = all_probabilities,
               original_interaction_p = interaction_p),
          file.path(output, "Table_2_results.rds"))
  writeLines(c(interaction_label, paste("Exact saved interaction P:", interaction_p),
               "", footnote), file.path(output, "analysis_notes.txt"))
  capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  print(final, row.names = FALSE)
  message(interaction_label, " (retained from input; not recomputed or forced to 0.12)")
  message("Outputs saved in: ", normalizePath(output, winslash = "/"))

  return(invisible(list(
    table = final,
    interaction_p = interaction_p,
    interaction_model = original_fit,
    overall_model = fit_all,
    output_dir = normalizePath(output, winslash = "/")
  )))
}

