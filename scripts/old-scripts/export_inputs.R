# Source this file to DEFINE export_inputs(); no analysis runs on source().
# Builds adjusted Table 2 objects, including adjusted probabilities, from df_pheno.
# Model and emmeans extraction follow the supplied revised Table 2 code.

export_inputs <- function(df_pheno, file = NULL) {
  required_packages <- c("dplyr", "geepack", "emmeans", "tidyr", "scales")
  absent <- required_packages[!vapply(required_packages, requireNamespace,
                                    logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install required packages: ", paste(absent, collapse = ", "))
  suppressPackageStartupMessages({
    library(dplyr)
    library(geepack)
    library(emmeans)
    library(scales)
  })
  if (!is.data.frame(df_pheno)) stop("df_pheno must be a data frame.")
  vars <- c("phenotype", "extub14", "treat", "GESTAGE", "MALE", "rand_id", "randcent")
  if (!all(vars %in% names(df_pheno))) {
    stop("df_pheno is missing columns: ", paste(setdiff(vars, names(df_pheno)), collapse = ", "))
  }
  pheno_order <- paste("Phenotype", 1:4)
  check_values <- function(x, allowed, label) {
    observed <- unique(as.character(x[!is.na(x)]))
    if (!all(observed %in% allowed)) {
      stop(label, " has unexpected values: ", paste(setdiff(observed, allowed), collapse = ", "))
    }
  }
  check_values(df_pheno$phenotype, pheno_order, "phenotype")
  check_values(df_pheno$treat, c("Placebo", "Hydrocortisone"), "treat")
  check_values(df_pheno$extub14, c("No", "Yes"), "extub14")
  check_values(df_pheno$MALE, c("0", "1"), "MALE")
  if (!is.numeric(df_pheno$GESTAGE)) stop("GESTAGE must be numeric.")
  if (any(!is.finite(df_pheno$GESTAGE[!is.na(df_pheno$GESTAGE)]))) {
    stop("GESTAGE contains infinite values.")
  }
  eligible <- complete.cases(df_pheno[vars])
  if (!any(eligible)) stop("No complete-case infants available for the adjusted GEE.")
  cells <- table(factor(df_pheno$phenotype[eligible], levels = pheno_order),
                 factor(df_pheno$treat[eligible], levels = c("Placebo", "Hydrocortisone")))
  if (any(cells == 0)) stop("Every phenotype must contain infants in both treatment arms.")
  message("Input rows: ", nrow(df_pheno), "; complete-case phenotype rows: ", sum(eligible))
  fig4_df <- df_pheno %>%
    filter(
      phenotype %in% pheno_order,
      !is.na(extub14),
      !is.na(treat),
      !is.na(GESTAGE),
      !is.na(MALE),
      !is.na(rand_id),
      !is.na(randcent)
    ) %>%
    mutate(
      phenotype = factor(phenotype, levels = pheno_order),
      treat = factor(treat, levels = c("Placebo", "Hydrocortisone")),
      male = factor(MALE, levels = c(0, 1), labels = c("Female", "Male")),
      extub14_num = as.numeric(extub14 == "Yes"),
      randcent = factor(randcent),
      rand_id = factor(rand_id)
    ) %>%
    arrange(rand_id)  # geeglm REQUIRES data sorted by cluster id

  cat("Adjusted analytic N:", nrow(fig4_df),
      "| unique families (rand_id):", n_distinct(fig4_df$rand_id),
      "| unique centers:", n_distinct(fig4_df$randcent), "\n")

  ## ------------------------------------------------------------
  ## 2. PRIMARY MODEL: GEE, clustered on rand_id, adjusted for GESTAGE,
  ##    sex, and study center (fixed effect), exchangeable correlation
  ## ------------------------------------------------------------

  fit_extub_gee <- geeglm(
    extub14_num ~ phenotype * treat + GESTAGE + male + randcent,
    id = rand_id,
    data = fig4_df,
    family = binomial,
    corstr = "exchangeable"
  )

  if (!is.null(fit_extub_gee$geese$error) && any(fit_extub_gee$geese$error != 0)) {
    stop("The adjusted interaction GEE did not converge.")
  }
  if (any(!is.finite(coef(fit_extub_gee))) || any(!is.finite(vcov(fit_extub_gee)))) {
    stop("The adjusted interaction GEE has non-finite coefficients or covariance.")
  }

  # Global interaction test (Wald, appropriate for GEE — GEE does not
  # support a likelihood-ratio test since it is not a full-likelihood method)
  anova_gee <- anova(fit_extub_gee)
  print(anova_gee)
  interaction_row <- match("phenotype:treat", trimws(rownames(anova_gee)))
  p_column <- intersect(c("P(>|Chi|)", "Pr(>Chi)", "Pr(>Chisq)"), names(anova_gee))
  if (is.na(interaction_row) || !length(p_column)) {
    stop("Could not identify the interaction test in anova(fit_extub_gee); inspect its output.")
  }
  p_interaction_gee <- as.numeric(anova_gee[interaction_row, p_column[1]])
  if (length(p_interaction_gee) != 1L || !is.finite(p_interaction_gee)) {
    stop("The interaction P-value was not estimable. Inspect the adjusted GEE.")
  }
  cat("Phenotype x Treatment interaction (GEE Wald test) p =",
      round(p_interaction_gee, 3), "\n")


  raw_counts <- fig4_df %>%
    group_by(phenotype, treat) %>%
    summarise(
      n_event = sum(extub14_num),
      N = n(),
      pct = n_event / N,
      .groups = "drop"
    )

  raw_wide <- raw_counts %>%
    select(phenotype, treat, n_event, N, pct) %>%
    tidyr::pivot_wider(
      names_from = treat,
      values_from = c(n_event, N, pct)
    ) %>%
    mutate(
      n_over_N = sprintf(
        "%d/%d vs %d/%d",
        n_event_Placebo, N_Placebo, n_event_Hydrocortisone, N_Hydrocortisone
      ),
      pct_placebo = percent(pct_Placebo, 1),
      pct_hc = percent(pct_Hydrocortisone, 1)
    ) %>%
    select(phenotype, n_over_N, pct_placebo, pct_hc)

  ## ------------------------------------------------------------
  ## 2. ODDS RATIO (native scale of the logit-link GEE model)
  ## ------------------------------------------------------------

  emm_link <- emmeans(fit_extub_gee, ~ treat | phenotype)  # link (logit) scale

  # NEW: display probabilities from the SAME grid used for RR and RD.
  # Back-transform AFTER logit-scale averaging; do not change the estimand.
  adjusted_probabilities <- as.data.frame(summary(emm_link, type = "response"))
  if (!all(c("phenotype", "treat", "prob") %in% names(adjusted_probabilities)) ||
      nrow(adjusted_probabilities) != 8L ||
      any(!is.finite(adjusted_probabilities$prob)) ||
      any(adjusted_probabilities$prob <= 0 | adjusted_probabilities$prob >= 1) ||
      anyDuplicated(adjusted_probabilities[c("phenotype", "treat")])) {
    stop("Expected one estimable adjusted probability per phenotype/treatment cell.")
  }
  adjusted_wide <- adjusted_probabilities %>%
    select(phenotype, treat, prob) %>%
    tidyr::pivot_wider(names_from = treat, values_from = prob) %>%
    mutate(
      phenotype = factor(phenotype, levels = pheno_order),
      adj_pct_placebo = sprintf("%.2f%%", 100 * Placebo),
      adj_pct_hc = sprintf("%.2f%%", 100 * Hydrocortisone)
    ) %>%
    select(phenotype, adj_pct_placebo, adj_pct_hc)
  cat("\nAdjusted treatment-specific probabilities (not observed proportions):\n")
  print(adjusted_probabilities, row.names = FALSE)
  or_contrast <- contrast(emm_link, method = "revpairwise", by = "phenotype", type = "response")
  or_summary <- as.data.frame(summary(or_contrast, infer = c(TRUE, TRUE)))

  lcl_col <- intersect(c("asymp.LCL", "lower.CL"), colnames(or_summary))[1]
  ucl_col <- intersect(c("asymp.UCL", "upper.CL"), colnames(or_summary))[1]

  or_df <- or_summary %>%
    mutate(
      phenotype = factor(phenotype, levels = pheno_order),
      OR_CI = sprintf("%.2f (%.2f-%.2f)", odds.ratio, .data[[lcl_col]], .data[[ucl_col]]),
      p_value = p.value
    ) %>%
    select(phenotype, OR_CI, p_value)

  ## ------------------------------------------------------------
  ## 3. RISK RATIO — regrid to the LOG-probability scale, then take a
  ##    ratio-style contrast. This is the standard emmeans approach for
  ##    getting risk ratios out of a model fit on the logit scale.
  ## ------------------------------------------------------------

  emm_log <- regrid(emm_link, transform = "log")   # now on log(probability) scale
  rr_contrast <- contrast(emm_log, method = "revpairwise", by = "phenotype", type = "response")
  rr_summary <- as.data.frame(summary(rr_contrast, infer = c(TRUE, TRUE)))

  lcl_col_rr <- intersect(c("asymp.LCL", "lower.CL"), colnames(rr_summary))[1]
  ucl_col_rr <- intersect(c("asymp.UCL", "upper.CL"), colnames(rr_summary))[1]
  ratio_col  <- intersect(c("ratio", "response"), colnames(rr_summary))[1]

  rr_df <- rr_summary %>%
    mutate(
      phenotype = factor(phenotype, levels = pheno_order),
      RR_CI = sprintf("%.2f (%.2f-%.2f)", .data[[ratio_col]], .data[[lcl_col_rr]], .data[[ucl_col_rr]])
    ) %>%
    select(phenotype, RR_CI)

  ## ------------------------------------------------------------
  ## 4. ABSOLUTE RISK DIFFERENCE — regrid to the RESPONSE (probability)
  ##    scale, then take a plain difference contrast.
  ## ------------------------------------------------------------

  emm_response <- regrid(emm_link, transform = "response")  # probability scale
  ard_contrast <- contrast(emm_response, method = "revpairwise", by = "phenotype")
  ard_summary <- as.data.frame(summary(ard_contrast, infer = c(TRUE, TRUE)))

  lcl_col_ard <- intersect(c("asymp.LCL", "lower.CL"), colnames(ard_summary))[1]
  ucl_col_ard <- intersect(c("asymp.UCL", "upper.CL"), colnames(ard_summary))[1]

  ard_df <- ard_summary %>%
    mutate(
      phenotype = factor(phenotype, levels = pheno_order),
      ARD_CI = sprintf(
        "%.1f (%.1f to %.1f)",
        estimate * 100, .data[[lcl_col_ard]] * 100, .data[[ucl_col_ard]] * 100
      )
    ) %>%
    select(phenotype, ARD_CI)

  ## ------------------------------------------------------------
  ## 5. GLOBAL INTERACTION P-VALUE (adjusted GEE Wald test — already
  ##    computed as p_interaction_gee above)
  ## ------------------------------------------------------------

  interaction_p_label <- sprintf("Global treatment-by-phenotype interaction P = %.2f",
                                 p_interaction_gee)
  # NOTE: This replaces the unadjusted P = 0.06 previously shown. State
  # explicitly in the caption that this is the covariate-adjusted,
  # family-clustered estimate (see Methods), and reference the unadjusted
  # comparison in the supplement rather than silently swapping the number.

  ## ------------------------------------------------------------
  ## 6. ASSEMBLE FINAL TABLE
  ## ------------------------------------------------------------

  table2_final <- raw_wide %>%
    left_join(adjusted_wide, by = "phenotype") %>%
    left_join(rr_df, by = "phenotype") %>%
    left_join(ard_df, by = "phenotype") %>%
    left_join(or_df, by = "phenotype") %>%
    mutate(
      p_fmt = case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.2f", p_value)
      )
    ) %>%
    select(
      phenotype, n_over_N, pct_placebo, pct_hc,
      adj_pct_placebo, adj_pct_hc, RR_CI, ARD_CI, OR_CI, p_fmt
    )

  print(as.data.frame(table2_final), row.names = FALSE)
  inputs <- list(
    fig4_df = fig4_df,
    fit_extub_gee = fit_extub_gee,
    table2_final = table2_final,
    adjusted_probabilities = adjusted_probabilities,
    p_interaction_gee = p_interaction_gee
  )
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(inputs, file = file)
    message("Saved adjusted inputs: ", normalizePath(file, winslash = "/"))
  }
  return(inputs)
}
