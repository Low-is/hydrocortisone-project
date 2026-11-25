# Sensitivity analyses for the hydrocortisone phenotyping rebuttal
# Run AFTER the primary pipeline has created df_pheno.
# Required packages: install.packages(c("sandwich", "mclust"))
# No results from the rebuttal are hard-coded or assumed.
# IMPORTANT: This specifies UNADJUSTED modified-Poisson RRs with HC1
# sandwich variance clustered by randomization unit. If the final primary
# analysis used a different adjustment/variance method, harmonize it here.
# Confirm that randcent + rand_id identify a randomization unit, including
# siblings, from your codebook. A combined key is ONE clustering dimension,
# not two-way clustering by center and randomization ID.

run_apgar_sensitivity <- function(df_pheno,
                                 output_dir = "apgar_sensitivity_results",
                                 unit_columns = c("randcent", "rand_id")) {
  for (pkg in c("sandwich", "mclust")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Install required package: ", pkg)
  }
  features <- c("GESTAGE", "BIRTHWTX", "APGAR5", "age_rndaysX", "map_fio2X")
  needed <- c("infant_row_id", "phenotype", "treatmnt", "extube14",
              features, unit_columns)
  missing_cols <- setdiff(needed, names(df_pheno))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  if (anyNA(df_pheno$infant_row_id) || anyDuplicated(df_pheno$infant_row_id))
    stop("infant_row_id must be nonmissing and unique per infant.")
  if (!all(vapply(df_pheno[features], is.numeric, logical(1))))
    stop("All clustering features must be numeric; verify codebook recoding.")
  valid_treatment <- !is.na(df_pheno$treatmnt) & df_pheno$treatmnt %in% c(1, 2)
  complete_features <- complete.cases(df_pheno[features])
  d <- as.data.frame(df_pheno[valid_treatment & complete_features, , drop = FALSE])
  if (nrow(d) < 8L) stop("Too few complete-case infants for this analysis.")
  if (any(!is.finite(as.matrix(d[features])))) stop("Nonfinite clustering values.")
  if (any(vapply(d[features], sd, numeric(1)) == 0)) stop("Constant clustering variable.")
  d$original <- as.character(d$phenotype)
  original_levels <- paste("Phenotype", 1:4)
  if (anyNA(d$original) || !setequal(unique(d$original), original_levels))
    stop("Expected original labels: Phenotype 1 through Phenotype 4, without missing assignments.")
  d$original <- factor(d$original, levels = original_levels)
  if (anyNA(d[unit_columns])) stop("Missing randomization-unit identifiers.")
  d$unit <- do.call(interaction, c(lapply(d[unit_columns], as.character),
                                  list(drop = TRUE, lex.order = TRUE)))
  d$hc <- as.integer(d$treatmnt == 2)
  if (!all(unique(na.omit(d$extube14)) %in% c(0, 1)))
    stop("extube14 must be raw numeric 0/1 (not No/Yes factor).")
  if (!is.numeric(d$extube14)) stop("extube14 must be numeric 0/1.")
  d$y <- d$extube14
  unit_treatments <- tapply(d$hc, d$unit, function(x) length(unique(x)))
  if (any(unit_treatments != 1L))
    stop("Conflicting treatments within a unit: check unit_columns/codebook.")

  # Verify that current features/rows reproduce the original partition.
  # ARI is invariant to arbitrary numerical cluster labels.
  cluster_four <- function(x) {
    cutree(hclust(dist(scale(as.matrix(x))), method = "ward.D2"), k = 4)
  }
  reproduced <- cluster_four(d[features])
  reproduction_ari <- mclust::adjustedRandIndex(d$original, reproduced)
  if (abs(reproduction_ari - 1) > 1e-10)
    stop("Five-variable clustering does not reproduce original assignments (ARI = ",
         signif(reproduction_ari, 4), "). Resolve row alignment/preprocessing first.")

  # Sensitivity 1: preserve original assignments and merge only groups 2/3.
  d$merged <- factor(ifelse(d$original %in% c("Phenotype 2", "Phenotype 3"),
                            "Phenotypes 2+3", as.character(d$original)),
                     levels = c("Phenotype 1", "Phenotypes 2+3", "Phenotype 4"))

  # Sensitivity 2: same infants, scaling, distance, linkage, and k; omit Apgar.
  no_apgar_features <- setdiff(features, "APGAR5")
  d$no_apgar <- factor(cluster_four(d[no_apgar_features]), levels = 1:4,
                       labels = paste("No-Apgar cluster", 1:4))
  ari <- mclust::adjustedRandIndex(d$original, d$no_apgar)
  # New cluster numbers have no inherent correspondence with original labels.
  membership <- as.data.frame.matrix(table(d$original, d$no_apgar))
  membership <- data.frame(original_phenotype = rownames(membership),
                           membership, row.names = NULL, check.names = FALSE)

  # Missing outcomes are excluded only from treatment estimation, not clustering.
  a <- d[!is.na(d$y), , drop = FALSE]
  for (gvar in c("original", "merged", "no_apgar")) {
    counts <- table(a[[gvar]], a$hc)
    if (ncol(counts) != 2L || any(counts == 0L))
      stop("Missing treatment arm in at least one ", gvar, " group.")
  }
  # Fail explicitly for boundary cells rather than silently emitting unstable RRs.
  check_events <- function(dat, group) {
    z <- split(dat$y, interaction(dat[[group]], dat$hc, drop = TRUE))
    if (any(vapply(z, function(x) sum(x) == 0L || sum(x) == length(x), logical(1))))
      stop("Zero events or zero nonevents in a treatment cell for ", group,
           "; prespecify a sparse-data method before interpreting Wald inference.")
  }
  lapply(c("original", "merged", "no_apgar"), function(g) check_events(a, g))

  fit_model <- function(formula, dat) {
    fit <- glm(formula, family = poisson(link = "log"), data = dat,
               na.action = na.fail)
    if (!fit$converged || any(!is.finite(coef(fit)))) stop("Model failed to converge.")
    if (nlevels(droplevels(dat$unit)) < 2L) stop("Insufficient independent units.")
    v <- sandwich::vcovCL(fit, cluster = dat$unit, type = "HC1", cadjust = TRUE)
    if (any(!is.finite(v))) stop("Invalid sandwich covariance matrix.")
    list(fit = fit, V = v)
  }
  contrast_rr <- function(model, L, label) {
    b <- coef(model$fit)
    L <- as.numeric(L)
    log_rr <- sum(L * b)
    variance <- as.numeric(t(L) %*% model$V %*% L)
    if (!is.finite(variance) || variance <= 0) stop("Nonpositive contrast variance.")
    se <- sqrt(variance)
    data.frame(group = label, RR = exp(log_rr),
               lower_95 = exp(log_rr - qnorm(.975) * se),
               upper_95 = exp(log_rr + qnorm(.975) * se),
               p_value = 2 * pnorm(-abs(log_rr / se)))
  }
  grouped_analysis <- function(dat, group_column) {
    dat$group <- droplevels(dat[[group_column]])
    model <- fit_model(y ~ group * hc, dat)
    lev <- levels(dat$group)
    contrasts <- lapply(lev, function(g) {
      nd <- data.frame(group = factor(c(g, g), levels = lev), hc = c(0, 1))
      X <- model.matrix(delete.response(terms(model$fit)), nd,
                        contrasts.arg = model$fit$contrasts)
      X[2, ] - X[1, ]
    })
    names(contrasts) <- lev
    effects <- do.call(rbind, lapply(lev, function(g)
      contrast_rr(model, contrasts[[g]], g)))
    counts <- do.call(rbind, lapply(lev, function(g) {
      z <- dat[dat$group == g, , drop = FALSE]
      data.frame(group = g,
                 placebo_n = sum(z$hc == 0),
                 placebo_events = sum(z$y[z$hc == 0]),
                 hydrocortisone_n = sum(z$hc == 1),
                 hydrocortisone_events = sum(z$y[z$hc == 1]))
    }))
    effects <- merge(counts, effects, by = "group", sort = FALSE)
    effects <- effects[match(lev, effects$group), ]
    effects$placebo_risk <- effects$placebo_events / effects$placebo_n
    effects$hydrocortisone_risk <- effects$hydrocortisone_events / effects$hydrocortisone_n
    # Joint test: all treatment-by-group terms equal zero on the log-RR scale.
    ix <- grep(":", names(coef(model$fit)), fixed = TRUE)
    b <- coef(model$fit)[ix]
    vv <- model$V[ix, ix, drop = FALSE]
    if (qr(vv)$rank < length(ix)) stop("Singular covariance for interaction test.")
    wald <- as.numeric(t(b) %*% solve(vv, b))
    interaction <- data.frame(analysis = group_column, chi_square = wald,
                              df = length(ix), p_value = pchisq(wald, length(ix), lower.tail = FALSE))
    list(effects = effects, interaction = interaction, model = model,
         contrasts = contrasts)
  }
  primary <- grouped_analysis(a, "original")
  merged <- grouped_analysis(a, "merged")
  no_apgar <- grouped_analysis(a, "no_apgar")
  pooled_data <- a[a$merged == "Phenotypes 2+3", , drop = FALSE]
  pooled_model <- fit_model(y ~ hc, pooled_data)
  pooled_L <- as.numeric(names(coef(pooled_model$fit)) == "hc")
  pooled <- contrast_rr(pooled_model, pooled_L, "Phenotypes 2+3 pooled")
  # Ratio of treatment RRs: original Phenotype 2 versus original Phenotype 3.
  compare_23 <- contrast_rr(primary$model,
                            primary$contrasts[["Phenotype 2"]] -
                              primary$contrasts[["Phenotype 3"]],
                            "Treatment RR in Phenotype 2 / treatment RR in Phenotype 3")
  names(compare_23)[names(compare_23) == "RR"] <- "ratio_of_RRs"

  # Baseline profiles for interpreting the new partition.
  profiles <- do.call(rbind, lapply(levels(d$no_apgar), function(g) {
    do.call(rbind, lapply(features, function(v) {
      x <- d[d$no_apgar == g, v]
      data.frame(cluster = g, variable = v, n = length(x),
                 median = median(x), q1 = unname(quantile(x, .25)),
                 q3 = unname(quantile(x, .75)))
    }))
  }))
  flow <- data.frame(
    measure = c("Input infants", "Randomized infants", "Same-infant clustering cohort",
                "Missing primary outcome in clustering cohort", "Treatment-analysis infants",
                "Treatment-analysis randomization units", "Five-variable reproduction ARI",
                "Original versus no-Apgar ARI"),
    value = c(nrow(df_pheno), sum(valid_treatment), nrow(d), sum(is.na(d$y)),
              nrow(a), nlevels(droplevels(a$unit)), reproduction_ari, ari))
  tables <- list(
    cohort_and_ARI = flow,
    primary_RRs = primary$effects,
    merged_phenotypes_RRs = merged$effects,
    pooled_phenotypes_2_3_RR = pooled,
    no_apgar_RRs = no_apgar$effects,
    interaction_tests = rbind(primary$interaction, merged$interaction, no_apgar$interaction),
    original_phenotype_2_vs_3_interaction = compare_23,
    cluster_cross_tabulation = membership,
    no_apgar_baseline_profiles = profiles)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(tables)) write.csv(tables[[nm]],
                                     file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
  # Avoid exporting infant IDs or raw participant-level data in reviewer tables.
  saveRDS(list(primary = primary$model$fit, merged = merged$model$fit,
               no_apgar = no_apgar$model$fit, pooled = pooled_model$fit),
          file.path(output_dir, "models_local_review_only.rds"))
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  writeLines(c(
    "Exploratory sensitivity analyses prompted by peer review.",
    "Clustering: same original five-variable complete-case cohort in both solutions;",
    "z-score scaling; Euclidean distance; Ward.D2; k=4; Apgar omitted only in alternative.",
    "Outcome: existing raw extube14 variable, coding 1=successful extubation.",
    "Treatment: treatmnt 1=placebo, 2=hydrocortisone, per supplied pipeline.",
    "RRs: unadjusted modified Poisson, HC1 sandwich covariance clustered by:",
    paste(unit_columns, collapse = " + "),
    "Confidence intervals use normal approximation; interaction uses joint robust Wald chi-square.",
    "Interaction is tested on the multiplicative RR scale; tests are not multiplicity-adjusted.",
    "The pooled significant RR alone does not establish treatment heterogeneity or distinct phenotypes.",
    "ARI quantifies partition agreement; it does not validate a biological phenotype.",
    "Do not assert diffuse effects, robustness, or a unique responder solely from significance patterns.",
    "Confirm the variance method and adjustment set against the final primary statistical analysis.",
    "RDS model objects may contain participant data; keep local, do not submit as reviewer tables."
  ), file.path(output_dir, "analysis_notes.txt"))
  print(flow, row.names = FALSE)
  cat("\nPooled Phenotypes 2+3 treatment effect:\n")
  print(pooled, row.names = FALSE)
  cat("\nNo-Apgar cluster treatment effects:\n")
  print(no_apgar$effects, row.names = FALSE)
  cat("\nJoint treatment-by-group interaction tests:\n")
  print(tables$interaction_tests, row.names = FALSE)
  cat("\nOriginal Phenotype 2 versus 3 interaction:\n")
  print(compare_23, row.names = FALSE)
  cat("\nSaved tables in: ", normalizePath(output_dir), "\n", sep = "")
  invisible(list(tables = tables, assignments = d[c("infant_row_id", "original", "merged", "no_apgar")]))
}

# Run after sourcing this file (uncomment):
# sensitivity_results <- run_apgar_sensitivity(df_pheno)
# If rand_id is globally unique by codebook, you may instead specify:
# sensitivity_results <- run_apgar_sensitivity(df_pheno, unit_columns = "rand_id")
