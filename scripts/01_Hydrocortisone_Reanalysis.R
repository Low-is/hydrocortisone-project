#####
#libraries-----
#####
# HydrocortisoneAJP2: complete, clean-session reanalysis.
# Research analysis only. Phenotype membership is not a clinical treatment rule.
# Read README.md before running. No pre-existing workspace objects are used.
# R >= 4.3 recommended. Run the ENTIRE file in a fresh R session.

required_packages <- c(
  "dplyr", "tidyr", "readr", "tibble", "purrr", "ggplot2", "scales",
  "patchwork", "cluster", "fpc", "mclust", "sandwich", "lmtest",
  "emmeans", "lme4", "lmerTest", "ordinal", "logistf", "ranger",
  "rpart", "rpart.plot", "pROC", "flextable", "officer", "clubSandwich",
  "brglm2"
)
INSTALL_MISSING_PACKAGES <- TRUE
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  if (!INSTALL_MISSING_PACKAGES) stop("Install: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package unavailable: ", p)
}))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(purrr); library(ggplot2); library(patchwork)
})
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE,
        contrasts = c("contr.treatment", "contr.poly"))
stopifnot(getRversion() >= "4.3.0", packageVersion("dplyr") >= "1.1.1")
set.seed(2024)

#####
#configuration and folders-----
#####
# Put the input files together and edit INPUT_DIR if necessary.
INPUT_DIR <- Sys.getenv("HYDRO_INPUT_DIR", unset = file.path(getwd(), "Data_Files_CSV"))
OUTPUT_ROOT <- file.path(INPUT_DIR, "Hydrocortisone_Reanalysis_Outputs")
RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT <- file.path(OUTPUT_ROOT, RUN_ID)
DIR <- as.list(setNames(file.path(OUT, c("figures", "tables", "models", "data", "logs")),
                        c("figures", "tables", "models", "data", "logs")))
invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(file.path(DIR$figures, "supplement"), showWarnings = FALSE)

ANALYSIS_FILE <- file.path(INPUT_DIR, "analysis_dash.csv")
LONGITUDINAL_FILE <- file.path(INPUT_DIR, "hc04_dash.csv")
PREVIOUS_DATASET <- file.path(INPUT_DIR, "df_pheno.rds")
# Optional genuine family linkage: randcent, rand_id, family_id.
# A family_id must link siblings, not reuse infant numbers from different centers.
FAMILY_FILE <- file.path(INPUT_DIR, "family_linkage.csv")
# Optional externally verified NICHD estimates. See README for required metadata.
NICHD_FILE <- file.path(INPUT_DIR, "nichd_verified_risk.csv")
BENCHMARK_FILE <- file.path(INPUT_DIR, "Independent_Primary_Benchmark.csv")
STRICT_ORIGINAL_COHORT <- TRUE
BOOTSTRAP_REPS <- 200L
SUBSAMPLE_REPS <- 100L
CONSENSUS_REPS <- 200L
RF_TREES <- 1000L
CV_FOLDS <- 5L
CONF_LEVEL <- 0.95
Z <- qnorm(1 - (1 - CONF_LEVEL) / 2)
PHENO <- paste("Phenotype", 1:4)
COLS <- setNames(c("#D1495B", "#2878B5", "#36966A", "#8B63A6"), PHENO)
CLUSTER_VARS <- c("gestage", "birthwtx", "apgar5", "age_rndaysx", "map_fio2x")
VAR_LABELS <- c(gestage = "Gestational age", birthwtx = "Birth weight",
                apgar5 = "5-minute Apgar", age_rndaysx = "Age at randomization",
                map_fio2x = "Respiratory index")
module_status <- tibble(module = character(), status = character(), detail = character())
warning_log <- character()

#####
#utility functions and output checks-----
#####
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
fmt_p <- function(p) ifelse(is.na(p), "Not estimable",
                           ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_ci <- function(x, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"), x, lo, hi)
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
record_status <- function(module, status, detail = "") {
  module_status <<- bind_rows(module_status, tibble(module, status, detail))
  write_csv(module_status, file.path(DIR$logs, "module_status.csv"))
}
run_module <- function(name, expr) {
  message("Running: ", name)
  tryCatch(withCallingHandlers({
    value <- force(expr)
    record_status(name, "COMPLETE")
    value
  }, warning = function(w) {
    warning_log <<- c(warning_log, paste(name, conditionMessage(w), sep = ": "))
  }), error = function(e) {
    record_status(name, "FAILED", conditionMessage(e))
    message("MODULE FAILED: ", name, ": ", conditionMessage(e))
    NULL
  })
}
save_model <- function(object, name) saveRDS(object, file.path(DIR$models, paste0(name, ".rds")))
save_table <- function(dat, name, title, note = "", landscape = FALSE) {
  dat <- as.data.frame(dat)
  assert(nrow(dat) > 0, paste("Empty table:", name))
  write_csv(dat, file.path(DIR$tables, paste0(name, ".csv")), na = "")
  ft <- flextable::flextable(dat)
  ft <- flextable::theme_vanilla(ft)
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::bg(ft, bg = "#E8EDF1", part = "header")
  ft <- flextable::color(ft, color = "black", part = "all")
  ft <- flextable::border_inner_h(ft, border = officer::fp_border(color = "#D9D9D9"))
  ft <- flextable::border_inner_v(ft, border = officer::fp_border(color = "#D9D9D9"))
  ft <- flextable::border_outer(ft, border = officer::fp_border(color = "#D9D9D9"))
  ft <- flextable::padding(ft, padding = 5, part = "all")
  ft <- flextable::autofit(ft)
  ft <- flextable::set_table_properties(ft, layout = "autofit", width = 1)
  ft <- flextable::set_caption(ft, caption = title)
  if (nzchar(note)) ft <- flextable::add_footer_lines(ft, values = note)
  doc <- officer::read_docx()
  doc <- flextable::body_add_flextable(doc, ft)
  sect <- officer::prop_section(
    page_size = officer::page_size(orient = if (landscape) "landscape" else "portrait"),
    page_margins = officer::page_mar(top = .65, bottom = .65, left = .65, right = .65))
  doc <- officer::body_set_default_section(doc, sect)
  print(doc, target = file.path(DIR$tables, paste0(name, ".docx")))
  invisible(ft)
}
save_fig <- function(plot, name, width = 10, height = 7, supplement = FALSE) {
  dest <- if (supplement) file.path(DIR$figures, "supplement") else DIR$figures
  ggplot2::ggsave(file.path(dest, paste0(name, ".pdf")), plot,
                 width = width, height = height, device = grDevices::pdf, useDingbats = FALSE)
  ggplot2::ggsave(file.path(dest, paste0(name, ".tiff")), plot,
                 width = width, height = height, dpi = 600, compression = "lzw", bg = "white")
  ggplot2::ggsave(file.path(dest, paste0(name, ".png")), plot,
                 width = width, height = height, dpi = 300, bg = "white")
  saveRDS(plot, file.path(dest, paste0(name, "_plot.rds")))
}
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#EEEEEE", linewidth = .3),
          plot.title = element_text(face = "bold", color = "black"),
          plot.subtitle = element_text(color = "#333333"),
          legend.position = "bottom", strip.background = element_rect(fill = "#F0F0F0"))
}
theme_set(theme_pub())
read_source <- function(path) {
  assert(file.exists(path), paste("Missing input:", path))
  d <- read_csv(path, na = c("", "NA", ".", ".Z", ".z"), show_col_types = FALSE,
                name_repair = "check_unique")
  assert(nrow(readr::problems(d)) == 0, paste("CSV parsing problems in", path))
  names(d) <- tolower(names(d))
  assert(!anyDuplicated(names(d)), paste("Column names collide after lowercasing:", path))
  d
}
make_id <- function(center, infant) paste(center, infant, sep = "_")
assert_unique <- function(d, keys, label) {
  assert(!anyNA(d[keys]), paste(label, "has missing keys"))
  assert(!anyDuplicated(d[keys]), paste(label, "has duplicate keys"))
}
wald_block <- function(b, V, pattern) {
  ii <- grep(pattern, names(b))
  assert(length(ii) > 0, paste("No coefficients match", pattern))
  a <- b[ii]; S <- V[ii, ii, drop = FALSE]
  assert(all(is.finite(a)) && all(is.finite(S)), "Nonfinite Wald input")
  assert(qr(S)$rank == length(ii), "Singular covariance for joint test")
  W <- as.numeric(crossprod(a, solve(S, a)))
  tibble(statistic = W, df = length(ii), p_value = pchisq(W, length(ii), lower.tail = FALSE))
}
check_glm <- function(fit) {
  assert(isTRUE(fit$converged), "GLM did not converge")
  assert(all(is.finite(coef(fit))), "GLM has aliased/nonfinite coefficients")
  assert(max(abs(coef(fit))) < 25, "Extreme GLM coefficients; inspect separation before inference")
  invisible(TRUE)
}
check_mer <- function(fit, label) {
  msgs <- fit@optinfo$conv$lme4$messages
  bad <- msgs[!grepl("boundary.*singular", msgs, ignore.case = TRUE)]
  assert(!length(bad), paste(label, paste(bad, collapse = "; ")))
  assert(all(fit@optinfo$conv$opt == 0), paste(label, "optimizer convergence code is not zero"))
  if (lme4::isSingular(fit)) warning(paste(label, "has a boundary/singular random-effects fit; inspect variance components"))
  invisible(TRUE)
}
robust_vcov <- function(fit, d) {
  assert(nobs(fit) == nrow(d), "Model row count differs from covariance data")
  # Composite infant key is unique at baseline. Only genuine supplied family
  # linkage triggers family-clustered inference.
  sandwich::vcovCL(fit, cluster = d$se_unit, type = "HC1", cadjust = TRUE)
}
coef_table <- function(fit, V, exponentiate = FALSE) {
  b <- coef(fit); s <- sqrt(diag(V)); lo <- b - Z * s; hi <- b + Z * s
  tibble(term = names(b), estimate = if (exponentiate) exp(b) else b,
         lower = if (exponentiate) exp(lo) else lo,
         upper = if (exponentiate) exp(hi) else hi,
         p_value = 2 * pnorm(abs(b / s), lower.tail = FALSE))
}
fit_binary <- function(d, outcome, group = "phenotype", adjusted = TRUE,
                       link = "logit", interaction = TRUE) {
  rhs <- if (interaction) paste0(group, " * hc") else "hc"
  if (adjusted) rhs <- paste(rhs, "+ gestage + male + center")
  form <- as.formula(paste(outcome, "~", rhs))
  keep <- unique(c(all.vars(form), "se_unit", "infant_id"))
  dd <- d[complete.cases(d[, keep]), , drop = FALSE]
  dd <- droplevels(dd)
  assert(nrow(dd) > 20 && length(unique(dd[[outcome]])) == 2,
         paste("Insufficient binary outcome data:", outcome))
  fam <- if (link == "log") poisson(link = "log") else binomial(link = "logit")
  fit <- glm(form, data = dd, family = fam, na.action = na.fail,
             control = glm.control(maxit = 100))
  check_glm(fit)
  list(fit = fit, V = robust_vcov(fit, dd), data = dd, outcome = outcome,
       group = group, adjusted = adjusted, link = link, interaction = interaction)
}
design <- function(fit, newdata) {
  model.matrix(delete.response(terms(fit)), newdata, contrasts.arg = fit$contrasts,
               xlev = fit$xlevels)[, names(coef(fit)), drop = FALSE]
}
# Empirical standardization WITHIN each phenotype. Overall row is standardized
# over all eligible infants from a separately fitted no-phenotype model.
# Delta-method uncertainty conditions on the observed covariate distribution.
binary_contrast <- function(obj, dd, label) {
  fit <- obj$fit; V <- obj$V; b <- coef(fit)
  nd0 <- dd; nd1 <- dd; nd0$hc <- 0; nd1$hc <- 1
  X0 <- design(fit, nd0); X1 <- design(fit, nd1)
  eta0 <- drop(X0 %*% b); eta1 <- drop(X1 %*% b)
  if (obj$link == "logit") {
    p0i <- plogis(eta0); p1i <- plogis(eta1)
    g0 <- colMeans(X0 * as.vector(p0i * (1 - p0i)))
    g1 <- colMeans(X1 * as.vector(p1i * (1 - p1i)))
  } else {
    p0i <- exp(eta0); p1i <- exp(eta1)
    g0 <- colMeans(X0 * as.vector(p0i)); g1 <- colMeans(X1 * as.vector(p1i))
  }
  p0 <- mean(p0i); p1 <- mean(p1i)
  assert(p0 > 0 && p1 > 0 && p0 < 1 && p1 < 1,
         "Standardized risk outside (0,1); inspect model")
  se <- function(g) sqrt(max(0, as.numeric(t(g) %*% V %*% g)))
  p0_ci <- plogis(qlogis(p0) + c(-1, 1) * Z * se(g0) / (p0 * (1 - p0)))
  p1_ci <- plogis(qlogis(p1) + c(-1, 1) * Z * se(g1) / (p1 * (1 - p1)))
  rr <- p1 / p0; grr <- g1 / p1 - g0 / p0; serr <- se(grr)
  rd <- p1 - p0; serd <- se(g1 - g0)
  # The conditional OR is constant within a group in this model. It is not
  # the odds ratio of the standardized risks (noncollapsibility).
  gor <- colMeans(X1 - X0); lor <- sum(gor * b); seor <- se(gor)
  out <- tibble(group = label, N = nrow(dd),
    n_placebo = sum(dd$hc == 0), n_hc = sum(dd$hc == 1),
    events_placebo = sum(dd[[obj$outcome]][dd$hc == 0]),
    events_hc = sum(dd[[obj$outcome]][dd$hc == 1]),
    adjusted_placebo = p0, placebo_lower = p0_ci[1], placebo_upper = p0_ci[2],
    adjusted_hc = p1, hc_lower = p1_ci[1], hc_upper = p1_ci[2],
    RR = rr, RR_lower = exp(log(rr) - Z * serr), RR_upper = exp(log(rr) + Z * serr),
    RR_p = 2 * pnorm(abs(log(rr) / serr), lower.tail = FALSE),
    RD = rd, RD_lower = rd - Z * serd, RD_upper = rd + Z * serd,
    conditional_OR = if (obj$link == "logit") exp(lor) else NA_real_,
    OR_lower = if (obj$link == "logit") exp(lor - Z * seor) else NA_real_,
    OR_upper = if (obj$link == "logit") exp(lor + Z * seor) else NA_real_,
    OR_p = if (obj$link == "logit") 2 * pnorm(abs(lor / seor), lower.tail = FALSE) else NA_real_)
  assert(abs(out$RR - out$adjusted_hc/out$adjusted_placebo) < 1e-10,
         "RR and standardized risks disagree")
  assert(abs(out$RD - (out$adjusted_hc-out$adjusted_placebo)) < 1e-10,
         "RD and standardized risks disagree")
  if (!obj$adjusted) {
    observed_rr <- (out$events_hc/out$n_hc)/(out$events_placebo/out$n_placebo)
    assert(abs(out$RR-observed_rr) < 1e-6, "Saturated unadjusted model failed raw-RR check")
  }
  out
}
all_binary_contrasts <- function(obj) {
  d <- obj$data
  if (!obj$interaction) return(binary_contrast(obj, d, "All infants"))
  lev <- levels(d[[obj$group]])
  bind_rows(lapply(lev, function(g) binary_contrast(obj, d[d[[obj$group]] == g, ], g)))
}
display_binary <- function(res) {
  res %>% transmute(Group = group,
    `Placebo n/N` = paste0(events_placebo, "/", n_placebo),
    `HC n/N` = paste0(events_hc, "/", n_hc),
    `Observed placebo %` = round(100 * events_placebo / n_placebo, 1),
    `Observed HC %` = round(100 * events_hc / n_hc, 1),
    `Adjusted placebo %` = round(100 * adjusted_placebo, 2),
    `Adjusted HC %` = round(100 * adjusted_hc, 2),
    `RR (95% CI)` = fmt_ci(RR, RR_lower, RR_upper),
    `RD percentage points (95% CI)` = fmt_ci(100 * RD, 100 * RD_lower, 100 * RD_upper, 1),
    `Conditional OR (95% CI)` = fmt_ci(conditional_OR, OR_lower, OR_upper),
    `OR P` = fmt_p(OR_p))
}
wilson <- function(events, total) {
  p <- events / total; den <- 1 + Z^2 / total
  center <- (p + Z^2 / (2 * total)) / den
  half <- Z * sqrt(p * (1 - p) / total + Z^2 / (4 * total^2)) / den
  tibble(prob = p, lower = center - half, upper = center + half)
}
observed_risk <- function(d, outcome) {
  a <- d %>% filter(!is.na(.data[[outcome]])) %>% group_by(phenotype) %>%
    summarise(events = sum(.data[[outcome]]), N = n(), .groups = "drop")
  bind_cols(a, wilson(a$events, a$N))
}
forest <- function(res, title, subtitle = NULL, effect = "RR") {
  d <- res %>% mutate(group = factor(group, levels = rev(unique(group))))
  if (effect == "RR") {
    ggplot(d, aes(RR, group)) + geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
      geom_segment(aes(x = RR_lower, xend = RR_upper, yend = group), linewidth = .7) +
      geom_point(size = 2.6, color = "#2878B5") + scale_x_log10() +
      labs(title = title, subtitle = subtitle, x = "Risk ratio (95% CI)", y = NULL)
  } else {
    ggplot(d, aes(100 * RD, group)) + geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      geom_segment(aes(x = 100 * RD_lower, xend = 100 * RD_upper, yend = group), linewidth = .7) +
      geom_point(size = 2.6, color = "#2878B5") +
      labs(title = title, subtitle = subtitle, x = "Risk difference, percentage points (95% CI)", y = NULL)
  }
}

#####
#load source data and validate identifiers-----
#####
raw <- read_source(ANALYSIS_FILE)
needed <- c(CLUSTER_VARS, "randcent", "rand_id", "male", "treatmnt", "extube14",
            "death36", "bpd_p_hc", "surv36_nobpd_p_hc", "bpd_jen36_hc",
            "days_vent36", "ndi3", "ndi_mod_sev", "ndi3cat", "exam5agecalcx",
            "death_fu", "race3x", "hisp", "multbrth", "antester", "obgmage", "sga")
assert(all(needed %in% names(raw)), paste("Missing columns:", paste(setdiff(needed, names(raw)), collapse = ", ")))
raw <- raw %>% mutate(across(all_of(needed), num),
                      source_row = row_number(), infant_id = make_id(randcent, rand_id))
assert_unique(raw, c("randcent", "rand_id"), "Analysis source")
assert(all(raw$treatmnt %in% 1:2), "Unknown/missing randomized treatment")
assert(all(raw$male %in% 0:1), "Unknown/missing sex code")
if (STRICT_ORIGINAL_COHORT) assert(nrow(raw) == 800L, "Expected 800 randomized infants")
write_csv(tibble(input = ANALYSIS_FILE, md5 = unname(tools::md5sum(ANALYSIS_FILE))),
          file.path(DIR$logs, "input_checksum.csv"))
id_audit <- tibble(rows = nrow(raw), unique_rand_id = n_distinct(raw$rand_id),
  unique_infants = n_distinct(raw$infant_id),
  rand_ids_reused_across_centers = sum(table(raw$rand_id) > 1))
write_csv(id_audit, file.path(DIR$logs, "identifier_audit.csv"))

#####
#outcome definitions and genuine family linkage-----
#####
df <- raw %>% mutate(
  hc = as.integer(treatmnt == 2),
  treatment = factor(hc, 0:1, c("Placebo", "Hydrocortisone")),
  center = factor(randcent), male = factor(male, 0:1, c("Female", "Male")),
  extub14 = extube14,
  survival_without_physiological_bpd = surv36_nobpd_p_hc,
  jensen_grade = ordered(bpd_jen36_hc, levels = 0:3,
                        labels = c("No BPD", "Grade 1", "Grade 2", "Grade 3")),
  ventilation_days_36 = days_vent36,
  # NDI3 is retained for comparability; the source category coding is audited below.
  ndi_source_binary = ndi3,
  ndi_categorized_modsev = ndi_mod_sev,
  death_or_grade3 = case_when(death36 == 1 ~ 1,
     death36 == 0 & !is.na(bpd_jen36_hc) ~ as.numeric(bpd_jen36_hc == 3),
     TRUE ~ NA_real_),
  race = factor(race3x, 1:3, c("Black", "White", "Other")),
  hispanic = factor(hisp, 0:1, c("Not Hispanic", "Hispanic")),
  respiratory_index = map_fio2x,
  ga_stratum = factor(ifelse(gestage < 27, "<27 weeks", ">=27 weeks")),
  se_unit = infant_id)
HAS_FAMILY_LINKAGE <- file.exists(FAMILY_FILE)
if (HAS_FAMILY_LINKAGE) {
  fam <- read_source(FAMILY_FILE)
  assert(all(c("randcent", "rand_id", "family_id") %in% names(fam)), "Incomplete family linkage columns")
  fam <- fam %>% transmute(randcent = num(randcent), rand_id = num(rand_id),
                          family_id = as.character(family_id))
  assert_unique(fam, c("randcent", "rand_id"), "Family linkage")
  assert(!anyNA(fam$family_id) && all(nzchar(fam$family_id)), "Missing genuine family IDs")
  df <- left_join(df, fam, by = c("randcent", "rand_id"), relationship = "one-to-one")
  assert(!anyNA(df$family_id), "Family linkage must cover every randomized infant")
  df <- df %>% mutate(se_unit = paste(randcent, family_id, sep = "_"))
  COVARIANCE_LABEL <- "HC1 sandwich, clustered by externally supplied genuine family ID"
} else {
  COVARIANCE_LABEL <- "HC1 sandwich with unique infants as independent units; sibling linkage unavailable"
}
# No fake vfd36 variable is created. This source endpoint is survivor ventilation duration.
assert(all(df$ventilation_days_36[!is.na(df$ventilation_days_36)] >= 0), "Negative source ventilation days")
for (v in c("extub14", "death36", "survival_without_physiological_bpd", "ndi_source_binary")) {
  assert(all(na.omit(df[[v]]) %in% 0:1), paste("Nonbinary values:", v))
}
assert(all(na.omit(df$bpd_jen36_hc) %in% 0:3), "Unknown Jensen grade")
assert(!any(df$death36 == 1 & !is.na(df$bpd_jen36_hc)), "Jensen grades assigned to deaths: investigate")
assert(all(df$survival_without_physiological_bpd ==
             as.numeric(df$death36 == 0 & !is.na(df$bpd_p_hc) & df$bpd_p_hc == 0)),
       "Physiological survival endpoint does not reconcile")
okcat <- !is.na(df$ndi3cat) & !is.na(df$ndi_source_binary)
assert(all(df$ndi_source_binary[okcat] == as.numeric(df$ndi3cat[okcat] %in% 2:3)),
       "NDI3 conflicts with observed NDI3CAT")
ndi_audit <- df %>% count(ndi3cat, ndi_source_binary, ndi_categorized_modsev, name = "N")
write_csv(ndi_audit, file.path(DIR$logs, "ndi_definition_crosswalk.csv"))
write_csv(df %>% count(death36, bpd_p_hc, bpd_jen36_hc,
                       survival_without_physiological_bpd, name = "N"),
          file.path(DIR$logs, "bpd_definition_crosswalk.csv"))

#####
#derive and preserve original phenotypes-----
#####
eligible <- complete.cases(df[, CLUSTER_VARS])
clustering_data <- df[eligible, ]
X <- scale(as.matrix(clustering_data[, CLUSTER_VARS]))
hc_original <- hclust(dist(X), method = "ward.D2")
assignments <- tibble(infant_id = clustering_data$infant_id,
                      phenotype = factor(paste("Phenotype", cutree(hc_original, 4)), levels = PHENO))
df_pheno <- left_join(df, assignments, by = "infant_id", relationship = "one-to-one")
assert(nrow(df_pheno) == nrow(raw), "Clustering merge changed row count")
analytic <- df_pheno %>% filter(!is.na(phenotype))
assert_unique(analytic, "infant_id", "Phenotyped cohort")
if (STRICT_ORIGINAL_COHORT) {
  assert(nrow(analytic) == 760L, "Expected 760 complete-case infants")
  assert(identical(as.integer(table(analytic$phenotype)), c(59L, 344L, 216L, 141L)),
         "Original phenotype sizes changed; STOP and investigate")
}
if (file.exists(PREVIOUS_DATASET)) {
  old <- as.data.frame(readRDS(PREVIOUS_DATASET)); names(old) <- tolower(names(old))
  assert(all(c("randcent", "rand_id", "phenotype") %in% names(old)), "Previous RDS lacks comparison fields")
  old <- old %>% transmute(infant_id = make_id(randcent, rand_id), previous = as.character(phenotype))
  assert_unique(old, "infant_id", "Previous working dataset")
  comp <- left_join(df_pheno %>% select(infant_id, phenotype), old,
                    by = "infant_id", relationship = "one-to-one")
  assert(nrow(comp) == nrow(df_pheno), "Previous-data comparison changed row count")
  same <- (is.na(comp$phenotype) & is.na(comp$previous)) |
    (!is.na(comp$phenotype) & !is.na(comp$previous) & as.character(comp$phenotype) == comp$previous)
  assert(all(same), "Saved phenotype assignments differ from reconstructed assignments")
  record_status("Original phenotype preservation", "COMPLETE", "Every saved assignment agrees")
}
saveRDS(df_pheno, file.path(DIR$data, "df_pheno_corrected.rds"))
saveRDS(analytic, file.path(DIR$data, "analytic_760.rds"))
write_csv(analytic, file.path(DIR$data, "analytic_760.csv"))
saveRDS(list(hclust = hc_original, scaled = X, center = attr(X, "scaled:center"),
             scale = attr(X, "scaled:scale"), assignments = assignments,
             variables = CLUSTER_VARS), file.path(DIR$models, "original_clustering.rds"))
record_status("Cohort and outcome validation", "COMPLETE", COVARIANCE_LABEL)

#####
#table 1 baseline characteristics-----
#####
continuous_summary <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return("Not available")
  q <- quantile(v, c(.25, .5, .75), names = FALSE)
  sprintf("%.1f (%.1f, %.1f); n=%d", q[2], q[1], q[3], length(v))
}
baseline_rows <- function(d) {
  vars <- c(CLUSTER_VARS, "obgmage")
  labs <- c(unname(VAR_LABELS), "Maternal age, years")
  r <- tibble(Characteristic = labs, Value = vapply(d[vars], continuous_summary, character(1)))
  for (v in c("male", "multbrth", "antester", "race", "hispanic", "sga")) {
    a <- d[[v]]; lev <- if (is.factor(a)) levels(a) else sort(unique(na.omit(a)))
    for (l in lev) r <- bind_rows(r, tibble(Characteristic = paste(v, l, sep = ": "),
      Value = sprintf("%d/%d (%.1f%%)", sum(a == l, na.rm = TRUE), sum(!is.na(a)),
                      100 * mean(a == l, na.rm = TRUE))))
    r <- bind_rows(r, tibble(Characteristic = paste(v, "missing"), Value = as.character(sum(is.na(a)))))
  }
  r
}
table1 <- baseline_rows(analytic) %>% rename(Overall = Value)
for (p in PHENO) {
  a <- baseline_rows(analytic[analytic$phenotype == p, ])
  names(a)[2] <- p; table1 <- left_join(table1, a, by = "Characteristic", relationship = "one-to-one")
}
save_table(table1, "Table1_Baseline", "Table 1 Baseline characteristics by phenotype",
  "Median (IQR) or n/N (%); missing values shown. No hypothesis tests for clustering variables. Respiratory index = MAP x FiO2 expressed as percent (equivalently MAP x fractional FiO2 x 100).",
  landscape = TRUE)
excluded_table <- bind_rows(
  baseline_rows(df_pheno[!eligible, ]) %>% rename(Excluded = Value),
  .id = NULL) %>% left_join(baseline_rows(analytic) %>% rename(Included = Value), by = "Characteristic")
save_table(excluded_table, "sTable1_Included_Excluded", "Supplementary Table 1 Included and excluded infants")

#####
#primary adjusted models and table 2-----
#####
primary <- fit_binary(analytic, "extub14")
primary_global <- wald_block(coef(primary$fit), primary$V, ":hc$|^hc:")
primary_by_pheno <- all_binary_contrasts(primary)
overall <- fit_binary(analytic, "extub14", interaction = FALSE)
table2_numeric <- bind_rows(all_binary_contrasts(overall), primary_by_pheno)
write_csv(table2_numeric, file.path(DIR$tables, "Table2_Unrounded.csv"))
TABLE2_NOTE <- paste(
  "Observed counts are unadjusted. Logistic model adjusted for gestational age, sex, and study center.",
  COVARIANCE_LABEL,
  "Risks are empirically standardized over observed covariates WITHIN each phenotype.",
  "The overall row uses a separate adjusted model without phenotype terms in the same analytic cohort.",
  "RR and RD derive from the same standardized risks; covariance is retained by the delta method.",
  "Conditional OR differs from the odds ratio of standardized risks. Row P values test conditional ORs.",
  "Global interaction is a robust joint Wald test on the log-odds scale: P", fmt_p(primary_global$p_value),
  "Intervals condition on observed covariates; no multiplicity adjustment. Results are exploratory.")
save_table(display_binary(table2_numeric), "Table2_Adjusted_Extubation",
           "Table 2 Hydrocortisone and extubation within 14 days", TABLE2_NOTE, landscape = TRUE)
save_model(primary, "primary_adjusted_logistic")
save_model(overall, "overall_adjusted_logistic")
save_model(list(global = primary_global, estimates = table2_numeric), "primary_results")
write_csv(coef_table(primary$fit, primary$V), file.path(DIR$tables, "Primary_All_Coefficients.csv"))
record_status("Primary model and Table 2", "COMPLETE", paste("Interaction P", fmt_p(primary_global$p_value)))
if (file.exists(BENCHMARK_FILE) && !HAS_FAMILY_LINKAGE &&
    identical(unname(tools::md5sum(ANALYSIS_FILE)), "f77e3d01afdb8fd0ebf790b467ea11f5")) {
  benchmark <- read_csv(BENCHMARK_FILE, show_col_types = FALSE)
  fields <- c("adjusted_placebo", "adjusted_hc", "RR", "RR_lower", "RR_upper", "RD", "RD_lower", "RD_upper")
  aligned <- match(primary_by_pheno$group, benchmark$group)
  assert(!anyNA(aligned), "Independent benchmark groups do not align")
  delta <- abs(as.matrix(primary_by_pheno[,fields]) - as.matrix(benchmark[aligned,fields]))
  tolerance <- 1e-7 + 1e-5 * pmax(abs(as.matrix(primary_by_pheno[,fields])),
                                 abs(as.matrix(benchmark[aligned,fields])))
  assert(all(delta <= tolerance), "R output differs from independently fitted primary benchmark")
  assert(abs(primary_global$p_value-benchmark$global_log_odds_interaction_p[1]) < 1e-5,
         "R primary interaction differs from independent benchmark")
  record_status("Independent numerical benchmark", "COMPLETE", paste("Maximum numeric discrepancy", max(delta)))
}

#####
#figure 1-----
#####
pca <- prcomp(X, center = FALSE, scale. = FALSE)
pcs <- as_tibble(pca$x[, 1:2]) %>% mutate(phenotype = analytic$phenotype)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
p_pca <- ggplot(pcs, aes(PC1, PC2, color = phenotype)) +
  geom_point(alpha = .45, size = 1.1) + stat_ellipse(type = "norm", level = .68, linewidth = .6) +
  scale_color_manual(values = COLS) +
  labs(title = "A  Baseline phenotype structure", subtitle = "Ellipses show 68% normal data contours",
       x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]), color = NULL)
centroids <- as_tibble(X) %>% mutate(phenotype = analytic$phenotype) %>%
  group_by(phenotype) %>% summarise(across(all_of(CLUSTER_VARS), mean)) %>%
  pivot_longer(-phenotype, names_to = "feature", values_to = "z")
profile_plot <- function(cen, title) {
  ggplot(cen, aes(factor(feature, levels = CLUSTER_VARS), phenotype, fill = z)) +
    geom_tile(color = "white", linewidth = .5) +
    geom_text(aes(label = sprintf("%.2f", z)), size = 3) +
    scale_fill_gradient2(low = "#6688BA", mid = "white", high = "#BE3948", midpoint = 0) +
    scale_x_discrete(labels = VAR_LABELS) + scale_y_discrete(limits = rev) +
    labs(title = title, x = NULL, y = NULL, fill = "Mean z score") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}
p_profiles <- profile_plot(centroids, "B  Standardized baseline features")
fig1 <- (p_pca | p_profiles) + plot_layout(widths = c(1.05, 1)) +
  plot_annotation(title = "Baseline phenotypes in ventilated preterm infants")
save_fig(fig1, "Figure1_Baseline_Phenotypes", 12, 5.8)
save_model(pca, "baseline_PCA")
write_csv(centroids, file.path(DIR$tables, "Baseline_Standardized_Centroids.csv"))

#####
#figure 2-----
#####
extub_obs <- observed_risk(analytic, "extub14")
ndi_obs <- observed_risk(analytic, "ndi_source_binary")
plot_observed <- function(d, title, ylab) {
  ggplot(d, aes(phenotype, prob, color = phenotype)) +
    geom_segment(aes(y = lower, yend = upper, xend = phenotype), linewidth = .7) +
    geom_point(size = 2.7) +
    geom_text(aes(label = paste0(events, "/", N)), nudge_x = .17, hjust = 0, size = 3, color = "black") +
    scale_color_manual(values = COLS) + scale_y_continuous(labels = scales::percent,
                                                        limits = c(0, 1)) +
    labs(title = title, x = NULL, y = ylab) + theme(legend.position = "none")
}
p2a <- plot_observed(extub_obs, "A  Successful extubation within 14 days", "Observed proportion (95% CI)")
bpd_dist <- analytic %>% filter(death36 == 0, !is.na(jensen_grade)) %>%
  count(phenotype, jensen_grade, .drop = FALSE) %>% group_by(phenotype) %>% mutate(prop = n / sum(n))
p2b <- ggplot(bpd_dist, aes(phenotype, prop, fill = jensen_grade)) +
  geom_col(width = .7) + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#DAE8EC", "#9CC5D4", "#4B95B4", "#164F6B")) +
  labs(title = "B  Jensen BPD severity among survivors", subtitle = paste("Observed severity, n =", sum(bpd_dist$n)),
       x = NULL, y = "Proportion", fill = NULL)
p2c <- ggplot(analytic %>% filter(death36 == 0, !is.na(ventilation_days_36)),
             aes(phenotype, ventilation_days_36, fill = phenotype)) +
  geom_boxplot(width = .6, outlier.alpha = .25) + scale_fill_manual(values = COLS) +
  labs(title = "C  Ventilation duration among survivors", x = NULL,
       y = "Mechanical ventilation days to 36 weeks PMA") + theme(legend.position = "none")
p2d <- plot_observed(ndi_obs, "D  Source-defined binary NDI at follow-up", "Observed proportion (95% CI)") +
  labs(subtitle = "Protocol target: 22-26 months corrected age; actual ages summarized separately")
fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_annotation(title = "Observed clinical outcomes by baseline phenotype",
                  caption = "Descriptive comparisons. BPD and NDI use available survivor follow-up data; denominators vary.")
save_fig(fig2, "Figure2_Clinical_Outcomes", 13, 10)

#####
#figure 3-----
#####
fig3 <- (forest(table2_numeric, "A  Successful extubation", effect = "RR") |
         forest(table2_numeric, "B  Absolute treatment contrast", effect = "RD")) +
  plot_annotation(title = "Hydrocortisone and successful extubation within 14 days",
    subtitle = paste("Adjusted log-odds treatment-by-phenotype interaction: P", fmt_p(primary_global$p_value)),
    caption = "Adjusted risks standardized within phenotype. Overall row uses a separate model. Exploratory subgroup estimates do not establish a treatment-selection rule.")
save_fig(fig3, "Figure3_Adjusted_Treatment_Contrasts", 12, 5.8)

#####
#supplementary binary outcomes-----
#####
secondary_names <- c(survival_without_physiological_bpd = "Survival without physiological moderate/severe BPD",
                     death36 = "Death before 36 weeks PMA",
                     death_or_grade3 = "Death or Jensen grade 3 BPD",
                     ndi_source_binary = "Source-defined binary NDI at follow-up",
                     ndi_categorized_modsev = "Categorized moderate/severe NDI at follow-up")
secondary <- list()
for (outcome in names(secondary_names)) {
  ans <- run_module(paste("Secondary binary", outcome), {
    obj <- fit_binary(analytic, outcome)
    est <- all_binary_contrasts(obj)
    test <- wald_block(coef(obj$fit), obj$V, ":hc$|^hc:")
    save_model(obj, paste0("secondary_", outcome))
    write_csv(est, file.path(DIR$tables, paste0("Secondary_", outcome, "_Unrounded.csv")))
    save_table(display_binary(est), paste0("sTable_", outcome), secondary_names[[outcome]],
      paste("Observed counts and proportions are unadjusted. Logistic models include treatment, phenotype, their interaction, continuous gestational age, sex, and fixed center effects.",
            COVARIANCE_LABEL,
            "Risks are standardized over observed covariates within each phenotype. RR and RD use the same risks; delta-method CIs condition on observed covariates.",
            "Conditional ORs and row P values refer to the fitted conditional treatment odds ratios, not the standardized risks.",
            "Outcome-specific complete cases. Global treatment-by-phenotype interaction is a joint Wald test on the log-odds scale: P", fmt_p(test$p_value),
            "Comparisons are exploratory and not adjusted for multiplicity. Genuine sibling dependence is not accounted for without verified linkage."), TRUE)
    list(model = obj, estimates = est, interaction = test)
  })
  secondary[[outcome]] <- ans
}
if (!is.null(secondary$survival_without_physiological_bpd)) {
  s <- secondary$survival_without_physiological_bpd
  save_fig(forest(s$estimates, "Survival without physiologically defined moderate/severe BPD",
     paste("Adjusted interaction P", fmt_p(s$interaction$p_value))),
     "sFigure_Physiological_BPD_Treatment", 9, 5, TRUE)
}

#####
#ventilation duration analysis not ventilator free days-----
#####
ventilation_results <- run_module("Ventilation duration", {
  vd <- analytic %>% filter(death36 == 0, !is.na(ventilation_days_36)) %>% droplevels()
  f <- lm(ventilation_days_36 ~ phenotype * hc + gestage + male + center, data = vd, na.action = na.fail)
  V <- robust_vcov(f, vd)
  test <- wald_block(coef(f), V, ":hc$|^hc:")
  est <- bind_rows(lapply(PHENO, function(p) {
    d0 <- vd[vd$phenotype == p, ]; d1 <- d0; d0$hc <- 0; d1$hc <- 1
    g <- colMeans(design(f, d1) - design(f, d0))
    b <- sum(g * coef(f)); s <- sqrt(as.numeric(t(g) %*% V %*% g))
    tibble(phenotype = p, N = nrow(d0), mean_difference_days = b,
           lower = b - Z * s, upper = b + Z * s,
           p_value = 2 * pnorm(abs(b / s), lower.tail = FALSE))
  }))
  save_model(list(fit = f, V = V, data = vd, estimates = est, interaction = test), "ventilation_duration")
  save_table(est, "sTable_Ventilation_Duration", "Ventilation duration to 36 weeks among survivors",
    paste("Adjusted mean difference HC minus placebo; negative values indicate shorter ventilation.",
          COVARIANCE_LABEL, "Survivor-only analysis; no death imputation or ventilator-free-days derivation.",
          "Robust joint interaction P", fmt_p(test$p_value)))
  est
})

#####
#ordinal BPD severity and threshold specific analyses-----
#####
ordinal_results <- run_module("Ordinal BPD severity", {
  od <- analytic %>% filter(death36 == 0, !is.na(jensen_grade)) %>% droplevels()
  fit <- ordinal::clm(jensen_grade ~ phenotype * hc + gestage + male + center,
                       data = od, link = "logit", model = TRUE)
  conv_code <- if (is.list(fit$convergence)) fit$convergence$code else fit$convergence
  writeLines(capture.output(str(list(convergence=fit$convergence,gradient=fit$maxGradient,
                                     Hessian=fit$Hessian))),
             file.path(DIR$logs,"ordinal_convergence_diagnostics.txt"))
  assert(all(is.finite(coef(fit))) && all(is.finite(vcov(fit))), "Invalid ordinal fit")
  assert(is.numeric(conv_code) && length(conv_code)>0 && all(is.finite(conv_code)) && all(conv_code==0),
         "Ordinal optimizer did not find a unique converged solution; inspect ordinal_convergence_diagnostics.txt")
  assert(is.finite(fit$maxGradient) && fit$maxGradient < 1e-4, "Ordinal model gradient too large")
  # clm fits a conditional cumulative logit. Positive location coefficients
  # imply a shift toward HIGHER severity; exponentiated terms are higher-category ORs.
  b <- coef(fit); V <- vcov(fit)
  test <- wald_block(b, V, ":hc$|^hc:")
  est <- bind_rows(lapply(PHENO, function(p) {
    g <- setNames(rep(0, length(b)), names(b)); g["hc"] <- 1
    if (p != PHENO[1]) {
      term <- paste0("phenotype", p, ":hc")
      assert(term %in% names(g), paste("Missing ordinal contrast", term)); g[term] <- 1
    }
    z <- sum(g * b); s <- sqrt(as.numeric(t(g) %*% V %*% g))
    tibble(phenotype = p, OR_higher_grade = exp(z), lower = exp(z - Z*s), upper = exp(z + Z*s),
           p_value = 2 * pnorm(abs(z/s), lower.tail = FALSE))
  }))
  save_table(est, "sTable_Ordinal_BPD", "Jensen severity among survivors",
    paste("Proportional-odds cumulative logit, adjusted for gestational age, sex, and center.",
      "Model-based covariance; this is not the previous ordinal GEE.",
      if (HAS_FAMILY_LINKAGE) "Family dependence is not accounted for in this CLM; see threshold-specific robust analyses." else "Sibling linkage unavailable.",
      "OR >1 indicates greater severity. Full phenotype-specific treatment contrasts, not interaction coefficients alone.",
      "Conditioning on survival can induce selection. Joint interaction P", fmt_p(test$p_value)))
  # Assess the proportional-odds assumption; a failed diagnostic is recorded.
  po_test <- tryCatch(ordinal::nominal_test(fit), error = function(e) {
    record_status("Ordinal proportional odds diagnostic", "FAILED", conditionMessage(e)); NULL
  })
  if (!is.null(po_test)) {
    pt <- as.data.frame(po_test); pt$term <- rownames(pt)
    write_csv(pt, file.path(DIR$tables, "sTable_Ordinal_Assumption_Diagnostic.csv"))
  }
  save_model(list(fit = fit, data = od, estimates = est, interaction = test, po_test = po_test), "ordinal_BPD")
  list(fit = fit, estimates = est, interaction = test)
})
threshold_results <- list()
for (k in 1:3) {
  res <- run_module(paste("BPD threshold", k), {
    od <- analytic %>% filter(death36 == 0) %>%
      mutate(threshold_outcome = ifelse(is.na(bpd_jen36_hc), NA, as.numeric(bpd_jen36_hc >= k)))
    obj <- fit_binary(od, "threshold_outcome")
    est <- all_binary_contrasts(obj) %>% mutate(threshold = paste0("Jensen grade >=", k))
    save_model(obj, paste0("BPD_threshold_", k)); est
  })
  threshold_results[[as.character(k)]] <- res
}
if (length(Filter(Negate(is.null), threshold_results))) {
  save_table(bind_rows(threshold_results), "sTable_BPD_Thresholds", "Exploratory Jensen severity threshold analyses",
    "Adjusted logistic models with the primary robust covariance. Among survivors with observed grade. These are overlapping exploratory outcomes without multiplicity correction.", TRUE)
}
save_table(analytic %>% filter(!is.na(jensen_grade)) %>%
             count(phenotype, treatment, jensen_grade, .drop = FALSE),
           "sTable_BPD_Grade_Counts", "Jensen grade counts by phenotype and treatment")

#####
#supplementary model sensitivity analyses-----
#####
run_module("Unadjusted and adjusted modified Poisson", {
  crude <- fit_binary(analytic, "extub14", adjusted = FALSE, link = "log")
  adj <- fit_binary(analytic, "extub14", adjusted = TRUE, link = "log")
  for (nm in c("crude", "adj")) {
    ob <- get(nm); save_model(ob, paste0("poisson_", nm))
    save_table(all_binary_contrasts(ob), paste0("sTable_Poisson_", nm),
      paste("Modified Poisson treatment contrasts", nm),
      paste(COVARIANCE_LABEL, "Joint log-RR interaction P",
            fmt_p(wald_block(coef(ob$fit), ob$V, ":hc$|^hc:")$p_value)), TRUE)
  }
  invisible(TRUE)
})
glmm_results <- run_module("Center random intercept sensitivity", {
  d <- primary$data
  random <- if (HAS_FAMILY_LINKAGE) "(1 | center) + (1 | se_unit)" else "(1 | center)"
  full <- as.formula(paste("extub14 ~ phenotype * hc + gestage + male +", random))
  reduced <- as.formula(paste("extub14 ~ phenotype + hc + gestage + male +", random))
  fit <- lme4::glmer(full, data = d, family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  fit0 <- lme4::glmer(reduced, data = d, family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  check_mer(fit, "GLMM full model")
  check_mer(fit0, "GLMM reduced model")
  em <- emmeans::emmeans(fit, ~ hc | phenotype, at = list(hc = c(0, 1)))
  ct <- emmeans::contrast(em, "revpairwise", by = "phenotype")
  est <- as.data.frame(summary(ct, infer = c(TRUE, TRUE), type = "response"))
  lrt <- as.data.frame(anova(fit0, fit)); lrt$model <- rownames(lrt)
  save_model(list(full = fit, reduced = fit0, contrasts = est, lrt = lrt), "GLMM_sensitivity")
  save_table(est, "sTable_GLMM", "Conditional GLMM odds ratios",
             "Center random intercept replaces fixed center effects. No artificial infant random intercept in one-row-per-infant data.")
  write_csv(lrt, file.path(DIR$tables, "GLMM_Interaction_LRT.csv"))
  write_csv(as.data.frame(lme4::VarCorr(fit)), file.path(DIR$tables, "GLMM_Variance_Components.csv"))
  est
})
run_module("Firth sparse-stratum sensitivity", {
  fits <- setNames(lapply(PHENO, function(p) {
    dd <- analytic[analytic$phenotype == p, ]
    logistf::logistf(extub14 ~ hc + gestage + male, data = dd)
  }), PHENO)
  res <- bind_rows(lapply(PHENO, function(p) {
    f <- fits[[p]]; i <- match("hc", names(coef(f)))
    tibble(phenotype = p, OR = exp(coef(f)[i]), lower = exp(f$ci.lower[i]),
           upper = exp(f$ci.upper[i]), p_value = f$prob[i])
  }))
  save_model(fits, "Firth_by_phenotype")
  save_table(res, "sTable_Firth", "Firth logistic sensitivity",
    "Adjusted for gestational age and sex. Site and sibling correlation not accounted for. Per-stratum sparse-data sensitivity, not the primary model.")
  res
})
run_module("Study center CR2 sensitivity", {
  Vcr2 <- clubSandwich::vcovCR(primary$fit, cluster = primary$data$center, type = "CR2")
  idx <- grep(":hc$|^hc:", names(coef(primary$fit)))
  C <- diag(length(coef(primary$fit)))[idx, , drop = FALSE]
  tst <- clubSandwich::Wald_test(primary$fit, constraints = C, vcov = Vcr2, test = "HTZ")
  save_table(as.data.frame(tst), "sTable_Center_CR2", "Center-clustered small-sample interaction sensitivity",
    "CR2 covariance and HTZ test. Fixed center effects remain in the mean model. Few centers can limit precision; this is a sensitivity analysis.")
  save_model(list(V = Vcr2, test = tst), "center_CR2_sensitivity")
  tst
})

#####
#Apgar exclusion merged phenotypes and supplementary table 13-----
#####
noap_vars <- setdiff(CLUSTER_VARS, "apgar5")
X_noap <- scale(as.matrix(analytic[, noap_vars]))
hc_noap <- hclust(dist(X_noap), method = "ward.D2")
noap_id <- cutree(hc_noap, 4)
analytic$noapgar <- factor(paste("No-Apgar cluster", noap_id), levels = paste("No-Apgar cluster", 1:4))
analytic$merged <- factor(ifelse(analytic$phenotype %in% PHENO[2:3], "Phenotypes 2+3",
                                 as.character(analytic$phenotype)),
                          levels = c(PHENO[1], "Phenotypes 2+3", PHENO[4]))
ari <- mclust::adjustedRandIndex(as.integer(analytic$phenotype), noap_id)
correspondence <- as.data.frame(table(Original = analytic$phenotype, Alternative = analytic$noapgar))
write_csv(correspondence, file.path(DIR$tables, "sTable13_Cluster_Correspondence_Long.csv"))
saveRDS(analytic, file.path(DIR$data, "analytic_with_sensitivity_partitions.rds"))
save_model(list(hclust = hc_noap, scaled = X_noap, assignments = analytic[, c("infant_id", "noapgar")],
                ARI = ari), "no_apgar_clustering")
s13 <- list(); s13_tests <- list(); sensitivity_models <- list()
for (partition in c("phenotype", "merged", "noapgar")) {
  for (adjusted in c(FALSE, TRUE)) {
    nm <- paste(partition, if (adjusted) "adjusted_logit" else "crude_logRR", sep = "_")
    ans <- run_module(nm, {
      obj <- fit_binary(analytic, "extub14", group = partition, adjusted = adjusted,
                        link = if (adjusted) "logit" else "log")
      est <- all_binary_contrasts(obj) %>% mutate(partition = partition,
                 model = if (adjusted) "Adjusted logistic standardized risks" else "Unadjusted modified Poisson")
      tst <- wald_block(coef(obj$fit), obj$V, ":hc$|^hc:") %>% mutate(analysis = nm)
      save_model(obj, paste0("sensitivity_", nm))
      list(obj = obj, est = est, test = tst)
    })
    if (!is.null(ans)) {
      sensitivity_models[[nm]] <- ans$obj; s13[[nm]] <- ans$est; s13_tests[[nm]] <- ans$test
    }
  }
}
pooled <- fit_binary(analytic %>% filter(phenotype %in% PHENO[2:3]), "extub14",
                     adjusted = FALSE, interaction = FALSE, link = "log")
pooled_est <- all_binary_contrasts(pooled) %>% mutate(group = "Original Phenotypes 2+3 pooled")
save_model(pooled, "pooled_original_2_3")
save_table(bind_rows(s13), "sTable13_Phenotype_Sensitivity", "Sensitivity to phenotype definition",
  paste("New labels are independent of original labels. Same baseline-complete cohort. ARI =", round(ari, 3),
    ". Adjusted log-odds and crude log-RR interaction tests are distinct. See interaction table.", COVARIANCE_LABEL), TRUE)
save_table(bind_rows(s13_tests), "sTable13_Interaction_Tests", "Interaction tests for each phenotype partition")
save_table(pooled_est, "sTable13_Pooled_2_3", "Pooled original Phenotypes 2 and 3",
  "Unadjusted treatment-only modified Poisson model; a pooled treatment association is not a heterogeneity test.", TRUE)
# Direct original Phenotype 2 versus 3 interaction contrast, rather than comparing
# one subgroup P value with another subgroup P value.
run_module("Direct Phenotype 2 versus 3 treatment contrast", {
  ob <- sensitivity_models$phenotype_crude_logRR
  b <- coef(ob$fit); g <- setNames(rep(0, length(b)), names(b))
  g[paste0("phenotype", PHENO[2], ":hc")] <- 1
  g[paste0("phenotype", PHENO[3], ":hc")] <- -1
  assert(length(g) == length(b), "Missing direct-contrast coefficient")
  a <- sum(g*b); se <- sqrt(as.numeric(t(g) %*% ob$V %*% g))
  out <- tibble(contrast = "Treatment RR in Phenotype 2 divided by treatment RR in Phenotype 3",
    ratio_of_RRs = exp(a), lower = exp(a-Z*se), upper = exp(a+Z*se),
    p_value = 2*pnorm(abs(a/se), lower.tail = FALSE))
  save_table(out, "sTable13_Direct_P2_P3", "Direct comparison of treatment associations")
  out
})
pca_noap <- prcomp(X_noap, center = FALSE, scale. = FALSE)
vn <- 100 * pca_noap$sdev^2/sum(pca_noap$sdev^2)
pn <- as_tibble(pca_noap$x[,1:2]) %>% mutate(cluster = analytic$noapgar)
noap_cols <- setNames(unname(COLS), levels(analytic$noapgar))
pn1 <- ggplot(pn, aes(PC1, PC2, color = cluster)) + geom_point(alpha = .4, size = 1) +
  stat_ellipse(type = "norm", level = .68) + scale_color_manual(values = noap_cols) +
  labs(title = "A  PCA without Apgar", x = sprintf("PC1 (%.1f%%)", vn[1]),
       y = sprintf("PC2 (%.1f%%)", vn[2]), color = NULL)
cn <- as_tibble(scale(as.matrix(analytic[,CLUSTER_VARS]))) %>% mutate(phenotype = analytic$noapgar) %>%
  group_by(phenotype) %>% summarise(across(all_of(CLUSTER_VARS), mean)) %>%
  pivot_longer(-phenotype, names_to = "feature", values_to = "z")
pn2 <- profile_plot(cn, "B  Features of no-Apgar clusters") +
  labs(subtitle = "Apgar displayed descriptively; excluded from clustering")
cr <- correspondence %>% group_by(Original) %>% mutate(row_pct = Freq/sum(Freq))
pn3 <- ggplot(cr, aes(Alternative, Original, fill = row_pct)) + geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Freq, 100*row_pct)), size = 3) +
  scale_fill_gradient(low = "white", high = "#6688BA", labels = scales::percent) +
  labs(title = "C  Correspondence with original phenotypes", subtitle = paste("Adjusted Rand index", round(ari, 3)),
       x = NULL, y = NULL, fill = "Row %") + theme(axis.text.x = element_text(angle = 30, hjust = 1))
noap_res <- s13$noapgar_crude_logRR
noap_test <- s13_tests$noapgar_crude_logRR
if (!is.null(noap_res)) {
  pn4 <- forest(bind_rows(noap_res, pooled_est), "D  Unadjusted treatment associations",
    paste("No-Apgar log-RR interaction P", fmt_p(noap_test$p_value)))
  save_fig((pn1 | pn2)/(pn3 | pn4) +
    plot_annotation(title = "Sensitivity to Apgar exclusion",
      caption = "Cluster labels do not imply correspondence. Pooled original 2+3 is a separate analysis. No within-group significance establishes a responder phenotype."),
    "sFigure_Apgar_Sensitivity", 14, 11, TRUE)
}

#####
#clustering robustness and supplemental figures-----
#####
run_module("Clustering robustness", {
  D <- dist(X); original_numeric <- as.integer(analytic$phenotype)
  sil <- cluster::silhouette(original_numeric, D)
  sil_d <- tibble(phenotype = analytic$phenotype, silhouette = sil[,3])
  k_stats <- bind_rows(lapply(2:6, function(k) {
    a <- cutree(hc_original, k); b <- cutree(hc_noap, k)
    bind_rows(
      tibble(solution = "Original", k = k, silhouette = mean(cluster::silhouette(a, D)[,3]),
        CH = fpc::cluster.stats(D, a)$ch),
      tibble(solution = "Without Apgar", k = k,
        silhouette = mean(cluster::silhouette(b, dist(X_noap))[,3]),
        CH = fpc::cluster.stats(dist(X_noap), b)$ch))
  }))
  set.seed(2024)
  km <- kmeans(X, centers = 4, nstart = 50)
  pam <- cluster::pam(X, k = 4)
  agreement <- tibble(comparison = c("Ward versus k-means", "Ward versus PAM"),
    ARI = c(mclust::adjustedRandIndex(original_numeric, km$cluster),
            mclust::adjustedRandIndex(original_numeric, pam$clustering)))
  boot <- fpc::clusterboot(as.matrix(X), B = BOOTSTRAP_REPS, clustermethod = fpc::hclustCBI,
                           method = "ward.D2", k = 4, seed = 2024)
  # Map clusterboot's cluster numbering to original labels using maximum overlap.
  boot_ref <- boot$result$partition
  assert(length(boot_ref) == nrow(X), "Unexpected clusterboot partition output")
  stability <- bind_rows(lapply(1:4, function(k) {
    overlap <- table(boot_ref[original_numeric == k]); lab <- as.integer(names(which.max(overlap)))
    tibble(phenotype = PHENO[k], bootstrap_cluster = lab, mean_Jaccard = boot$bootmean[lab])
  }))
  set.seed(2025)
  subset_ari <- replicate(SUBSAMPLE_REPS, {
    keep <- sort(sample(seq_len(nrow(X)), floor(.9*nrow(X))))
    alt <- cutree(hclust(dist(X[keep,,drop=FALSE]), "ward.D2"), 4)
    mclust::adjustedRandIndex(original_numeric[keep], alt)
  })
  save_table(k_stats, "sTable_Cluster_Number", "Clustering validity across candidate cluster counts")
  save_table(stability, "sTable_Bootstrap_Jaccard", "Bootstrap cluster stability")
  save_table(agreement, "sTable_Cross_Algorithm", "Cross-algorithm cluster agreement")
  write_csv(tibble(iteration = seq_along(subset_ari), ARI = subset_ari),
             file.path(DIR$tables, "Subsampling_ARI.csv"))
  ps <- ggplot(sil_d, aes(silhouette, fill = phenotype)) + geom_histogram(bins = 30) +
    geom_vline(xintercept = 0, linetype = 2) + facet_wrap(~phenotype, ncol = 2) +
    scale_fill_manual(values = COLS) + theme(legend.position = "none") +
    labs(title = "Silhouette widths", x = "Silhouette width", y = "Infants")
  pk <- k_stats %>% pivot_longer(c(silhouette, CH), names_to = "metric", values_to = "value") %>%
    ggplot(aes(k, value, color = solution)) + geom_line() + geom_point() +
    facet_wrap(~metric, scales = "free_y") + labs(x = "Number of clusters", y = NULL, color = NULL)
  pj <- ggplot(stability, aes(phenotype, mean_Jaccard, fill = phenotype)) + geom_col(width = .6) +
    scale_fill_manual(values = COLS) + coord_cartesian(ylim = c(0,1)) +
    labs(title = "Bootstrap stability", x = NULL, y = "Mean Jaccard") + theme(legend.position = "none")
  pa <- ggplot(tibble(ARI = subset_ari), aes(ARI)) + geom_histogram(bins = 20, fill = "#6688BA") +
    labs(title = "90% subsample stability", subtitle = "Not an imputation or missingness-mechanism analysis", y = "Resamples")
  save_fig(ps, "sFigure_Silhouette", 10, 7, TRUE)
  save_fig(pk, "sFigure_Cluster_Number", 10, 5, TRUE)
  save_fig(pj | pa, "sFigure_Cluster_Stability", 11, 5, TRUE)
  save_model(list(bootstrap = boot, kmeans = km, pam = pam, k_stats = k_stats,
                  silhouette = sil, subsample_ARI = subset_ari), "clustering_robustness")
  invisible(TRUE)
})
run_module("Consensus co-clustering", {
  set.seed(2026); n <- nrow(X)
  together <- matrix(0, n, n); sampled <- matrix(0, n, n)
  for (b in seq_len(CONSENSUS_REPS)) {
    idx <- sample.int(n, n, replace = TRUE)
    cl <- cutree(hclust(dist(X[idx,,drop=FALSE]), "ward.D2"), 4)
    u <- sort(unique(idx)); lab <- cl[match(u, idx)]
    sampled[u,u] <- sampled[u,u] + 1
    together[u,u] <- together[u,u] + outer(lab, lab, "==")
  }
  consensus <- together/pmax(sampled, 1); diag(consensus) <- 1
  ord <- order(analytic$phenotype); M <- consensus[ord,ord]
  grid <- expand.grid(i = seq_len(n), j = seq_len(n)); grid$p <- as.vector(M)
  p <- ggplot(grid, aes(i,j,fill=p)) + geom_raster() + coord_equal() +
    scale_fill_gradient(low = "white", high = "#2878B5", limits = c(0,1)) +
    labs(title = "Bootstrap co-clustering", subtitle = "Infants ordered by original phenotype",
         x = "Infant order", y = "Infant order", fill = "Probability") +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  save_fig(p, "sFigure_Consensus", 7, 7, TRUE)
  save_model(list(consensus = consensus, sampled = sampled, order = ord), "consensus_clustering")
  invisible(TRUE)
})

#####
#simpler treatment effect modifiers and added predictive value-----
#####
simple_results <- run_module("Simple baseline treatment modifiers", {
  d <- analytic %>% mutate(ga_z = as.numeric(scale(gestage)), ri_z = as.numeric(scale(map_fio2x)))
  forms <- list(Gestational_age = extub14 ~ hc * ga_z + male + center,
                Respiratory_index = extub14 ~ hc * ri_z + gestage + male + center,
                GA_and_RI = extub14 ~ hc * (ga_z + ri_z) + male + center,
                GA_strata = extub14 ~ hc * ga_stratum + male + center)
  out <- bind_rows(lapply(names(forms), function(nm) {
    f <- glm(forms[[nm]], data = d, family = binomial, na.action = na.fail)
    check_glm(f); V <- robust_vcov(f, d)
    save_model(list(fit = f, V = V), paste0("simple_modifier_", nm))
    wald_block(coef(f), V, ":hc$|^hc:") %>% mutate(model = nm)
  }))
  save_table(out, "sTable_Simple_Modifiers", "Exploratory simpler baseline treatment modifiers",
    "Continuous predictors standardized without selecting cutoffs from treatment results. Robust log-odds interaction tests. These are not NICHD estimator analyses.")
  out
})
# Optional NICHD comparison: predictions must come from a verified external
# implementation. This script DOES NOT invent or independently implement a
# clinical calculator or substitute an internally fitted risk score for NICHD.
if (file.exists(NICHD_FILE)) {
  run_module("Externally verified NICHD risk comparison", {
    nr <- read_source(NICHD_FILE)
    req <- c("randcent", "rand_id", "risk", "model_version", "risk_endpoint", "assessment_day", "verified")
    assert(all(req %in% names(nr)), "NICHD file missing required provenance columns")
    assert_unique(nr, c("randcent", "rand_id"), "NICHD prediction file")
    assert(all(tolower(as.character(nr$verified)) %in% c("true", "1", "yes")),
           "NICHD estimates are not marked externally verified")
    assert(!anyNA(nr[,c("model_version", "risk_endpoint", "assessment_day")]),
           "NICHD model provenance contains missing fields")
    assert(all(num(nr$risk) >= 0 & num(nr$risk) <= 1), "Invalid supplied NICHD risk")
    assert(n_distinct(nr$model_version) == 1 && n_distinct(nr$risk_endpoint) == 1 &&
             n_distinct(nr$assessment_day) == 1, "Do not pool different NICHD models, endpoints, or assessment days")
    assert(all(nzchar(as.character(nr$model_version))) && all(nzchar(as.character(nr$risk_endpoint))),
           "NICHD model provenance is blank")
    nr <- nr %>% mutate(randcent = num(randcent), rand_id = num(rand_id), risk = num(risk))
    d <- left_join(analytic, nr, by = c("randcent", "rand_id"), relationship = "one-to-one") %>% filter(!is.na(risk))
    assert(nrow(d) > 100, "Insufficient verified NICHD predictions")
    assert(all(num(d$assessment_day) <= d$age_rndaysx), "NICHD assessment occurs after treatment randomization")
    f <- glm(extub14 ~ hc * risk + gestage + male + center, data = d, family = binomial)
    check_glm(f); V <- robust_vcov(f,d)
    tst <- wald_block(coef(f), V, ":hc$|^hc:")
    save_model(list(fit=f,V=V,data=d,test=tst), "externally_verified_NICHD_modifier")
    save_table(tst, "sTable_NICHD_Modifier", "Externally verified NICHD baseline risk modifier",
      paste("Continuous risk interaction; no invented risk-group thresholds. Model:", unique(d$model_version),
             "Endpoint:", unique(d$risk_endpoint), "Assessment day:", unique(d$assessment_day)))
    tst
  })
} else record_status("NICHD estimator comparison", "NOT RUN - INPUT UNAVAILABLE",
                      "Requires externally verified predictions with model, endpoint, and pretreatment assessment-day provenance; no substitute calculated")

#####
#observed extubation failures not causal nonresponders-----
#####
run_module("Observed outcomes among hydrocortisone-treated infants", {
  dd <- analytic %>% filter(hc == 1)
  a <- baseline_rows(dd %>% filter(extub14 == 0)) %>% rename(`Not extubated within 14 days` = Value)
  b <- baseline_rows(dd %>% filter(extub14 == 1)) %>% rename(`Extubated within 14 days` = Value)
  tab <- left_join(a,b,by="Characteristic",relationship="one-to-one")
  save_table(tab, "sTable_HC_Observed_Extubation", "Baseline features by observed extubation outcome in the hydrocortisone arm",
    "Observed success/failure is not individual corticosteroid response. No untreated counterfactual is observed. Descriptive; no causal responder/nonresponder classification.")
  invisible(TRUE)
})

#####
#classification helpers and split objects-----
#####
metric_values <- function(truth, prob, threshold = .5, positive = "Phenotype2") {
  truth <- as.integer(as.character(truth) == positive)
  assert(length(prob) == length(truth) && all(is.finite(prob)), "Invalid prediction vector")
  assert(all(prob >= 0 & prob <= 1), "Predictions outside [0,1]")
  pred <- as.integer(prob >= threshold)
  TP <- sum(pred==1 & truth==1); FN <- sum(pred==0 & truth==1)
  FP <- sum(pred==1 & truth==0); TN <- sum(pred==0 & truth==0)
  div <- function(a,b) if (b == 0) NA_real_ else a/b
  auc <- if (length(unique(truth))==2) as.numeric(pROC::auc(pROC::roc(truth,prob,
                     levels=c(0,1),direction="<",quiet=TRUE))) else NA_real_
  tibble(N=length(truth),TP=TP,FN=FN,FP=FP,TN=TN,AUC=auc,
    Accuracy=div(TP+TN,length(truth)),Sensitivity=div(TP,TP+FN),
    Specificity=div(TN,TN+FP),PPV=div(TP,TP+FP),NPV=div(TN,TN+FN),
    F1=div(2*TP,2*TP+FP+FN),Brier=mean((truth-prob)^2),
    Log_loss=-mean(truth*log(pmax(prob,1e-10))+(1-truth)*log(pmax(1-prob,1e-10))))
}
group_split <- function(d, prop=.7, seed=42) {
  set.seed(seed)
  units <- d %>% group_by(se_unit) %>% summarise(label=as.character(target[1]))
  train_units <- unlist(lapply(split(units$se_unit, units$label), function(v) {
    sample(v, max(1L,min(length(v)-1L,floor(prop*length(v)))))
  }), use.names=FALSE)
  tr <- which(d$se_unit %in% train_units); te <- setdiff(seq_len(nrow(d)),tr)
  assert(!length(intersect(d$se_unit[tr], d$se_unit[te])), "Group leakage across split")
  list(train=tr,test=te)
}
fold_assignment <- function(d, K=5, seed=43) {
  set.seed(seed)
  u <- d %>% group_by(se_unit) %>% summarise(label=as.character(target[1]))
  u$fold <- NA_integer_
  for (lab in unique(u$label)) {
    ii <- which(u$label==lab); u$fold[ii] <- sample(rep(seq_len(K),length.out=length(ii)))
  }
  u$fold[match(d$se_unit,u$se_unit)]
}
classifier_fit <- function(d, seed) {
  f <- as.formula(paste("target ~",paste(CLUSTER_VARS,collapse=" + ")))
  # Mean bias reduction prevents complete separation in the membership classifier.
  # This is declared and differs from ordinary logistic regression in earlier code.
  logistic <- glm(f,data=d,family=binomial(),method=brglm2::brglmFit,
                   control=brglm2::brglmControl(type="AS_mean"))
  assert(isTRUE(logistic$converged) && all(is.finite(coef(logistic))),
         "Bias-reduced membership classifier did not converge")
  tree <- rpart::rpart(f,data=d,method="class",model=TRUE,x=TRUE,y=TRUE,
    control=rpart::rpart.control(maxdepth=3,minbucket=20,cp=.01,xval=0))
  rf <- ranger::ranger(f,data=as.data.frame(d),probability=TRUE,num.trees=RF_TREES,
    mtry=2,min.node.size=10,importance="permutation",seed=seed,num.threads=1)
  list(Logistic=logistic,Decision_tree=tree,Random_forest=rf)
}
classifier_predict <- function(fits,d) {
  tibble(Logistic=as.numeric(predict(fits$Logistic,newdata=d,type="response")),
         Decision_tree=as.numeric(predict(fits$Decision_tree,newdata=d,type="prob")[,"Phenotype2"]),
         Random_forest=as.numeric(predict(fits$Random_forest,data=as.data.frame(d))$predictions[,"Phenotype2"]))
}

#####
#phenotype 2 classification internal validation and importance-----
#####
classification <- run_module("Phenotype 2 internal classification", {
  md <- analytic %>% mutate(target=factor(ifelse(phenotype==PHENO[2],"Phenotype2","Other"),
                                         levels=c("Other","Phenotype2")))
  sp <- group_split(md); train <- md[sp$train,]; test <- md[sp$test,]
  saveRDS(list(train=train,test=test,split=sp),file.path(DIR$data,"classifier_train_test.rds"))
  folds <- fold_assignment(train,CV_FOLDS)
  cv_predictions <- list(); cv_models <- list()
  for (k in seq_len(CV_FOLDS)) {
    tr <- train[folds!=k,]; va <- train[folds==k,]
    assert(length(unique(tr$target))==2 && length(unique(va$target))==2,"CV fold lacks a target class")
    assert(!length(intersect(tr$se_unit,va$se_unit)),"Group leakage in CV")
    fits <- classifier_fit(tr,100+k)
    pred <- classifier_predict(fits,va)
    cv_predictions[[k]] <- bind_cols(va %>% select(infant_id,target),pred) %>% mutate(fold=k)
    cv_models[[k]] <- fits
  }
  cv <- bind_rows(cv_predictions)
  final_fits <- classifier_fit(train,42)
  preds <- bind_cols(test %>% select(infant_id,target,race),classifier_predict(final_fits,test))
  metrics <- bind_rows(lapply(c("Logistic","Decision_tree","Random_forest"),function(m) {
    bind_rows(metric_values(cv$target,cv[[m]]) %>% mutate(Sample="Training out-of-fold"),
              metric_values(preds$target,preds[[m]]) %>% mutate(Sample="Held-out test")) %>% mutate(Model=m)
  }))
  save_table(metrics,"sTable_Classifier_Performance","Binary Phenotype 2 membership classification",
    paste(CV_FOLDS,"-fold CV occurs only within the training split. Fixed hyperparameters, threshold 0.5.",
      "Logistic classifier uses mean bias reduction. Positive class: Phenotype 2.",
      "Labels come from clustering the full analytic cohort: validation is conditional internal label reproduction, not independent phenotyping validation or treatment-benefit prediction."),TRUE)
  write_csv(preds,file.path(DIR$tables,"Classifier_Heldout_Probabilities.csv"))
  write_csv(cv,file.path(DIR$tables,"Classifier_OutOfFold_Probabilities.csv"))
  write_csv(tibble(infant_id=train$infant_id,fold=folds),file.path(DIR$tables,"Classifier_CV_Folds.csv"))
  rf_metric <- metric_values(preds$target,preds$Random_forest)
  write_csv(rf_metric,file.path(DIR$tables,"Random_Forest_Unrounded_Metrics.csv"))
  cm <- table(Predicted=factor(ifelse(preds$Random_forest>=.5,"Phenotype2","Other"),levels=c("Other","Phenotype2")),
              Observed=preds$target)
  assert(sum(cm)==nrow(test) && cm[2,2]==rf_metric$TP && cm[1,1]==rf_metric$TN,
         "Confusion matrix and exported metrics disagree")
  write_csv(as.data.frame.matrix(cm) %>% rownames_to_column("Predicted"),
            file.path(DIR$tables,"Random_Forest_Confusion_Matrix.csv"))
  imp <- tibble(feature=names(final_fits$Random_forest$variable.importance),
                importance=as.numeric(final_fits$Random_forest$variable.importance)) %>% arrange(desc(importance))
  save_table(imp,"sTable_RF_Importance","Random forest training out-of-bag permutation importance",
    "Importance concerns phenotype membership, not a causal treatment-response mechanism. No hard-coded ranking.")
  race_metrics <- bind_rows(lapply(unique(as.character(preds$race)),function(r) {
    dd <- if (is.na(r)) preds[is.na(preds$race),] else preds[!is.na(preds$race)&preds$race==r,]
    if (!nrow(dd)) return(NULL)
    metric_values(dd$target,dd$Random_forest) %>% mutate(Race=ifelse(is.na(r),"Missing",r))
  }))
  save_table(race_metrics,"sTable_Classifier_By_Race","Descriptive held-out classifier performance by race",
    "Small subgroup denominators can produce unstable estimates. Undefined class-specific metrics are left missing. Does not establish fairness of a treatment-selection rule.",TRUE)
  roc_data <- bind_rows(lapply(c("Logistic","Decision_tree","Random_forest"),function(m) {
    rr <- pROC::roc(preds$target,preds[[m]],levels=c("Other","Phenotype2"),direction="<",quiet=TRUE)
    tibble(FPR=1-rr$specificities,TPR=rr$sensitivities,Model=m)
  }))
  pr <- ggplot(roc_data,aes(FPR,TPR,color=Model))+geom_path(linewidth=.8)+
    geom_abline(slope=1,intercept=0,linetype=2,color="grey60")+coord_equal()+
    labs(title="A  Held-out phenotype membership ROC",x="1 - specificity",y="Sensitivity")
  calib <- preds %>% mutate(bin=ntile(Random_forest,5),y=as.numeric(target=="Phenotype2")) %>%
    group_by(bin) %>% summarise(predicted=mean(Random_forest),observed=mean(y),N=n())
  pc <- ggplot(calib,aes(predicted,observed,size=N))+geom_point(color="#2878B5")+
    geom_abline(slope=1,intercept=0,linetype=2)+coord_equal(xlim=c(0,1),ylim=c(0,1))+
    labs(title="B  Held-out random forest calibration",x="Predicted membership probability",y="Observed Phenotype 2 fraction")
  pi <- ggplot(imp,aes(importance,reorder(feature,importance)))+geom_col(fill="#2878B5")+
    scale_y_discrete(labels=VAR_LABELS)+labs(title="C  Random forest feature importance",x="OOB permutation importance",y=NULL)
  save_fig((pr|pc)/pi,"sFigure_Classifier_Validation",12,10,TRUE)
  pdf(file.path(DIR$figures,"supplement","sFigure_Decision_Tree.pdf"),width=10,height=7,useDingbats=FALSE)
  tryCatch(rpart.plot::rpart.plot(final_fits$Decision_tree,extra=104,fallen.leaves=TRUE,
      roundint=FALSE,main="Internal classifier of Phenotype 2 membership"),finally=dev.off())
  save_model(final_fits,"classifier_models_training_only")
  save_model(cv_models,"classifier_cross_validation_models")
  save_model(list(metrics=metrics,predictions=preds,cv_predictions=cv,importance=imp),"classifier_results")
  list(metrics=metrics,importance=imp)
})

#####
#incremental internal outcome prediction beyond simpler variables-----
#####
run_module("Incremental internal prediction comparison", {
  d <- analytic %>% mutate(target=factor(ifelse(extub14==1,"Extubated","Not_extubated"),
                                         levels=c("Not_extubated","Extubated")))
  folds <- fold_assignment(d,CV_FOLDS,77)
  forms <- list(Simple_GA_RI = extub14 ~ hc + gestage + map_fio2x,
    Five_baseline_variables = extub14 ~ hc + gestage + birthwtx + apgar5 + age_rndaysx + map_fio2x,
    Five_plus_phenotype = extub14 ~ hc + gestage + birthwtx + apgar5 + age_rndaysx + map_fio2x + phenotype,
    Five_plus_phenotype_interaction = extub14 ~ hc * phenotype + gestage + birthwtx + apgar5 + age_rndaysx + map_fio2x)
  all_pred <- list()
  for (k in seq_len(CV_FOLDS)) for (nm in names(forms)) {
    tr <- d[folds!=k,]; va <- d[folds==k,]
    f <- glm(forms[[nm]],data=tr,family=binomial(),method=brglm2::brglmFit,
              control=brglm2::brglmControl(type="AS_mean"))
    prob <- as.numeric(predict(f,newdata=va,type="response"))
    all_pred[[paste(k,nm)]] <- tibble(infant_id=va$infant_id,truth=va$target,probability=prob,model=nm,fold=k)
  }
  pp <- bind_rows(all_pred)
  out <- bind_rows(lapply(names(forms),function(nm) {
    a <- pp[pp$model==nm,]; metric_values(a$truth,a$probability,positive="Extubated") %>% mutate(Model=nm)
  }))
  save_table(out,"sTable_Incremental_Prediction","Conditional internal prediction of observed extubation",
    "Common folds and bias-reduced logistic models. Scores describe observed extubation prediction, not treatment-benefit prediction. Center omitted from these transport-oriented predictive comparisons. Phenotypes were discovered in the full cohort, so this is conditional internal assessment, not external or full-pipeline validation. No superiority claim is inferred from small AUC differences.",TRUE)
  write_csv(pp,file.path(DIR$tables,"Outcome_Prediction_OutOfFold.csv"))
  out
})

#####
#longitudinal respiratory data and supplementary figure-----
#####
if (file.exists(LONGITUDINAL_FILE)) {
  run_module("Longitudinal respiratory trajectories", {
    long <- read_source(LONGITUDINAL_FILE)
    req <- c("randcent","rand_id","hc4day","hc4rmodex","hc4rfio2","hc4rpresx","hc4rpco2x")
    assert(all(req %in% names(long)),"Longitudinal columns missing")
    long <- long %>% mutate(across(all_of(req),num),infant_id=make_id(randcent,rand_id))
    assert_unique(long,c("infant_id","hc4day"),"Respiratory snapshots")
    long <- inner_join(long,analytic %>% select(infant_id,phenotype,hc,center,se_unit),
                       by="infant_id",relationship="many-to-one")
    long <- long %>% filter(hc4day>=1,hc4day<=14)
    assert(all(na.omit(long$hc4rfio2)>=0 & na.omit(long$hc4rfio2)<=1),"FiO2 is not fractional")
    avail <- long %>% group_by(hc4day,phenotype) %>% summarise(
      snapshots=n(),invasive=sum(hc4rmodex %in% c(1,2)),
      observed_FiO2=sum(hc4rmodex %in% c(1,2)&!is.na(hc4rfio2)))
    save_table(avail,"sTable_Respiratory_Availability","Respiratory snapshot availability",
      "HC4DAY is study day, with no day 0 observation in these files. HFV=1, conventional ventilation=2. No values imputed after extubation or death.")
    inv <- long %>% filter(hc4rmodex %in% c(1,2)) %>%
      mutate(day=factor(hc4day,levels=sort(unique(long$hc4day))))
    long_plots <- list(); long_models <- list()
    measures <- c(hc4rfio2="Fraction of inspired oxygen",hc4rpresx="Mean airway pressure (cm H2O)",
                  hc4rpco2x="Highest pCO2 (mmHg)")
    for (v in names(measures)) {
      dd <- inv %>% filter(!is.na(.data[[v]])) %>% droplevels()
      ff <- as.formula(paste(v,"~ phenotype * day + hc + (1 | center) + (1 | infant_id)"))
      fit <- lmerTest::lmer(ff,data=dd,REML=TRUE,
        control=lme4::lmerControl(optimizer="bobyqa"))
      check_mer(fit, paste("Trajectory model", v))
      em <- as.data.frame(emmeans::emmeans(fit,~day|phenotype,at=list(hc=.5)))
      low <- intersect(c("lower.CL","asymp.LCL"),names(em))[1]
      high <- intersect(c("upper.CL","asymp.UCL"),names(em))[1]
      assert(!is.na(low)&&!is.na(high),"Trajectory interval columns unavailable")
      em$lower <- em[[low]]; em$upper <- em[[high]]; em$study_day <- num(em$day)
      long_plots[[v]] <- ggplot(em,aes(study_day,emmean,color=phenotype,fill=phenotype))+
        geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.12,color=NA)+geom_line()+geom_point()+
        scale_color_manual(values=COLS)+scale_fill_manual(values=COLS)+
        scale_x_continuous(breaks=sort(unique(dd$hc4day)))+
        labs(x="Study day",y=measures[[v]],color=NULL,fill=NULL)
      write_csv(em,file.path(DIR$tables,paste0("Trajectory_",v,"_Predictions.csv")))
      long_models[[v]] <- fit
    }
    save_fig(wrap_plots(long_plots,ncol=1,guides="collect")+
      plot_annotation(title="Respiratory snapshots while invasively ventilated",
        caption="Available HFV/CV snapshots only. Study day modeled categorically; no baseline extrapolation. Adjusted for treatment. Changing availability after extubation/death limits interpretation; these are not treatment-response trajectories."),
      "sFigure_Respiratory_Trajectories",10,12,TRUE)
    save_model(long_models,"respiratory_snapshot_mixed_models")
    saveRDS(long,file.path(DIR$data,"respiratory_snapshots_linked.rds"))
    invisible(TRUE)
  })
} else record_status("Longitudinal respiratory trajectories","NOT RUN - INPUT UNAVAILABLE","hc04_dash.csv not found")

#####
#follow up timing race reporting and cohort flow-----
#####
age_summary <- analytic %>% filter(!is.na(ndi_source_binary)) %>% group_by(phenotype) %>%
  summarise(N_with_NDI=n(),N_with_recorded_age=sum(!is.na(exam5agecalcx)),
    median_months=median(exam5agecalcx,na.rm=TRUE),
    minimum_months=min(exam5agecalcx,na.rm=TRUE),maximum_months=max(exam5agecalcx,na.rm=TRUE),
    recorded_outside_protocol_22_26=sum(!is.na(exam5agecalcx)&!between(exam5agecalcx,22,26)))
save_table(age_summary,"sTable_Followup_Ages","Recorded corrected ages at neurodevelopmental assessment",
  "Protocol target 22-26 months; observed ages summarized without retrospective exclusion. Deidentification can set outlying ages to missing. Missing recorded age is not proof of absent assessment.")
page <- ggplot(analytic %>% filter(!is.na(ndi_source_binary),!is.na(exam5agecalcx)),aes(exam5agecalcx))+
  annotate("rect",xmin=22,xmax=26,ymin=-Inf,ymax=Inf,fill="#D6E7EF",alpha=.7)+
  geom_histogram(binwidth=1,boundary=.5,fill="#2878B5",color="white")+
  labs(title="Recorded corrected age at follow-up",subtitle="Shading: protocol target window",
       x="Corrected age (months)",y="Infants")
save_fig(page,"sFigure_Followup_Ages",8,5,TRUE)
race_tab <- analytic %>% mutate(race_display=ifelse(is.na(race),"Missing",as.character(race))) %>%
  count(phenotype,race_display) %>% group_by(phenotype) %>% mutate(percent_of_all=100*n/sum(n))
save_table(race_tab,"sTable_Race_Completeness","Available race categories and missingness",
  "No Other-coded observations are invented. Missing is not recoded as Other. Categories reflect the supplied analytic source.")
flow <- tibble(x=c(0,0,1.8),y=c(2,0,1),label=c(
  paste0("Randomized infants\nn = ",nrow(df_pheno)),
  paste0("Baseline-complete phenotype cohort\nn = ",nrow(analytic)),
  paste0("Missing clustering variables\nn = ",sum(!eligible))))
pflow <- ggplot(flow,aes(x,y))+geom_label(aes(label=label),size=4,label.padding=grid::unit(.6,"lines"))+
  annotate("segment",x=0,xend=0,y=1.65,yend=.4,arrow=grid::arrow(length=grid::unit(.15,"inches")))+
  annotate("segment",x=.1,xend=1.2,y=1,yend=1,arrow=grid::arrow(length=grid::unit(.15,"inches")))+
  coord_cartesian(xlim=c(-1,3),ylim=c(-.6,2.7),clip="off")+theme_void()+
  labs(title="Analytic cohort derivation")
save_fig(pflow,"sFigure_Cohort_Flow",9,5,TRUE)

#####
#graphical abstract-----
#####
profile_text <- analytic %>% group_by(phenotype) %>% summarise(n=n(),
  ga=median(gestage),bw=median(birthwtx),apgar=median(apgar5),ri=median(map_fio2x)) %>%
  mutate(label=paste0(phenotype,"  (n=",n,")\n",
     "Median GA ",ga," wk | BW ",round(bw)," g\n",
     "Apgar ",apgar," | Respiratory index ",round(ri)))
pcards <- ggplot(profile_text,aes(1,phenotype,fill=phenotype))+
  geom_tile(width=.98,height=.85,alpha=.16,color="white")+
  geom_text(aes(label=label),size=3.2,lineheight=1.1)+
  scale_fill_manual(values=COLS)+scale_y_discrete(limits=rev)+
  scale_x_continuous(limits=c(.45,1.55))+theme_void()+theme(legend.position="none")+
  labs(title="Baseline phenotypes")
overview <- ggplot()+annotate("text",x=0,y=3.8,label=paste0(nrow(analytic)," infants"),size=8,fontface="bold")+
  annotate("text",x=0,y=2.9,label="Ventilated preterm infants\nSecondary randomized-trial analysis",size=4.5)+
  annotate("text",x=0,y=1.7,label="Five pretreatment variables\nStandardized Ward clustering\nFour exploratory phenotypes",size=4.3)+
  annotate("text",x=0,y=.3,label="Outcome\nSuccessful extubation within 14 days",size=4.3)+
  coord_cartesian(xlim=c(-1,1),ylim=c(-.5,4.4))+theme_void()
pg <- forest(table2_numeric,"Hydrocortisone versus placebo",
  paste("Adjusted log-odds interaction P",fmt_p(primary_global$p_value)))
ga <- (overview|pcards|pg)+plot_layout(widths=c(.9,1.1,1.2))+
  plot_annotation(title="Baseline phenotypes and short-term respiratory outcomes after hydrocortisone",
    caption="Exploratory subgroup analysis. Phenotype-specific estimates and stability require validation. This analysis does not establish a phenotype-guided treatment rule.")
save_fig(ga,"Graphical_Abstract",15,6.5)

#####
#methods output model checks and reproducibility manifest-----
#####
methods_text <- c(
  "ANALYSIS METHODS FOR MANUSCRIPT REVISION - CHECK AGAINST MODULE STATUS",
  paste("Randomized source cohort:",nrow(df_pheno),"; baseline-complete phenotype cohort:",nrow(analytic)),
  "Phenotypes: z-score standardization of five pretreatment variables, Euclidean distance, Ward.D2 linkage, fixed k=4. No treatment or outcome used to derive clusters.",
  "Primary outcome: source successful-extubation-within-14-days variable. No endpoint reconstruction from daily snapshots.",
  "Primary adjusted analysis: logistic regression with treatment-by-phenotype interaction, continuous gestational age, sex, and fixed study-center effects.",
  paste("Covariance:",COVARIANCE_LABEL),
  "Infant identity is the composite of center and infant number. Reused infant numbers are never interpreted as families.",
  "Adjusted treatment-specific risks are standardized over observed covariates within each phenotype. Overall estimate uses a separate adjusted model without phenotype. Delta-method intervals condition on the observed covariate distribution.",
  "RR and RD are derived from the same standardized risks. Conditional ORs are reported separately. Primary global interaction is a robust joint Wald test on the log-odds scale. Modified Poisson interaction is a separate log-RR sensitivity analysis.",
  "BPD endpoints: survival without physiologically defined moderate/severe BPD is distinct from Jensen severity grades among survivors. Death is not a Jensen severity category.",
  "Ventilation endpoint: source duration of mechanical ventilation to 36 weeks among survivors; no 36-minus-days ventilator-free-days calculation.",
  "NDI endpoint: source NDI3 flag with NDI3CAT coding audit; NDI_MOD_SEV is analyzed separately. Protocol assessment target is 22-26 months corrected age; actual recorded ages and missingness are reported.",
  "Ordinal analysis: exploratory proportional-odds CLM with model-based covariance; assumptions diagnosed, with robust threshold-specific analyses also reported. This is a changed analysis and must not be described as the old ordinal GEE.",
  "Apgar sensitivity: identical original cohort, four remaining standardized features, independent alternative labels, ARI/cross-tabulation; merged original Phenotypes 2+3 also evaluated.",
  "Clustering robustness: silhouette, candidate cluster counts, bootstrap Jaccard, cross-algorithm agreement, subsampling, and consensus. Subsampling is not a missing-data-mechanism analysis.",
  "Classification: binary Phenotype 2 membership, fixed-parameter models, 70/30 group-aware split, training-only five-fold CV, untouched held-out set, saved probabilities, real RF permutation importance. Logistic membership model uses mean bias reduction. Threshold 0.5.",
  "Classifier labels derive from the full-cohort clustering. Performance reflects conditional internal membership reproduction, not independent phenotyping validation or individualized treatment benefit.",
  "Longitudinal snapshots: observed invasive-ventilation modes on recorded study days only, categorical time, treatment adjustment, center and unique-infant random intercepts. No day-0 extrapolation or postextubation imputation. Availability is informative.",
  "All new analyses are exploratory, conducted in response to peer review, with no multiplicity adjustment. Do not claim a prospective analysis plan without documentation.",
  "NIH/NICHD comparator is performed only with externally verified compatible pretreatment risk predictions; absence is explicitly logged.",
  "All manuscript numeric results, table footnotes, figure legends, references to GEE, and response-letter claims must be reconciled with these regenerated outputs."
)
writeLines(methods_text,file.path(OUT,"Methods_and_Interpretation_Notes.txt"))
writeLines(capture.output(sessionInfo()),file.path(DIR$logs,"sessionInfo.txt"))
write_csv(tibble(package=required_packages,
                 version=vapply(required_packages,function(p)as.character(packageVersion(p)),character(1))),
          file.path(DIR$logs,"package_versions.csv"))
writeLines(unique(warning_log),file.path(DIR$logs,"warnings.txt"))
saveRDS(list(primary=primary,overall=overall,secondary=secondary,
             sensitivity_models=sensitivity_models,classification=classification,
             ordinal=ordinal_results,ventilation=ventilation_results,
             module_status=module_status,configuration=list(seed=2024,run_id=RUN_ID,
               family_linkage=HAS_FAMILY_LINKAGE,covariance=COVARIANCE_LABEL,
               bootstrap=BOOTSTRAP_REPS,consensus=CONSENSUS_REPS,cv_folds=CV_FOLDS)),
        file.path(DIR$models,"analysis_bundle.rds"))
failed <- module_status %>% filter(status=="FAILED")
notrun <- module_status %>% filter(grepl("NOT RUN",status))
writeLines(c(
  if (nrow(failed)) "COMPLETED WITH FAILED MODULES - REVIEW REQUIRED" else "CORE ANALYSIS COMPLETED - REVIEW OUTPUTS BEFORE SUBMISSION",
  paste("Failed modules:",nrow(failed)),paste("Modules requiring unavailable inputs:",nrow(notrun)),
  "Figures and tables require visual/scientific review. A script completing does not certify manuscript accuracy.",
  if (!HAS_FAMILY_LINKAGE) "Sibling dependence could not be modeled because no genuine family linkage was supplied." else "Genuine supplied family linkage used where documented.",
  paste("Output folder:",OUT)),file.path(OUT,"RUN_STATUS.txt"))
paths <- list.files(OUT,recursive=TRUE,full.names=TRUE)
manifest <- tibble(file=substring(paths,nchar(OUT)+2),bytes=file.info(paths)$size,
                   md5=unname(tools::md5sum(paths)))
write_csv(manifest,file.path(OUT,"output_manifest.csv"))
message("Output folder: ",OUT)
message("Failed modules: ",nrow(failed),"; unavailable-input modules: ",nrow(notrun))
message("Read RUN_STATUS.txt, logs/module_status.csv, and logs/warnings.txt before using outputs.")
if (nrow(failed)) warning("Some modules failed; results are incomplete. See logs/module_status.csv.")
