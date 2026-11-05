# Create the primary five-variable phenotype dataset for the sensitivity pipeline.
# Base R only; no additional packages required.
# Raw coding from supplied script: treatmnt 1=placebo, 2=hydrocortisone;
# extube14 0=no, 1=yes. Numeric special missing-value codes must be recoded to
# NA according to the codebook BEFORE calling this function.
# Retains all original columns, adds/updates treat, extub14, phenotype, and
# infant_row_id. Does not derive secondary clinical outcomes.
make_df_pheno <- function(analysis_dash, expected_n = NULL,
                          unit_columns = c("randcent", "rand_id")) {
  features <- c("GESTAGE", "BIRTHWTX", "APGAR5", "age_rndaysX", "map_fio2X")
  required <- unique(c("treatmnt", "extube14", features, unit_columns))
  missing_columns <- setdiff(required, names(analysis_dash))
  if (length(missing_columns))
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  d <- as.data.frame(analysis_dash)
  numeric_columns <- c("treatmnt", "extube14", features)
  if (!all(vapply(d[numeric_columns], is.numeric, logical(1))))
    stop("Treatment, outcome and clustering variables must be numeric; check coding before import.")
  # Save source-row IDs before filtering. Never use shared randomization IDs
  # to join infant-level phenotype assignments.
  if (!"infant_row_id" %in% names(d)) d$infant_row_id <- seq_len(nrow(d))
  if (anyNA(d$infant_row_id) || anyDuplicated(d$infant_row_id))
    stop("infant_row_id must be nonmissing and unique per infant.")
  d <- d[!is.na(d$treatmnt) & d$treatmnt %in% c(1, 2), , drop = FALSE]
  rownames(d) <- NULL
  if (!nrow(d)) stop("No randomized infants with treatmnt 1 or 2.")
  if (!is.null(expected_n) && nrow(d) != expected_n)
    stop("Expected ", expected_n, " randomized infants; found ", nrow(d), ".")
  if (!all(na.omit(d$extube14) %in% c(0, 1)))
    stop("extube14 must contain 0, 1, or NA.")
  if (length(unit_columns) && anyNA(d[unit_columns]))
    stop("Missing randomization identifiers; resolve these before treatment inference.")
  if (length(unit_columns)) {
    unit <- do.call(interaction, c(lapply(d[unit_columns], as.character),
                                   list(drop = TRUE, lex.order = TRUE)))
    if (any(vapply(split(d$treatmnt, unit), function(x) length(unique(x)) != 1L,
                   logical(1))))
      stop("Conflicting treatment codes within a randomization unit; check unit_columns.")
  }
  d$treat <- factor(d$treatmnt, levels = c(1, 2),
                     labels = c("Placebo", "Hydrocortisone"))
  d$extub14 <- factor(d$extube14, levels = c(0, 1), labels = c("No", "Yes"))
  eligible <- complete.cases(d[features])
  if (sum(eligible) <= 4L) stop("Too few complete cases for four clusters.")
  x <- as.matrix(d[eligible, features, drop = FALSE])
  if (any(!is.finite(x))) stop("Clustering values must be finite or NA.")
  feature_sd <- apply(x, 2, sd)
  if (any(!is.finite(feature_sd) | feature_sd == 0))
    stop("At least one clustering feature has no variation.")
  # Same scaling, feature order, distance, linkage and k as supplied primary code.
  # hclust is deterministic; a random seed is unnecessary.
  hc <- hclust(dist(scale(x), method = "euclidean"), method = "ward.D2")
  assignments <- rep(NA_integer_, nrow(d))
  assignments[eligible] <- cutree(hc, k = 4)
  d$phenotype <- factor(assignments, levels = 1:4,
                         labels = paste("Phenotype", 1:4))
  message("Randomized infants retained: ", nrow(d),
          "; assigned phenotypes: ", sum(eligible),
          "; missing baseline features: ", sum(!eligible), ".")
  print(table(d$phenotype, useNA = "ifany"))
  d
}

# Example after importing your raw analysis_dash dataset:
# source("make_df_pheno.R")
# df_pheno <- make_df_pheno(analysis_dash, expected_n = 800)
#
# Your original figure showed counts 59, 344, 216, 141. Verify these against
# the output before calling groups 'original Phenotypes 1-4'. Input row order
# and cohort/preprocessing changes can affect cluster numbering or membership.
#
# Continue using the previous supplementary pipeline:
# source("apgar_sensitivity_analysis.R")
# source("export_supplement.R")
# source("plot_reference_style.R")
# results <- run_apgar_sensitivity(df_pheno)
# export_apgar_supplement(results, df_pheno = df_pheno,
#                         table_number = "S13", figure_number = "SXX")
