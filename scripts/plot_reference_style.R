# Four-panel supplementary figure styled to match the supplied baseline figure.
# Actual values are calculated from the analysis results and original data.
# Requires ggplot2 and patchwork. No synthetic results or fixed P values.
plot_apgar_reference_style <- function(results, df_pheno,
                                       output_dir = "apgar_sensitivity_results",
                                       figure_number = "SXX") {
  for (pkg in c("ggplot2", "patchwork", "mclust"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
  ids <- results$assignments
  if (anyDuplicated(df_pheno$infant_row_id) || anyDuplicated(ids$infant_row_id))
    stop("Infant IDs must be unique.")
  ii <- match(ids$infant_row_id, df_pheno$infant_row_id)
  if (anyNA(ii)) stop("Some analysis infants are absent from df_pheno.")
  d <- as.data.frame(df_pheno[ii, , drop = FALSE])
  features <- c("GESTAGE", "BIRTHWTX", "APGAR5", "age_rndaysX", "map_fio2X")
  if (!all(features %in% names(d))) stop("Missing baseline variables.")
  if (anyNA(d[features])) stop("Missing baseline values in matched analysis cohort.")
  d$original <- ids$original
  d$no_apgar <- ids$no_apgar
  new_levels <- paste("No-Apgar cluster", 1:4)
  old_levels <- paste("Phenotype", 1:4)
  colors <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), new_levels)
  four_features <- setdiff(features, "APGAR5")
  reproduced <- cutree(hclust(dist(scale(as.matrix(d[four_features]))),
                              method = "ward.D2"), k = 4)
  if (abs(mclust::adjustedRandIndex(reproduced, d$no_apgar) - 1) > 1e-10)
    stop("Input data do not reproduce the supplied no-Apgar assignments.")
  n_by <- table(factor(d$no_apgar, levels = new_levels))
  theme_pub <- function(base_size = 11) {
    ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = base_size + 1,
                                           margin = ggplot2::margin(b = 8)),
        plot.subtitle = ggplot2::element_text(size = base_size - 2, colour = "grey30"),
        panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = .3),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "grey35", linewidth = .3),
        axis.text = ggplot2::element_text(colour = "grey25"),
        plot.background = ggplot2::element_rect(fill = "grey98", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        legend.background = ggplot2::element_rect(fill = "grey98", colour = NA),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(12, 14, 12, 12))
  }
  # A: PCA uses only the FOUR variables actually used in alternative clustering.
  pc <- prcomp(d[four_features], center = TRUE, scale. = TRUE)
  variance <- summary(pc)$importance[2, 1:2] * 100
  a <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], cluster = d$no_apgar)
  # Draw ellipses only when a nonsingular within-cluster covariance is available.
  ellipse_groups <- new_levels[vapply(new_levels, function(g) {
    x <- a[a$cluster == g, c("PC1", "PC2"), drop = FALSE]
    nrow(x) >= 4 && all(is.finite(cov(x))) && det(cov(x)) > 1e-12
  }, logical(1))]
  pa <- ggplot2::ggplot(a, ggplot2::aes(PC1, PC2, colour = cluster)) +
    ggplot2::geom_point(alpha = .55, size = 1.5) +
    ggplot2::stat_ellipse(data = a[a$cluster %in% ellipse_groups, , drop = FALSE],
                           type = "norm", level = .95, linewidth = .8, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = colors, breaks = new_levels,
      labels = paste0(new_levels, " (n=", as.integer(n_by), ")"), drop = FALSE) +
    ggplot2::labs(title = "A. Principal component analysis without Apgar",
      subtitle = "Projection of the four clustering variables",
      x = sprintf("PC1 (%.1f%% variance)", variance[1]),
      y = sprintf("PC2 (%.1f%% variance)", variance[2])) +
    theme_pub() + ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  # B: include Apgar descriptively, explicitly labelled as excluded from clustering.
  zz <- scale(as.matrix(d[features]))
  labs <- c(GESTAGE = "Gestational age", BIRTHWTX = "Birth weight",
            APGAR5 = "Apgar 5 min\n(not used)", age_rndaysX = "Age at treatment",
            map_fio2X = "MAP x FiO2")
  hh <- do.call(rbind, lapply(new_levels, function(g) {
    m <- colMeans(zz[d$no_apgar == g, , drop = FALSE])
    data.frame(cluster = g, feature = names(m), z = unname(m))
  }))
  hh$cluster <- factor(hh$cluster, levels = rev(new_levels))
  hh$feature <- factor(hh$feature, levels = c("age_rndaysX", "APGAR5", "BIRTHWTX", "GESTAGE", "map_fio2X"))
  pb <- ggplot2::ggplot(hh, ggplot2::aes(feature, cluster, fill = z)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", z),
                       colour = abs(z) > 1), fontface = "bold", size = 3.1) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = "white"), guide = "none") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, name = "Z-score") +
    ggplot2::scale_x_discrete(labels = labs) +
    ggplot2::labs(title = "B. Standardized baseline features",
      subtitle = "Apgar is shown descriptively; it was excluded from clustering", x = NULL, y = NULL) +
    theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                                 panel.grid = ggplot2::element_blank(),
                                 legend.title = ggplot2::element_text(size = 9))

  # C: original-to-alternative membership, normalized WITHIN original phenotype.
  tab <- table(factor(d$original, levels = old_levels),
               factor(d$no_apgar, levels = new_levels))
  ct <- as.data.frame(tab)
  names(ct) <- c("original", "alternative", "n")
  ct$pct <- 100 * as.numeric(prop.table(tab, margin = 1))
  ct$original <- factor(ct$original, levels = rev(old_levels))
  ct$label <- sprintf("%d\n(%.1f%%)", ct$n, ct$pct)
  ari <- mclust::adjustedRandIndex(d$original, d$no_apgar)
  pcross <- ggplot2::ggplot(ct, ggplot2::aes(alternative, original, fill = pct)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = pct > 55), size = 3.2) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "white"), guide = "none") +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", limits = c(0, 100),
                                 name = "Row %") +
    ggplot2::scale_x_discrete(labels = setNames(paste("Cluster", 1:4), new_levels)) +
    ggplot2::labs(title = "C. Correspondence with original phenotypes",
      subtitle = sprintf("Adjusted Rand index = %.3f; cells show n (row %%)", ari),
      x = "Clustering without Apgar", y = NULL) +
    theme_pub() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                 legend.title = ggplot2::element_text(size = 9))

  # D: four no-Apgar treatment effects + original groups 2/3 pooled, in grey.
  f <- results$tables$no_apgar_RRs[c("group", "RR", "lower_95", "upper_95")]
  pooled <- results$tables$pooled_phenotypes_2_3_RR[names(f)]
  pooled$group <- "Original phenotypes 2+3 (pooled)"
  f <- rbind(f, pooled)
  f$group <- factor(f$group, levels = rev(as.character(f$group)))
  f$label <- sprintf("%.2f (%.2f–%.2f)", f$RR, f$lower_95, f$upper_95)
  if (any(!is.finite(as.matrix(f[c("RR", "lower_95", "upper_95")]))))
    stop("Nonfinite treatment estimates; inspect models before plotting.")
  if (any(f$lower_95 <= 0)) stop("Nonpositive lower bound for log RR plot.")
  colors_f <- c(colors, `Original phenotypes 2+3 (pooled)` = "#555555")
  interaction_p <- results$tables$interaction_tests$p_value[
    results$tables$interaction_tests$analysis == "no_apgar"]
  ptext <- if (interaction_p < .001) "<0.001" else sprintf("= %.3f", interaction_p)
  max_ci <- max(f$upper_95)
  annotation_x <- max_ci * 1.35
  pd <- ggplot2::ggplot(f, ggplot2::aes(RR, group, colour = group)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(x = lower_95, xend = upper_95, yend = group), linewidth = .8) +
    ggplot2::geom_point(size = 3, shape = 18) +
    ggplot2::geom_text(ggplot2::aes(x = annotation_x, label = label),
                       hjust = 0, size = 2.8, colour = "grey25") +
    ggplot2::scale_colour_manual(values = colors_f, guide = "none") +
    ggplot2::scale_y_discrete(labels = function(x) sub("Original phenotypes 2\\+3 \\(pooled\\)",
                                                     "Original 2+3\n(pooled)", x)) +
    ggplot2::scale_x_log10(limits = c(min(.8, min(f$lower_95) * .85), annotation_x * 4)) +
    ggplot2::labs(title = "D. Hydrocortisone treatment effects",
      subtitle = paste0("No-Apgar treatment-by-cluster interaction P ", ptext),
      x = "Risk ratio (95% CI) for successful extubation within 14 days", y = NULL) +
    theme_pub() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  figure <- patchwork::wrap_plots(pa, pb, pcross, pd, ncol = 2) +
    patchwork::plot_annotation(
      title = paste0("Supplementary Figure ", figure_number,
                     ". Sensitivity of phenotypes and treatment effects to Apgar exclusion"),
      caption = "Alternative cluster numbers and colors do not imply correspondence with original phenotype labels.",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.caption = ggplot2::element_text(size = 9, hjust = 0),
        plot.background = ggplot2::element_rect(fill = "grey98", colour = NA)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(output_dir, "Supplementary_Figure.pdf"), figure,
    width = 15, height = 11.5, units = "in", device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(file.path(output_dir, "Supplementary_Figure.tiff"), figure,
    width = 15, height = 11.5, units = "in", dpi = 600, compression = "lzw", bg = "grey98")
  ggplot2::ggsave(file.path(output_dir, "Supplementary_Figure.png"), figure,
    width = 15, height = 11.5, units = "in", dpi = 300, bg = "grey98")
  writeLines(c(
    paste0("Supplementary Figure ", figure_number,
      ". Sensitivity of phenotypes and treatment effects to exclusion of 5-minute Apgar score."),
    "A, Principal component analysis of standardized gestational age, birth weight, age at randomization, and baseline respiratory severity. Colors denote independently derived no-Apgar clusters; ellipses are 95% normal-data ellipses, not confidence intervals for means or cluster boundaries.",
    "B, Mean standardized features by no-Apgar cluster. Apgar score is displayed descriptively and was not used to derive these clusters. Z-scores are standardized across the common clustering cohort.",
    sprintf("C, Original versus no-Apgar memberships. Cells show infant counts and percentages within each original phenotype. Adjusted Rand index = %.3f.", ari),
    "D, Hydrocortisone versus placebo risk ratios for successful extubation within 14 days. Diamonds denote RRs; horizontal lines denote 95% CIs. The grey estimate combines original Phenotypes 2 and 3. Estimates use unadjusted modified Poisson models with HC1 sandwich variance clustered by the specified randomization unit. The interaction P value is the joint robust Wald test across the four no-Apgar clusters on the log-RR scale; it excludes the separately displayed pooled group.",
    "Both partitions use the same complete-case infants, standardized inputs, Euclidean distance, Ward.D2 linkage and four clusters. Outcome analyses exclude missing extubation outcomes. Tests are exploratory and not multiplicity-adjusted. New cluster numbers and colors do not establish correspondence with original phenotype labels."
  ), file.path(output_dir, "Supplementary_Figure_legend.txt"))
  invisible(figure)
}
