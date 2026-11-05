# Source after apgar_sensitivity_analysis.R; requires flextable, officer, ggplot2.
# Journal-neutral formatting; inspect the actual Word/PDF output before submission.
export_apgar_supplement <- function(results, output_dir = "apgar_sensitivity_results",
                                    table_number = "SXX", figure_number = "SXX",
                                    make_figure = TRUE, df_pheno = NULL) {
  pkgs <- c("flextable", "officer", if (make_figure) "ggplot2")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)
  t <- results$tables
  required <- c("primary_RRs", "merged_phenotypes_RRs", "pooled_phenotypes_2_3_RR",
                "no_apgar_RRs", "interaction_tests", "cohort_and_ARI", "cluster_cross_tabulation")
  if (!all(required %in% names(t))) stop("Use results from run_apgar_sensitivity().")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  p_fmt <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))
  add_analysis <- function(x, label) { x$analysis <- label; x }
  # Use the separately fitted pooled model for the merged 2+3 row.
  merged <- t$merged_phenotypes_RRs
  j <- match("Phenotypes 2+3", merged$group)
  for (v in c("RR", "lower_95", "upper_95", "p_value"))
    merged[j, v] <- t$pooled_phenotypes_2_3_RR[[v]][1]
  z <- rbind(add_analysis(t$primary_RRs, "A. Original phenotypes"),
             add_analysis(merged, "B. Phenotypes 2 and 3 merged"),
             add_analysis(t$no_apgar_RRs, "C. Clustering without Apgar"))
  # Interaction P belongs to the whole partition, never to a single cluster.
  ip <- setNames(t$interaction_tests$p_value, t$interaction_tests$analysis)
  z$interaction_p <- ""
  for (key in c("original", "merged", "no_apgar")) {
    label <- c(original = "A. Original phenotypes", merged = "B. Phenotypes 2 and 3 merged",
               no_apgar = "C. Clustering without Apgar")[[key]]
    z$interaction_p[which(z$analysis == label)[1]] <- p_fmt(ip[[key]])
  }
  pub <- data.frame(
    Analysis = z$analysis,
    Group = z$group,
    Placebo = sprintf("%d/%d (%.1f)", z$placebo_events, z$placebo_n, 100 * z$placebo_risk),
    Hydrocortisone = sprintf("%d/%d (%.1f)", z$hydrocortisone_events,
                             z$hydrocortisone_n, 100 * z$hydrocortisone_risk),
    RR_CI = sprintf("%.2f (%.2f–%.2f)", z$RR, z$lower_95, z$upper_95),
    P = p_fmt(z$p_value), Interaction_P = z$interaction_p,
    check.names = FALSE)
  # Keep unrounded numerical estimates in a separate machine-readable table.
  write.csv(z, file.path(output_dir, "Supplement_Treatment_Effects_numeric.csv"), row.names = FALSE)
  write.csv(pub, file.path(output_dir, "Supplement_Treatment_Effects_formatted.csv"), row.names = FALSE)
  ari <- t$cohort_and_ARI$value[t$cohort_and_ARI$measure == "Original versus no-Apgar ARI"]
  cohort_n <- t$cohort_and_ARI$value[t$cohort_and_ARI$measure == "Same-infant clustering cohort"]
  analysis_n <- t$cohort_and_ARI$value[t$cohort_and_ARI$measure == "Treatment-analysis infants"]
  notes <- c(
    "Values are successful extubations within 14 days / infants with observed outcome (%).",
    "RR, risk ratio; CI, confidence interval. Hydrocortisone is compared with placebo.",
    "RRs and 95% CIs use unadjusted modified Poisson regression with HC1 sandwich variance clustered by the specified randomization unit; P values use two-sided normal-approximation tests.",
    "Interaction P values are joint robust Wald tests of treatment-by-group interaction on the log-RR scale, one per partition. They are shown once per analysis and do not test treatment within an individual cluster. P values are not adjusted for multiple testing.",
    "The merged 2+3 row uses a separate pooled treatment model; other group rows use partition-specific interaction models.",
    sprintf("Both clustering solutions include the same %d infants; %d have observed primary outcomes. Four clusters were derived using standardized variables, Euclidean distance, and Ward.D2 linkage.", cohort_n, analysis_n),
    sprintf("Adjusted Rand index comparing original and no-Apgar partitions: %.3f. New cluster numbers do not imply correspondence with original phenotype numbers.", ari),
    "These sensitivity analyses were performed in response to peer review. Statistical significance within one group alone does not establish treatment-effect heterogeneity."
  )
  style_table <- function(ft, widths) {
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::padding(ft, padding = 4, part = "all")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::align(ft, align = "left", part = "all")
    for (i in seq_along(widths)) ft <- flextable::width(ft, j = i, width = widths[i])
    flextable::set_table_properties(ft, layout = "fixed", opts_word = list(split = FALSE))
  }
  ft <- flextable::flextable(pub)
  ft <- flextable::set_header_labels(ft, Analysis = "Analysis", Group = "Group",
    Placebo = "Placebo\nn/N (%)", Hydrocortisone = "Hydrocortisone\nn/N (%)",
    RR_CI = "RR (95% CI)", P = "Treatment P", Interaction_P = "Interaction P")
  ft <- style_table(ft, c(1.35, 1.25, 1.1, 1.2, 1.35, .75, .85))
  ft <- flextable::merge_v(ft, j = "Analysis")
  ft <- flextable::valign(ft, j = "Analysis", valign = "top")
  doc <- officer::read_docx()
  doc <- officer::body_set_default_section(doc, officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = .55, bottom = .55, left = .55, right = .55)))
  doc <- officer::body_add_par(doc, paste0("Supplementary Table ", table_number,
    ". Sensitivity of treatment effects to phenotype definition"), style = "heading 1")
  doc <- flextable::body_add_flextable(doc, ft)
  note_fp <- officer::fp_text(font.size = 9, font.family = "Times New Roman")
  for (note in notes) doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(note, note_fp)))
  doc <- officer::body_add_break(doc)
  doc <- officer::body_add_par(doc, paste0("Supplementary Table ", table_number,
    " continued. Correspondence of original and no-Apgar clusters"), style = "heading 1")
  cross <- t$cluster_cross_tabulation
  names(cross)[1] <- "Original phenotype"
  counts <- as.matrix(cross[-1])
  cross$Total <- rowSums(counts)
  total <- data.frame("Original phenotype" = "Total", as.list(colSums(counts)),
                       Total = sum(counts), check.names = FALSE)
  names(total) <- names(cross)
  cross <- rbind(cross, total)
  ft_cross <- style_table(flextable::flextable(cross), c(1.6, rep(1.25, 4), 1.25))
  doc <- flextable::body_add_flextable(doc, ft_cross)
  doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(
    sprintf("Counts include all %d infants in the clustering cohort, regardless of primary-outcome availability. Adjusted Rand index = %.3f. Cluster labels are arbitrary.", cohort_n, ari), note_fp)))
  print(doc, target = file.path(output_dir, "Supplementary_Table.docx"))
  writeLines(notes, file.path(output_dir, "Supplementary_Table_notes.txt"))

  if (make_figure) {
    if (is.null(df_pheno)) stop("Supply df_pheno for the four-panel reference-style figure.")
    if (!exists("plot_apgar_reference_style", mode = "function"))
      stop("Source plot_reference_style.R before exporting the figure.")
    plot_apgar_reference_style(results, df_pheno, output_dir, figure_number)
  }
  message("Supplement exports saved to ", normalizePath(output_dir),
          ". Inspect Word pagination and figure labels before submission.")
  invisible(pub)
}
