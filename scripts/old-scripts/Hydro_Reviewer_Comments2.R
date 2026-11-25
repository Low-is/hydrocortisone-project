############################################################
# FINAL ANALYSIS OUTPUTS (REVIEWER-COMPLIANT)
############################################################

library(sandwich)
library(lmtest)
library(tidyverse)
library(gtsummary)
library(flextable)
library(sandwich)
library(lmtest)
library(ggplot2)



############################################################
# TABLE 1 — BASELINE (NO P-VALUES)
############################################################

table1 <- df_pheno %>%
  select(
    phenotype,
    GESTAGE,
    BIRTHWTX,
    male,
    sga,
    ANTESTER,
    MULTBRTH,
    APGAR5,
    age_rndaysX,
    map_fio2X
  ) %>%
  tbl_summary(
    by = phenotype,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    missing = "no"
  ) %>%
  bold_labels() %>%
  modify_caption(
    "**Table 1. Baseline characteristics by phenotype (descriptive only; no hypothesis testing performed).**"
  )

table1 %>%
  as_flex_table() %>%
  save_as_docx(path = "Table_1.docx")

############################################################
# TABLE 2 — TREATMENT EFFECT (COUNTS + RR + RD + OR + RD CI)
############################################################

library(tidyverse)
library(flextable)
library(sandwich)
library(lmtest)

df_tx <- df_pheno %>%
  filter(!is.na(phenotype), !is.na(treat), !is.na(extub14)) %>%
  mutate(extub_bin = as.numeric(extub14 == "Yes"))

############################################################
# COUNTS + RATES
############################################################

counts <- df_tx %>%
  group_by(phenotype, treat) %>%
  summarise(events = sum(extub_bin), total = n(), .groups = "drop") %>%
  pivot_wider(names_from = treat, values_from = c(events, total)) %>%
  mutate(
    pct_placebo = 100 * events_Placebo / total_Placebo,
    pct_hc = 100 * events_Hydrocortisone / total_Hydrocortisone
  )

############################################################
# RR (Poisson robust)
############################################################

rr <- df_tx %>%
  group_by(phenotype) %>%
  group_modify(~{
    model <- glm(extub_bin ~ treat, family = poisson(link = "log"), data = .x)
    est <- coeftest(model, vcov = sandwich)[2,]
    
    tibble(
      RR = exp(est[1]),
      RR_LCL = exp(est[1] - 1.96 * est[2]),
      RR_UCL = exp(est[1] + 1.96 * est[2])
    )
  })

############################################################
# OR (logistic)
############################################################

or <- df_tx %>%
  group_by(phenotype) %>%
  group_modify(~{
    model <- glm(extub_bin ~ treat, family = binomial, data = .x)
    est <- summary(model)$coef[2,]
    
    tibble(
      OR = exp(est[1]),
      OR_LCL = exp(est[1] - 1.96 * est[2]),
      OR_UCL = exp(est[1] + 1.96 * est[2]),
      p_value = est[4]
    )
  })

############################################################
# RD (ABSOLUTE RISK DIFFERENCE + 95% CI)
############################################################

rd <- df_tx %>%
  group_by(phenotype) %>%
  group_modify(~{
    p_placebo <- mean(.x$extub_bin[.x$treat == "Placebo"])
    p_hc <- mean(.x$extub_bin[.x$treat == "Hydrocortisone"])
    
    rd <- p_hc - p_placebo
    
    # Standard error for difference in proportions
    n1 <- sum(.x$treat == "Placebo")
    n2 <- sum(.x$treat == "Hydrocortisone")
    
    se <- sqrt(
      (p_placebo * (1 - p_placebo) / n1) +
        (p_hc * (1 - p_hc) / n2)
    )
    
    tibble(
      RD = rd,
      RD_LCL = rd - 1.96 * se,
      RD_UCL = rd + 1.96 * se
    )
  })

############################################################
# COMBINE ALL
############################################################

table2 <- counts %>%
  left_join(rr, by = "phenotype") %>%
  left_join(or, by = "phenotype") %>%
  left_join(rd, by = "phenotype") %>%
  mutate(
    n_placebo = sprintf("%d/%d", events_Placebo, total_Placebo),
    n_hc = sprintf("%d/%d", events_Hydrocortisone, total_Hydrocortisone),
    
    pct_placebo_fmt = sprintf("%.0f%%", pct_placebo),
    pct_hc_fmt = sprintf("%.0f%%", pct_hc),
    
    RR_ci = sprintf("%.2f (%.2f–%.2f)", RR, RR_LCL, RR_UCL),
    
    RD_ci = sprintf("%.1f%% (%.1f to %.1f)",
                    100 * RD,
                    100 * RD_LCL,
                    100 * RD_UCL),
    
    OR_ci = sprintf("%.2f (%.2f–%.2f)", OR, OR_LCL, OR_UCL),
    
    p_fmt = ifelse(p_value < 0.001, "<0.001", sprintf("%.2f", p_value))
  )

############################################################
# INTERACTION P
############################################################

interaction_model <- glm(
  extub_bin ~ phenotype * treat,
  family = binomial,
  data = df_tx
)

interaction_p <- anova(interaction_model, test = "LRT")$`Pr(>Chi)`[3]

############################################################
# FINAL TABLE FORMAT (JOURNAL STYLE)
############################################################

final_table <- table2 %>%
  select(
    phenotype,
    n_placebo,
    n_hc,
    pct_placebo_fmt,
    pct_hc_fmt,
    RR_ci,
    RD_ci,
    OR_ci,
    p_fmt
  ) %>%
  setNames(c(
    "Phenotype",
    "n/N Placebo",
    "n/N HC",
    "Extubation % (Placebo)",
    "Extubation % (HC)",
    "Risk Ratio (95% CI)",
    "Absolute Risk Difference (95% CI)",
    "Odds Ratio (95% CI)",
    "p-value"
  ))

############################################################
# EXPORT
############################################################

final_table %>%
  flextable() %>%
  set_caption(
    paste0(
      "Table 2. Hydrocortisone vs placebo by baseline phenotype. ",
      "Global treatment-by-phenotype interaction P = ",
      round(interaction_p, 2)
    )
  ) %>%
  autofit() %>%
  save_as_docx(path = "Table_2_Final.docx")


############################################################
# FIGURE 3 — FOREST PLOT (RR)
############################################################

p_forest <- ggplot(table2, aes(x = RR, y = phenotype)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(xmin = RR_LCL, xmax = RR_UCL),
    height = 0.2,
    orientation = "y"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_classic() +
  labs(
    title = paste0(
      "Hydrocortisone Effect by Phenotype (Interaction P = ",
      round(interaction_p, 2), ")"
    ),
    x = "Risk Ratio (log scale)",
    y = NULL
  )

ggsave("Figure_3.png", p_forest, width = 8, height = 6)

############################################################
# FLOW DIAGRAM
############################################################

library(ggplot2)
library(grid)
library(tibble)

flow_df <- tibble(
  label = c(
    "Randomized\n(n = 800)",
    "Complete baseline data\n(n = 760)",
    "Phenotype derivation\n(n = 760)",
    "Included in phenotype-based analyses\n(n = 760)",
    
    "Excluded from phenotype derivation\n(n = 40)",
    
    "Phenotype-stratified analyses",
    "Extubation\n(n = 760)",
    "BPD\n(n = 704)",
    "Neurodevelopment\n(n = 606)"
  ),
  
  x = c(
    1, 1, 1, 1,
    1.9,     # excluded box
    1,
    0.4, 1, 1.6
  ),
  
  y = c(
    10, 9, 8, 7,
    8.6,
    6,
    5, 5, 5
  )
)

p_flow <- ggplot(flow_df, aes(x = x, y = y)) +
  
  geom_label(
    aes(label = label),
    size = 3.6,
    label.size = 0.4,
    label.padding = unit(0.35, "lines")
  ) +
  
  # 🔥 vertical arrows (fixed)
  annotate("segment", x=1, xend=1, y=9.7, yend=9.5,
           arrow=arrow(length=unit(0.2,"cm"))) +
  annotate("segment", x=1, xend=1, y=8.6, yend=8.2,
           arrow=arrow(length=unit(0.2,"cm"))) +
  annotate("segment", x=1, xend=1, y=7.6, yend=7.2,
           arrow=arrow(length=unit(0.2,"cm"))) +
  annotate("segment", x=1, xend=1, y=6.7, yend=6.2,
           arrow=arrow(length=unit(0.2,"cm"))) +
  
  # 🔥 proper fork to excluded (clean diagonal)
  annotate("segment", x=1, xend=1.9, y=9, yend=8.6,
           arrow=arrow(length=unit(0.2,"cm"))) +
  
  # 🔥 branching arrows
  annotate("segment", x=1, xend=0.4, y=6, yend=5.3,
           arrow=arrow(length=unit(0.2,"cm"))) +
  annotate("segment", x=1, xend=1, y=6, yend=5.3,
           arrow=arrow(length=unit(0.2,"cm"))) +
  annotate("segment", x=1, xend=1.6, y=6, yend=5.3,
           arrow=arrow(length=unit(0.2,"cm"))) +
  
  theme_void() +
  xlim(0, 2.2) +
  ylim(4.5, 10.5) +
  theme(
    plot.background = element_rect(fill="white", color=NA),
    panel.background = element_rect(fill="white", color=NA)
  )

ggsave(
  "Figure_Flow_Final.png",
  p_flow,
  width = 7,
  height = 9,
  dpi = 300,
  bg = "white"
)

############################################################
# SUPPLEMENTAL — SURVIVAL WITHOUT BPD
############################################################

# 🔥 Step 1: build correctly (REMOVE non-phenotype patients)
supp_bpd <- df_pheno %>%
  filter(!is.na(phenotype)) %>%   # <-- CRITICAL FIX
  group_by(phenotype, treat) %>%
  summarise(
    events = sum(surv36_no_bpd == "Yes"),
    total = n(),
    .groups = "drop"
  )

# 🔥 Step 2: format exactly how you want
supp_bpd_clean <- supp_bpd %>%
  mutate(
    pct = 100 * events / total,
    label = sprintf("%d / %d (%.1f%%)", events, total, pct)
  ) %>%
  select(phenotype, treat, label) %>%
  pivot_wider(
    names_from = treat,
    values_from = label
  ) %>%
  arrange(phenotype) %>%   # keeps order clean
  rename(
    "Phenotype" = phenotype,
    "Placebo n/N (%)" = Placebo,
    "Hydrocortisone n/N (%)" = Hydrocortisone
  )

# 🔥 Step 3: export clean table
supp_bpd_clean %>%
  flextable() %>%
  set_caption("Supplemental Table. Survival without BPD by phenotype and treatment.") %>%
  save_as_docx(path = "Supplement_Table_Survival_without_BPD.docx")


############################################################
# SUPPLEMENTAL — NDI FOLLOW-UP
############################################################
# 🔥 Step 1: restrict to phenotype cohort + NDI available
supp_ndi <- df_pheno %>%
  filter(!is.na(phenotype), !is.na(ndi_any)) %>%   # CRITICAL
  group_by(phenotype, treat) %>%
  summarise(
    events = sum(ndi_any == "Yes"),
    total = n(),
    .groups = "drop"
  )

# 🔥 Step 2: format like BPD table
supp_ndi_clean <- supp_ndi %>%
  mutate(
    pct = 100 * events / total,
    label = sprintf("%d / %d (%.1f%%)", events, total, pct)
  ) %>%
  select(phenotype, treat, label) %>%
  pivot_wider(
    names_from = treat,
    values_from = label
  ) %>%
  arrange(phenotype) %>%
  rename(
    "Phenotype" = phenotype,
    "Placebo n/N (%)" = Placebo,
    "Hydrocortisone n/N (%)" = Hydrocortisone
  )

# 🔥 Step 3: export
supp_ndi_clean %>%
  flextable() %>%
  set_caption("Supplemental Table 2. Neurodevelopmental impairment by phenotype and treatment.") %>%
  save_as_docx(path = "Supplement_Table_2.docx")


############################################################
# SUPPLEMENTAL FIGURE — CLUSTER ROBUSTNESS
# Appends to the existing pipeline. Assumes the following objects
# already exist in the environment from your script:
#   cluster_df_complete  (rand_id + 5 raw clustering vars, drop_na'd, n~760)
#   cluster_scaled        (cluster_df_complete with vars scaled, rand_id kept)
#   dist_mat               (dist() on cluster_scaled minus rand_id)
#   hc                      (hclust(dist_mat, method="ward.D2"))
#   df_pheno                (800-row analytic dataset with phenotype merged)
#   pheno_cols, pheno_order (color/level definitions from Figure 1/2/4)
############################################################

library(tidyverse)
library(cluster)
library(fpc)
library(patchwork)
library(pheatmap)
library(mclust)      # adjustedRandIndex

set.seed(2024)

pheno_order <- c("Phenotype 1", "Phenotype 2", "Phenotype 3", "Phenotype 4")

pheno_cols <- c(
  "Phenotype 1" = "#E41A1C",
  "Phenotype 2" = "#377EB8",
  "Phenotype 3" = "#4DAF4A",
  "Phenotype 4" = "#984EA3"
)

# Numeric matrix used for all clustering diagnostics (matches what hc was built on)
clust_vars <- c("GESTAGE", "BIRTHWTX", "APGAR5", "age_rndaysX", "map_fio2X")
X <- cluster_scaled %>% select(all_of(clust_vars)) %>% as.matrix()
rownames(X) <- cluster_df_complete$rand_id

phen_num <- cutree(hc, k = 4)  # integer 1:4, same order as cluster_df_complete/X

############################################################
# A. SILHOUETTE DISTRIBUTION BY PHENOTYPE
############################################################

sil <- silhouette(phen_num, dist_mat)
sil_df <- tibble(
  silhouette = sil[, 3],
  phenotype = factor(paste("Phenotype", sil[, 1]), levels = pheno_order)
)

sil_means <- sil_df %>%
  group_by(phenotype) %>%
  summarise(mean_sil = round(mean(silhouette), 2), n = n())

p_sil <- ggplot(sil_df, aes(x = silhouette, fill = phenotype)) +
  geom_histogram(bins = 35, alpha = 0.85, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(
    data = sil_means,
    aes(x = 0.3, y = Inf, label = paste0("mean = ", mean_sil, "\nn = ", n)),
    inherit.aes = FALSE, vjust = 1.5, size = 3
  ) +
  facet_wrap(~phenotype, ncol = 1) +
  scale_fill_manual(values = pheno_cols) +
  theme_classic() +
  labs(title = "A. Silhouette Distribution by Phenotype",
       x = "Silhouette width", y = "Count") +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

cat("Overall mean silhouette width:", round(mean(sil[, 3]), 3), "\n")

############################################################
# B. BOOTSTRAP STABILITY (Jaccard) WITH REFERENCE THRESHOLD
############################################################

boot <- clusterboot(
  data = X,
  B = 200,
  clustermethod = hclustCBI,
  method = "ward.D2",
  k = 4,
  seed = 2024
)

stab_df <- tibble(
  phenotype = factor(pheno_order, levels = pheno_order),
  stability = boot$bootmean,
  dissolved = boot$bootbrd  # times cluster "dissolved" across resamples
)

p_stab <- ggplot(stab_df, aes(x = phenotype, y = stability, fill = phenotype)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.7, y = 0.63, label = "moderate stability (0.6)",
           size = 3, hjust = 0, color = "grey30") +
  geom_text(aes(label = round(stability, 2)), vjust = -0.6, size = 3.5) +
  scale_fill_manual(values = pheno_cols) +
  ylim(0, 1) +
  theme_classic() +
  labs(title = "B. Cluster Stability (200 bootstrap resamples, Jaccard mean)",
       x = NULL, y = "Stability") +
  theme(legend.position = "none")

print(stab_df)
# NOTE: report these numbers honestly in text. Values <0.6 (Hennig's rule
# of thumb for "moderate" stability) should be flagged as a limitation,
# e.g. "Phenotype 1 showed lower bootstrap stability (mean Jaccard = X.XX),
# consistent with its less distinct silhouette profile."

############################################################
# C. PCA PROJECTION (already have this — add % variance + ellipse level note)
############################################################
# Your existing p_pca/p_fig1a already does this well. Only addition:
# report cumulative variance explained by PC1+PC2 in the caption, since
# ellipse overlap in 2D doesn't necessarily mean overlap in full 5D space.

pca_res_check <- prcomp(X, center = TRUE, scale. = FALSE)  # X already scaled
cum_var <- summary(pca_res_check)$importance["Cumulative Proportion", 2]
cat("Cumulative variance explained by PC1+PC2:", round(cum_var * 100, 1), "%\n")

############################################################
# D. CONSENSUS / CO-CLUSTERING HEATMAP
############################################################

n <- nrow(X)
B <- 200
co_matrix <- matrix(0, n, n)
count_matrix <- matrix(0, n, n)

for (b in 1:B) {
  idx <- sample(1:n, n, replace = TRUE)
  d_b <- dist(X[idx, , drop = FALSE])
  hc_b <- hclust(d_b, method = "ward.D2")
  cl_b <- cutree(hc_b, k = 4)
  
  uidx <- unique(idx)
  for (i in uidx) {
    pos_i <- which(idx == i)[1]
    for (j in uidx) {
      pos_j <- which(idx == j)[1]
      count_matrix[i, j] <- count_matrix[i, j] + 1
      if (cl_b[pos_i] == cl_b[pos_j]) co_matrix[i, j] <- co_matrix[i, j] + 1
    }
  }
}
# For n ~ 760 this loop is slow in base R (minutes, not seconds).
# If runtime is a problem, swap for the diceR package:
#   diceR::consensus_cluster(X, nk = 4, reps = 200, algorithms = "hc",
#                             hc.method = "ward.D2")
# which vectorizes this and returns the same consensus matrix concept.

consensus <- co_matrix / pmax(count_matrix, 1)
diag(consensus) <- 1

ord <- order(phen_num)
consensus_ord <- consensus[ord, ord]
ann <- data.frame(Phenotype = factor(paste("Phenotype", phen_num[ord]), levels = pheno_order))
rownames(ann) <- rownames(consensus_ord) <- colnames(consensus_ord) <- 1:n

pheatmap(
  consensus_ord,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  annotation_row = ann, annotation_col = ann,
  annotation_colors = list(Phenotype = pheno_cols),
  color = colorRampPalette(c("white", "#377EB8"))(100),
  main = "D. Consensus Co-Clustering Matrix (200 bootstrap resamples)",
  filename = "Supplement_Consensus_Heatmap.png",
  width = 6.5, height = 6
)

############################################################
# E. CLUSTER NUMBER (k) JUSTIFICATION
############################################################

k_range <- 2:6
diag_df <- map_dfr(k_range, function(k) {
  cl_k <- cutree(hc, k = k)
  cs <- cluster.stats(dist_mat, cl_k)
  tibble(
    k = k,
    mean_silhouette = mean(silhouette(cl_k, dist_mat)[, 3]),
    ch_index = cs$ch
  )
})

p_k1 <- ggplot(diag_df, aes(k, mean_silhouette)) +
  geom_line() + geom_point(size = 2) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = k_range) +
  theme_classic() +
  labs(title = "E1. Mean Silhouette by k", x = "k", y = "Mean silhouette width")

p_k2 <- ggplot(diag_df, aes(k, ch_index)) +
  geom_line() + geom_point(size = 2) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = k_range) +
  theme_classic() +
  labs(title = "E2. Calinski-Harabasz by k", x = "k", y = "CH index")

p_kselect <- p_k1 | p_k2
print(diag_df)
# Be candid in the manuscript if k=4 isn't the sharpest peak on these
# indices — justify k=4 on clinical grounds (e.g. maps onto known GA/
# respiratory-severity strata, or matches BPD severity categories used
# elsewhere in the trial) rather than implying pure data-driven selection.

############################################################
# F. CROSS-ALGORITHM AGREEMENT (Adjusted Rand Index)
############################################################

km <- kmeans(X, centers = 4, nstart = 50)
pam_fit <- pam(X, k = 4)

ari_kmeans <- adjustedRandIndex(phen_num, km$cluster)
ari_pam    <- adjustedRandIndex(phen_num, pam_fit$clustering)

agreement_df <- tibble(
  Comparison = c("Ward hierarchical vs k-means", "Ward hierarchical vs PAM"),
  ARI = c(ari_kmeans, ari_pam)
)
print(agreement_df)

p_ari <- ggplot(agreement_df, aes(x = Comparison, y = ARI)) +
  geom_col(fill = "grey40", width = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ylim(0, 1) +
  theme_classic() +
  labs(title = "F. Cross-Algorithm Agreement", y = "Adjusted Rand Index", x = NULL) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
# ARI > 0.5 = strong agreement; report directly, e.g.:
# "Phenotype assignments showed [strong/moderate] agreement with an
#  independent k-means solution (ARI = 0.XX) and PAM solution (ARI = 0.XX),
#  indicating the structure is not an artifact of the hierarchical/Ward
#  linkage choice."

############################################################
# G. FEATURE-LEVEL CENTROID HEATMAP (clinically interpretable complement to PCA)
############################################################
# Reuses Figure 1D-style logic but placed here explicitly as a robustness/
# structure panel. If you already have p_fig1d from Figure 1, you can skip
# regenerating this and just cross-reference it.

centroid_df <- df_pheno %>%
  filter(phenotype %in% pheno_order) %>%
  select(phenotype, all_of(clust_vars)) %>%
  mutate(across(all_of(clust_vars), scale)) %>%
  group_by(phenotype) %>%
  summarise(across(all_of(clust_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-phenotype, names_to = "Variable", values_to = "Zscore")

p_centroid <- ggplot(centroid_df, aes(x = Variable, y = phenotype, fill = Zscore)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Zscore, 2)), size = 3) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0) +
  theme_minimal() +
  labs(title = "G. Standardized Feature Centroids by Phenotype", x = NULL, y = NULL, fill = "Mean\n(z-score)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

############################################################
# H. SENSITIVITY TO MISSINGNESS
# You already do complete-case phenotyping (cluster_df_complete, n~760
# out of 800). This checks whether the ~40 dropped infants' exclusion
# meaningfully changes structure — i.e., a leave-out-missing sensitivity
# check rather than a true imputation comparison, since no imputed
# version currently exists in the pipeline.
############################################################

# Re-cluster on a random 90% subsample repeatedly and check agreement
# with the full complete-case solution, as a proxy for sensitivity to
# which subjects happen to be complete cases.
set.seed(2024)
n_reps <- 50
ari_subsample <- numeric(n_reps)
for (r in 1:n_reps) {
  keep <- sample(1:n, size = floor(0.9 * n))
  d_sub <- dist(X[keep, , drop = FALSE])
  hc_sub <- hclust(d_sub, method = "ward.D2")
  cl_sub <- cutree(hc_sub, k = 4)
  ari_subsample[r] <- adjustedRandIndex(phen_num[keep], cl_sub)
}
cat("Subsample sensitivity — mean ARI (90% resamples vs full solution):",
    round(mean(ari_subsample), 3), "SD:", round(sd(ari_subsample), 3), "\n")

# If you later build a formally imputed dataset (e.g. mice) for the full
# n=800, the cleaner comparison is:
#   cl_imputed <- cutree(hclust(dist(scale(imputed_vars)), "ward.D2"), 4)
#   adjustedRandIndex(cl_imputed[complete_ids], phen_num)

############################################################
# COMBINE INTO SUPPLEMENT PAGES
############################################################

p_row1 <- (p_sil | p_stab)
p_row2 <- p_kselect
p_row3 <- (p_ari | p_centroid)

ggsave("Supplement_p1_Silhouette_Stability.png", p_row1, width = 11, height = 9, dpi = 300, bg = "white")
ggsave("Supplement_p2_k_Justification.png", p_row2, width = 10, height = 5, dpi = 300, bg = "white")
ggsave("Supplement_p3_Agreement_Centroids.png", p_row3, width = 11, height = 5, dpi = 300, bg = "white")
# Supplement_Consensus_Heatmap.png already saved directly by pheatmap() above.







############################################################
# REVISED: TREATMENT-EFFECT HETEROGENEITY MODEL
# Addresses reviewer comment: "The treatment-effect heterogeneity model
# needs more detail, including covariate adjustment, handling of center
# and multiple births, and whether random effects were used."
#
# Replaces the original fit_extub <- glm(extub14 ~ phenotype * treat, ...)
# in the Figure 4 / Table 2 block, which had NO covariate adjustment,
# NO clustering handling, and NO random effects despite the Methods text
# claiming GEE with robust SEs accounting for center and multiple births.
############################################################
library(geepack)   # GEE with robust (sandwich) SEs, clustered on family
library(lme4)       # sensitivity: GLMM with crossed/nested random effects
library(emmeans)
library(broom)
library(broom.mixed)

## ------------------------------------------------------------
## 0. Rationale for clustering choice
## ------------------------------------------------------------
# rand_id = randomization unit. From your own sanity_check(), multiple
# births SHARE a rand_id, so clustering on rand_id directly captures
# within-family (twin/triplet) correlation. randcent = study center.
# Center is a separate, crossed source of correlation (infants at the
# same center are more similar to each other than to infants elsewhere,
# independent of family).
#
# GEE (geeglm) supports exactly ONE clustering/id variable. To satisfy
# both the "center" and "multiple births" clauses in your Methods text,
# report BOTH of the following as primary + sensitivity, rather than
# implying a single model does everything:
#
#   PRIMARY:     GEE clustered on rand_id (multiple births), with study
#                center (randcent) included as a fixed-effect covariate
#                to adjust for between-center differences in outcome risk.
#   SENSITIVITY: GLMM (glmer) with crossed/nested random intercepts for
#                BOTH randcent and rand_id simultaneously, to confirm the
#                treatment-by-phenotype interaction is not sensitive to
#                the choice of clustering strategy.
#
# If you intended the original Methods sentence to describe only the GEE
# approach, revise the Methods text to drop the "random effects" framing
# for the primary model and reserve it for the sensitivity analysis (see
# suggested Methods wording at the bottom of this file).

## ------------------------------------------------------------
## 1. Build adjusted analytic dataset (adds covariates + clustering IDs)
## ------------------------------------------------------------

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

summary(fit_extub_gee)

# Global interaction test (Wald, appropriate for GEE — GEE does not
# support a likelihood-ratio test since it is not a full-likelihood method)
anova_gee <- anova(fit_extub_gee)
print(anova_gee)
p_interaction_gee <- anova_gee["phenotype:treat", "P(>|Chi|)"]
cat("Phenotype x Treatment interaction (GEE Wald test) p =",
    round(p_interaction_gee, 3), "\n")

# Phenotype-specific treatment ORs (marginal means from the adjusted model)

## ------------------------------------------------------------
## FIX: emmeans column names differ for geeglm objects
## ------------------------------------------------------------

or_summary_gee <- as.data.frame(summary(or_contrast_gee, infer = c(TRUE, TRUE)))
print(colnames(or_summary_gee))   # inspect once — likely "lower.CL"/"upper.CL" or "asymp.LCL"/"asymp.UCL"

# Robust version that works regardless of which naming emmeans used
lcl_col <- intersect(c("asymp.LCL", "lower.CL"), colnames(or_summary_gee))[1]
ucl_col <- intersect(c("asymp.UCL", "upper.CL"), colnames(or_summary_gee))[1]

or_df_gee <- or_summary_gee %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    OR = odds.ratio,
    LCL = .data[[lcl_col]],
    UCL = .data[[ucl_col]],
    label = sprintf("%.2f (%.2f–%.2f)", OR, LCL, UCL)
  )

print(or_df_gee %>% select(phenotype, OR, LCL, UCL, label, p.value))

# These adjusted ORs (GESTAGE-, sex-, and center-adjusted, family-clustered)
# are what should replace the unadjusted fit_extub estimates in Figure 4A
# and Table 2.

## ------------------------------------------------------------
## 3. SENSITIVITY MODEL: GLMM with crossed random intercepts for center
##    and family (rand_id), to directly test whether results depend on
##    modeling correlation via random effects instead of GEE/robust SE
## ------------------------------------------------------------

fit_extub_glmm <- glmer(
  extub14_num ~ phenotype * treat + GESTAGE + male +
    (1 | randcent) + (1 | rand_id),
  data = fig4_df,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(fit_extub_glmm)

# LRT for interaction is valid here since GLMM is a full-likelihood model
fit_extub_glmm_noint <- glmer(
  extub14_num ~ phenotype + treat + GESTAGE + male +
    (1 | randcent) + (1 | rand_id),
  data = fig4_df,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

lrt_glmm <- anova(fit_extub_glmm_noint, fit_extub_glmm)
print(lrt_glmm)
p_interaction_glmm <- lrt_glmm$`Pr(>Chisq)`[2]
cat("Phenotype x Treatment interaction (GLMM LRT) p =",
    round(p_interaction_glmm, 3), "\n")

# Report the variance components — this is what tells you whether center
# or family clustering matters much in practice (small variance = little
# excess correlation, supporting that simpler models aren't badly biased)
print(VarCorr(fit_extub_glmm))

emm_or_glmm <- emmeans(fit_extub_glmm, ~ treat | phenotype, type = "response")
or_contrast_glmm <- contrast(emm_or_glmm, method = "revpairwise", by = "phenotype")
or_df_glmm <- summary(or_contrast_glmm, infer = c(TRUE, TRUE)) %>%
  as.data.frame() %>%
  mutate(label = sprintf("%.2f (%.2f–%.2f)", odds.ratio, asymp.LCL, asymp.UCL))
print(or_df_glmm)

## ------------------------------------------------------------
## 4. SIDE-BY-SIDE COMPARISON TABLE (for reviewer response / supplement)
## ------------------------------------------------------------

comparison_df <- or_df_gee %>%
  select(phenotype, GEE_OR_CI = label) %>%
  left_join(
    or_df_glmm %>% mutate(phenotype = factor(phenotype, levels = pheno_order)) %>%
      select(phenotype, GLMM_OR_CI = label),
    by = "phenotype"
  )

print(comparison_df)
# If GEE and GLMM estimates are similar (which is expected when random-
# effect variances are small), this table is strong evidence for the
# reviewer that phenotype-specific treatment effects are not an artifact
# of how clustering/correlation was modeled.

## ------------------------------------------------------------
## 5. UPDATE TABLE 3 (Firth models) — add covariate adjustment,
##    note clustering limitation explicitly
## ------------------------------------------------------------
# logistf does not support clustered/robust SEs, so this remains an
# UNCLUSTERED sparse-data sensitivity analysis. Now adjusted for GESTAGE
# and sex to align with the primary model's covariate set; center is
# omitted here to preserve degrees of freedom given small per-phenotype N.

library(logistf)

df_tx <- fig4_df %>%
  filter(!is.na(phenotype), !is.na(treat), !is.na(extub14))

firth_results <- df_tx %>%
  group_by(phenotype) %>%
  group_modify(~ {
    fit <- logistf(extub14_num ~ treat + GESTAGE + male, data = .x)
    tibble(
      OR  = exp(coef(fit)["treatHydrocortisone"]),
      LCL = exp(confint(fit)["treatHydrocortisone", 1]),
      UCL = exp(confint(fit)["treatHydrocortisone", 2]),
      p   = fit$prob["treatHydrocortisone"]
    )
  }) %>%
  ungroup()

print(firth_results)

# NOTE FOR MANUSCRIPT: state explicitly that Table 3 Firth estimates are
# adjusted for gestational age and sex but, unlike the primary GEE model,
# do not account for center or multiple-birth clustering, since penalized
# likelihood methods for sparse binary data with clustered SEs are not
# well established in standard software. Frame Table 3 as a small-sample
# bias-correction sensitivity analysis, not the primary inferential model.

############################################################
# SUGGESTED METHODS TEXT REVISION
# (Replace the "Treatment Effect Heterogeneity" paragraph)
############################################################

# Treatment Effect Heterogeneity
#
# To evaluate potential heterogeneity of treatment effect across
# phenotypes, generalized estimating equations (GEE) with an exchangeable
# working correlation structure and robust (sandwich) standard errors were
# used, with infants clustered by randomization unit (rand_id) to account
# for correlation among multiple births. Models included treatment group,
# phenotype, and a treatment-by-phenotype interaction term, and were
# adjusted for gestational age, sex, and study center. Global interaction
# was assessed using a Wald test appropriate for the GEE framework. As a
# sensitivity analysis, generalized linear mixed models with crossed
# random intercepts for study center and randomization unit were fit to
# confirm that treatment-by-phenotype estimates were not sensitive to the
# choice of correlation-modeling strategy (Supplemental Table SX). Odds
# ratios and absolute risk differences were estimated using marginal means
# derived from the adjusted GEE model. When sparse events were present
# within phenotype strata, penalized logistic regression with Firth
# correction, adjusted for gestational age and sex, was applied as an
# additional sensitivity analysis; because clustered standard errors are
# not well established for penalized likelihood methods, these estimates
# do not account for center or multiple-birth correlation and are
# presented as a small-sample bias-correction check rather than the
# primary analysis. Given the exploratory nature of phenotype derivation,
# interaction results were interpreted cautiously.




############################################################
# TABLE 2 (REVISED): Hydrocortisone vs Placebo by Baseline Phenotype
# Adds Risk Ratio + Absolute Risk Difference alongside Odds Ratio,
# using the ADJUSTED GEE model (fit_extub_gee) as the primary analysis
# per the reviewer's request re: covariate adjustment / clustering.
############################################################

library(emmeans)
library(dplyr)
library(scales)
library(flextable)

## ------------------------------------------------------------
## 1. RAW OBSERVED COUNTS (unadjusted — for the descriptive columns)
## ------------------------------------------------------------

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
      "%.1f%% (%.1f to %.1f)",
      estimate * 100, .data[[lcl_col_ard]] * 100, .data[[ucl_col_ard]] * 100
    )
  ) %>%
  select(phenotype, ARD_CI)

## ------------------------------------------------------------
## 5. GLOBAL INTERACTION P-VALUE (adjusted GEE Wald test — already
##    computed as p_interaction_gee earlier: 0.121)
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
    RR_CI, ARD_CI, OR_CI, p_fmt
  )

print(table2_final)

## ------------------------------------------------------------
## 7. FORMATTED EXPORT (matches the reviewer's requested layout)
## ------------------------------------------------------------

ft <- table2_final %>%
  flextable() %>%
  set_header_labels(
    phenotype = "Phenotype",
    n_over_N = "n/N (Placebo vs HC)",
    pct_placebo = "Extubation %\n(Placebo)",
    pct_hc = "Extubation %\n(HC)",
    RR_CI = "Risk Ratio\n(95% CI)",
    ARD_CI = "Absolute Risk\nDifference (95% CI)",
    OR_CI = "Odds Ratio\n(95% CI)",
    p_fmt = "p-value"
  ) %>%
  add_header_lines(values = interaction_p_label) %>%
  add_header_lines(values = "Table 2. Hydrocortisone vs placebo by baseline phenotype (adjusted).") %>%
  bold(part = "header", i = 1) %>%
  italic(part = "header", i = 2) %>%
  autofit() %>%
  add_footer_lines(
    "Values represent observed counts and proportions (unadjusted). Risk ratios, absolute risk differences, and odds ratios are model-based estimates from a generalized estimating equation adjusted for gestational age, sex, and study center, with infants clustered by randomization unit (rand_id) to account for correlation among multiple births. Analyses restricted to infants with assigned phenotype (n = 760). Unadjusted estimates and a sensitivity analysis using mixed-effects models with random intercepts for center and randomization unit are provided in the Supplement."
  )

ft

save_as_docx(ft, path = "Table2_Adjusted_RR_ARD_OR_by_Phenotype.docx")

############################################################
# SUPPLEMENTAL TABLE: UNADJUSTED vs ADJUSTED COMPARISON
# Put this in the supplement to directly show the reviewer the effect
# of adjustment — this is the transparency piece Table 2 alone can't show.
############################################################

# Unadjusted OR (original fit_extub, no covariates/clustering) --
# rerun quickly here for the comparison table, or reuse if still in memory
fit_extub_unadj <- glm(extub14_num ~ phenotype * treat, data = fig4_df, family = binomial)
emm_unadj <- emmeans(fit_extub_unadj, ~ treat | phenotype, type = "response")
or_contrast_unadj <- contrast(emm_unadj, method = "revpairwise", by = "phenotype")
or_unadj_df <- as.data.frame(summary(or_contrast_unadj, infer = c(TRUE, TRUE))) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    Unadjusted_OR_CI = sprintf("%.2f (%.2f-%.2f)", odds.ratio, asymp.LCL, asymp.UCL)
  ) %>%
  select(phenotype, Unadjusted_OR_CI)

supp_compare <- or_unadj_df %>%
  left_join(or_df %>% rename(Adjusted_OR_CI = OR_CI), by = "phenotype")

print(supp_compare)

supp_compare %>%
  flextable() %>%
  set_header_labels(
    phenotype = "Phenotype",
    Unadjusted_OR_CI = "Unadjusted OR (95% CI)",
    Adjusted_OR_CI = "Adjusted OR (95% CI)\nGESTAGE + sex + center, family-clustered"
  ) %>%
  add_footer_lines(
    sprintf(
      "Global interaction P = 0.06 (unadjusted) vs P = %.2f (adjusted GEE). See Supplemental Table SX2 for GLMM and Firth sensitivity analyses.",
      p_interaction_gee
    )
  ) %>%
  autofit() %>%
  save_as_docx(path = "Supplement_Table_Unadjusted_vs_Adjusted.docx")







############################################################
# SUPPLEMENTAL TABLE SX2: CROSS-METHOD AGREEMENT
# Treatment effect on extubation <=14 days by phenotype,
# GEE (primary) vs GLMM vs Firth penalized logistic regression.
#
# Assumes the following objects already exist from prior scripts:
#   or_df_gee   (from fit_extub_gee, via regrid/contrast — has phenotype, OR/label/p.value)
#   or_df_glmm  (from fit_extub_glmm)
#   firth_results (from logistf loop)
#   p_interaction_gee   (GEE Wald p for phenotype:treat, = 0.121)
#   p_interaction_glmm  (GLMM LRT p for phenotype:treat, = 0.083)
############################################################

library(dplyr)
library(flextable)

## ------------------------------------------------------------
## 1. Harmonize each source into phenotype | estimate_CI | p_fmt
## ------------------------------------------------------------

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.2f", p))

gee_col <- or_df_gee %>%
  transmute(
    phenotype = factor(phenotype, levels = pheno_order),
    GEE = sprintf("%s, P=%s", label, fmt_p(p.value))
  )

glmm_col <- or_df_glmm %>%
  transmute(
    phenotype = factor(phenotype, levels = pheno_order),
    GLMM = sprintf("%s, P=%s", label, fmt_p(p.value))
  )

firth_col <- firth_results %>%
  transmute(
    phenotype = factor(phenotype, levels = pheno_order),
    Firth = sprintf("%.2f (%.2f-%.2f), P=%s", OR, LCL, UCL, fmt_p(p))
  )

## ------------------------------------------------------------
## 2. Combine into one table
## ------------------------------------------------------------

cross_method_tbl <- gee_col %>%
  left_join(glmm_col, by = "phenotype") %>%
  left_join(firth_col, by = "phenotype") %>%
  arrange(phenotype)

print(cross_method_tbl)

## ------------------------------------------------------------
## 3. Global interaction test row (Firth has no pooled interaction
##    test across strata without a separate pooled model, so leave
##    that cell blank rather than force a number)
## ------------------------------------------------------------

interaction_row <- tibble(
  phenotype = "Global interaction test",
  GEE   = sprintf("Wald P = %.2f", p_interaction_gee),
  GLMM  = sprintf("LRT P = %.2f", p_interaction_glmm),
  Firth = "—"
)

cross_method_tbl_final <- bind_rows(
  cross_method_tbl %>% mutate(phenotype = as.character(phenotype)),
  interaction_row
)

## ------------------------------------------------------------
## 4. Formatted export
## ------------------------------------------------------------

ft_cross <- cross_method_tbl_final %>%
  flextable() %>%
  set_header_labels(
    phenotype = "Phenotype",
    GEE   = "GEEᵃ\nOR (95% CI), P",
    GLMM  = "GLMMᵇ\nOR (95% CI), P",
    Firth = "Firth penalized logisticᶜ\nOR (95% CI), P"
  ) %>%
  add_header_lines(
    values = "Supplemental Table SX2. Cross-method agreement for phenotype-specific hydrocortisone effect on extubation ≤14 days."
  ) %>%
  bold(part = "header", i = 1) %>%
  bold(i = ~ phenotype == "Global interaction test") %>%
  border(i = which(cross_method_tbl_final$phenotype == "Global interaction test") - 1,
         border.bottom = officer::fp_border(width = 1.5), part = "body") %>%
  autofit() %>%
  add_footer_lines(
    c(
      "a. GEE: exchangeable working correlation, clustered on randomization unit (rand_id) to account for correlation among multiple births; adjusted for gestational age, sex, and study center. Primary model (Table 2).",
      "b. GLMM: crossed random intercepts for study center and randomization unit. Family (rand_id) variance was estimated at the boundary (~0), consistent with negligible within-family correlation in this cohort; center variance was non-negligible (Supplemental Table SX1).",
      "c. Firth penalized logistic regression, adjusted for gestational age and sex only. Does not account for center or family clustering, as clustered standard errors are not well established for penalized-likelihood methods; included as a small-sample bias-correction sensitivity analysis given sparse events, particularly in Phenotype 1.",
      "Global interaction test not computed for the Firth approach, as penalized per-stratum models do not yield a single pooled interaction test."
    )
  )

ft_cross

save_as_docx(ft_cross, path = "Supplement_Table_SX2_Cross_Method_Agreement.docx")

## ------------------------------------------------------------
## 5. Optional companion: variance components table (SX1), referenced
##    in footnote b above — small standalone table showing why family
##    clustering turned out to matter little vs. center.
## ------------------------------------------------------------

varcomp_tbl <- as.data.frame(VarCorr(fit_extub_glmm)) %>%
  transmute(
    `Random effect` = grp,
    Variance = round(vcov, 3),
    `Std. Dev.` = round(sdcor, 3)
  )

varcomp_tbl %>%
  flextable() %>%
  add_header_lines(values = "Supplemental Table SX1. Random-effect variance components (GLMM sensitivity model).") %>%
  bold(part = "header", i = 1) %>%
  autofit() %>%
  add_footer_lines(
    "Near-zero variance for randomization unit (rand_id) indicates negligible residual correlation among multiple births after accounting for fixed effects; nonzero center variance indicates meaningful between-center heterogeneity, consistent with the significant center term in the GEE model (Wald P < 0.001)."
  ) %>%
  save_as_docx(path = "Supplement_Table_SX1_Variance_Components.docx")




## =========================================================
## FIGURE 2. Clinical outcomes stratified by baseline phenotype
## FIX: p-values were displaying as "p = 0.000" in Panels A and C,
## which implies an impossible exact-zero p-value. Both panels now
## use a shared fmt_p() helper that switches to "<0.001" below that
## threshold, consistent with formatting used in Table 2 and Results.
## FIX 2: Panel D's "χ²" label was getting clipped on the left edge
## (unicode glyph rendering issue in annotate()). Now built as a
## plotmath expression (chi^2) with parse = TRUE for reliable rendering.
## =========================================================

library(tidyverse)
library(patchwork)
library(scales)
library(ggtext)
library(emmeans)
library(broom)

# -------------------------
# Constants
# -------------------------
pheno_order <- c("Phenotype 1","Phenotype 2","Phenotype 3","Phenotype 4")

pheno_cols <- c(
  "Phenotype 1" = "#E41A1C",
  "Phenotype 2" = "#377EB8",
  "Phenotype 3" = "#4DAF4A",
  "Phenotype 4" = "#984EA3"
)

# -------------------------
# Shared p-value formatter
# -------------------------
fmt_p <- function(p) {
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# -------------------------
# Publication theme
# -------------------------
theme_publication <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 8)),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey30"),
      axis.ticks = element_line(color = "grey30"),
      axis.text = element_text(size = 11, color = "grey20"),
      axis.title = element_text(size = 11),
      plot.background = element_rect(fill = "grey98", color = NA),
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

## =========================================================
## PANEL A — Successful extubation ≤14 days (GLM + EMMEANS)
## =========================================================

df_extub <- df_pheno %>%
  filter(!is.na(extube14), phenotype %in% pheno_order) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    extube14 = as.numeric(extube14 == 1)
  )

glm_extub <- glm(extube14 ~ phenotype, data = df_extub, family = binomial)
lrt_p <- anova(glm_extub, test = "LRT")$`Pr(>Chi)`[2]

emm_extub <- emmeans(glm_extub, ~ phenotype, type = "response")

emm_df <- as_tibble(emm_extub) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    label = sprintf(
      "%d/%d (%s)",
      tapply(df_extub$extube14, df_extub$phenotype, sum)[phenotype],
      tapply(df_extub$extube14, df_extub$phenotype, length)[phenotype],
      percent(prob, 0.1)
    )
  ) %>%
  rename(rate = prob, lcl = asymp.LCL, ucl = asymp.UCL) %>%
  mutate(phenotype = fct_rev(phenotype))

p_fig2a <- ggplot(emm_df, aes(x = rate, y = phenotype, color = phenotype)) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.25, linewidth = 1.2) +
  geom_point(size = 5) +
  geom_text(aes(label = label), nudge_x = 0.12, size = 3.8,
            fontface = "bold", color = "grey10") +
  annotate(
    "text",
    x = 0.02, y = length(pheno_order) + 0.4,
    label = sprintf("Logistic regression: overall p %s",
                    ifelse(lrt_p < 0.001, "< 0.001", sprintf("= %.3f", lrt_p))),
    hjust = 0, size = 3.5, color = "grey30"
  ) +
  scale_x_continuous(labels = percent, limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = pheno_cols) +
  theme_publication() +
  labs(title = "A. Successful extubation ≤14 days",
       x = "Model-estimated probability (95% CI)", y = NULL) +
  guides(color = "none")

## =========================================================
## PANEL B — BPD severity at 36 weeks (Jensen) — unchanged, no p-value shown
## =========================================================
bpd_df <- df_pheno %>%
  filter(!is.na(bpd_jen36_hc), phenotype %in% pheno_order) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    jensen = factor(
      bpd_jen36_hc, levels = 0:3,
      labels = c("No support","≤2L NC",">2L/NIPPV","IPPV")
    )
  ) %>%
  dplyr::count(phenotype, jensen) %>%
  group_by(phenotype) %>%
  mutate(
    prop = n / sum(n),
    total_n = sum(n),
    label = if_else(prop >= 0.08, percent(prop, 0), "")
  )

bpd_n <- bpd_df %>%
  distinct(phenotype, total_n) %>%
  mutate(label = paste0(phenotype, "\n(n=", total_n, ")"))

p_fig2b <- ggplot(bpd_df, aes(x = phenotype, y = prop, fill = jensen)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.9) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 3.8, fontface = "bold", color = "white") +
  scale_x_discrete(labels = setNames(bpd_n$label, bpd_n$phenotype)) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_fill_manual(values = c("#F0F0F0","#BDBDBD","#737373","#252525")) +
  theme_publication() +
  labs(title = "B. BPD severity at 36 weeks (Jensen)",
       x = NULL, y = "Distribution", fill = "Severity")

## =========================================================
## PANEL C — Duration of mechanical ventilation
## =========================================================
vent_df <- df_pheno %>%
  filter(!is.na(days_VENT36), phenotype %in% pheno_order) %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order), vent_days = days_VENT36)

vent_sum <- vent_df %>%
  group_by(phenotype) %>%
  summarise(
    n = n(), median = median(vent_days), q25 = quantile(vent_days, 0.25),
    q75 = quantile(vent_days, 0.75), mean = mean(vent_days),
    label = sprintf("Median %.0f\n(IQR %.0f–%.0f)", median, q25, q75),
    .groups = "drop"
  )

kw_p <- kruskal.test(vent_days ~ phenotype, data = vent_df)$p.value

p_fig2c <- ggplot(vent_df, aes(x = phenotype, y = vent_days, fill = phenotype)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = phenotype), width = 0.15, alpha = 0.4, size = 1.3) +
  geom_point(data = vent_sum, aes(y = mean), size = 4.5, shape = 23, color = "white") +
  geom_text(data = vent_sum, aes(y = q75 + 6, label = label), size = 3.2, fontface = "bold") +
  annotate(
    "text",
    x = 2.5, y = max(vent_df$vent_days) * 0.95,
    label = sprintf("Kruskal–Wallis p %s",
                    ifelse(kw_p < 0.001, "< 0.001", sprintf("= %.3f", kw_p))),
    size = 3.5, color = "grey30"
  ) +
  scale_fill_manual(values = pheno_cols) +
  scale_color_manual(values = pheno_cols) +
  theme_publication() +
  labs(title = "C. Duration of mechanical ventilation", x = NULL, y = "Days on ventilator") +
  guides(color = "none", fill = "none")

## =========================================================
## PANEL D — Any neurodevelopmental impairment
## FIX: chi-squared symbol now rendered via plotmath (parse = TRUE)
## instead of a raw unicode "χ²" glyph, which was getting clipped.
## =========================================================
ndi_indiv <- df_pheno %>%
  filter(!is.na(NDI3), phenotype %in% pheno_order) %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order), y = as.numeric(NDI3 == 1))

ndi_sum <- ndi_indiv %>%
  group_by(phenotype) %>%
  summarise(
    n = n(), n_ndi = sum(y), rate = mean(y),
    se = sqrt(rate * (1 - rate) / n),
    lcl = pmax(0, rate - 1.96 * se), ucl = pmin(1, rate + 1.96 * se),
    label = sprintf("%d/%d\n%s", n_ndi, n, percent(rate, 0.1)),
    .groups = "drop"
  )

tbl <- table(ndi_indiv$phenotype, ndi_indiv$y)
ndi_p <- if (any(tbl < 5)) fisher.test(tbl)$p.value else chisq.test(tbl)$p.value

# Plotmath-safe test label: chi^2 for chi-squared, quoted string for Fisher's
ndi_test_lab <- if (any(tbl < 5)) "\"Fisher's exact\"" else "chi^2"
ndi_p_lab <- ifelse(ndi_p < 0.001, "\"< 0.001\"", sprintf("\"= %.3f\"", ndi_p))

p_fig2d <- ggplot() +
  geom_jitter(data = ndi_indiv, aes(x = phenotype, y = y, color = phenotype),
              width = 0.12, height = 0.02, alpha = 0.35, size = 1.4) +
  geom_linerange(data = ndi_sum, aes(x = phenotype, ymin = lcl, ymax = ucl, color = phenotype),
                 linewidth = 1.6) +
  geom_point(data = ndi_sum, aes(x = phenotype, y = rate, fill = phenotype),
             size = 5.5, shape = 21, color = "white") +
  geom_text(data = ndi_sum, aes(x = phenotype, y = rate, label = label),
            nudge_y = 0.14, size = 3.2, fontface = "bold") +
  annotate(
    "text",
    x = 1.5, y = 0.88,
    label = sprintf("%s~p~%s", ndi_test_lab, ndi_p_lab),
    parse = TRUE,
    size = 3, color = "grey30"
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
  scale_color_manual(values = pheno_cols) +
  scale_fill_manual(values = pheno_cols) +
  theme_publication() +
  labs(title = "D. Any neurodevelopmental impairment at 18–24 months",
       x = NULL, y = "Probability") +
  guides(color = "none", fill = "none")

## =========================================================
## FINAL LAYOUT
## =========================================================
fig2_final <-
  (p_fig2a | p_fig2b) /
  (p_fig2c | p_fig2d) +
  plot_annotation(
    title = "Figure 2. Clinical outcomes stratified by baseline phenotype",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )

fig2_final

ggsave("Figure2_enhanced.png", fig2_final,
       width = 12, height = 10, dpi = 300, bg = "white")

ggsave("Figure2_enhanced.pdf", fig2_final,
       width = 12, height = 10, dpi = 300, bg = "white")

##########################
#figure 3 adjusted----
##########################

library(tidyverse )
library(geepack)
library(emmeans)
library(patchwork)
library(scales)

## =========================
## Publication Theme
## =========================
theme_publication <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "grey40", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "grey20", size = base_size - 1),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      strip.background = element_blank(),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA)
    )
}

## =========================
## Colors
## =========================
pheno_cols <- c(
  "Phenotype 1" = "#E57373",
  "Phenotype 2" = "#64B5F6",
  "Phenotype 3" = "#81C784",
  "Phenotype 4" = "#BA68C8"
)

treat_cols <- c(
  "Placebo" = "grey70",
  "Hydrocortisone" = "#4C72B0"
)

## =========================
## Prepare dataset (adds covariates + clustering IDs, matches fig4_df
## from the Table 2 / treatment-heterogeneity revision)
## =========================
fig3_df <- df_pheno %>%
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
  arrange(rand_id)  # required for geeglm

## =========================
## PRIMARY MODEL — adjusted GEE (same as Table 2 / fit_extub_gee)
## =========================
fit_extub_gee <- geeglm(
  extub14_num ~ phenotype * treat + GESTAGE + male + randcent,
  id = rand_id,
  data = fig3_df,
  family = binomial,
  corstr = "exchangeable"
)

anova_gee <- anova(fit_extub_gee)
p_interaction_gee <- anova_gee["phenotype:treat", "P(>|Chi|)"]

## =========================
## PANEL A — Forest plot (RISK RATIO, not OR)
## =========================

emm_link <- emmeans(fit_extub_gee, ~ treat | phenotype)  # logit scale
emm_log  <- regrid(emm_link, transform = "log")           # log-probability scale

rr_contrast <- contrast(emm_log, method = "revpairwise", by = "phenotype", type = "response")
rr_summary  <- as.data.frame(summary(rr_contrast, infer = c(TRUE, TRUE)))

lcl_col <- intersect(c("asymp.LCL", "lower.CL"), colnames(rr_summary))[1]
ucl_col <- intersect(c("asymp.UCL", "upper.CL"), colnames(rr_summary))[1]
ratio_col <- intersect(c("ratio", "response"), colnames(rr_summary))[1]

rr_df <- rr_summary %>%
  mutate(
    phenotype = factor(phenotype, levels = rev(pheno_order)),
    RR = .data[[ratio_col]],
    lower = .data[[lcl_col]],
    upper = .data[[ucl_col]],
    sig = p.value < 0.05,
    label = sprintf("%.2f (%.2f-%.2f)", RR, lower, upper)
  )

p_fig3a <- ggplot(rr_df, aes(x = RR, y = phenotype)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -Inf, ymax = Inf, fill = "grey95") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25, linewidth = 1.1, color = "#4C72B0") +
  geom_point(size = 4.5, shape = 21, fill = "white", color = "#4C72B0", stroke = 1) +
  geom_text(aes(label = label), x = 6, hjust = 0, size = 3.3) +
  annotate(
    "text", x = 0.45, y = 4.4,
    label = sprintf("Interaction P = %.2f", p_interaction_gee),
    hjust = 0, size = 3.2, color = "grey30", fontface = "italic"
  ) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), limits = c(0.4, 9)) +
  theme_publication() +
  labs(
    title = "A. Hydrocortisone effect on successful extubation ≤14 days",
    subtitle = "Risk ratios with 95% confidence intervals, adjusted for gestational age, sex, and center",
    x = "Risk ratio (Hydrocortisone vs Placebo)",
    y = NULL
  )

## =========================
## PANEL B — Model-estimated probability of extubation (adjusted)
## =========================
prob_df <- as.data.frame(emmeans(fit_extub_gee, ~ treat | phenotype, type = "response")) %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order))

lcl_col_b <- intersect(c("asymp.LCL", "lower.CL"), colnames(prob_df))[1]
ucl_col_b <- intersect(c("asymp.UCL", "upper.CL"), colnames(prob_df))[1]

or_contrast_for_p <- contrast(emm_link, method = "revpairwise", by = "phenotype")  # for p-values only
pvals_extub <- as.data.frame(summary(or_contrast_for_p)) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    label = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.05  ~ sprintf("p = %.3f", p.value),
      TRUE            ~ "NS"
    )
  )

annot_b <- prob_df %>%
  group_by(phenotype) %>%
  summarise(y = max(.data[[ucl_col_b]]) + 0.08, .groups = "drop") %>%
  left_join(pvals_extub, by = "phenotype")

p_fig3b <- ggplot(prob_df, aes(x = treat, y = prob, fill = treat)) +
  geom_col(width = 0.6, color = "grey40") +
  geom_errorbar(aes(ymin = .data[[lcl_col_b]], ymax = .data[[ucl_col_b]]), width = 0.15) +
  geom_text(aes(label = percent(prob, 0.1)), position = position_stack(vjust = 0.5),
            fontface = "bold", size = 4) +
  geom_text(data = annot_b, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 4) +
  facet_wrap(~ phenotype, nrow = 2) +
  scale_fill_manual(values = treat_cols) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  theme_publication() +
  labs(
    title = "B. Model-estimated probability of extubation ≤14 days",
    subtitle = "Adjusted for gestational age, sex, and center",
    x = NULL, y = "Probability"
  )

## =========================
## PANEL C — Survival without BPD (NOW GEE-adjusted, same covariate
## set and clustering as the primary extubation model, for internal
## consistency within the figure)
## =========================
fig3_bpd <- fig3_df %>%
  filter(!is.na(surv36_no_bpd)) %>%
  mutate(surv36_no_bpd_num = as.numeric(surv36_no_bpd == "Yes")) %>%
  arrange(rand_id)

fit_bpd_gee <- geeglm(
  surv36_no_bpd_num ~ phenotype * treat + GESTAGE + male + randcent,
  id = rand_id,
  data = fig3_bpd,
  family = binomial,
  corstr = "exchangeable"
)

emm_bpd_link <- emmeans(fit_bpd_gee, ~ treat | phenotype)
bpd_df <- as.data.frame(emmeans(fit_bpd_gee, ~ treat | phenotype, type = "response")) %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order))

lcl_col_c <- intersect(c("asymp.LCL", "lower.CL"), colnames(bpd_df))[1]
ucl_col_c <- intersect(c("asymp.UCL", "upper.CL"), colnames(bpd_df))[1]

bpd_contrast <- contrast(emm_bpd_link, method = "revpairwise", by = "phenotype")
bpd_p <- as.data.frame(summary(bpd_contrast)) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    label = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.05  ~ sprintf("p = %.3f", p.value),
      TRUE            ~ "NS"
    )
  )

annot_c <- bpd_df %>%
  group_by(phenotype) %>%
  summarise(y = max(.data[[ucl_col_c]]) + 0.06, .groups = "drop") %>%
  left_join(bpd_p, by = "phenotype")

p_fig3c <- ggplot(bpd_df, aes(x = treat, y = prob, fill = treat)) +
  geom_col(width = 0.6, color = "grey40") +
  geom_errorbar(aes(ymin = .data[[lcl_col_c]], ymax = .data[[ucl_col_c]]), width = 0.15) +
  geom_text(aes(label = percent(prob, 0.1)), position = position_stack(vjust = 0.5),
            fontface = "bold", size = 4) +
  geom_text(data = annot_c, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 4) +
  facet_wrap(~ phenotype, nrow = 2) +
  scale_fill_manual(values = treat_cols) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  theme_publication() +
  labs(
    title = "C. Survival without BPD at 36 weeks",
    subtitle = "Adjusted for gestational age, sex, and center",
    x = NULL, y = "Probability"
  )
# NOTE: Panel C's model uses the same covariate set and clustering as
# the primary extubation model for internal consistency, but this GEE
# was newly added for this revision and is not otherwise described in
# Methods -> add a sentence to Methods/Table footnotes noting that
# BPD-survival interaction estimates use the same adjusted approach.

## =========================
## PANEL D — Absolute risk difference (model-based, via regrid to the
## response/probability scale — replaces the manual SE-summation
## approximation used previously)
## =========================

emm_response <- regrid(emm_link, transform = "response")
ard_contrast <- contrast(emm_response, method = "revpairwise", by = "phenotype")
ard_summary <- as.data.frame(summary(ard_contrast, infer = c(TRUE, TRUE)))

lcl_col_d <- intersect(c("asymp.LCL", "lower.CL"), colnames(ard_summary))[1]
ucl_col_d <- intersect(c("asymp.UCL", "upper.CL"), colnames(ard_summary))[1]

rd_df <- ard_summary %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    risk_diff = estimate,
    lower = .data[[lcl_col_d]],
    upper = .data[[ucl_col_d]],
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

p_fig3d <- ggplot(rd_df, aes(x = phenotype, y = risk_diff, fill = phenotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_col(width = 0.65, color = "grey40", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.8) +
  geom_text(
    aes(label = percent(risk_diff, 0.1), vjust = ifelse(risk_diff > 0, -0.7, 1.6)),
    fontface = "bold", size = 4
  ) +
  geom_text(aes(y = upper + 0.04, label = stars), size = 5, fontface = "bold") +
  scale_fill_manual(values = pheno_cols) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0.2, 0.2))) +
  theme_publication() +
  labs(
    title = "D. Absolute risk difference (Hydrocortisone \u2212 Placebo)",
    subtitle = "*** p < 0.001   Adjusted for gestational age, sex, and center",
    x = NULL, y = "Risk difference"
  )

## =========================
## FINAL FIGURE (renamed Figure 3)
## =========================
fig3_final <- (p_fig3a | p_fig3b) / (p_fig3c | p_fig3d) 

fig3_final

ggsave("Figure3_enhanced.pdf", fig3_final,
       width = 12, height = 10, dpi = 300, bg = "white")

ggsave(
  "Figure3_Phenotype_Treatment_Adjusted.png",
  fig3_final,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)



############################################################
# EXPLORATORY ANALYSES TO RECOVER POWER ON BPD-RELATED OUTCOMES
#
# The binary "survival without BPD" GEE showed no significant treatment-
# by-phenotype interaction for any stratum. Two extensions that preserve
# more information from the same underlying data may recover signal that
# dichotomization discarded:
#
#   (1) Ordinal GEE on the full 4-level Jensen BPD severity scale
#       (already visualized in Figure 2B) instead of collapsing to a
#       binary "no BPD" outcome.
#   (2) Continuous GEE on ventilator-free days to 36 weeks (vfd36),
#       already derived in your pipeline but not yet used in the
#       treatment-heterogeneity analysis.
#
# Both use the same adjustment set and clustering strategy as the
# primary extubation model (GESTAGE, sex, center; clustered on rand_id)
# for internal consistency across all treatment-heterogeneity analyses.
############################################################

library(tidyverse)
library(multgee)   # ordinal GEE (proportional-odds, GEE-based)
library(ordinal)   # clmm — GLMM-style sensitivity check for ordinal outcome
library(geepack)   # continuous-outcome GEE
library(lme4)       # LMM sensitivity check for continuous outcome
library(lmerTest)
library(emmeans)

## ------------------------------------------------------------
## 0. Base analytic dataset (same covariate/clustering setup as
##    fit_extub_gee — reuse fig3_df if already in the environment,
##    otherwise rebuild it here)
## ------------------------------------------------------------

fig3_df <- df_pheno %>%
  filter(
    phenotype %in% pheno_order,
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
    randcent = factor(randcent),
    rand_id_num = as.numeric(factor(rand_id))  # multgee needs a numeric/integer cluster id
  )

############################################################
# (1) ORDINAL GEE — BPD severity (Jensen, 4 levels) by phenotype x treat
############################################################

## ------------------------------------------------------------
## 1a. Build the ordinal outcome, ordered worst-to-best or best-to-worst
##     consistently with how Figure 2B displayed it (No support < 2L NC
##     < >2L/NIPPV < IPPV, i.e., increasing severity)
## ------------------------------------------------------------

bpd_ord_df <- fig3_df %>%
  filter(!is.na(bpd_jen36_hc)) %>%
  mutate(
    jensen_ord = factor(
      bpd_jen36_hc,
      levels = 0:3,
      labels = c("No support", "<=2L NC", ">2L/NIPPV", "IPPV"),
      ordered = TRUE
    ),
    jensen_num = as.integer(jensen_ord)  # 1..4, required by multgee
  ) %>%
  arrange(rand_id_num)

cat("Ordinal BPD analytic N:", nrow(bpd_ord_df), "\n")
table(bpd_ord_df$jensen_ord, bpd_ord_df$phenotype)

## ------------------------------------------------------------
## 1b. Primary ordinal GEE (proportional-odds cumulative-link model,
##     GEE estimation, clustered on family). multgee::ordLORgee()
##     handles clustered ordinal responses directly.
## ------------------------------------------------------------

bpd_ord_df <- bpd_ord_df %>%
  group_by(rand_id_num) %>%
  mutate(within_id_index = row_number()) %>%
  ungroup() %>%
  arrange(rand_id_num, within_id_index)

fit_bpd_ordgee <- ordLORgee(
  formula = jensen_num ~ phenotype * treat + GESTAGE + male + randcent,
  data = as.data.frame(bpd_ord_df),
  id = rand_id_num,
  repeated = within_id_index,   # 1, 2, (3, 4) — position within the family cluster,
  # NOT the outcome. Singletons all get index 1;
  # twins get 1 and 2; etc.
  LORstr = "independence"
)

summary(fit_bpd_ordgee)

# NOTE ON INTERACTION TESTING: multgee's summary() gives per-coefficient
# Wald tests but not a built-in multi-df joint Wald test for the full
# phenotype:treat interaction block. Build one manually from the
# coefficient vector and robust covariance matrix:

coefs <- coef(fit_bpd_ordgee)
vcov_mat <- fit_bpd_ordgee$robust.variance  # confirm this slot name exists;
# if it errors, run names(fit_bpd_ordgee)
# to find the correct covariance matrix slot

int_terms <- c(
  "phenotypePhenotype 2:treatHydrocortisone",
  "phenotypePhenotype 3:treatHydrocortisone",
  "phenotypePhenotype 4:treatHydrocortisone"
)

beta_int <- coefs[int_terms]
V_int <- vcov_mat[int_terms, int_terms]

wald_stat <- as.numeric(t(beta_int) %*% solve(V_int) %*% beta_int)
p_interaction_bpd_ord <- pchisq(wald_stat, df = 3, lower.tail = FALSE)

cat("Ordinal BPD severity: phenotype x treatment joint Wald p =",
    round(p_interaction_bpd_ord, 3), "\n")

## ------------------------------------------------------------
## 1c. Sensitivity check — cumulative link MIXED model (clmm) with
##     random intercepts for center and family, mirroring the GLMM
##     sensitivity approach already used for the primary binary outcome
## ------------------------------------------------------------

fit_bpd_clmm <- clmm(
  jensen_ord ~ phenotype * treat + GESTAGE + male + (1 | randcent) + (1 | rand_id),
  data = bpd_ord_df
)

fit_bpd_clmm_noint <- clmm(
  jensen_ord ~ phenotype + treat + GESTAGE + male + (1 | randcent) + (1 | rand_id),
  data = bpd_ord_df
)

lrt_bpd_clmm <- anova(fit_bpd_clmm_noint, fit_bpd_clmm)
print(lrt_bpd_clmm)
p_interaction_bpd_clmm <- lrt_bpd_clmm$`Pr(>Chisq)`[2]
cat("Ordinal BPD severity (clmm sensitivity): interaction LRT p =",
    round(p_interaction_bpd_clmm, 3), "\n")

############################################################
# (2) CONTINUOUS GEE — Ventilator-free days to 36 weeks (vfd36)
############################################################

## ------------------------------------------------------------
## 2a. Analytic dataset
## ------------------------------------------------------------

vfd_df <- fig3_df %>%
  filter(!is.na(vfd36)) %>%
  arrange(rand_id)

cat("VFD analytic N:", nrow(vfd_df), "\n")

## ------------------------------------------------------------
## 2b. Primary continuous GEE — identity link, gaussian family,
##     same adjustment/clustering as the primary extubation model
## ------------------------------------------------------------

fit_vfd_gee <- geeglm(
  vfd36 ~ phenotype * treat + GESTAGE + male + randcent,
  id = rand_id,
  data = vfd_df,
  family = gaussian,
  corstr = "exchangeable"
)

summary(fit_vfd_gee)

anova_vfd <- anova(fit_vfd_gee)
print(anova_vfd)
p_interaction_vfd <- anova_vfd["phenotype:treat", "P(>|Chi|)"]
cat("VFD36: phenotype x treatment interaction (GEE Wald) p =",
    round(p_interaction_vfd, 3), "\n")

## Phenotype-specific treatment effect (mean difference in VFD, adjusted)
emm_vfd <- emmeans(fit_vfd_gee, ~ treat | phenotype)
vfd_contrast <- contrast(emm_vfd, method = "revpairwise", by = "phenotype")
vfd_summary <- as.data.frame(summary(vfd_contrast, infer = c(TRUE, TRUE)))

lcl_col_v <- intersect(c("asymp.LCL", "lower.CL"), colnames(vfd_summary))[1]
ucl_col_v <- intersect(c("asymp.UCL", "upper.CL"), colnames(vfd_summary))[1]

vfd_df_results <- vfd_summary %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    label = sprintf("%.1f days (%.1f to %.1f)", estimate, .data[[lcl_col_v]], .data[[ucl_col_v]])
  ) %>%
  select(phenotype, label, p.value)

print(vfd_df_results)

## ------------------------------------------------------------
## 2c. Sensitivity check — linear mixed model with crossed random
##     intercepts for center and family
## ------------------------------------------------------------

fit_vfd_lmm <- lmer(
  vfd36 ~ phenotype * treat + GESTAGE + male + (1 | randcent) + (1 | rand_id),
  data = vfd_df,
  REML = FALSE
)

fit_vfd_lmm_noint <- lmer(
  vfd36 ~ phenotype + treat + GESTAGE + male + (1 | randcent) + (1 | rand_id),
  data = vfd_df,
  REML = FALSE
)

lrt_vfd_lmm <- anova(fit_vfd_lmm_noint, fit_vfd_lmm)
print(lrt_vfd_lmm)
p_interaction_vfd_lmm <- lrt_vfd_lmm$`Pr(>Chisq)`[2]
cat("VFD36 (LMM sensitivity): interaction LRT p =",
    round(p_interaction_vfd_lmm, 3), "\n")

print(VarCorr(fit_vfd_lmm))  # check whether family/center variance is nonzero here too

############################################################
# SUMMARY TABLE FOR BOTH ANALYSES
############################################################

summary_new_analyses <- tibble(
  Outcome = c("BPD severity (ordinal, 4-level Jensen)", "Ventilator-free days to 36wk"),
  `Primary model` = c("Ordinal GEE (multgee)", "Continuous GEE (geepack)"),
  `Interaction P (primary)` = c(round(p_interaction_bpd_ord, 3), round(p_interaction_vfd, 3)),
  `Sensitivity model` = c("clmm (cumulative-link mixed)", "lmer (linear mixed)"),
  `Interaction P (sensitivity)` = c(round(p_interaction_bpd_clmm, 3), round(p_interaction_vfd_lmm, 3))
)

print(summary_new_analyses)

# If either interaction test is now significant (or notably lower than
# the binary BPD-survival P-value), this supports the interpretation
# that dichotomizing BPD severity / restricting to a binary composite
# discarded real signal. If both remain non-significant and similar in
# magnitude to the binary result, that strengthens the case that the
# null BPD finding reflects a genuine absence of a detectable phenotype-
# specific effect on this endpoint rather than a power/measurement
# artifact -- also a useful, honest thing to be able to say.


library(cluster)
library(fpc)

# dist_mat and hc should already exist from your clustering pipeline
# (dist_mat = dist() object used for hclust; hc = hclust(dist_mat, method="ward.D2"))

k_range <- 2:6
ch_values <- numeric(length(k_range))
silhouette_values <- numeric(length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  cl_k <- cutree(hc, k = k)
  
  # Calinski-Harabasz index
  ch_values[i] <- cluster.stats(dist_mat, cl_k)$ch
  
  # Mean silhouette width (for cross-check against your existing values)
  silhouette_values[i] <- mean(silhouette(cl_k, dist_mat)[, 3])
}

ch_table <- data.frame(
  k = k_range,
  Mean_Silhouette_Width = round(silhouette_values, 3),
  Calinski_Harabasz = round(ch_values, 1)
)

print(ch_table)


############################################################
# SENSITIVITY ANALYSES — Reviewer 2, Comment 1
# "Additional work to demonstrate/understand if phenotype 2 is
#  distinct from phenotype 3, and whether treatment heterogeneity
#  results might differ if there is minimal distinction."
#
# (1) PRIMARY: re-cluster WITHOUT 5-minute Apgar score, compare
#     structure to the original 4-variable solution, and re-fit the
#     adjusted treatment-heterogeneity model on the new assignments.
# (2) SECONDARY: force-merge original Phenotypes 2 and 3, re-fit the
#     adjusted model on a 3-level phenotype variable, and check
#     whether the Phenotype 2 treatment signal survives/attenuates.
#
# Assumes cluster_df_complete, cluster_scaled, df_pheno, pheno_cols,
# pheno_order, fig3_df, and fit_extub_gee already exist from prior
# scripts in this analysis (see earlier sections of the pipeline).
############################################################

library(tidyverse)
library(cluster)
library(fpc)
library(mclust)     # adjustedRandIndex
library(geepack)
library(emmeans)
library(patchwork)

set.seed(2024)

############################################################
# PART 1 (PRIMARY): RE-CLUSTER WITHOUT APGAR
############################################################

## ------------------------------------------------------------
## 1a. Build the reduced variable set and re-cluster
## ------------------------------------------------------------

clust_vars_noapgar <- c("GESTAGE", "BIRTHWTX", "age_rndaysX", "map_fio2X")

cluster_noapgar <- cluster_df_complete %>%
  select(rand_id, all_of(clust_vars_noapgar)) %>%
  mutate(across(-rand_id, scale))

dist_noapgar <- dist(cluster_noapgar %>% select(-rand_id))
hc_noapgar <- hclust(dist_noapgar, method = "ward.D2")

## ------------------------------------------------------------
## 1b. Check whether k=4 is still supported without Apgar, using the
##     same internal-validity indices as the primary analysis
## ------------------------------------------------------------

k_range <- 2:6
diag_noapgar <- map_dfr(k_range, function(k) {
  cl_k <- cutree(hc_noapgar, k = k)
  tibble(
    k = k,
    mean_silhouette = mean(silhouette(cl_k, dist_noapgar)[, 3]),
    ch_index = cluster.stats(dist_noapgar, cl_k)$ch
  )
})
print(diag_noapgar)

## Retain k=4 for direct comparability to the primary analysis
phen_noapgar <- cutree(hc_noapgar, k = 4)

## ------------------------------------------------------------
## 1c. Compare structure to the ORIGINAL (Apgar-included) solution
## ------------------------------------------------------------

phen_original <- cutree(hc, k = 4)  # from the primary pipeline (hc = original hclust object)

ari_noapgar_vs_original <- adjustedRandIndex(phen_original, phen_noapgar)
cat("ARI: original (with Apgar) vs. no-Apgar re-clustering =",
    round(ari_noapgar_vs_original, 3), "\n")
# High ARI (>0.5-0.6) => structure is robust to dropping Apgar.
# Low ARI => Apgar was a major driver of the original partition,
# which would itself be an important, honest finding to report.

## Cross-tabulation for the reviewer response (mirrors sTable5's format)
crosstab_noapgar <- table(
  Original = phen_original,
  NoApgar = phen_noapgar
)
print(crosstab_noapgar)

## ------------------------------------------------------------
## 1d. Characterize the new phenotypes (centroid profile), to see
##     whether the P2/P3-like split persists or collapses
## ------------------------------------------------------------

noapgar_assignments <- tibble(
  rand_id = cluster_df_complete$rand_id,
  phenotype_noapgar = factor(phen_noapgar)
)

centroid_noapgar <- cluster_df_complete %>%
  left_join(noapgar_assignments, by = "rand_id") %>%
  select(phenotype_noapgar, all_of(clust_vars_noapgar)) %>%
  mutate(across(-phenotype_noapgar, scale)) %>%
  group_by(phenotype_noapgar) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-phenotype_noapgar, names_to = "Variable", values_to = "Zscore")

p_centroid_noapgar <- ggplot(centroid_noapgar,
                             aes(x = Variable, y = phenotype_noapgar, fill = Zscore)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Zscore, 2)), size = 3) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0) +
  theme_minimal() +
  labs(title = "Standardized Feature Centroids: No-Apgar Re-Clustering (k=4)",
       x = NULL, y = "Cluster", fill = "Mean\n(z-score)")

p_centroid_noapgar

## ------------------------------------------------------------
## 1e. Re-fit the adjusted treatment-heterogeneity model on the
##     NEW (no-Apgar) phenotype assignments
## ------------------------------------------------------------

df_pheno_noapgar <- df_pheno %>%
  left_join(
    noapgar_assignments %>% rename(rand_id_join = rand_id),
    by = c("rand_id" = "rand_id_join")
  )

fig_noapgar_df <- df_pheno_noapgar %>%
  filter(
    !is.na(phenotype_noapgar),
    !is.na(treat),
    !is.na(GESTAGE),
    !is.na(MALE),
    !is.na(rand_id),
    !is.na(randcent)
  ) %>%
  mutate(
    treat = factor(treat, levels = c("Placebo", "Hydrocortisone")),
    male = factor(MALE, levels = c(0, 1), labels = c("Female", "Male")),
    extub14_num = as.numeric(extub14 == "Yes"),
    randcent = factor(randcent),
    rand_id = factor(rand_id)
  ) %>%
  arrange(rand_id)

fit_extub_noapgar_gee <- geeglm(
  extub14_num ~ phenotype_noapgar * treat + GESTAGE + male + randcent,
  id = rand_id,
  data = fig_noapgar_df,
  family = binomial,
  corstr = "exchangeable"
)

summary(fit_extub_noapgar_gee)

anova_noapgar <- anova(fit_extub_noapgar_gee)
print(anova_noapgar)
p_interaction_noapgar <- anova_noapgar["phenotype_noapgar:treat", "P(>|Chi|)"]
cat("No-Apgar re-clustering: phenotype x treatment interaction (GEE Wald) p =",
    round(p_interaction_noapgar, 3), "\n")

## Phenotype-specific treatment effects (RR + OR) under the new clustering
emm_link_noapgar <- emmeans(fit_extub_noapgar_gee, ~ treat | phenotype_noapgar)

# Risk ratios
emm_log_noapgar <- regrid(emm_link_noapgar, transform = "log")
rr_contrast_noapgar <- contrast(emm_log_noapgar, method = "revpairwise",
                                by = "phenotype_noapgar", type = "response")
rr_noapgar_df <- as.data.frame(summary(rr_contrast_noapgar, infer = c(TRUE, TRUE)))
print(rr_noapgar_df)

# Odds ratios (native scale)
or_contrast_noapgar <- contrast(emm_link_noapgar, method = "revpairwise",
                                by = "phenotype_noapgar", type = "response")
or_noapgar_df <- as.data.frame(summary(or_contrast_noapgar, infer = c(TRUE, TRUE)))
print(or_noapgar_df)

# NOTE: whichever new cluster most closely resembles original Phenotype 2
# (identify via the crosstab above -- likely the no-Apgar cluster with the
# largest overlap with phen_original == 2) is the one to report as the
# "does the signal survive without Apgar" answer for the reviewer.

############################################################
# PART 2 (SECONDARY): FORCE-MERGE ORIGINAL PHENOTYPES 2 AND 3
############################################################

## ------------------------------------------------------------
## 2a. Collapse Phenotype 2 and Phenotype 3 into a single stratum
## ------------------------------------------------------------

df_pheno_merged <- df_pheno %>%
  mutate(
    phenotype_merged = case_when(
      phenotype %in% c("Phenotype 2", "Phenotype 3") ~ "Phenotype 2+3 (merged)",
      TRUE ~ as.character(phenotype)
    ),
    phenotype_merged = factor(
      phenotype_merged,
      levels = c("Phenotype 1", "Phenotype 2+3 (merged)", "Phenotype 4")
    )
  )

fig_merged_df <- df_pheno_merged %>%
  filter(
    !is.na(phenotype_merged),
    !is.na(treat),
    !is.na(GESTAGE),
    !is.na(MALE),
    !is.na(rand_id),
    !is.na(randcent)
  ) %>%
  mutate(
    treat = factor(treat, levels = c("Placebo", "Hydrocortisone")),
    male = factor(MALE, levels = c(0, 1), labels = c("Female", "Male")),
    extub14_num = as.numeric(extub14 == "Yes"),
    randcent = factor(randcent),
    rand_id = factor(rand_id)
  ) %>%
  arrange(rand_id)

cat("Merged-phenotype analytic N by group:\n")
print(table(fig_merged_df$phenotype_merged))

## ------------------------------------------------------------
## 2b. Re-fit the adjusted model with 3-level (merged) phenotype
## ------------------------------------------------------------

fit_extub_merged_gee <- geeglm(
  extub14_num ~ phenotype_merged * treat + GESTAGE + male + randcent,
  id = rand_id,
  data = fig_merged_df,
  family = binomial,
  corstr = "exchangeable"
)

summary(fit_extub_merged_gee)

anova_merged <- anova(fit_extub_merged_gee)
print(anova_merged)
p_interaction_merged <- anova_merged["phenotype_merged:treat", "P(>|Chi|)"]
cat("Merged Phenotype 2+3: phenotype x treatment interaction (GEE Wald) p =",
    round(p_interaction_merged, 3), "\n")

emm_link_merged <- emmeans(fit_extub_merged_gee, ~ treat | phenotype_merged)

emm_log_merged <- regrid(emm_link_merged, transform = "log")
rr_contrast_merged <- contrast(emm_log_merged, method = "revpairwise",
                               by = "phenotype_merged", type = "response")
rr_merged_df <- as.data.frame(summary(rr_contrast_merged, infer = c(TRUE, TRUE)))
print(rr_merged_df)

# INTERPRETATION GUIDE:
# - If the merged "Phenotype 2+3" stratum shows an attenuated, non-significant
#   treatment effect relative to original Phenotype 2 alone (RR 1.83, P<0.001),
#   this indicates the treatment signal was driven specifically by the
#   original Phenotype 2 subgroup and is diluted by combining with Phenotype 3
#   -- i.e., the split between them IS clinically/statistically meaningful
#   for treatment heterogeneity, supporting retention of 4 phenotypes.
# - If the merged stratum RETAINS a similar significant effect size, this
#   suggests the signal is robust to the P2/P3 boundary and less sensitive
#   to exact cluster granularity at k=4 vs k=3 in this region of the data.

############################################################
# SUMMARY TABLE FOR REVIEWER RESPONSE
############################################################

summary_sensitivity <- tibble(
  Analysis = c("Original (4 phenotypes, Apgar included)",
               "No-Apgar re-clustering (4 clusters)",
               "Merged Phenotype 2+3 (3 groups)"),
  `ARI vs. original` = c(1.00, round(ari_noapgar_vs_original, 3), NA),
  `Global interaction P` = c(0.12, round(p_interaction_noapgar, 3), round(p_interaction_merged, 3))
)

print(summary_sensitivity)
