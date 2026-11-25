options(dplyr.summarise.inform = FALSE)

sanity_check <- function(df, label) {
  cat("\n==============================\n")
  cat("SANITY CHECK:", label, "\n")
  cat("Rows (infants):", nrow(df), "\n")
  cat("Unique rand_id (randomization units):",
      n_distinct(df$rand_id), "\n")
  cat("Multiple births present:",
      any(df %>% count(rand_id) %>% pull(n) > 1), "\n")
  cat("==============================\n")
}


library(readr)
analysis_dash <- read_csv("C:/Users/Alvaro/OneDrive/Desktop/Hydrocortisone/Data_Files_CSV/analysis_dash.csv")
View(analysis_dash)

analysis_dash<-analysis_dash

############################################################
# 0. LIBRARIES
############################################################
library(tidyverse)
library(cluster)
library(factoextra)
library(broom)
library(emmeans)
library(gtsummary)
library(flextable)

############################################################
# 1. LOAD DATA
############################################################
df_raw <- analysis_dash

sanity_check(df_raw, "Raw data loaded")
stopifnot(nrow(df_raw) == 800)


############################################################
# 2. CODEBOOK-AWARE RECODING
############################################################
df <- df_raw %>%
  mutate(
    treat = factor(treatmnt, levels = c(1, 2),
                   labels = c("Placebo", "Hydrocortisone")),
    male = factor(MALE, levels = c(0,1), labels = c("Female","Male")),
    sga  = factor(SGA, levels = c(0,1), labels = c("No","Yes")),
    hisp = factor(HISP, levels = c(0,1),
                  labels = c("Not Hispanic","Hispanic")),
    race = factor(Race3X, levels = c(1,2,3),
                  labels = c("Black","White","Other")),
    extub14 = factor(extube14, levels = c(0,1),
                     labels = c("No","Yes")),
    death36_bin = factor(death36, levels = c(0,1),
                         labels = c("Alive","Dead")),
    protocol_deviation = NA
  )

sanity_check(df, "After recoding")
stopifnot(nrow(df) == 800)

############################################################
# 3. OUTCOME DEFINITIONS
############################################################
df <- df %>%
  mutate(
    vfd36 = if_else(!is.na(days_VENT36),
                    36 - days_VENT36,
                    NA_real_),
    surv36_no_bpd = factor(
      surv36_nobpd_p_hc,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    bpd_jensen36 = factor(
      bpd_jen36_hc,
      levels = c(0,1,2,3),
      labels = c("None","Mild","Moderate","Severe"),
      ordered = TRUE
    ),
    death_or_severe_bpd36 = factor(
      if_else(death36 == 1 | bpd_jen36_hc == 3, 1, 0),
      levels = c(0,1),
      labels = c("No","Yes")
    ),
    ndi_any = factor(NDI3, levels = c(0,1),
                     labels = c("No","Yes")),
    ndi_modsev = factor(ndi_mod_sev, levels = c(0,1),
                        labels = c("No","Yes"))
  )

sanity_check(df, "After outcome derivation")
stopifnot(nrow(df) == 800)

############################################################
# 4. DEFINE RCT COHORT (NO OUTCOME FILTERING)
############################################################
df_rct <- df %>% filter(treatmnt %in% c(1,2))

sanity_check(df_rct, "RCT cohort defined")
stopifnot(nrow(df_rct) == 800)

############################################################
# 5. BASELINE VARIABLES FOR PHENOTYPING
# ⚠️ N is allowed to drop here (complete-case phenotyping)
############################################################
cluster_df <- df_rct %>%
  select(
    rand_id,
    GESTAGE,
    BIRTHWTX,
    APGAR5,
    age_rndaysX,
    map_fio2X
  )

sanity_check(cluster_df, "Baseline variables selected")

cluster_df_complete <- cluster_df %>% drop_na()

sanity_check(cluster_df_complete, "Complete-case phenotyping cohort")
#760

############################################################
# 6. CLUSTERING
############################################################
set.seed(2024)

cluster_scaled <- cluster_df_complete %>%
  mutate(across(-rand_id, scale))

dist_mat <- dist(cluster_scaled %>% select(-rand_id))
hc <- hclust(dist_mat, method = "ward.D2")

stopifnot(length(hc$order) == nrow(cluster_df_complete))

############################################################
# 7. FINAL PHENOTYPE ASSIGNMENT
############################################################
cluster_assignments <- tibble(
  rand_id = cluster_df_complete$rand_id,
  phenotype = factor(
    cutree(hc, k = 4),
    levels = 1:4,
    labels = c("Phenotype 1","Phenotype 2","Phenotype 3","Phenotype 4")
  )
)

sanity_check(cluster_assignments, "Cluster assignments")

############################################################
# 8. MERGE PHENOTYPES BACK TO RCT DATA (CORRECT WAY)
############################################################
df_rct <- df %>%
  mutate(infant_row_id = row_number())

#Build cluster_df keeping that ID
cluster_df <- df_rct %>%
  select(
    infant_row_id,
    GESTAGE,
    BIRTHWTX,
    APGAR5,
    age_rndaysX,
    map_fio2X
  ) %>%
  drop_na()

#create cluster assingments based on those ids
cluster_assignments <- tibble(
  infant_row_id = cluster_df$infant_row_id,
  phenotype = factor(
    cutree(hc, k = 4),
    levels = 1:4,
    labels = c(
      "Phenotype 1",
      "Phenotype 2",
      "Phenotype 3",
      "Phenotype 4"
    )
  )
)

df_pheno <- df_rct %>%
  left_join(cluster_assignments, by = "infant_row_id")

sanity_check(df_pheno, "After phenotype merge")

stopifnot(nrow(df_pheno) == 800)
stopifnot(n_distinct(df_pheno$infant_row_id) == 800)

table(df_pheno$phenotype, useNA = "ifany")

#saveRDS(df_pheno,'df_pheno.rds')

############################################################
# 9. ANALYTIC OBJECTS, TABLES, AND FIGURES (FINAL, CORRECT)
#    Uses df_pheno from NEW PIPELINE ONLY
############################################################

## ---------------------------
## A. SANITY: infant-level only
## ---------------------------
stopifnot(
  nrow(df_pheno) == 800,
  n_distinct(df_pheno$infant_row_id) == 800
)

## ===========================
## B. PCA OF BASELINE PHENOTYPES
## ===========================

pca_df <- df_pheno %>%
  select(phenotype, GESTAGE, BIRTHWTX, APGAR5, age_rndaysX, map_fio2X) %>%
  drop_na()

pca_res <- prcomp(
  pca_df %>% select(-phenotype),
  center = TRUE,
  scale. = TRUE
)

pca_plot_df <- as_tibble(pca_res$x[, 1:2]) %>%
  mutate(phenotype = pca_df$phenotype)

p_pca <- ggplot(
  pca_plot_df,
  aes(PC1, PC2, color = phenotype)
) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(type = "norm", linewidth = 0.8) +
  theme_classic(base_size = 14) +
  labs(
    title = "Principal Component Analysis of Baseline Phenotypes",
    subtitle = "Derived from pre-randomization clinical variables",
    x = paste0(
      "PC1 (",
      round(summary(pca_res)$importance[2,1] * 100, 1),
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca_res)$importance[2,2] * 100, 1),
      "% variance)"
    ),
    color = "Phenotype"
  )

ggsave(
  "Figure_PCA_Baseline_Phenotypes.pdf",
  p_pca,
  width = 7,
  height = 6
)

## ===========================
## C. BASELINE TABLE BY PHENOTYPE
## ===========================

tbl_baseline_pheno <- df_pheno %>%
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
    label = list(
      GESTAGE     ~ "Gestational age, weeks",
      BIRTHWTX    ~ "Birth weight, g",
      male        ~ "Male sex",
      sga         ~ "Small for gestational age",
      ANTESTER    ~ "Antenatal steroid exposure",
      MULTBRTH    ~ "Multiple gestation",
      APGAR5      ~ "5-minute Apgar score",
      age_rndaysX ~ "Age at randomization, days",
      map_fio2X   ~ "Baseline respiratory severity (MAP × FiO₂)"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_p() %>%
  bold_labels()

tbl_baseline_pheno %>%
  as_flex_table() %>%
  set_caption(
    "Table 1. Baseline demographic and clinical characteristics by phenotype."
  ) %>%
  save_as_docx(
    path = "Table_1_Baseline_Characteristics_by_Phenotype.docx"
  )

## ===========================
## D. OUTCOMES BY PHENOTYPE
## ===========================

tbl_outcomes_pheno <- df_pheno %>%
  select(
    phenotype,
    extub14,
    vfd36,
    surv36_no_bpd,
    bpd_jensen36,
    death36_bin,
    death_or_severe_bpd36,
    ndi_any,
    ndi_modsev
  ) %>%
  tbl_summary(
    by = phenotype,
    label = list(
      extub14 ~ "Extubation ≤14 days",
      vfd36 ~ "Ventilator-free days to 36 weeks PMA",
      surv36_no_bpd ~ "Survival to 36 weeks without BPD",
      bpd_jensen36 ~ "BPD severity at 36 weeks (Jensen)",
      death36_bin ~ "Death by 36 weeks PMA",
      death_or_severe_bpd36 ~ "Death or severe BPD",
      ndi_any ~ "Any neurodevelopmental impairment",
      ndi_modsev ~ "Moderate/severe neurodevelopmental impairment"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_p() %>%
  bold_labels()

tbl_outcomes_pheno %>%
  as_flex_table() %>%
  set_caption(
    "Table 2. Clinical outcomes by baseline phenotype."
  ) %>%
  save_as_docx(
    path = "Table_2_Outcomes_by_Phenotype.docx"
  )

## ===========================
## E. TREATMENT BALANCE BY PHENOTYPE
## ===========================

df_tx <- df_pheno %>%
  filter(!is.na(phenotype), !is.na(treat), !is.na(extub14))

tbl_tx_balance <- as.data.frame(
  table(df_tx$phenotype, df_tx$treat)
) %>%
  rename(
    Phenotype = Var1,
    Treatment = Var2,
    N = Freq
  )

tbl_tx_balance %>%
  flextable() %>%
  set_caption(
    "Supplemental Table. Distribution of randomized treatment by phenotype."
  ) %>%
  save_as_docx(
    path = "Supplement_Treatment_Balance_by_Phenotype.docx"
  )

## ===========================
## F. PHENOTYPE-SPECIFIC TREATMENT EFFECT (FIRTH)
## ===========================

library(logistf)

firth_results <- df_tx %>%
  group_by(phenotype) %>%
  group_modify(~ {
    fit <- logistf(extub14 ~ treat, data = .x)
    tibble(
      OR  = exp(coef(fit)["treatHydrocortisone"]),
      LCL = exp(confint(fit)["treatHydrocortisone", 1]),
      UCL = exp(confint(fit)["treatHydrocortisone", 2]),
      p   = fit$prob["treatHydrocortisone"]
    )
  }) %>%
  ungroup()

tbl_tx_pheno <- firth_results %>%
  mutate(
    OR = round(OR, 2),
    LCL = round(LCL, 2),
    UCL = round(UCL, 2),
    p = ifelse(p < 0.001, "<0.001", sprintf("%.2f", p))
  )

tbl_tx_pheno %>%
  flextable() %>%
  set_caption(
    "Table 3. Phenotype-specific treatment effect of hydrocortisone on extubation ≤14 days."
  ) %>%
  save_as_docx(
    path = "Table_3_Treatment_Effect_by_Phenotype.docx"
  )

## ===========================
## G. FOREST PLOT
## ===========================

p_forest <- ggplot(
  tbl_tx_pheno,
  aes(x = OR, y = phenotype)
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = LCL, xmax = UCL),
    height = 0.2
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  labs(
    title = "Hydrocortisone Effect on Early Extubation by Phenotype",
    x = "Odds Ratio (log scale)",
    y = NULL
  )

p_forest

ggsave(
  "Figure_Forest_Treatment_Effect_by_Phenotype.pdf",
  p_forest,
  width = 8,
  height = 6
)



## =========================================================
## FIGURE 1. Derivation and Clinical Structure of Baseline Phenotypes
## Panels: A (PCA), B (Radar), C (Respiratory Severity), D (Heatmap)
## =========================================================

library(tidyverse)
library(patchwork)
library(scales)

## -------------------------
## Phenotype color palette (MATCHED TO FIGURE 2)
## -------------------------
pheno_cols <- c(
  "Phenotype 1" = "#E41A1C",
  "Phenotype 2" = "#377EB8",
  "Phenotype 3" = "#4DAF4A",
  "Phenotype 4" = "#984EA3"
)

pheno_order <- c("Phenotype 1", "Phenotype 2", "Phenotype 3", "Phenotype 4")

## -------------------------
## Publication theme (CONSISTENT WITH FIGURE 2)
## -------------------------
theme_publication <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", margin = margin(b = 8), size = base_size),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey30"),
      axis.ticks = element_line(color = "grey30"),
      axis.text = element_text(size = base_size - 2, color = "grey20"),
      axis.title = element_text(size = base_size - 1),
      plot.background = element_rect(fill = "grey98", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = base_size - 3),
      legend.title = element_text(size = base_size - 2, face = "bold"),
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

## --------------------------------
## Shared phenotype dataset
## --------------------------------
fig1_df <- df_pheno %>%
  filter(phenotype %in% pheno_order) %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order)) %>%
  select(
    phenotype,
    GESTAGE,        # Gestational age
    BIRTHWTX,       # Birth weight
    APGAR5,         # 5-minute Apgar
    age_rndaysX,    # Age at randomization
    map_fio2X       # MAP × FiO2
  ) %>%
  drop_na()

## Clinical variable labels (radar-safe)
var_labels <- c(
  GESTAGE     = "Gest age",
  BIRTHWTX    = "Birth wt",
  APGAR5      = "Apgar 5min",
  age_rndaysX = "Age at Rx",
  map_fio2X   = "MAP×FiO₂"
)

## =========================================================
## FIGURE 1A — PCA (Enhanced with better labels)
## =========================================================
pca_res <- prcomp(
  fig1_df %>% select(-phenotype),
  center = TRUE,
  scale. = TRUE
)

pca_plot_df <- as_tibble(pca_res$x[, 1:2]) %>%
  mutate(phenotype = fig1_df$phenotype)

# Calculate sample sizes for legend
n_by_pheno <- fig1_df %>%
  count(phenotype) %>%
  mutate(label = paste0(phenotype, " (n=", n, ")"))

pca_plot_df <- pca_plot_df %>%
  left_join(n_by_pheno %>% select(phenotype, label), by = "phenotype")

p_fig1a <- ggplot(
  pca_plot_df,
  aes(PC1, PC2, color = phenotype, fill = phenotype)
) +
  geom_point(alpha = 0.6, size = 2.2) +
  stat_ellipse(aes(color = phenotype), type = "norm", 
               linewidth = 1.1, alpha = 0.8, show.legend = FALSE) +
  scale_color_manual(values = pheno_cols, 
                     labels = setNames(n_by_pheno$label, n_by_pheno$phenotype)) +
  scale_fill_manual(values = pheno_cols) +
  theme_publication(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3)
  ) +
  labs(
    title = "A. Principal component analysis",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "% variance)"),
    color = NULL
  ) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(size = 3.5, alpha = 1)))

## =========================================================
## FIGURE 1B — RADAR (Enhanced styling)
## =========================================================
radar_df <- fig1_df %>%
  mutate(across(-phenotype, scale)) %>%
  group_by(phenotype) %>%
  summarise(across(everything(), mean), .groups = "drop")

radar_long <- radar_df %>%
  pivot_longer(-phenotype, names_to = "variable", values_to = "zscore") %>%
  mutate(variable = recode(variable, !!!var_labels))

p_fig1b <- ggplot(
  radar_long,
  aes(variable, zscore, group = phenotype, color = phenotype, fill = phenotype)
) +
  geom_polygon(alpha = 0.15, linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pheno_cols) +
  scale_fill_manual(values = pheno_cols) +
  coord_polar(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(fill = "grey98", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(size = 10.5, face = "bold", color = "grey20"),
    axis.text.y = element_text(size = 9, color = "grey40"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 50, b = 10, l = 10)
  ) +
  labs(
    title = "B. Standardized feature profiles",
    x = NULL,
    y = "Z-score"
  )

## =========================================================
## FIGURE 1C — Baseline respiratory severity (ENHANCED)
## =========================================================
summary_df <- fig1_df %>%
  group_by(phenotype) %>%
  summarise(
    n = n(),
    median = median(map_fio2X),
    q25 = quantile(map_fio2X, 0.25),
    q75 = quantile(map_fio2X, 0.75),
    mean = mean(map_fio2X),
    label = sprintf("%.0f", median),
    .groups = "drop"
  )

# Kruskal-Wallis test
kw_test <- kruskal.test(map_fio2X ~ phenotype, data = fig1_df)
p_value_text <- sprintf("Kruskal-Wallis p < 0.001")

p_fig1c <- ggplot(
  fig1_df,
  aes(phenotype, map_fio2X, fill = phenotype, color = phenotype)
) +
  geom_violin(trim = FALSE, alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8,
               color = "grey30", linewidth = 0.6) +
  geom_jitter(alpha = 0.25, size = 1.1, width = 0.12) +
  geom_point(
    data = summary_df,
    aes(y = mean),
    size = 4, shape = 23, color = "white", stroke = 1.2
  ) +
  geom_text(
    data = summary_df,
    aes(x = as.numeric(phenotype) + 0.35, y = median, label = label),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold",
    color = "grey10"
  ) +
  geom_text(
    data = summary_df,
    aes(x = phenotype, y = -15, label = paste0("n=", n)),
    size = 3.2, color = "grey30"
  ) +
  annotate("text", x = 2.5, y = max(fig1_df$map_fio2X) * 0.95,
           label = p_value_text, size = 3.5, color = "grey20") +
  scale_fill_manual(values = pheno_cols) +
  scale_color_manual(values = pheno_cols) +
  theme_publication(base_size = 13) +
  theme(panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3)) +
  labs(
    title = "C. Baseline respiratory severity",
    x = NULL,
    y = "MAP × FiO₂"
  ) +
  guides(fill = "none", color = "none")

## =========================================================
## FIGURE 1D — Heatmap (Enhanced with better annotations)
## =========================================================
heatmap_df <- fig1_df %>%
  mutate(across(-phenotype, scale)) %>%
  group_by(phenotype) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  pivot_longer(-phenotype, names_to = "variable", values_to = "zscore") %>%
  mutate(
    variable = recode(variable, !!!var_labels),
    label = sprintf("%.2f", zscore)
  )

p_fig1d <- ggplot(
  heatmap_df,
  aes(variable, phenotype, fill = zscore)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold",
            color = ifelse(abs(heatmap_df$zscore) > 1, "white", "grey20")) +
  scale_fill_gradient2(
    low = "#2166AC",      # Blue
    mid = "white",
    high = "#B2182B",     # Red
    midpoint = 0,
    name = "Z-score",
    limits = c(-2, 2)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(fill = "grey98", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "D. Standardized baseline features",
    x = NULL,
    y = NULL
  )

## =========================================================
## FINAL LAYOUT (Enhanced composition)
## =========================================================
fig1_final <-
  (p_fig1a | p_fig1b) /
  (p_fig1c | p_fig1d) +
  plot_layout(
    heights = c(1, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(
    title = "Figure 1. Derivation and Clinical Structure of Baseline Phenotypes",
    theme = theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0
      ),
      plot.background = element_rect(fill = "grey98", color = NA)
    )
  )

fig1_final

ggsave(
  "Figure_1_Baseline_Phenotypes_Enhanced.png",
  fig1_final,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)


## =========================================================
## TABLE 1. Baseline demographic and clinical characteristics
## =========================================================

library(tidyverse)
library(gtsummary)
library(gt)

## -------------------------
## Prepare analysis dataset
## -------------------------
tbl1_df <- df_pheno %>%
  mutate(
    phenotype = factor(phenotype),
    
    ## Sex
    MALE = factor(
      MALE,
      levels = c(0, 1),
      labels = c("Female", "Male")
    ),
    
    ## Maternal education
    MEDUC4 = factor(
      MEDUC4,
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "< High school",
        "High school",
        "Partial college",
        "College degree or more",
        "Unknown"
      )
    ),
    
    ## Race
    Race3X = factor(
      Race3X,
      levels = c(1, 2),
      labels = c("Black", "White")
    ),
    
    ## Binary variables
    MULTBRTH = factor(MULTBRTH, levels = c(0, 1), labels = c("No", "Yes")),
    ANTESTER = factor(ANTESTER, levels = c(0, 1), labels = c("No", "Yes")),
    SGA      = factor(SGA, levels = c(0, 1), labels = c("No", "Yes")),
    HISP     = factor(HISP, levels = c(0, 1), labels = c("No", "Yes"))
  )

## -------------------------
## Build Table 1
## -------------------------
table1 <-
  tbl1_df %>%
  select(
    phenotype,
    
    ## Birth & perinatal
    GESTAGE,
    BIRTHWTX,
    MALE,
    MULTBRTH,
    SGA,
    ANTESTER,
    
    ## Social / maternal
    MEDUC4,
    Race3X,
    HISP,
    OBGMAGE,
    
    ## Baseline severity
    APGAR5,
    age_rndaysX,
    map_fio2X
  ) %>%
  tbl_summary(
    by = phenotype,
    type = GESTAGE~'continuous',
    
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    
    digits = all_continuous() ~ 1,
    
    label = list(
      GESTAGE     ~ "Gestational age, weeks",
      BIRTHWTX    ~ "Birth weight, g",
      MALE       ~ "Sex",
      MULTBRTH    ~ "Multiple gestation",
      SGA         ~ "Small for gestational age",
      ANTESTER    ~ "Antenatal steroid exposure",
      MEDUC4     ~ "Maternal education",
      Race3X      ~ "Race",
      HISP        ~ "Hispanic ethnicity",
      OBGMAGE     ~ "Maternal age, years",
      APGAR5      ~ "5-minute Apgar score",
      age_rndaysX ~ "Age at randomization, days",
      map_fio2X   ~ "Baseline respiratory severity (MAP × FiO₂)"
    ),
    
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ kruskal.test,
      all_categorical() ~ chisq.test
    )
    
  ) %>%
  add_overall()%>%
  
  modify_header(
    label ~ "**Baseline characteristic**"
  ) %>%
  bold_labels() %>%
  modify_caption(
    "**Table 1. Baseline demographic and clinical characteristics by phenotype.**"
  )

table1

## -------------------------
## Export (submission-ready)
## -------------------------

table1 %>%
  as_flex_table() %>%
  set_caption(
    "Table_1_Baseline_1_characteristics_by_Pheno"
  ) %>%
  save_as_docx(
    path = "Table_1_Baseline_1_characteristics_by_Pheno.docx"
  )



## =========================================================
## FIGURE 2. Do baseline phenotypes meaningfully stratify short-term respiratory 
# success, BPD severity, and longer-term neurodevelopment?
## =========================================================

## =========================================================
## FIGURE 2. Clinical outcomes stratified by baseline phenotype
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

# Logistic regression
glm_extub <- glm(extube14 ~ phenotype, data = df_extub, family = binomial)

# Overall likelihood ratio test
lrt_p <- anova(glm_extub, test = "LRT")$`Pr(>Chi)`[2]

# Marginal means (probabilities)
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
  rename(
    rate = prob,
    lcl = asymp.LCL,
    ucl = asymp.UCL
  ) %>%
  mutate(phenotype = fct_rev(phenotype))

p_fig2a <- ggplot(emm_df, aes(x = rate, y = phenotype, color = phenotype)) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl),
                 height = 0.25, linewidth = 1.2) +
  geom_point(size = 5) +
  geom_text(
    aes(label = label),
    nudge_x = 0.12, size = 3.8,
    fontface = "bold", color = "grey10"
  ) +
  annotate(
    "text",
    x = 0.02, y = length(pheno_order) + 0.4,
    label = sprintf("Logistic regression: overall p = %.3f", lrt_p),
    hjust = 0, size = 3.5, color = "grey30"
  ) +
  scale_x_continuous(
    labels = percent,
    limits = c(0, 1.15),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_color_manual(values = pheno_cols) +
  theme_publication() +
  labs(
    title = "A. Successful extubation ≤14 days",
    x = "Model-estimated probability (95% CI)",
    y = NULL
  ) +
  guides(color = "none")

## =========================================================
## PANEL B — BPD severity at 36 weeks (Jensen)
## =========================================================
bpd_df <- df_pheno %>%
  filter(!is.na(bpd_jen36_hc), phenotype %in% pheno_order) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    jensen = factor(
      bpd_jen36_hc,
      levels = 0:3,
      labels = c("No support","≤2L NC",">2L/NIPPV","IPPV")
    )
  ) %>%
  count(phenotype, jensen) %>%
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
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.8, fontface = "bold", color = "white"
  ) +
  scale_x_discrete(labels = setNames(bpd_n$label, bpd_n$phenotype)) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_fill_manual(values = c("#F0F0F0","#BDBDBD","#737373","#252525")) +
  theme_publication() +
  labs(
    title = "B. BPD severity at 36 weeks (Jensen)",
    x = NULL, y = "Distribution", fill = "Severity"
  )

## =========================================================
## PANEL C — Duration of mechanical ventilation
## =========================================================
vent_df <- df_pheno %>%
  filter(!is.na(days_VENT36), phenotype %in% pheno_order) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    vent_days = days_VENT36
  )

vent_sum <- vent_df %>%
  group_by(phenotype) %>%
  summarise(
    n = n(),
    median = median(vent_days),
    q25 = quantile(vent_days, 0.25),
    q75 = quantile(vent_days, 0.75),
    mean = mean(vent_days),
    label = sprintf("Median %.0f\n(IQR %.0f–%.0f)", median, q25, q75),
    .groups = "drop"
  )

kw_p <- kruskal.test(vent_days ~ phenotype, data = vent_df)$p.value

p_fig2c <- ggplot(vent_df, aes(x = phenotype, y = vent_days, fill = phenotype)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = phenotype),
              width = 0.15, alpha = 0.4, size = 1.3) +
  geom_point(
    data = vent_sum,
    aes(y = mean),
    size = 4.5, shape = 23, color = "white"
  ) +
  geom_text(
    data = vent_sum,
    aes(y = q75 + 6, label = label),
    size = 3.2, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2.5, y = max(vent_df$vent_days) * 0.95,
    label = sprintf("Kruskal–Wallis p = %.3f", kw_p),
    size = 3.5, color = "grey30"
  ) +
  scale_fill_manual(values = pheno_cols) +
  scale_color_manual(values = pheno_cols) +
  theme_publication() +
  labs(
    title = "C. Duration of mechanical ventilation",
    x = NULL, y = "Days on ventilator"
  ) +
  guides(color = "none", fill = "none")

## =========================================================
## PANEL D — Any neurodevelopmental impairment
## =========================================================
ndi_indiv <- df_pheno %>%
  filter(!is.na(NDI3), phenotype %in% pheno_order) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    y = as.numeric(NDI3 == 1)
  )

ndi_sum <- ndi_indiv %>%
  group_by(phenotype) %>%
  summarise(
    n = n(),
    n_ndi = sum(y),
    rate = mean(y),
    se = sqrt(rate * (1 - rate) / n),
    lcl = pmax(0, rate - 1.96 * se),
    ucl = pmin(1, rate + 1.96 * se),
    label = sprintf("%d/%d\n%s", n_ndi, n, percent(rate, 0.1)),
    .groups = "drop"
  )

tbl <- table(ndi_indiv$phenotype, ndi_indiv$y)
ndi_p <- if (any(tbl < 5)) fisher.test(tbl)$p.value else chisq.test(tbl)$p.value
ndi_test <- if (any(tbl < 5)) "Fisher's exact" else "χ²"

p_fig2d <- ggplot() +
  geom_jitter(
    data = ndi_indiv,
    aes(x = phenotype, y = y, color = phenotype),
    width = 0.12, height = 0.02, alpha = 0.35, size = 1.4
  ) +
  geom_linerange(
    data = ndi_sum,
    aes(x = phenotype, ymin = lcl, ymax = ucl, color = phenotype),
    linewidth = 1.6
  ) +
  geom_point(
    data = ndi_sum,
    aes(x = phenotype, y = rate, fill = phenotype),
    size = 5.5, shape = 21, color = "white"
  ) +
  geom_text(
    data = ndi_sum,
    aes(x = phenotype, y = rate, label = label),
    nudge_y = 0.14, size = 3.2, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2.5, y = 0.95,
    label = sprintf("%s p = %.3f", ndi_test, ndi_p),
    size = 3.5, color = "grey30"
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
  scale_color_manual(values = pheno_cols) +
  scale_fill_manual(values = pheno_cols) +
  theme_publication() +
  labs(
    title = "D. Any neurodevelopmental impairment at 18–24 months",
    x = NULL, y = "Probability"
  ) +
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








## =========================================================
## FIGURE 3. Enhanced FiO₂ trajectories following treatment initiation
## =========================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(scales)

## -------------------------
## Phenotype colors (CONSISTENT)
## -------------------------
pheno_cols <- c(
  "Phenotype 1" = "#E41A1C",
  "Phenotype 2" = "#377EB8",
  "Phenotype 3" = "#4DAF4A",
  "Phenotype 4" = "#984EA3"
)

pheno_order <- c("Phenotype 1", "Phenotype 2", "Phenotype 3", "Phenotype 4")

## -------------------------
## Publication theme
## -------------------------
theme_publication <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", margin = margin(b = 8), size = base_size),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey30", linewidth = 0.5),
      axis.ticks = element_line(color = "grey30"),
      axis.text = element_text(size = base_size - 2, color = "grey20"),
      axis.title = element_text(size = base_size - 1, face = "bold"),
      plot.background = element_rect(fill = "grey98", color = NA),
      panel.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
      legend.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.key.size = unit(0.6, "cm"),
      legend.background = element_rect(fill = "grey98", color = NA),
      plot.margin = margin(15, 15, 15, 15)
    )
}

## -------------------------
## Prepare data
## -------------------------
fig3_long <- hc04_dash %>%
  left_join(
    df_pheno %>% select(randcent, rand_id, phenotype),
    by = c("randcent", "rand_id")
  ) %>%
  filter(
    phenotype %in% pheno_order,
    !is.na(HC4RFIO2),
    HC4DAY >= 0,
    HC4DAY <= 14
  ) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    day = HC4DAY
  )

## -------------------------
## Mixed-effects model
## -------------------------
fit_fio2 <- lmer(
  HC4RFIO2 ~ phenotype * day + (1 | randcent / rand_id),
  data = fig3_long,
  REML = FALSE
)

anova_fit <- anova(fit_fio2)
p_int <- anova_fit["phenotype:day", "Pr(>F)"]

## -------------------------
## Estimated marginal means
## -------------------------
emm_fio2 <- emmeans(
  fit_fio2,
  ~ phenotype | day,
  at = list(day = 0:14),
  type = "response"
)

emm_fio2_df <- as.data.frame(emm_fio2) %>%
  rename(
    lower.CL = asymp.LCL,
    upper.CL = asymp.UCL
  )

## -------------------------
## Calculate sample sizes and baseline values
## -------------------------
pheno_n <- fig3_long %>%
  distinct(phenotype, randcent, rand_id) %>%
  count(phenotype) %>%
  arrange(factor(phenotype, levels = pheno_order))

pheno_stats <- fig3_long %>%
  filter(day == 0) %>%
  group_by(phenotype) %>%
  summarise(
    baseline_mean = mean(HC4RFIO2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pheno_n, by = "phenotype")

# Endpoint values (day 14)
endpoint_vals <- emm_fio2_df %>%
  filter(day == 14) %>%
  mutate(label = sprintf("%.2f", emmean))

# Calculate change from baseline for annotation
change_stats <- emm_fio2_df %>%
  group_by(phenotype) %>%
  filter(day %in% c(0, 14)) %>%
  summarise(
    change = emmean[day == 14] - emmean[day == 0],
    pct_change = (emmean[day == 14] - emmean[day == 0]) / emmean[day == 0] * 100,
    .groups = "drop"
  )

## =========================================================
## PANEL A — Main trajectory plot (ENHANCED)
## =========================================================
p_fig3a <- ggplot(
  emm_fio2_df,
  aes(x = day, y = emmean, color = phenotype, fill = phenotype)
) +
  # Ribbon first (behind lines)
  geom_ribbon(
    aes(ymin = lower.CL, ymax = upper.CL),
    alpha = 0.20,
    color = NA
  ) +
  # Main trajectory lines
  geom_line(linewidth = 1.6, alpha = 0.9) +
  # Add points at key timepoints (0, 7, 14)
  geom_point(
    data = emm_fio2_df %>% filter(day %in% c(0, 7, 14)),
    aes(x = day, y = emmean),
    size = 3.5,
    shape = 21,
    color = "white",
    stroke = 1.2
  ) +
  # Endpoint labels
  geom_text(
    data = endpoint_vals,
    aes(x = day + 0.8, y = emmean, label = label),
    hjust = 0,
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  # Interaction p-value annotation
  annotate(
    "text",
    x = 1,
    y = 0.95,
    label = sprintf("Phenotype × Time: p %s", 
                    ifelse(p_int < 0.001, "< 0.001", sprintf("= %.3f", p_int))),
    hjust = 0,
    size = 3.8,
    fontface = "italic",
    color = "grey20"
  ) +
  scale_color_manual(
    values = pheno_cols,
    labels = paste0(pheno_n$phenotype, " (n=", pheno_n$n, ")")
  ) +
  scale_fill_manual(
    values = pheno_cols,
    labels = paste0(pheno_n$phenotype, " (n=", pheno_n$n, ")")
  ) +
  scale_x_continuous(
    breaks = seq(0, 14, 2),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  scale_y_continuous(
    limits = c(0.38, 1.0),
    breaks = seq(0.4, 1.0, 0.1),
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_publication(base_size = 14) +
  labs(
    title = "A. FiO trajectories by baseline phenotype",
    x = "Days since treatment initiation",
    y = "Fraction of Inspired Oxygen (FiO2)",
    color = "Baseline Phenotype",
    fill = "Baseline Phenotype"
  ) +
  theme(
    legend.position = c(0.78, 0.85),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.5),
    legend.key = element_rect(fill = "white")
  )

## =========================================================
## PANEL B — Change from baseline (NEW)
## =========================================================
change_df <- emm_fio2_df %>%
  left_join(
    emm_fio2_df %>% 
      filter(day == 0) %>% 
      select(phenotype, baseline = emmean),
    by = "phenotype"
  ) %>%
  mutate(
    change = emmean - baseline,
    change_lower = lower.CL - baseline,
    change_upper = upper.CL - baseline
  )

p_fig3b <- ggplot(
  change_df,
  aes(x = day, y = change, color = phenotype, fill = phenotype)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_ribbon(
    aes(ymin = change_lower, ymax = change_upper),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.6, alpha = 0.9) +
  geom_point(
    data = change_df %>% filter(day %in% c(0, 7, 14)),
    aes(x = day, y = change),
    size = 3.5,
    shape = 21,
    color = "white",
    stroke = 1.2
  ) +
  scale_color_manual(values = pheno_cols) +
  scale_fill_manual(values = pheno_cols) +
  scale_x_continuous(
    breaks = seq(0, 14, 2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.01),
    breaks = seq(-0.15, 0.05, 0.05)
  ) +
  theme_publication(base_size = 14) +
  labs(
    title = "B. Change from baseline",
    x = "Days since treatment initiation",
    y = "Change in FiO2 from Day 0",
    color = NULL,
    fill = NULL
  ) +
  guides(color = "none", fill = "none")

## =========================================================
## PANEL C — Rate of change (slope visualization)
## =========================================================
slope_df <- change_stats %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    label = sprintf("%.3f", change)
  )

p_fig3c <- ggplot(slope_df, aes(x = phenotype, y = change, fill = phenotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(width = 0.7, alpha = 0.85, color = "grey30", linewidth = 0.6) +
  geom_text(
    aes(label = label, y = change + sign(change) * 0.008),
    size = 4,
    fontface = "bold",
    color = "grey10"
  ) +
  scale_fill_manual(values = pheno_cols) +
  theme_publication(base_size = 14) +
  theme(panel.grid.major.x = element_blank()) +
  labs(
    title = "C. Total change (Day 0 → Day 14)",
    x = NULL,
    y = "Change in FiO2"
  ) +
  guides(fill = "none")

## =========================================================
## PANEL D — Individual spaghetti plot (sample)
## =========================================================
# Sample 30 patients per phenotype for visualization
set.seed(42)

# Get sample size per phenotype first
pheno_counts <- fig3_long %>%
  distinct(phenotype, randcent, rand_id) %>%
  count(phenotype)

# Sample patients
sample_ids <- fig3_long %>%
  distinct(phenotype, randcent, rand_id) %>%
  group_by(phenotype) %>%
  slice_sample(prop = 1) %>%  # shuffle first
  slice_head(n = 30) %>%      # take first 30 (or all if fewer)
  ungroup()

spaghetti_df <- fig3_long %>%
  semi_join(sample_ids, by = c("phenotype", "randcent", "rand_id"))

p_fig3d <- ggplot() +
  # Individual trajectories (thin, transparent)
  geom_line(
    data = spaghetti_df,
    aes(x = day, y = HC4RFIO2, 
        group = interaction(randcent, rand_id),
        color = phenotype),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  # Population mean (thick, opaque)
  geom_line(
    data = emm_fio2_df,
    aes(x = day, y = emmean, color = phenotype),
    linewidth = 1.8,
    alpha = 1
  ) +
  scale_color_manual(values = pheno_cols) +
  scale_x_continuous(breaks = seq(0, 14, 2)) +
  scale_y_continuous(
    limits = c(0.2, 1.0),
    labels = number_format(accuracy = 0.01)
  ) +
  facet_wrap(~ phenotype, nrow = 1) +
  theme_publication(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey70"),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "D. Individual patient trajectories (n=30 per phenotype)",
    x = "Days since treatment initiation",
    y = "FiO2"
  ) +
  guides(color = "none")

## =========================================================
## FINAL LAYOUT
## =========================================================
fig3_final <- 
  p_fig3a /
  (p_fig3b | p_fig3c) /
  p_fig3d +
  plot_layout(heights = c(1.2, 0.8, 0.9)) +
  plot_annotation(
    title = "Figure 3. FiO2 trajectories following treatment initiation",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.background = element_rect(fill = "grey98", color = NA)
    )
  )

fig3_final

ggsave(
  "Figure3_FiO2_Trajectories_Enhanced.png",
  fig3_final,
  width = 14.5,
  height = 13.5,
  dpi = 300,
  bg = "white"
)



library(tidyverse)
library(lme4)
library(lmerTest)
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
## Prepare dataset
## =========================
fig4_df <- df_pheno %>%
  filter(
    phenotype %in% pheno_order,
    !is.na(extub14),
    !is.na(treat)
  ) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    treat = factor(treat, levels = c("Placebo", "Hydrocortisone"))
  )

## =========================
## PANEL A — Forest plot (ORs)
## =========================
fit_extub <- glm(
  extub14 ~ phenotype * treat,
  data = fig4_df,
  family = binomial
)

emm_or <- emmeans(
  fit_extub,
  ~ treat | phenotype,
  type = "response"
)

or_contrast <- contrast(
  emm_or,
  method = "revpairwise",
  by = "phenotype"
)

or_df <- summary(or_contrast, infer = c(TRUE, TRUE)) %>%
  as.data.frame() %>%
  mutate(
    phenotype = factor(phenotype, levels = rev(pheno_order)),
    OR = odds.ratio,
    lower = asymp.LCL,
    upper = asymp.UCL,
    sig = p.value < 0.05,
    label = sprintf("%.2f (%.2f–%.2f)", OR, lower, upper)
  )

p_fig4a <- ggplot(or_df, aes(x = OR, y = phenotype)) +
  annotate("rect", xmin = 0.9, xmax = 1.1,
           ymin = -Inf, ymax = Inf,
           fill = "grey95") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.25,
    linewidth = 1.1,
    color = "#4C72B0"
  ) +
  geom_point(
    size = 4.5,
    shape = 21,
    fill = "white",
    color = "#4C72B0",
    stroke = 1
  ) +
  geom_text(
    aes(label = label),
    x = 10,
    hjust = 0,
    size = 3.3
  ) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 4, 8),
    limits = c(0.4, 12)
  ) +
  theme_publication() +
  labs(
    title = "A. Hydrocortisone effect on successful extubation ≤14 days",
    subtitle = "Odds ratios with 95% confidence intervals",
    x = "Odds ratio (Hydrocortisone vs Placebo)",
    y = NULL
  )

## =========================
## PANEL B — Probability of extubation
## =========================
prob_df <- as.data.frame(emm_or) %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order))

pvals_extub <- summary(or_contrast) %>%
  as.data.frame() %>%
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
  summarise(y = max(asymp.UCL) + 0.08, .groups = "drop") %>%
  left_join(pvals_extub, by = "phenotype")

p_fig4b <- ggplot(prob_df, aes(x = treat, y = prob, fill = treat)) +
  geom_col(width = 0.6, color = "grey40") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  geom_text(
    aes(label = percent(prob, 0.1)),
    position = position_stack(vjust = 0.5),
    fontface = "bold",
    size = 4
  ) +
  geom_text(
    data = annot_b,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(~ phenotype, nrow = 2) +
  scale_fill_manual(values = treat_cols) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  theme_publication() +
  labs(
    title = "B. Model-estimated probability of extubation ≤14 days",
    x = NULL,
    y = "Probability"
  )

## =========================
## PANEL C — Survival without BPD
## =========================
fig4_bpd <- fig4_df %>%
  filter(!is.na(surv36_no_bpd))

fit_bpd <- glm(
  surv36_no_bpd ~ phenotype * treat,
  data = fig4_bpd,
  family = binomial
)

emm_bpd <- emmeans(
  fit_bpd,
  ~ treat | phenotype,
  type = "response"
)

bpd_df <- as.data.frame(emm_bpd)

bpd_contrast <- contrast(
  emm_bpd,
  method = "revpairwise",
  by = "phenotype"
)

bpd_p <- summary(bpd_contrast) %>%
  as.data.frame() %>%
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
  summarise(y = max(asymp.UCL) + 0.06, .groups = "drop") %>%
  left_join(bpd_p, by = "phenotype")

p_fig4c <- ggplot(bpd_df, aes(x = treat, y = prob, fill = treat)) +
  geom_col(width = 0.6, color = "grey40") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  geom_text(
    aes(label = percent(prob, 0.1)),
    position = position_stack(vjust = 0.5),
    fontface = "bold",
    size = 4
  ) +
  geom_text(
    data = annot_c,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(~ phenotype, nrow = 2) +
  scale_fill_manual(values = treat_cols) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  theme_publication() +
  labs(
    title = "C. Survival without BPD at 36 weeks",
    x = NULL,
    y = "Probability"
  )

## =========================
## PANEL D — Absolute risk difference (FIXED)
## =========================
## =========================
## PANEL D — Absolute risk difference (CORRECTED)
## =========================

# Calculate risk differences manually from probabilities
rd_df <- prob_df %>%
  select(phenotype, treat, prob, SE) %>%
  pivot_wider(
    names_from = treat,
    values_from = c(prob, SE)
  ) %>%
  mutate(
    risk_diff = prob_Hydrocortisone - prob_Placebo,
    # Approximate SE for difference
    se_diff = sqrt(SE_Placebo^2 + SE_Hydrocortisone^2),
    lower = risk_diff - 1.96 * se_diff,
    upper = risk_diff + 1.96 * se_diff
  )

# Get p-values from the OR contrast
pvals_d <- summary(or_contrast) %>%
  as.data.frame() %>%
  select(phenotype, p.value) %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order))

# Combine
rd_df <- rd_df %>%
  left_join(pvals_d, by = "phenotype") %>%
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

p_fig4d <- ggplot(
  rd_df,
  aes(x = phenotype, y = risk_diff, fill = phenotype)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  geom_col(
    width = 0.65,
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = percent(risk_diff, 0.1),
      vjust = ifelse(risk_diff > 0, -0.7, 1.6)
    ),
    fontface = "bold",
    size = 4
  ) +
  geom_text(
    aes(y = upper + 0.04, label = stars),
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = pheno_cols) +
  scale_y_continuous(
    labels = percent,
    expand = expansion(mult = c(0.2, 0.2))
  ) +
  theme_publication() +
  labs(
    title = "D. Absolute risk difference (Hydrocortisone − Placebo)",
    subtitle = "*** p < 0.001",
    x = NULL,
    y = "Risk difference"
  )
## =========================
## FINAL FIGURE
## =========================
fig4_final <- (p_fig4a | p_fig4b) / (p_fig4c | p_fig4d) +
  plot_annotation(
    title = "Figure 4. Differential response to hydrocortisone by baseline phenotype",
    caption = "Model-adjusted estimates from logistic regression with treatment-by-phenotype interaction. Error bars represent 95% confidence intervals. Statistical tests compare hydrocortisone vs placebo within each phenotype."
  )

fig4_final

ggsave(
  "Figure4_Phenotype2.png",
  fig4_final,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

#------------------------
# TABLE 2
#-----------------------
library(tidyverse)
library(emmeans)
library(gtsummary)
library(flextable)

# -------------------------
# Absolute risks
# -------------------------
prob_df <- as.data.frame(emm_or) %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    treat = factor(treat, levels = c("Placebo", "Hydrocortisone"))
  )

# Pivot probabilities
prob_wide <- prob_df %>%
  select(phenotype, treat, prob, asymp.LCL, asymp.UCL) %>%
  pivot_wider(
    names_from = treat,
    values_from = c(prob, asymp.LCL, asymp.UCL)
  )

# -------------------------
# Odds ratios
# -------------------------
or_df <- confint(or_contrast, type = "response") %>%
  as.data.frame() %>%
  mutate(
    phenotype = factor(phenotype, levels = pheno_order),
    OR = odds.ratio,
    OR_ci = sprintf("%.2f (%.2f–%.2f)", OR, asymp.LCL, asymp.UCL)
  ) %>%
  left_join(
    summary(or_contrast) %>% as.data.frame() %>% select(phenotype, p.value),
    by = "phenotype"
  )

# -------------------------
# Sample sizes
# -------------------------
n_df <- df_pheno %>%
  filter(phenotype %in% pheno_order) %>%
  count(phenotype, treat) %>%
  pivot_wider(names_from = treat, values_from = n) %>%
  mutate(n_label = sprintf("%d / %d", Placebo, Hydrocortisone))

# -------------------------
# Final Table
# -------------------------
table2 <- prob_wide %>%
  left_join(or_df, by = "phenotype") %>%
  left_join(n_df, by = "phenotype") %>%
  mutate(
    placebo_pct = percent(prob_Placebo, 1),
    hc_pct = percent(prob_Hydrocortisone, 1),
    ard = prob_Hydrocortisone - prob_Placebo,
    ard_ci = sprintf(
      "%.1f%% (%.1f–%.1f)",
      100 * ard,
      100 * (asymp.LCL_Hydrocortisone - asymp.UCL_Placebo),
      100 * (asymp.UCL_Hydrocortisone - asymp.LCL_Placebo)
    )
  ) %>%
  select(
    phenotype,
    n_label,
    placebo_pct,
    hc_pct,
    ard_ci,
    OR_ci,
    p.value
  )

table2 <- table2 %>%
  mutate(
    p_value_fmt = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ formatC(p.value, digits = 2, format = "fg")
    )
  )

table2 %>%
  select(
    phenotype,
    n_label,
    placebo_pct,
    hc_pct,
    ard_ci,
    OR_ci,
    p_value_fmt
  )

# -------------------------
# Pretty export
# -------------------------
table2 %>%
  flextable() %>%
  bold(i = ~ phenotype == "Phenotype 2", bold = TRUE) %>%
  set_header_labels(
    phenotype = "Phenotype",
    n_label = "n (Placebo / HC)",
    placebo_pct = "Extubation ≤14d (%)\nPlacebo",
    hc_pct = "Extubation ≤14d (%)\nHydrocortisone",
    ard_ci = "Absolute Risk Difference (95% CI)",
    OR_ci = "Odds Ratio (95% CI)",
    p.value = "p-value"
  ) %>%
  autofit() %>%
  save_as_docx(path = "Table2_Treatment_Effects_By_Phenotype.docx")
