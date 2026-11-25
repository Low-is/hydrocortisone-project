# ===============================================================
# 14. SAVE ALL MODELS AND OUTPUTS
# ===============================================================
# saveRDS(tree_fit, "phenotype2_decision_tree.rds")
# saveRDS(logistic_fit, "phenotype2_logistic_model.rds")
# saveRDS(rf_fit, "phenotype2_random_forest.rds")

# cat("\n\n========================================")
# cat("\n✓ ALL OUTPUTS SAVED:")
# cat("\n========================================")
# cat("\n\nCSV FILES:")
# cat("\n  • model_performance_comparison.csv")
# cat("\n  • phenotype2_risk_score.csv")
# cat("\n  • baseline_characteristics_by_responder.csv")
# cat("\n  • pvalues_responder_comparison.csv")
# cat("\n  • comprehensive_summary.csv")
# cat("\n\nMODEL FILES (.rds):")
# cat("\n  • phenotype2_decision_tree.rds")
# cat("\n  • phenotype2_logistic_model.rds")
# cat("\n  • phenotype2_random_forest.rds")
# cat("\n\nIMAGES:")
# cat("\n  • decision_tree_phenotype2.png")
# cat("\n\n========================================")
# cat("\n\nCLINICAL INTERPRETATION:")
# cat("\n========================================")
# cat("\n\nKEY FINDINGS:")
# cat("\n1. APGAR5 is the STRONGEST predictor (↑ Apgar = ↑ Phenotype 2)")
# cat("\n2. Lower GA, lower BW, higher MAP×FiO2 = LESS likely Phenotype 2")
# cat("\n3. Random Forest has best performance (AUC ~0.97, Accuracy ~90%)")
# cat("\n\nDECISION RULE (from tree):")
# cat("\n• High Apgar (≥3.5) + Higher weight (≥860g) + Higher illness (MAP×FiO2 ≥948)")
# cat("\n  → 57% chance of Phenotype 2 → CONSIDER HYDROCORTISONE")
# cat("\n\nRISK SCORE:")
# cat("\n• Calculate: (APGAR5 × 36) - (GESTAGE × 23) - (MAP×FiO2 × 16) - (BW × 13) + (age × 12)")
# cat("\n• Score > 0 → Likely Phenotype 2 responder")
# cat("\n========================================\n")

# ===============================================================
# CLINICAL DECISION SUPPORT: Hydrocortisone Responder Prediction
# Goal: Identify babies likely to benefit from hydrocortisone (Phenotype 2)
# Based on: GESTAGE, BIRTHWTX, APGAR5, age_rndaysX, map_fio2X
# ===============================================================

library(tidymodels)
library(rpart.plot)
library(vip)
library(ggplot2)
library(dplyr)

set.seed(42)

# ===============================================================
# 1. Setup - Define Responder vs Non-Responder
# ===============================================================

df_pheno <- readRDS("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/R_objs/df_pheno.rds")

df <- df_pheno %>%
  filter(!is.na(phenotype)) %>%
  mutate(
    # Phenotype 2 = Responders (benefit from hydrocortisone)
    # Phenotypes 1, 3, 4 = Non-responders
    responder_status = factor(
      ifelse(phenotype == "Phenotype 2", "Responder", "Non_Responder"),
      levels = c("Non_Responder", "Responder")
    )
  ) %>%
  select(
    responder_status,
    GESTAGE,
    BIRTHWTX,
    APGAR5,
    age_rndaysX,
    map_fio2X
  ) %>%
  drop_na()

cat("\n=== COHORT DISTRIBUTION ===\n")
cat("Total babies with complete data:", nrow(df), "\n")
print(table(df$responder_status))
cat("\nPhenotype 2 (Responder) rate:", 
    round(mean(df$responder_status == "Responder") * 100, 1), "%\n")

# Split
split <- initial_split(df, prop = 0.7, strata = responder_status)
train <- training(split)
test <- testing(split)

# ===============================================================
# 2. Exploratory: How do responders differ at baseline?
# ===============================================================
cat("\n=== BASELINE CHARACTERISTICS BY RESPONDER STATUS ===\n")

comparison <- train %>%
  group_by(responder_status) %>%
  summarise(
    n = n(),
    mean_GA = round(mean(GESTAGE), 1),
    mean_BW = round(mean(BIRTHWTX), 0),
    mean_Apgar5 = round(mean(APGAR5), 1),
    mean_age_rand = round(mean(age_rndaysX), 1),
    mean_MAP_FiO2 = round(mean(map_fio2X), 1),
    .groups = "drop"
  )

print(comparison)

# T-tests for continuous variables
cat("\n=== STATISTICAL TESTS ===\n")
vars_to_test <- c("GESTAGE", "BIRTHWTX", "APGAR5", "age_rndaysX", "map_fio2X")

for(var in vars_to_test) {
  test_result <- t.test(
    train[[var]][train$responder_status == "Responder"],
    train[[var]][train$responder_status == "Non_Responder"]
  )
  cat(sprintf("%s: p = %.4f\n", var, test_result$p.value))
}

# ===============================================================
# 3. MODEL 1: Simple Decision Tree (Bedside Tool)
# ===============================================================
cat("\n\n=== MODEL 1: CLINICAL DECISION TREE ===\n")

tree_spec <- decision_tree(
  tree_depth = 3,
  min_n = 20
) %>%
  set_engine("rpart") %>%
  set_mode("classification")

rec_tree <- recipe(responder_status ~ ., data = train) %>%
  step_impute_median(all_numeric_predictors())

tree_wf <- workflow() %>%
  add_recipe(rec_tree) %>%
  add_model(tree_spec)

tree_fit <- fit(tree_wf, train)

# Visualize
cat("\nDecision tree for identifying Phenotype 2 responders:\n")
rpart.plot(
  extract_fit_engine(tree_fit),
  main = "Hydrocortisone Responder Decision Tree",
  extra = 101,
  roundint = FALSE,
  fallen.leaves = TRUE,
  cex = 0.8
)

# ===============================================================
# 4. MODEL 2: Logistic Regression (Risk Score)
# ===============================================================
cat("\n\n=== MODEL 2: LOGISTIC REGRESSION RISK SCORE ===\n")

logistic_spec <- logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification")

rec_logistic <- recipe(responder_status ~ ., data = train) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors())

logistic_wf <- workflow() %>%
  add_recipe(rec_logistic) %>%
  add_model(logistic_spec)

logistic_fit <- fit(logistic_wf, train)

# Get coefficients
coefs <- extract_fit_parsnip(logistic_fit) %>%
  tidy() %>%
  arrange(desc(abs(estimate)))

cat("\nLogistic Regression Coefficients:\n")
print(coefs)

# ===============================================================
# 5. MODEL 3: Random Forest (Best Performance)
# ===============================================================
cat("\n\n=== MODEL 3: RANDOM FOREST ===\n")

rf_spec <- rand_forest(trees = 500, min_n = 10) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

rec_rf <- recipe(responder_status ~ ., data = train) %>%
  step_impute_median(all_numeric_predictors())

rf_wf <- workflow() %>%
  add_recipe(rec_rf) %>%
  add_model(rf_spec)

rf_fit <- fit(rf_wf, train)

# Variable importance
cat("\nVariable Importance:\n")
vip_plot <- vip(extract_fit_parsnip(rf_fit), num_features = 5) +
  ggtitle("Top Predictors of Hydrocortisone Response") +
  theme_minimal()
print(vip_plot)

# ===============================================================
# 6. PERFORMANCE COMPARISON
# ===============================================================
cat("\n\n=== MODEL PERFORMANCE ON TEST SET ===\n")

# Predictions from all models
tree_pred <- predict(tree_fit, test, type = "prob") %>%
  bind_cols(predict(tree_fit, test)) %>%
  bind_cols(test %>% select(responder_status))

logistic_pred <- predict(logistic_fit, test, type = "prob") %>%
  bind_cols(predict(logistic_fit, test)) %>%
  bind_cols(test %>% select(responder_status))

rf_pred <- predict(rf_fit, test, type = "prob") %>%
  bind_cols(predict(rf_fit, test)) %>%
  bind_cols(test %>% select(responder_status))

# Compare (FIXED AUC calculation)
models <- list(
  "Decision Tree" = tree_pred,
  "Logistic Regression" = logistic_pred,
  "Random Forest" = rf_pred
)

results <- map_dfr(names(models), function(name) {
  pred <- models[[name]]
  
  # Calculate metrics properly
  cm <- conf_mat(pred, truth = responder_status, estimate = .pred_class)
  
  tibble(
    Model = name,
    AUC = roc_auc(pred, truth = responder_status, .pred_Responder, 
                  event_level = "second")$.estimate,
    Accuracy = accuracy(pred, truth = responder_status, estimate = .pred_class)$.estimate,
    Sensitivity = sensitivity(pred, truth = responder_status, estimate = .pred_class,
                              event_level = "second")$.estimate,
    Specificity = specificity(pred, truth = responder_status, estimate = .pred_class,
                              event_level = "second")$.estimate,
    PPV = precision(pred, truth = responder_status, estimate = .pred_class,
                    event_level = "second")$.estimate,
    NPV = npv(pred, truth = responder_status, estimate = .pred_class,
              event_level = "second")$.estimate
  )
})

cat("\n=== MODEL PERFORMANCE (CORRECTED) ===\n")
print(results)

# Save as formatted table
results_formatted <- results %>%
  mutate(across(where(is.numeric), ~round(., 3)))

write.csv(results_formatted, "model_performance_comparison.csv", row.names = FALSE)
cat("\n✓ Performance table saved: model_performance_comparison.csv\n")

# ===============================================================
# 7. CONFUSION MATRICES
# ===============================================================
# cat("\n\n=== CONFUSION MATRICES ===\n")
# 
# cat("\nDecision Tree:\n")
# print(conf_mat(tree_pred, truth = responder_status, estimate = .pred_class))
# 
# cat("\nLogistic Regression:\n")
# print(conf_mat(logistic_pred, truth = responder_status, estimate = .pred_class))
# 
# cat("\nRandom Forest:\n")
# print(conf_mat(rf_pred, truth = responder_status, estimate = .pred_class))

# Supports either the original or updated outcome labels.
truth_col <- if ("phenotype2_status" %in% names(rf_pred)) {
  "phenotype2_status"
} else {
  "responder_status"
}

class_labels <- c(
  "Other_Phenotypes" = "Others (Phenotypes 1, 3, 4)",
  "Non_Responder"    = "Others (Phenotypes 1, 3, 4)",
  "Phenotype_2"      = "Phenotype 2",
  "Responder"       = "Phenotype 2"
)

class_order <- c("Others (Phenotypes 1, 3, 4)", "Phenotype 2")

rf_eval <- data.frame(
  Actual = factor(
    unname(class_labels[as.character(rf_pred[[truth_col]])]),
    levels = class_order
  ),
  Predicted = factor(
    unname(class_labels[as.character(rf_pred$.pred_class)]),
    levels = class_order
  )
)

if (nrow(rf_eval) == 0L || anyNA(rf_eval)) {
  stop("Missing predictions or unexpected class labels in rf_pred.")
}

# Confusion matrix: rows = predicted; columns = actual.
rf_cm <- yardstick::conf_mat(
  rf_eval, truth = Actual, estimate = Predicted
)

cat("\nRANDOM FOREST BINARY CONFUSION MATRIX\n")
cat("Rows = predicted; columns = actual\n")
print(rf_cm)

# One combined performance table.
# Each row treats the named class as positive.
# Overall accuracy applies to the entire model and repeats in both rows.
rf_metrics_table <- dplyr::bind_rows(
  lapply(c(2L, 1L), function(i) {
    j <- if (i == 2L) 1L else 2L
    event <- if (i == 2L) "second" else "first"
    
    data.frame(
      Class = class_order[i],
      Actual_N = sum(rf_eval$Actual == class_order[i]),
      TP = unname(rf_cm$table[i, i]),
      FN = unname(rf_cm$table[j, i]),
      FP = unname(rf_cm$table[i, j]),
      TN = unname(rf_cm$table[j, j]),
      Sensitivity = yardstick::sens(
        rf_eval, Actual, Predicted, event_level = event
      )$.estimate,
      Specificity = yardstick::spec(
        rf_eval, Actual, Predicted, event_level = event
      )$.estimate,
      PPV = yardstick::ppv(
        rf_eval, Actual, Predicted, event_level = event
      )$.estimate,
      NPV = yardstick::npv(
        rf_eval, Actual, Predicted, event_level = event
      )$.estimate,
      F1 = yardstick::f_meas(
        rf_eval, Actual, Predicted, event_level = event
      )$.estimate,
      Overall_accuracy = yardstick::accuracy(
        rf_eval, Actual, Predicted
      )$.estimate
    )
  })
)

# Format metrics as percentages for display.
rf_metrics_display <- rf_metrics_table %>%
  dplyr::mutate(
    dplyr::across(
      c(Sensitivity, Specificity, PPV, NPV, F1, Overall_accuracy),
      ~ ifelse(is.na(.x), "Not estimable", sprintf("%.1f%%", 100 * .x))
    )
  )

cat("\nRANDOM FOREST BINARY CLASSIFICATION METRICS\n")
print(as.data.frame(rf_metrics_display), row.names = FALSE)
cat("\nOverall accuracy is a model-wide metric, repeated for reference.\n")

# Export the combined table, including counts and performance metrics.
save_dir <- if (exists("output_dir")) output_dir else "."
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  rf_metrics_display,
  file.path(save_dir, "random_forest_binary_metrics_table.csv"),
  row.names = FALSE
)

# Preserve unrounded estimates for verification.
write.csv(
  rf_metrics_table,
  file.path(save_dir, "random_forest_binary_metrics_unrounded.csv"),
  row.names = FALSE
)

# ===============================================================
# 8. CLINICAL RISK SCORE CALCULATOR (SAVE TO TABLE)
# ===============================================================
cat("\n\n=== CLINICAL RISK SCORE FOR PHENOTYPE 2 (RESPONDERS) ===\n")

# Create simplified risk score from logistic regression
risk_score <- coefs %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # Scale to 0-100 point system
    points = round((estimate / sum(abs(estimate))) * 100, 0),
    direction = ifelse(estimate > 0, "↑ increases", "↓ decreases")
  ) %>%
  select(Feature = term, Coefficient = estimate, Points = points, 
         Direction = direction, `P-value` = p.value) %>%
  arrange(desc(abs(Coefficient)))

print(risk_score)

# Save risk score table
write.csv(risk_score, "phenotype2_risk_score.csv", row.names = FALSE)
cat("\n✓ Risk score table saved: phenotype2_risk_score.csv\n")

cat("\n=== INTERPRETATION ===")
cat("\nPositive points = Higher likelihood of being Phenotype 2 (Responder)")
cat("\nNegative points = Lower likelihood of being Phenotype 2")
cat("\n\nCLINICAL USE:")
cat("\n1. Calculate total score for each patient")
cat("\n2. Score > 0: Consider Phenotype 2 → Likely to benefit from hydrocortisone")
cat("\n3. Score < 0: Consider other phenotypes → Less likely to benefit")

# ===============================================================
# 8B. BASELINE CHARACTERISTICS TABLE
# ===============================================================
cat("\n\n=== BASELINE CHARACTERISTICS BY RESPONDER STATUS ===\n")

baseline_table <- train %>%
  group_by(responder_status) %>%
  summarise(
    N = n(),
    `GA (weeks)` = sprintf("%.1f ± %.1f", mean(GESTAGE), sd(GESTAGE)),
    `BW (grams)` = sprintf("%.0f ± %.0f", mean(BIRTHWTX), sd(BIRTHWTX)),
    `Apgar 5min` = sprintf("%.1f ± %.1f", mean(APGAR5), sd(APGAR5)),
    `Age at Rx (days)` = sprintf("%.1f ± %.1f", mean(age_rndaysX), sd(age_rndaysX)),
    `MAP×FiO2` = sprintf("%.1f ± %.1f", mean(map_fio2X), sd(map_fio2X)),
    .groups = "drop"
  )

print(baseline_table)

write.csv(baseline_table, "baseline_characteristics_by_responder.csv", row.names = FALSE)
cat("\n✓ Baseline table saved: baseline_characteristics_by_responder.csv\n")

# Add p-values
pvalue_table <- tibble(
  Variable = c("Gestational age", "Birth weight", "Apgar 5min", 
               "Age at randomization", "MAP×FiO2"),
  `P-value` = sapply(vars_to_test, function(var) {
    t.test(
      train[[var]][train$responder_status == "Responder"],
      train[[var]][train$responder_status == "Non_Responder"]
    )$p.value
  })
) %>%
  mutate(`P-value` = sprintf("%.4f", `P-value`))

print(pvalue_table)
write.csv(pvalue_table, "pvalues_responder_comparison.csv", row.names = FALSE)
cat("\n✓ P-value table saved: pvalues_responder_comparison.csv\n")

# ===============================================================
# 9. CALIBRATION PLOT (Are predicted probabilities accurate?)
# ===============================================================
cat("\n\n=== CALIBRATION ASSESSMENT ===\n")

# Use logistic model for calibration
calib_df <- logistic_pred %>%
  mutate(
    pred_bin = cut(.pred_Responder, 
                   breaks = seq(0, 1, 0.1),
                   include.lowest = TRUE)
  ) %>%
  group_by(pred_bin) %>%
  summarise(
    n = n(),
    predicted = mean(.pred_Responder),
    observed = mean(responder_status == "Responder"),
    .groups = "drop"
  ) %>%
  filter(n >= 5)  # Only show bins with ≥5 patients

p_calib <- ggplot(calib_df, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Calibration Plot - Logistic Model",
    subtitle = "Perfect calibration = points on diagonal line",
    x = "Predicted Probability of Phenotype 2",
    y = "Observed Proportion of Phenotype 2",
    size = "N patients"
  ) +
  theme_minimal()

print(p_calib)

# ===============================================================
# 10. PRACTICAL EXAMPLE
# ===============================================================
cat("\n\n=== EXAMPLE PATIENT CALCULATION ===\n")

example_patient <- data.frame(
  GESTAGE = 25.0,
  BIRTHWTX = 750,
  APGAR5 = 6,
  age_rndaysX = 3,
  map_fio2X = 45
)

cat("\nExample Patient:\n")
cat("  Gestational age:", example_patient$GESTAGE, "weeks\n")
cat("  Birth weight:", example_patient$BIRTHWTX, "g\n")
cat("  5-min Apgar:", example_patient$APGAR5, "\n")
cat("  Age at randomization:", example_patient$age_rndaysX, "days\n")
cat("  MAP × FiO2:", example_patient$map_fio2X, "\n\n")

# Predict using all models
tree_prob <- predict(tree_fit, example_patient, type = "prob")$.pred_Responder
logistic_prob <- predict(logistic_fit, example_patient, type = "prob")$.pred_Responder
rf_prob <- predict(rf_fit, example_patient, type = "prob")$.pred_Responder

cat("PREDICTIONS:\n")
cat("  Decision Tree:       ", round(tree_prob * 100, 1), "% probability of Phenotype 2\n")
cat("  Logistic Regression: ", round(logistic_prob * 100, 1), "% probability of Phenotype 2\n")
cat("  Random Forest:       ", round(rf_prob * 100, 1), "% probability of Phenotype 2\n")
cat("\n  Average prediction:  ", round(mean(c(tree_prob, logistic_prob, rf_prob)) * 100, 1), "%\n")

if(mean(c(tree_prob, logistic_prob, rf_prob)) > 0.5) {
  cat("\n✓ RECOMMENDATION: Likely Phenotype 2 - Consider HYDROCORTISONE\n")
} else {
  cat("\n✗ RECOMMENDATION: Unlikely Phenotype 2 - Alternative treatment may be preferred\n")
}

# ===============================================================
# 12. CREATE COMPREHENSIVE SUMMARY TABLE
# ===============================================================
cat("\n\n=== COMPREHENSIVE RESULTS SUMMARY ===\n")

# Key findings summary
summary_table <- tibble(
  Metric = c(
    "Total patients",
    "Phenotype 2 (Responders)",
    "Other phenotypes",
    "Responder rate",
    "",
    "Best model (AUC)",
    "Best model accuracy",
    "Best model sensitivity",
    "Best model specificity",
    "",
    "Most important feature",
    "2nd most important",
    "3rd most important"
  ),
  Value = c(
    as.character(nrow(df)),
    as.character(sum(df$responder_status == "Responder")),
    as.character(sum(df$responder_status == "Non_Responder")),
    sprintf("%.1f%%", mean(df$responder_status == "Responder") * 100),
    "",
    results$Model[which.max(results$AUC)],
    sprintf("%.1f%%", max(results$Accuracy) * 100),
    sprintf("%.1f%%", results$Sensitivity[which.max(results$AUC)] * 100),
    sprintf("%.1f%%", results$Specificity[which.max(results$AUC)] * 100),
    "",
    risk_score$Feature[1],
    risk_score$Feature[2],
    risk_score$Feature[3]
  )
)

print(summary_table)
write.csv(summary_table, "comprehensive_summary.csv", row.names = FALSE)
cat("\n✓ Summary table saved: comprehensive_summary.csv\n")

# ===============================================================
# 13. DECISION TREE RULES (TEXT FORMAT)
# ===============================================================
cat("\n\n=== DECISION TREE RULES (FOR CLINICAL USE) ===\n")

tree_rules <- extract_fit_engine(tree_fit) %>%
  rpart.plot::rpart.rules(., style = "tall", roundint = FALSE) # changed rpart.rules:: -> rpart.plot

cat("\nSimplified decision rules:\n")
cat("1. If APGAR5 < 3.5 → Non-responder (20% Phenotype 2)\n")
cat("2. If APGAR5 ≥ 3.5:\n")
cat("   - If BIRTHWTX < 860g → Non-responder (17% Phenotype 2)\n")
cat("   - If BIRTHWTX ≥ 860g:\n")
cat("     • If MAP×FiO2 < 948 → Non-responder (7% Phenotype 2)\n")
cat("     • If MAP×FiO2 ≥ 948 → RESPONDER (57% Phenotype 2)\n")

# Save tree as image
png("decision_tree_phenotype2.png", width = 10, height = 8, units = "in", res = 300)
rpart.plot(
  extract_fit_engine(tree_fit),
  main = "Hydrocortisone Responder Decision Tree",
  extra = 101,
  roundint = FALSE,
  fallen.leaves = TRUE,
  cex = 0.8
)
dev.off()
cat("\n✓ Decision tree saved: decision_tree_phenotype2.png\n")

saveRDS(training_data,'training_data.rds')
saveRDS(testing_data,'testing_data.rds')
