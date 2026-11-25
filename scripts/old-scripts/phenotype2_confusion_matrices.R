# REPLACEMENT FOR SECTION 7 OF YOUR SCRIPT
# Run sections 1-6 first, then run this file in the same R session.
# Reuses existing test-set predictions; does not refit models or resplit data.
# The legacy labels map to phenotype membership, not proven treatment response:
#   Responder     = Phenotype 2 (positive)
#   Non_Responder = Phenotypes 1, 3, and 4 combined (negative)
# These outputs do not establish individual hydrocortisone benefit.

if (!requireNamespace("yardstick", quietly = TRUE)) {
  stop("Install the yardstick package before running this section.")
}

prediction_objects <- c("tree_pred", "logistic_pred", "rf_pred")
if (!all(vapply(prediction_objects, exists, logical(1), inherits = TRUE))) {
  stop("Run sections 1-6 first to create tree_pred, logistic_pred, and rf_pred.")
}

binary_predictions <- list(
  decision_tree = tree_pred,
  logistic_regression = logistic_pred,
  random_forest = rf_pred
)

output_dir <- "phenotype2_classification_outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phenotype_levels <- c("Other phenotypes (1, 3, 4)", "Phenotype 2")

relabel_membership <- function(x) {
  x <- as.character(x)
  if (anyNA(x) || any(!x %in% c("Non_Responder", "Responder"))) {
    stop("Expected nonmissing legacy labels Non_Responder and Responder.")
  }
  factor(x, levels = c("Non_Responder", "Responder"),
         labels = phenotype_levels)
}

confusion_matrices <- list()
binary_metric_rows <- list()

for (model_name in names(binary_predictions)) {
  pred <- binary_predictions[[model_name]]
  if (!all(c("responder_status", ".pred_class") %in% names(pred)) ||
      nrow(pred) == 0L) {
    stop("Missing prediction columns or empty predictions for ", model_name)
  }

  evaluated <- data.frame(
    Actual = relabel_membership(pred$responder_status),
    Predicted = relabel_membership(pred$.pred_class)
  )
  if (any(table(evaluated$Actual) == 0L)) {
    stop("Both actual classes must be present for ", model_name)
  }

  cm <- yardstick::conf_mat(evaluated, truth = Actual, estimate = Predicted)
  confusion_matrices[[model_name]] <- cm
  cat("\n", model_name, ": Phenotype 2 vs other phenotypes\n", sep = "")
  cat("Rows = predicted; columns = actual. Positive class = Phenotype 2.\n")
  print(cm)

  cm_csv <- data.frame(
    Predicted = rownames(cm$table),
    as.data.frame.matrix(cm$table),
    check.names = FALSE, row.names = NULL
  )
  write.csv(cm_csv, file.path(output_dir, paste0(model_name, "_confusion_matrix.csv")),
            row.names = FALSE)

  model_metrics <- rbind(
    yardstick::accuracy(evaluated, truth = Actual, estimate = Predicted),
    yardstick::sens(evaluated, truth = Actual, estimate = Predicted,
                    event_level = "second"),
    yardstick::spec(evaluated, truth = Actual, estimate = Predicted,
                    event_level = "second")
  )
  binary_metric_rows[[model_name]] <- data.frame(
    Model = model_name, as.data.frame(model_metrics), row.names = NULL
  )
}

binary_metrics <- do.call(rbind, binary_metric_rows)
rownames(binary_metrics) <- NULL
cat("\nTEST-SET METRICS (proportions; displayed with 10 significant digits)\n")
print(binary_metrics, digits = 10, row.names = FALSE)
write.csv(binary_metrics, file.path(output_dir, "binary_classification_metrics.csv"),
          row.names = FALSE)

rf_confusion_matrix <- confusion_matrices[["random_forest"]]
rf_binary_metrics <- binary_metrics[binary_metrics$Model == "random_forest", ]

cat("\nRANDOM FOREST COUNTS FOR REVIEWER\n")
print(data.frame(
  TP = unname(rf_confusion_matrix$table[2, 2]),
  FN = unname(rf_confusion_matrix$table[1, 2]),
  FP = unname(rf_confusion_matrix$table[2, 1]),
  TN = unname(rf_confusion_matrix$table[1, 1])
), row.names = FALSE)
cat("\nRANDOM FOREST METRICS\n")
print(rf_binary_metrics, digits = 10, row.names = FALSE)
cat("\nCSV files saved in: ", normalizePath(output_dir), "\n", sep = "")
