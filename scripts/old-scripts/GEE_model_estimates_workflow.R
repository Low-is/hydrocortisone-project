library(dplyr)
library(sandwich)
library(lmtest)
library(tidyverse)
library(gtsummary)
library(flextable)
library(sandwich)
library(lmtest)
library(ggplot2)
library(grid)
library(tibble)
library(geepack)
library(emmeans)
library(patchwork)
library(scales)
library(scales)
library(broom)


project_dir <- "C:/Users/lrando11/Documents/projects/bpd_hydrocortisone"
scripts_dir <- file.path(project_dir, "scripts")
objects_dir <- file.path(project_dir, "R_objs")
output_dir <- file.path(project_dir, "table2_output")


source("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/scripts/export_inputs.R")
source("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/scripts/table2_all_infants.R")

df_pheno <- readRDS("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/R_objs/df_pheno.rds")

table2_inputs <- export_inputs(
  df_pheno,
  file = file.path(objects_dir, "table2_inputs.rds")
)


table2_results <- run_table2(
  table2_inputs,
  output_dir = output_dir,
  install_packages = FALSE
)
GEE
print(table2_results$table)

