# install.packages(c("sandwich", "mclust"))

library(sandwich)
library(mclust)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

source("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/scripts/apgar_sensitivity_analysis.R")
source("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/scripts/export_supplement.R")
source("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/scripts/plot_reference_style.R")
source("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/scripts/make_df_pheno.R")

output_dir <- "C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/supplemental"
analysis_dash <- read_csv("C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/JKim_Project/Data_Files_CSV/analysis_dash.csv")

df_pheno <- make_df_pheno(
  analysis_dash,
  expected_n = 800
)

output_dir_pheno_df <- "C:/Users/lrando11/Documents/projects/bpd_hydrocortisone/R_objs/df_pheno.rds"
saveRDS(df_pheno, file = output_dir_pheno_df)

results <- run_apgar_sensitivity(df_pheno,
                                 output_dir = output_dir)

export_apgar_supplement(
  results,
  df_pheno = df_pheno,
  output_dir = output_dir,
  table_number = "S13",   # Change to your actual table number
  figure_number = "S1", # Change if including the figure
  make_figure = TRUE     # FALSE for tables only
)