# Supporting analysis code

These are the supplied analysis workflow and targeted follow-up scripts supporting the revised analyses. Raw infant-level data, model objects, individual predictions and linkage files are deliberately excluded.

## Running the analysis

Use a fresh R session. Set HYDRO_INPUT_DIR to the folder containing analysis_dash.csv, hc04_dash.csv and df_pheno.rds. The default is Data_Files_CSV under the R working directory. Run 01_Hydrocortisone_Reanalysis.R. Each run creates dated figures, tables, models, data and logs folders. Examine module_status.csv, RUN_STATUS.txt, warnings and model diagnostics; completion does not establish scientific validity.

Set HYDRO_RUN_DIR to the resulting dated output directory before running 02_Repair_Ordinal_Tree.R and then 03_BPD_Flexibility.R. The latter creates a new dated BPD_Flexibility subfolder and provides the retained threshold-specific Firth sensitivity analyses. The unsuccessful center-flexible ordinal model and nonretained outputs must not be substituted for these analyses.

The primary and secondary results in the manuscript were reconciled with the supplied final numerical exports and independent numerical checks. This environment has no R executable; these archived scripts have not been rerun end to end here. Presentation-only updates in this archive correct the shared secondary-table footnote and make input/run paths configurable. No model specification or numerical result was changed by those updates.

## Figures and reporting

The submission figures preserve the original artwork and panel structures where valid. They are the files in the Main Figures, Supplement and Graphical Abstract folders; the default plotting sections in this analysis script do not reproduce all final artwork. Do not replace the reviewed submission figures with default regenerated candidate figures. The supplementary document states which analyses were retained and the relevant assumptions and limitations.

The code supports optional verified family linkage and NICHD predictions. Neither was available for the reported analyses. No sibling-adjusted or direct NICHD-estimator comparison is claimed. Any new optional inputs change the analysis and require renewed reconciliation.

