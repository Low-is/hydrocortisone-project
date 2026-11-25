#####
#libraries-----
#####
# Targeted repair for the completed run uploaded on November 24, 2026.
# Does NOT rerun clustering, primary treatment models, classifier fitting,
# bootstrap analyses, or the other completed modules.
packages <- c("dplyr", "readr", "tibble", "ordinal", "officer", "flextable", "rpart.plot")
for (p in packages) if (!requireNamespace(p, quietly=TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble)})

#####
#existing output folder-----
#####
OUT <- Sys.getenv("HYDRO_RUN_DIR", unset = "SET_COMPLETED_RUN_DIRECTORY_HERE")
DIR <- as.list(setNames(file.path(OUT,c("figures","tables","models","data","logs")),
                        c("figures","tables","models","data","logs")))
if (!file.exists(file.path(DIR$data,"analytic_760.rds"))) {
  stop("Set OUT to the completed run directory containing data/analytic_760.rds")
}
PHENO <- paste("Phenotype",1:4)
Z <- qnorm(.975)
backup <- file.path(OUT,paste0("repair_backup_",format(Sys.time(),"%Y%m%d_%H%M%S")))
dir.create(backup,recursive=TRUE,showWarnings=FALSE)
protect <- c(file.path(DIR$logs,c("module_status.csv","warnings.txt")),
             file.path(OUT,c("RUN_STATUS.txt","output_manifest.csv")),
             file.path(DIR$models,c("analysis_bundle.rds","ordinal_BPD.rds")),
             file.path(DIR$figures,"supplement","sFigure_Decision_Tree.pdf"))
for (f in protect[file.exists(protect)]) {
  destination <- file.path(backup,substring(f,nchar(OUT)+2))
  dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
  if (!file.copy(f,destination,overwrite=FALSE)) stop("Could not back up: ",f)
}
primary_table <- file.path(DIR$tables,"Table2_Unrounded.csv")
primary_md5 <- unname(tools::md5sum(primary_table))
status <- read_csv(file.path(DIR$logs,"module_status.csv"),show_col_types=FALSE)
repair_warnings <- character()
assert <- function(ok,msg) if (!isTRUE(ok)) stop(msg,call.=FALSE)
update_status <- function(module,state,detail="") {
  status <<- status %>% filter(.data$module != !!module) %>%
    bind_rows(tibble(module=module,status=state,detail=detail))
}
export_table <- function(d,name,title,note="") {
  write_csv(as.data.frame(d),file.path(DIR$tables,paste0(name,".csv")),na="")
  ft <- flextable::flextable(as.data.frame(d))
  ft <- flextable::theme_vanilla(ft)
  ft <- flextable::fontsize(ft,size=10,part="all")
  ft <- flextable::font(ft,fontname="Arial",part="all")
  ft <- flextable::bg(ft,bg="#E8EDF1",part="header")
  ft <- flextable::padding(ft,padding=5,part="all")
  ft <- flextable::autofit(ft)
  ft <- flextable::set_table_properties(ft,layout="autofit",width=1)
  ft <- flextable::set_caption(ft,caption=title)
  if (nzchar(note)) ft <- flextable::add_footer_lines(ft,values=note)
  doc <- officer::read_docx()
  doc <- flextable::body_add_flextable(doc,ft)
  print(doc,target=file.path(DIR$tables,paste0(name,".docx")))
}

#####
#ordinal BPD repair-----
#####
ordinal_result <- tryCatch(withCallingHandlers({
  d <- readRDS(file.path(DIR$data,"analytic_760.rds"))
  od <- d %>% filter(death36==0,!is.na(jensen_grade)) %>% droplevels()
  assert(nrow(od)==704,"Unexpected ordinal cohort size; investigate before proceeding")
  fit <- ordinal::clm(jensen_grade ~ phenotype * hc + gestage + male + center,
                      data=od,link="logit",model=TRUE)
  # Current ordinal versions store convergence as a LIST with $code.
  # Comparing the entire list to zero is invalid and can create coercion NAs.
  code <- if (is.list(fit$convergence)) fit$convergence$code else fit$convergence
  H <- fit$Hessian
  ev <- if (all(is.finite(H))) eigen(H,symmetric=TRUE,only.values=TRUE)$values else NA_real_
  diagnostics <- list(N=nrow(od),convergence=fit$convergence,max_gradient=fit$maxGradient,
                       Hessian_eigenvalues=ev,condition_number=fit$cond.H)
  saveRDS(list(fit=fit,diagnostics=diagnostics),file.path(DIR$models,"ordinal_BPD_raw_diagnostic.rds"))
  writeLines(capture.output(str(diagnostics)),file.path(DIR$logs,"ordinal_convergence_diagnostics.txt"))
  assert(is.numeric(code)&&length(code)>0&&all(is.finite(code))&&all(code==0),
         "Ordinal convergence code is not zero; inspect ordinal_convergence_diagnostics.txt")
  assert(is.finite(fit$maxGradient)&&fit$maxGradient<1e-4,"Ordinal gradient check failed")
  assert(all(is.finite(ev))&&all(ev>0),"Ordinal Hessian is not positive definite")
  b <- coef(fit); V <- vcov(fit)
  assert(all(is.finite(b))&&all(is.finite(V)),"Nonfinite ordinal coefficients/covariance")
  ii <- grep(":hc$|^hc:",names(b))
  assert(length(ii)==3,"Expected three phenotype-by-treatment interaction terms")
  B <- b[ii]; S <- V[ii,ii,drop=FALSE]
  assert(qr(S)$rank==3,"Singular interaction covariance")
  W <- as.numeric(t(B)%*%solve(S,B))
  interaction_test <- tibble(statistic=W,df=3L,p_value=pchisq(W,3,lower.tail=FALSE))
  estimates <- bind_rows(lapply(PHENO,function(p) {
    g <- setNames(rep(0,length(b)),names(b)); g["hc"] <- 1
    if (p!=PHENO[1]) {
      nm <- paste0("phenotype",p,":hc")
      assert(nm %in% names(g),paste("Missing contrast:",nm));g[nm] <- 1
    }
    z <- sum(g*b); se <- sqrt(as.numeric(t(g)%*%V%*%g))
    tibble(phenotype=p,OR_higher_grade=exp(z),lower=exp(z-Z*se),upper=exp(z+Z*se),
           p_value=2*pnorm(abs(z/se),lower.tail=FALSE))
  }))
  export_table(estimates,"sTable_Ordinal_BPD","Jensen severity among survivors",
    paste("Proportional-odds cumulative logit, adjusted for gestational age, sex, and center.",
      "Model-based covariance; sibling dependence not modeled. OR >1 indicates higher severity.",
      "Full phenotype-specific treatment contrasts, not interaction coefficients alone.",
      "Survivor-only exploratory analysis. Joint log-odds interaction P =",
      format.pval(interaction_test$p_value,digits=4)))
  write_csv(interaction_test,file.path(DIR$tables,"Ordinal_BPD_Interaction.csv"))
  po_test <- tryCatch(ordinal::nominal_test(fit),error=function(e) {
    update_status("Ordinal proportional odds diagnostic","FAILED",conditionMessage(e));NULL
  })
  if (!is.null(po_test)) {
    pt <- as.data.frame(po_test);pt$term <- rownames(pt)
    write_csv(pt,file.path(DIR$tables,"sTable_Ordinal_Assumption_Diagnostic.csv"))
    update_status("Ordinal proportional odds diagnostic","COMPLETE",
                   "Review term-specific diagnostic results; successful execution does not prove proportional odds")
  }
  result <- list(fit=fit,data=od,estimates=estimates,interaction=interaction_test,
                  po_test=po_test,diagnostics=diagnostics)
  saveRDS(result,file.path(DIR$models,"ordinal_BPD.rds"))
  update_status("Ordinal BPD severity","COMPLETE","Convergence list inspected through $code; model, gradient, and Hessian checks passed")
  print(estimates)
  print(interaction_test)
  result
},warning=function(w) {
  repair_warnings <<- c(repair_warnings,paste0("Ordinal BPD severity: ",conditionMessage(w)))
}),error=function(e) {
  update_status("Ordinal BPD severity","FAILED",conditionMessage(e))
  message("Ordinal repair not complete: ",conditionMessage(e));NULL
})

#####
#decision tree plot repair without refitting models-----
#####
tree_ok <- tryCatch(withCallingHandlers({
  fits <- readRDS(file.path(DIR$models,"classifier_models_training_only.rds"))
  tree <- fits$Decision_tree
  destination <- file.path(DIR$figures,"supplement","sFigure_Decision_Tree.pdf")
  tmp <- file.path(dirname(destination),"sFigure_Decision_Tree_repair_tmp.pdf")
  grDevices::pdf(tmp,width=10,height=7,useDingbats=FALSE)
  tryCatch(rpart.plot::rpart.plot(tree,extra=104,fallen.leaves=TRUE,roundint=FALSE,
    main="Internal classifier of Phenotype 2 membership"),finally=grDevices::dev.off())
  assert(file.exists(tmp)&&file.info(tmp)$size>1000,"Decision-tree PDF was not generated")
  assert(file.copy(tmp,destination,overwrite=TRUE),"Could not save repaired tree PDF")
  unlink(tmp)
  update_status("Decision tree figure repair","COMPLETE","Plot regenerated without refitting any classifier")
  TRUE
},warning=function(w) {
  repair_warnings <<- c(repair_warnings,paste0("Decision tree figure repair: ",conditionMessage(w)))
}),error=function(e) {
  update_status("Decision tree figure repair","FAILED",conditionMessage(e))
  message("Tree figure repair not complete: ",conditionMessage(e));FALSE
})

#####
#update logs and bundle while preserving the primary analysis-----
#####
write_csv(status,file.path(DIR$logs,"module_status.csv"))
original_warnings <- if(file.exists(file.path(DIR$logs,"warnings.txt")))
  readLines(file.path(DIR$logs,"warnings.txt"),warn=FALSE) else character()
# Archive exists above. Remove only the specifically resolved warning entries
# from the active warning log, retaining all other warnings and new diagnostics.
keep <- rep(TRUE,length(original_warnings));skip_tree <- FALSE
prefixes <- paste0(status$module,": ")
for (i in seq_along(original_warnings)) {
  line <- original_warnings[i]
  is_header <- any(startsWith(line,prefixes))
  if (is_header) skip_tree <- FALSE
  if (!is.null(ordinal_result)&&startsWith(line,"Ordinal BPD severity: NAs introduced by coercion")) keep[i] <- FALSE
  if (tree_ok&&startsWith(line,"Phenotype 2 internal classification: Cannot retrieve the data used to build the model")) skip_tree <- TRUE
  if (skip_tree) keep[i] <- FALSE
}
writeLines(unique(c(original_warnings[keep],repair_warnings)),file.path(DIR$logs,"warnings.txt"))
writeLines(capture.output(sessionInfo()),file.path(DIR$logs,"repair_sessionInfo.txt"))
bundle_path <- file.path(DIR$models,"analysis_bundle.rds")
if (file.exists(bundle_path)) {
  bundle <- readRDS(bundle_path)
  if (!is.null(ordinal_result)) bundle$ordinal <- ordinal_result
  bundle$module_status <- status
  saveRDS(bundle,bundle_path)
}
assert(identical(primary_md5,unname(tools::md5sum(primary_table))),"Primary Table 2 changed unexpectedly")
failed <- sum(status$status=="FAILED")
not_run <- sum(grepl("NOT RUN",status$status))
writeLines(c(
  if(failed>0) "COMPLETED WITH FAILED MODULES - REVIEW REQUIRED" else "CORE ANALYSIS COMPLETED - REVIEW OUTPUTS BEFORE SUBMISSION",
  paste("Failed modules:",failed),paste("Modules requiring unavailable inputs:",not_run),
  "Primary Table 2 was not changed by this repair.",
  "Sibling dependence remains unmodeled without genuine family linkage.",
  "Review ordinal convergence diagnostics, proportional-odds diagnostics, warnings, and final figures.",
  paste("Backup folder:",backup),paste("Output folder:",OUT)),file.path(OUT,"RUN_STATUS.txt"))
paths <- list.files(OUT,recursive=TRUE,full.names=TRUE)
paths <- paths[basename(paths)!="output_manifest.csv"]
write_csv(tibble(file=substring(paths,nchar(OUT)+2),bytes=file.info(paths)$size,
                 md5=unname(tools::md5sum(paths))),file.path(OUT,"output_manifest.csv"))
message("Repair finished. Failed modules: ",failed,"; unavailable-input modules: ",not_run)
message("Read logs/ordinal_convergence_diagnostics.txt and tables/sTable_Ordinal_Assumption_Diagnostic.csv")
message("Send the updated module_status.csv, warnings.txt, and ordinal results for review.")
