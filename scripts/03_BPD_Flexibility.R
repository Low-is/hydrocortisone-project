#####
#libraries-----
#####
# Targeted exploratory BPD analyses following proportional-odds diagnostics.
# Uses the existing 704-infant survivor cohort. No clustering or primary
# extubation model is refitted. No existing output files are overwritten.
pkgs <- c("dplyr","readr","tibble","ordinal","logistf","ggplot2","officer","flextable")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) stop("Missing installed package: ",p)
suppressPackageStartupMessages({library(dplyr);library(readr);library(tibble);library(ggplot2)})
options(contrasts=c("contr.treatment","contr.poly"))

#####
#existing run and new output folders-----
#####
RUN_DIR <- Sys.getenv("HYDRO_RUN_DIR", unset = "SET_COMPLETED_RUN_DIRECTORY_HERE")
data_file <- file.path(RUN_DIR,"data","analytic_760.rds")
if (!file.exists(data_file)) stop("Set RUN_DIR to the completed run containing data/analytic_760.rds")
OUT <- file.path(RUN_DIR,paste0("BPD_Flexibility_",format(Sys.time(),"%Y%m%d_%H%M%S")))
DIR <- as.list(setNames(file.path(OUT,c("tables","models","figures","logs")),
                        c("tables","models","figures","logs")))
invisible(lapply(DIR,dir.create,recursive=TRUE,showWarnings=FALSE))
PHENO <- paste("Phenotype",1:4)
Z <- qnorm(.975)
assert <- function(ok,msg) if (!isTRUE(ok)) stop(msg,call.=FALSE)
fmt_p <- function(p) ifelse(is.na(p),"Not estimable",ifelse(p<.001,"<0.001",sprintf("%.3f",p)))
status <- tibble(analysis=character(),status=character(),detail=character())
warnings_found <- character()
record <- function(analysis,state,detail="") {
  status <<- bind_rows(status,tibble(analysis=analysis,status=state,detail=detail))
  write_csv(status,file.path(DIR$logs,"BPD_Analysis_Status.csv"))
}
capture_analysis <- function(name,expr) {
  tryCatch(withCallingHandlers({
    value <- force(expr)
    record(name,"COMPLETE")
    value
  },warning=function(w) {
    warnings_found <<- c(warnings_found,paste0(name,": ",conditionMessage(w)))
  }),error=function(e) {
    record(name,"NOT ESTIMABLE / REVIEW REQUIRED",conditionMessage(e))
    message(name,": ",conditionMessage(e));NULL
  })
}
export <- function(d,name,title,note="",digits=NULL) {
  write_csv(as.data.frame(d),file.path(DIR$tables,paste0(name,".csv")),na="")
  display <- as.data.frame(d)
  if (!is.null(digits)) for (v in names(display)) {
    if (is.numeric(display[[v]])) display[[v]] <- round(display[[v]],digits)
  }
  ft <- flextable::flextable(display)
  ft <- flextable::theme_vanilla(ft)
  ft <- flextable::fontsize(ft,size=10,part="all")
  ft <- flextable::padding(ft,padding=5,part="all")
  ft <- flextable::autofit(ft)
  ft <- flextable::set_table_properties(ft,layout="autofit",width=1)
  ft <- flextable::set_caption(ft,caption=title)
  if(nzchar(note)) ft <- flextable::add_footer_lines(ft,values=note)
  doc <- officer::read_docx()
  doc <- flextable::body_add_flextable(doc,ft)
  print(doc,target=file.path(DIR$tables,paste0(name,".docx")))
}
joint_wald <- function(b,V) {
  ii <- grep(":hc$|^hc:",names(b))
  assert(length(ii)==3,"Expected three interaction coefficients")
  a <- b[ii];S <- V[ii,ii,drop=FALSE]
  assert(all(is.finite(S))&&qr(S)$rank==3,"Invalid interaction covariance")
  W <- as.numeric(t(a)%*%solve(S,a))
  tibble(statistic=W,df=3L,p_value=pchisq(W,3,lower.tail=FALSE))
}
pheno_ORs <- function(b,V) {
  assert(all(is.finite(b))&&all(is.finite(V)),"Nonfinite coefficients/covariance")
  bind_rows(lapply(PHENO,function(p) {
    g <- setNames(rep(0,length(b)),names(b));assert("hc" %in% names(g),"Treatment coefficient missing");g["hc"] <- 1
    if(p!=PHENO[1]) {
      nm <- paste0("phenotype",p,":hc")
      assert(nm %in% names(g),paste("Contrast coefficient missing:",nm));g[nm] <- 1
    }
    a <- sum(g*b);se <- sqrt(as.numeric(t(g)%*%V%*%g))
    tibble(phenotype=p,OR_higher_severity=exp(a),lower=exp(a-Z*se),upper=exp(a+Z*se),
           Wald_P=2*pnorm(abs(a/se),lower.tail=FALSE))
  }))
}

#####
#load exactly the same BPD cohort-----
#####
d <- readRDS(data_file)
required <- c("infant_id","phenotype","hc","gestage","male","center","death36","bpd_jen36_hc","jensen_grade")
assert(all(required %in% names(d)),"The saved analytic dataset lacks required columns")
assert(!anyDuplicated(d$infant_id),"Duplicate infant keys in saved analytic cohort")
od <- d %>% filter(death36==0,!is.na(bpd_jen36_hc)) %>% droplevels()
od$phenotype <- factor(od$phenotype,levels=PHENO)
od$jensen_grade <- ordered(od$bpd_jen36_hc,levels=0:3,labels=c("No BPD","Grade 1","Grade 2","Grade 3"))
# Centering changes the intercept parameterization, not the gestational-age effect.
od$ga_c <- od$gestage-mean(od$gestage)
assert(nrow(od)==704,"Expected the same 704 infants with observed survivor BPD grade")
assert(!anyNA(od[,required]),"Missing required model fields in the ordinal cohort")
counts <- od %>% count(center,jensen_grade,.drop=FALSE)
export(counts,"BPD_Center_Grade_Counts","BPD severity by center",
       "Zero category counts are retained; these can cause boundary estimates when center-specific thresholds are fitted.")
saveRDS(od,file.path(DIR$models,"BPD_704_Analysis_Data.rds"))
primary_table <- file.path(RUN_DIR,"tables","Table2_Unrounded.csv")
primary_md5 <- unname(tools::md5sum(primary_table))

#####
#partial proportional odds with center specific thresholds-----
#####
# Center is moved to nominal; it must NOT also appear in the location formula.
# Do not accept an unstable fit merely because it produces a P value.
ppo <- capture_analysis("Center-specific-threshold ordinal model",{
  f <- ordinal::clm(jensen_grade ~ phenotype*hc + ga_c + male,
                    nominal=~center,data=od,link="logit",model=TRUE,
                    control=ordinal::clm.control(maxIter=200))
  code <- if(is.list(f$convergence)) f$convergence$code else f$convergence
  eigenvalues <- if(all(is.finite(f$Hessian)))
    eigen(f$Hessian,symmetric=TRUE,only.values=TRUE)$values else NA_real_
  diag <- list(convergence=f$convergence,max_gradient=f$maxGradient,
                Hessian_eigenvalues=eigenvalues,condition_number=f$cond.H,
                coefficients=coef(f),zero_center_grade_cells=counts[counts$n==0,])
  # Save the attempted fit and diagnostics even when it is not usable.
  saveRDS(list(fit=f,diagnostics=diag),file.path(DIR$models,"BPD_Center_PPO_Attempt.rds"))
  writeLines(capture.output(str(diag)),file.path(DIR$logs,"BPD_Center_PPO_Diagnostics.txt"))
  assert(is.numeric(code)&&length(code)>0&&all(is.finite(code))&&all(code==0),
    "Flexible ordinal model did not pass convergence/identifiability checks; inspect its diagnostics, do not report its ORs")
  assert(is.finite(f$maxGradient)&&f$maxGradient<1e-4,"Flexible ordinal model gradient check failed")
  assert(all(is.finite(eigenvalues))&&all(eigenvalues>0),"Flexible ordinal Hessian is not positive definite")
  b <- coef(f);V <- vcov(f)
  est <- pheno_ORs(b,V);test <- joint_wald(b,V)
  export(est,"BPD_Center_PPO","BPD severity allowing center-specific thresholds",
    "Exploratory partial proportional-odds model. Adjusted for gestational age and sex; center effects differ by severity threshold. OR >1 indicates higher severity. Model-based covariance; sibling dependence not accounted for. Survivor-only analysis.",4)
  export(test,"BPD_Center_PPO_Interaction","Ordinal treatment-by-phenotype interaction",digits=6)
  # These are checks, not proof that the remaining proportionality assumptions hold.
  nominal <- tryCatch(ordinal::nominal_test(f),error=function(e) {
    record("Remaining ordinal proportionality diagnostics","UNAVAILABLE",conditionMessage(e));NULL
  })
  if(!is.null(nominal)) {
    nt <- as.data.frame(nominal);nt$term <- rownames(nt)
    write_csv(nt,file.path(DIR$tables,"BPD_Center_PPO_Assumptions.csv"))
    if(any(is.na(nt[["Pr(>Chi)"]][nt$term!="<none>"]))) {
      record("Remaining ordinal proportionality diagnostics","PARTLY UNAVAILABLE",
              "At least one term has no diagnostic P value; do not count NA as a passed check")
    } else record("Remaining ordinal proportionality diagnostics","REVIEW RESULTS",
                   "Read term-specific diagnostics before deciding which ordinal model to retain")
  }
  list(fit=f,estimates=est,interaction=test,nominal_test=nominal)
})

#####
#threshold specific Firth models without proportional odds restriction-----
#####
# Three prespecified thresholds, all fitted regardless of significance:
# grade >=1 versus 0; grade >=2 versus 0-1; grade 3 versus 0-2.
# Firth handles the center with no grade-0 infants without dropping that center.
# These are separate binary analyses, not one joint ordinal interaction test.
threshold_models <- list();threshold_tables <- list();threshold_tests <- list()
for(k in 1:3) {
  name <- paste0("Firth Jensen grade >=",k)
  answer <- capture_analysis(name,{
    td <- od;td$higher_grade <- as.integer(td$bpd_jen36_hc>=k)
    ctrl <- logistf::logistf.control(maxit=300)
    f <- logistf::logistf(higher_grade ~ phenotype*hc + ga_c + male + center,
                          data=td,pl=FALSE,control=ctrl,model=TRUE)
    diagnostics <- list(convergence=f$conv,iterations=f$iter,control=f$control)
    saveRDS(list(fit=f,data=td,diagnostics=diagnostics),
              file.path(DIR$models,paste0("BPD_Threshold_",k,"_Firth_Attempt.rds")))
    writeLines(capture.output(str(diagnostics)),
                 file.path(DIR$logs,paste0("BPD_Threshold_",k,"_Diagnostics.txt")))
    assert(length(f$conv)==3&&all(is.finite(f$conv)),"Firth convergence diagnostics unavailable")
    tolerances <- c(f$control$lconv,f$control$gconv,f$control$xconv)
    assert(length(tolerances)==3&&all(abs(f$conv)<=tolerances),"Firth convergence criteria were not all satisfied")
    b <- coef(f);V <- vcov(f)
    est <- pheno_ORs(b,V) %>% mutate(threshold=paste0("Jensen grade >=",k))
    ii <- grep(":hc$|^hc:",names(b));assert(length(ii)==3,"Wrong interaction block")
    # Penalized likelihood ratio test, with all three interactions fixed to zero.
    # Numeric coefficient positions are validated above.
    critical_warnings <- character()
    lr <- withCallingHandlers(logistf::logistftest(f,test=ii,values=rep(0,3),control=ctrl),
      warning=function(w) {
        critical_warnings <<- c(critical_warnings,conditionMessage(w))
      })
    assert(!length(critical_warnings),paste("Penalized LR test warning:",paste(critical_warnings,collapse="; ")))
    assert(length(lr$prob)==1&&is.finite(lr$prob),"Invalid penalized LR P value")
    test <- tibble(threshold=paste0("Jensen grade >=",k),df=lr$df,
                   interaction_P=as.numeric(lr$prob),test="Firth penalized likelihood ratio")
    n <- td %>% group_by(phenotype,hc) %>% summarise(events=sum(higher_grade),N=n())
    export(n,paste0("BPD_Threshold_",k,"_Counts"),paste("Observed",test$threshold,"counts"))
    saveRDS(list(fit=f,data=td,estimates=est,interaction=test,LR_test=lr),
                file.path(DIR$models,paste0("BPD_Threshold_",k,"_Firth.rds")))
    list(estimates=est,test=test,fit=f)
  })
  if(!is.null(answer)) {
    threshold_models[[as.character(k)]] <- answer$fit
    threshold_tables[[as.character(k)]] <- answer$estimates
    threshold_tests[[as.character(k)]] <- answer$test
  }
}
if(length(threshold_tables)>0) {
  estimates <- bind_rows(threshold_tables)
  tests <- bind_rows(threshold_tests)
  # Adjust across the defined family of THREE threshold interaction tests.
  # n=3 remains fixed even if a threshold could not be fitted.
  tests$interaction_P_Holm <- p.adjust(tests$interaction_P,method="holm",n=3)
  export(estimates,"BPD_Threshold_Firth","Treatment contrasts at each BPD threshold",
    "Firth penalized logistic models adjusted for gestational age, sex, and center, separately at each threshold. Model-based Wald 95% CIs and row P values. OR >1 indicates higher severity. Row tests are unadjusted for multiplicity. Sibling dependence is not modeled. Survivor-only exploratory associations.",4)
  export(tests,"BPD_Threshold_Interactions","Treatment-by-phenotype interactions at three BPD thresholds",
    "Raw penalized likelihood-ratio P values and Holm adjustment across the three threshold interaction tests. No combined ordinal P value is calculated from these correlated outcomes.",6)
  pdat <- estimates %>% mutate(phenotype=factor(phenotype,levels=rev(PHENO)))
  p <- ggplot(pdat,aes(OR_higher_severity,phenotype))+
    geom_vline(xintercept=1,linetype=2,color="grey50")+
    geom_segment(aes(x=lower,xend=upper,yend=phenotype),linewidth=.7)+
    geom_point(size=2.5,color="#2878B5")+scale_x_log10()+
    facet_wrap(~threshold,nrow=1)+theme_bw(base_size=11)+
    labs(title="Hydrocortisone and BPD severity thresholds among survivors",
      subtitle="Separate center-adjusted Firth models; no common effect assumed across severity thresholds",
      x="Odds ratio for higher severity (95% Wald CI)",y=NULL,
      caption="Exploratory associations. No phenotype-guided treatment rule established. See table for threshold-specific global interaction tests.")
  ggsave(file.path(DIR$figures,"BPD_Threshold_Firth.pdf"),p,width=12,height=5,
          device=grDevices::pdf,useDingbats=FALSE)
  ggsave(file.path(DIR$figures,"BPD_Threshold_Firth.png"),p,width=12,height=5,dpi=300,bg="white")
  saveRDS(p,file.path(DIR$figures,"BPD_Threshold_Firth_plot.rds"))
}

#####
#save review notes and confirm primary results untouched-----
#####
assert(identical(primary_md5,unname(tools::md5sum(primary_table))),"Primary Table 2 changed unexpectedly")
saveRDS(list(partial_proportional_odds=ppo,threshold_models=threshold_models,
             threshold_results=threshold_tables,threshold_tests=threshold_tests,
             status=status,source_data=data_file),file.path(DIR$models,"BPD_Flexibility_Bundle.rds"))
writeLines(unique(warnings_found),file.path(DIR$logs,"BPD_Warnings.txt"))
writeLines(capture.output(sessionInfo()),file.path(DIR$logs,"sessionInfo.txt"))
writeLines(c(
  "EXPLORATORY BPD MODEL CHECK - RESULTS REQUIRE REVIEW",
  "The source cohort and primary extubation analysis were not changed.",
  "Center-specific thresholds address the detected center proportional-odds restriction, if this larger model is estimable.",
  "A center with no observed grade-0 infants can cause boundary/identifiability problems; such a model is not accepted just because it returns estimates.",
  "Three separate Firth logistic threshold models provide analyses without a proportional-odds restriction. They are not a single ordinal model.",
  "Firth row CIs/tests are Wald and model-based; global interaction tests are penalized likelihood-ratio tests, with Holm adjustment across three thresholds.",
  "No family linkage is available; sibling dependence is not accounted for in these models.",
  "All analyses condition on survival and available BPD severity, and may be affected by survivor selection.",
  "Do not equate nonsignificance with no effect or choose models based on favorable treatment P values.",
  "Do not finalize the old ordinal P=0.5906 before reviewing these diagnostics and results.",
  if(is.null(ppo)) "Center-specific ordinal model: NOT ESTIMABLE / REVIEW REQUIRED." else "Center-specific ordinal model: numeric checks passed; remaining assumptions require review.",
  paste("Threshold models completed:",length(threshold_tables),"of 3."),
  paste("Output folder:",OUT)),file.path(OUT,"BPD_Review_Notes.txt"))
message("New BPD outputs: ",OUT)
message("Primary Table 2 unchanged. Threshold models completed: ",length(threshold_tables),"/3.")
message("Send BPD_Analysis_Status.csv, BPD_Warnings.txt, BPD_Center_PPO_Diagnostics.txt, BPD_Threshold_Firth.csv, and BPD_Threshold_Interactions.csv.")
