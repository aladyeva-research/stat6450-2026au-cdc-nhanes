# NHANES: DATA CLEANING AND LASSO FEATURE SELECTION
# Plain R script, exported November 15, 2026.
#
# Use the data/raw/ and inputs/ folders from the earlier project handoff ZIP.
# Put these .R files in that extracted project folder, set it as your working
# directory, then run:
#   source("setup_nhanes_R.R")
#   source("NHANES_Cleaning_LASSO.R")
# No Rmd, Knit, knitr, rmarkdown or Pandoc is required.
# Results go to results/cleaning_lasso_R/ by default.
# Optional overrides, set before sourcing this script:
#   nhanes_params <- list(repeated_cv_repeats=20,
#                         output_dir="results/cleaning_lasso_R")
#   options(nhanes.project_dir="C:/path/to/extracted/project")
#
# Outcome: existing self-reported diagnosed diabetes, DIQ010 1 versus 2.
# Starting groups: BMI, age, sex, race/ethnicity and education.
# Logistic LASSO tunes lambda using weighted log loss, not Gaussian SSE.
# Selection frequencies are descriptive and do not establish significance.
# Method references: doi:10.18637/jss.v033.i01 and doi:10.1186/1471-2105-7-91.

# ---- setup ----
options(stringsAsFactors=FALSE,width=100)
# Locate the existing project by its input folders, without requiring an Rmd.
document_dir<-function(file) {
  if(is.character(file)&&length(file)==1L&&!is.na(file)&&nzchar(file)) dirname(file) else character()
}
script_args<-grep("^--file=",commandArgs(FALSE),value=TRUE)
source_files<-unlist(lapply(sys.frames(),function(frame) {
  if(exists("ofile",envir=frame,inherits=FALSE))get("ofile",envir=frame,inherits=FALSE)else NULL
}),use.names=FALSE)
script_candidates<-c(sub("^--file=","",script_args),source_files)
script_candidates<-script_candidates[!is.na(script_candidates)&nzchar(script_candidates)]
script_candidates<-script_candidates[file.exists(script_candidates)]
script_file<-script_candidates[basename(script_candidates)=="NHANES_Cleaning_LASSO.R"]
script_file<-if(length(script_file))normalizePath(tail(script_file,1),winslash="/")else character()
project_candidates<-unique(c(getOption("nhanes.project_dir"),getwd(),
  file.path(getwd(),"nhanes-project","diabetes"),unlist(lapply(script_candidates,document_dir))))
project_candidates<-project_candidates[
  file.exists(file.path(project_candidates,"data/raw/P_DEMO.xpt")) &
  dir.exists(file.path(project_candidates,"inputs/cleaning_lasso_reference"))]
if(!length(project_candidates))stop(
  "Set the working directory to your existing NHANES project folder containing data/raw/ and inputs/cleaning_lasso_reference/. These input folders are in the earlier handoff ZIP.",
  call.=FALSE)
project<-normalizePath(project_candidates[1],winslash="/",mustWork=TRUE)

task_library<-file.path(project,".r-library")
if(dir.exists(task_library)).libPaths(c(task_library,.libPaths()))
required<-c("foreign","glmnet","jsonlite")
missing_packages<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages))stop("Run source(\"setup_nhanes_R.R\") first. Missing packages: ",
  paste(missing_packages,collapse=", "),call.=FALSE)
default_params<-list(output_dir="results/cleaning_lasso_R",repeated_cv_repeats=20,
  repeated_cv_seed=64512000)
analysis_params<-utils::modifyList(default_params,
  if(exists("nhanes_params",inherits=TRUE)&&is.list(nhanes_params))nhanes_params else list())
out<-file.path(project,analysis_params$output_dir)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ref<-file.path(project,"inputs/cleaning_lasso_reference")
save_csv<-function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na="")
read_ref<-function(name)read.csv(file.path(ref,name),check.names=FALSE)
show_table<-function(x,digits=4) {
  display<-x
  numeric_columns<-vapply(display,is.numeric,logical(1))
  display[numeric_columns]<-lapply(display[numeric_columns],round,digits=digits)
  print(display,row.names=FALSE)
  invisible(x)
}
verification<-list()
FACTORS<-c("BMI_5","age_10","female","race","education")
FORM<-y~BMI_5+age_10+female+race+education
LAMBDAS<-10^seq(-4,0,length.out=13)
set.seed(6450L)
canonical<-function(x) {
  x<-sub("^\\(Intercept\\)$","Intercept",x)
  x<-sub("^female1$","female_1.0",x)
  x<-sub("^race([0-9]+)$","race_\\1.0",x)
  sub("^education([0-9]+)$","education_\\1.0",x)
}

# ---- cleaning ----
components<-c("DEMO","DIQ","BMX")
paths<-file.path(project,"data/raw",paste0("P_",components,".xpt"))
stopifnot(all(file.exists(paths)))
hashes<-data.frame(file=basename(paths),md5=unname(tools::md5sum(paths)))
stopifnot(identical(hashes,read_ref("raw_file_md5.csv")))
raw<-setNames(lapply(paths,foreign::read.xport),components)
audit<-do.call(rbind,lapply(components,function(component) {
  x<-raw[[component]]
  stopifnot(!anyNA(x$SEQN),!anyDuplicated(x$SEQN))
  data.frame(component=component,n=nrow(x),duplicate_IDs=sum(duplicated(x$SEQN)))
}))
# foreign decodes genuine XPORT zeros correctly. Never replace arbitrary
# small measurements with zero or alter the original XPT files.
for(v in c("RIDAGEYR","RIDAGEMN","WTMECPRP","INDFMPIR")) {
  x<-raw$DEMO[[v]]
  stopifnot(!any(!is.na(x)&x!=0&abs(x)<1e-60))
}
d<-raw$DEMO
for(component in c("DIQ","BMX")) {
  stopifnot(!length(intersect(setdiff(names(d),"SEQN"),setdiff(names(raw[[component]]),"SEQN"))))
  added<-raw[[component]][match(d$SEQN,raw[[component]]$SEQN),setdiff(names(raw[[component]]),"SEQN"),drop=FALSE]
  rownames(added)<-NULL
  d<-cbind(d,added)
}
stopifnot(nrow(d)==nrow(raw$DEMO),!anyDuplicated(d$SEQN))
d$WTMECPRP[!(d$RIDSTATR %in% 2)]<-0
d$design<-d$RIDSTATR %in% 2 & !is.na(d$WTMECPRP)&d$WTMECPRP>0 &
  complete.cases(d[c("SDMVSTRA","SDMVPSU")])
d$BMI_invalid<-!is.na(d$BMXBMI)&(!is.finite(d$BMXBMI)|d$BMXBMI<=0)
d$BMI_clean<-ifelse(d$BMI_invalid,NA_real_,d$BMXBMI)
d$BMI_5<-(d$BMI_clean-25)/5
d$age_10<-(d$RIDAGEYR-50)/10
d$female<-ifelse(d$RIAGENDR %in% 1:2,d$RIAGENDR-1,NA_real_)
d$race<-ifelse(d$RIDRETH3 %in% c(1,2,3,4,6,7),d$RIDRETH3,NA_real_)
d$education<-ifelse(d$DMDEDUC2 %in% 1:5,d$DMDEDUC2,NA_real_)
d$y<-ifelse(d$DIQ010 %in% 1:2,as.numeric(d$DIQ010==1),NA_real_)
d$PSU_ID<-paste(d$SDMVSTRA,d$SDMVPSU,sep="_")
d$eligible_adult_domain<-d$design & !is.na(d$RIDAGEYR)&d$RIDAGEYR>=20 & !(d$RIDEXPRG %in% 1)
d$eligible_known_outcome<-d$eligible_adult_domain & !is.na(d$y)
d$eligible<-d$eligible_known_outcome & complete.cases(d[FACTORS])
d$BMI_missing<-is.na(d$BMI_5);d$education_missing<-is.na(d$education)
d$age_topcoded_80plus<-d$RIDAGEYR %in% 80
d$BMI_above_50_review_only<-!is.na(d$BMXBMI)&d$BMXBMI>50

# Verify the rebuilt analysis variables against the previous cleaned data.
# The reference is used for verification only, not as a substitute for cleaning.
old<-read_ref("analysis_reference.csv")
stopifnot(identical(as.numeric(d$SEQN),as.numeric(old$SEQN)))
errors<-do.call(rbind,lapply(setdiff(names(old),"SEQN"),function(v) {
  x<-d[[v]];z<-old[[v]];stopifnot(identical(is.na(x),is.na(z)))
  ok<-!is.na(x)
  if(is.numeric(x)) {
    error<-max(abs(x[ok]-z[ok]),0);stopifnot(error<1e-9)
  } else {
    stopifnot(identical(tolower(as.character(x[ok])),tolower(as.character(z[ok]))));error<-0
  }
  data.frame(variable=v,maximum_error=error)
}))
flow<-data.frame(step=c("Demographic records","MEC with positive weight","Adults 20+, not known pregnant",
  "Definite yes/no outcome","Complete five predictors","Diabetes cases in complete cases"),
  n=c(nrow(d),sum(d$design),sum(d$eligible_adult_domain),sum(d$eligible_known_outcome),
    sum(d$eligible),sum(d$y[d$eligible])))
stopifnot(identical(as.numeric(flow$n),c(15560,14300,8457,8206,8038,1278)))
save_csv(audit,"raw_file_audit.csv");save_csv(flow,"cohort_flow.csv")
save_csv(errors,"cleaning_comparison_with_previous.csv")
save_csv(d,"cleaned_master_with_flags.csv")
save_csv(d[d$eligible_known_outcome,],"cleaned_diabetes_candidates.csv")
save_csv(d[d$eligible,],"cleaned_diabetes_complete_cases.csv")
verification$cleaning<-list(status="passed",same_inputs_as_previous=TRUE,
  candidate_n=8206,complete_case_n=8038,complete_case_cases=1278)
show_table(flow,0)

# ---- cleaning-checks ----
a<-d[d$eligible_known_outcome,]
missingness<-do.call(rbind,lapply(FACTORS,function(v)data.frame(variable=v,n=nrow(a),
  missing=sum(is.na(a[[v]])),missing_percent=100*mean(is.na(a[[v]])))))
category_counts<-do.call(rbind,lapply(c("female","race","education"),function(v) {
  do.call(rbind,lapply(sort(unique(a[[v]][!is.na(a[[v]])])),function(level) {
    known<-!is.na(a[[v]])&a[[v]]==level
    cc<-known&a$eligible
    data.frame(variable=v,code=level,candidate_n=sum(known),complete_case_n=sum(cc),
      complete_case_cases=sum(a$y[cc]),complete_case_noncases=sum(1-a$y[cc]))
  }))
}))
cleaning_decisions<-data.frame(
  check=c("Duplicate IDs","Unexpected sex/race/education/diabetes codes among adults",
    "Nonpositive/nonfinite BMI among adults","Missing BMI in known-outcome adults",
    "Missing/invalid derived education in known-outcome adults","BMI >50 in complete cases",
    "Age released as 80+ in complete cases","Borderline diabetes among eligible adults"),
  n=c(sum(audit$duplicate_IDs),
    sum(d$eligible_adult_domain & ((!is.na(d$RIAGENDR)&!(d$RIAGENDR %in% 1:2)) |
      (!is.na(d$RIDRETH3)&!(d$RIDRETH3 %in% c(1,2,3,4,6,7))) |
      (!is.na(d$DMDEDUC2)&!(d$DMDEDUC2 %in% c(1:5,7,9))) |
      (!is.na(d$DIQ010)&!(d$DIQ010 %in% c(1,2,3,7,9))))),
    sum(d$eligible_adult_domain&d$BMI_invalid),sum(is.na(a$BMXBMI)),sum(is.na(a$education)),
    sum(d$eligible&d$BMI_above_50_review_only),sum(d$eligible&d$age_topcoded_80plus),
    sum(d$eligible_adult_domain&d$DIQ010 %in% 3)),
  decision=c("No rows removed","Validate against documented codes; preserve raw codes",
    "Set invalid derived BMI missing; no invalid adult values found",
    "Keep in candidate data; exclude from complete-case LASSO",
    "7/9 are missing, not education levels; preserve raw codes",
    "Keep and flag; high BMI alone is not evidence of error",
    "Keep as 80+; do not invent precise ages","Exclude from the primary yes/no outcome"))
stopifnot(sum(is.na(a$BMXBMI))==158,sum(is.na(a$education))==12,
  sum(is.na(a$BMXBMI)&is.na(a$education))==2,
  all(category_counts$complete_case_cases>=5),all(category_counts$complete_case_noncases>=5),
  all(vapply(d[d$eligible,FACTORS],function(x)length(unique(x))>1,logical(1))))
save_csv(missingness,"missingness.csv");save_csv(category_counts,"category_counts.csv")
save_csv(cleaning_decisions,"cleaning_decisions.csv")
show_table(missingness)
show_table(category_counts,0)
show_table(cleaning_decisions,0)

# ---- missingness-comparison ----
# Compare who is excluded; missing-data handling is a scientific decision,
# not a way to remove participants until apparent model fit improves.
included_excluded<-do.call(rbind,lapply(c("Included_complete","Excluded_missing_predictor"),function(group) {
  mask<-if(group=="Included_complete")d$eligible else d$eligible_known_outcome & !d$eligible
  z<-d[mask,]
  data.frame(group=group,n=nrow(z),cases=sum(z$y),
    weighted_mean_age=weighted.mean(z$RIDAGEYR,z$WTMECPRP),
    weighted_diabetes_percent=100*weighted.mean(z$y,z$WTMECPRP))
}))
save_csv(included_excluded,"included_excluded_comparison.csv")
show_table(included_excluded,2)

# ---- predictor-coding ----
d$female<-factor(d$female,levels=0:1)
d$race<-factor(d$race,levels=c(1,2,3,4,6,7))
d$education<-factor(d$education,levels=1:5)
g<-d[d$eligible,];rownames(g)<-NULL
coding<-data.frame(variable=c("BMI_5","age_10","female","race","education"),
  definition=c("(BMI - 25)/5","(Age - 50)/10","0 male (reference); 1 female",
    "1 Mexican American (reference); 2 other Hispanic; 3 non-Hispanic White; 4 non-Hispanic Black; 6 non-Hispanic Asian; 7 other/multiracial",
    "1 <9th grade (reference); 2 9-11th grade; 3 high school/GED; 4 some college/AA; 5 college graduate+"))
x_matrix<-function(frame)model.matrix(FORM,data=frame)[,-1,drop=FALSE]
x<-x_matrix(g)
stopifnot(ncol(x)==12,all(is.finite(x)),qr(cbind(1,x))$rank==13)
owners<-attr(terms(FORM),"term.labels")[attr(model.matrix(FORM,g),"assign")[-1]]
save_csv(coding,"variable_coding.csv")
show_table(coding)

# ---- LASSO-functions ----
roc_points<-function(y,p,w) {
  stopifnot(all(is.finite(p)),all(w>0),setequal(unique(y),0:1))
  z<-aggregate(cbind(case=w*y,control=w*(1-y)),list(cut=p),sum)
  z<-z[order(z$cut,decreasing=TRUE),]
  data.frame(cut=c(Inf,z$cut),TPR=c(0,cumsum(z$case)/sum(z$case)),
    FPR=c(0,cumsum(z$control)/sum(z$control)))
}
scores<-function(y,p,w) {
  r<-roc_points(y,p,w)
  auc<-sum(diff(r$FPR)*(head(r$TPR,-1)+tail(r$TPR,-1))/2)
  p<-pmin(pmax(p,1e-12),1-1e-12)
  c(AUC=auc,log_loss=-weighted.mean(y*log(p)+(1-y)*log1p(-p),w),
    Brier=weighted.mean((y-p)^2,w))
}
lasso_path<-function(train,test) {
  x<-x_matrix(train);xt<-x_matrix(test)
  stopifnot(identical(colnames(x),colnames(xt)))
  w<-train$WTMECPRP/sum(train$WTMECPRP)
  mu<-colSums(x*w);sd<-sqrt(colSums(sweep(x,2,mu,"-")^2*w));sd[sd<1e-12]<-1
  xs<-sweep(sweep(x,2,mu,"-"),2,sd,"/")
  xv<-sweep(sweep(xt,2,mu,"-"),2,sd,"/")
  args<-list(x=xs,y=train$y,family="binomial",alpha=1,weights=w*nrow(train),
    lambda=sort(LAMBDAS,decreasing=TRUE),standardize=FALSE,intercept=TRUE)
  controls<-list(thresh=1e-14,maxit=1000000)
  if("control" %in% names(formals(glmnet::glmnet)))args$control<-controls else args<-c(args,controls)
  fit<-do.call(glmnet::glmnet,args)
  stopifnot(fit$jerr==0)
  p<-as.matrix(predict(fit,newx=xv,s=LAMBDAS,type="response"))
  b<-as.matrix(coef(fit,s=LAMBDAS))
  eta<-cbind(1,xs)%*%b
  gradient<-crossprod(xs,(plogis(eta)-train$y)*w)
  intercept_error<-abs(colSums((plogis(eta)-train$y)*w))
  kkt<-vapply(seq_along(LAMBDAS),function(j) {
    bj<-b[-1,j];gr<-gradient[,j];lam<-LAMBDAS[j]
    err<-ifelse(abs(bj)>1e-6,abs(gr+lam*sign(bj)),pmax(abs(gr)-lam,0))
    max(err,intercept_error[j])
  },numeric(1))
  stopifnot(max(kkt)<2e-5,all(is.finite(p)),all(p>=0&p<=1))
  raw<-b;raw[-1,]<-sweep(b[-1,,drop=FALSE],1,sd,"/")
  raw[1,]<-b[1,]-colSums(raw[-1,,drop=FALSE]*mu)
  list(p=p,beta=raw,KKT=kkt,training_mean=mu,training_SD=sd)
}
assert_folds<-function(frame,fold) {
  stopifnot(length(fold)==nrow(frame),!anyNA(fold),setequal(unique(fold),1:5),
    all(vapply(split(fold,frame$PSU_ID),function(z)length(unique(z))==1,logical(1))))
  for(k in 1:5)stopifnot(!length(intersect(frame$PSU_ID[fold==k],frame$PSU_ID[fold!=k])),
    length(unique(frame$y[fold==k]))==2,length(unique(frame$y[fold!=k]))==2)
  invisible(TRUE)
}
tune_lasso<-function(train,fold,audit=FALSE) {
  assert_folds(train,fold)
  loss<-matrix(NA_real_,5,length(LAMBDAS))
  preprocessing<-list()
  for(k in 1:5) {
    tr<-train[fold!=k,];te<-train[fold==k,];path<-lasso_path(tr,te)
    # Tuning needs only log loss. This is the same clipped binomial loss as
    # scores(), without computing an unused ROC curve for every lambda.
    for(j in seq_along(LAMBDAS)) {
      pp<-pmin(pmax(path$p[,j],1e-12),1-1e-12)
      loss[k,j]<- -weighted.mean(te$y*log(pp)+(1-te$y)*log1p(-pp),te$WTMECPRP)
    }
    if(audit)preprocessing[[k]]<-data.frame(inner_fold=k,term=canonical(names(path$training_mean)),
      training_mean=unname(path$training_mean),training_SD=unname(path$training_SD),
      training_n=nrow(tr),validation_n=nrow(te))
  }
  avg<-colMeans(loss);se<-apply(loss,2,sd)/sqrt(5)
  j<-which.min(avg);one<-max(which(avg<=avg[j]+se[j]))
  list(choices=c(minimum=LAMBDAS[j],one_se=LAMBDAS[one]),
    table=data.frame(lambda=LAMBDAS,mean_log_loss=avg,SE=se,minimum=seq_along(LAMBDAS)==j,
      one_se=seq_along(LAMBDAS)==one),fold_loss=loss,
    preprocessing=if(audit)do.call(rbind,preprocessing)else NULL)
}
outer_map<-read_ref("outer_PSU_folds.csv")
inner_map<-read_ref("inner_PSU_folds.csv")
outer<-outer_map$fold[match(g$PSU_ID,outer_map$PSU_ID)]
assert_folds(g,outer)

# ---- nested-CV ----
predictions<-list();tuning<-list();fold_choices<-list();fold_coefficients<-list()
for(k in 1:5) {
  tr<-g[outer!=k,];te<-g[outer==k,]
  map<-inner_map[inner_map$outer_fold==k,]
  inner<-map$inner_fold[match(tr$PSU_ID,map$PSU_ID)]
  tune<-tune_lasso(tr,inner);path<-lasso_path(tr,te)
  tuning[[k]]<-cbind(outer_fold=k,tune$table)
  for(rule in c("minimum","one_se")) {
    j<-match(tune$choices[rule],LAMBDAS)
    predictions[[length(predictions)+1]]<-data.frame(fold=k,method=paste0("LASSO_",rule),
      SEQN=te$SEQN,y=te$y,weight=te$WTMECPRP,probability=path$p[,j])
    selected<-FACTORS[FACTORS %in% owners[abs(path$beta[-1,j])>1e-6]]
    fold_choices[[length(fold_choices)+1]]<-data.frame(fold=k,rule=rule,lambda=LAMBDAS[j],
      selected_factors=paste(selected,collapse=" | "),KKT_max=path$KKT[j])
    fold_coefficients[[length(fold_coefficients)+1]]<-data.frame(fold=k,rule=rule,
      term=canonical(rownames(path$beta)),coefficient=path$beta[,j])
  }
}
predictions<-do.call(rbind,predictions)
CV<-do.call(rbind,lapply(unique(predictions$method),function(method) {
  z<-predictions[predictions$method==method,]
  stopifnot(nrow(z)==8038,!anyDuplicated(z$SEQN),setequal(z$SEQN,g$SEQN))
  cbind(data.frame(method=method,n=nrow(z)),as.data.frame(as.list(scores(z$y,z$probability,z$weight))))
}))
old<-read_ref("previous_CV.csv");at<-match(CV$method,old$method)
comparison<-CV
for(v in c("AUC","log_loss","Brier"))comparison[[paste0("delta_",v)]]<-CV[[v]]-old[[v]][at]
stopifnot(!anyNA(at),max(abs(comparison$delta_AUC))<1e-4,
  max(abs(comparison$delta_log_loss))<1e-5,max(abs(comparison$delta_Brier))<1e-5)
save_csv(predictions,"LASSO_heldout_predictions.csv")
save_csv(CV,"LASSO_CV_performance.csv")
save_csv(comparison,"LASSO_comparison_with_previous.csv")
save_csv(do.call(rbind,tuning),"LASSO_inner_tuning.csv")
save_csv(do.call(rbind,fold_choices),"LASSO_outer_fold_choices.csv")
save_csv(do.call(rbind,fold_coefficients),"LASSO_outer_coefficients.csv")
verification$nested_CV<-list(status="passed",outer_folds=5,inner_folds=5,PSU_overlap=0,
  max_AUC_difference=max(abs(comparison$delta_AUC)),
  max_log_loss_difference=max(abs(comparison$delta_log_loss)),
  max_Brier_difference=max(abs(comparison$delta_Brier)))
show_table(CV,6)
show_table(do.call(rbind,fold_choices),6)

# ---- final-LASSO ----
full_tune<-tune_lasso(g,outer)
full_path<-lasso_path(g,g)
coefficient_tables<-list();selection_tables<-list()
for(rule in c("minimum","one_se")) {
  j<-match(full_tune$choices[rule],LAMBDAS)
  b<-full_path$beta[,j]
  coefficient_tables[[rule]]<-data.frame(method=paste0("LASSO_",rule),lambda=LAMBDAS[j],
    term=canonical(names(b)),factor=c("Intercept",owners),coefficient=unname(b),
    selected=c(NA,abs(b[-1])>1e-6))
  selection_tables[[rule]]<-data.frame(method=paste0("LASSO_",rule),lambda=LAMBDAS[j],
    selected_factors=paste(FACTORS[FACTORS %in% owners[abs(b[-1])>1e-6]],collapse=" | "),
    nonzero_columns=sum(abs(b[-1])>1e-6),KKT_max=full_path$KKT[j])
}
coefficients<-do.call(rbind,coefficient_tables)
selection<-do.call(rbind,selection_tables)
old<-read_ref("previous_coefficients.csv")
at<-match(paste(coefficients$method,coefficients$term),paste(old$method,old$term))
coefficient_error<-max(abs(coefficients$coefficient-old$coefficient[at]))
old_selection<-read_ref("previous_selection.csv")
ats<-match(selection$method,old_selection$method)
stopifnot(!anyNA(at),coefficient_error<1e-4,max(abs(selection$lambda-old_selection$penalty[ats]))<1e-12,
  all(selection$selected_factors==old_selection$features[ats]))
save_csv(full_tune$table,"LASSO_full_data_tuning.csv")
save_csv(coefficients,"LASSO_final_coefficients.csv")
save_csv(selection,"LASSO_selected_predictors.csv")
verification$final_LASSO<-list(status="passed",maximum_coefficient_difference=coefficient_error,
  same_selected_factors=TRUE,same_lambda_choices=TRUE)
model_export<-list(outcome="Existing self-reported diagnosed diabetes; DIQ010 1 versus 2",
  predictor_formula=paste(deparse(FORM),collapse=" "),n=nrow(g),
  transformations=list(BMI_5="(BMI-25)/5",age_10="(age-50)/10"),
  references=list(sex="Male",race="Mexican American",education="Less than 9th grade"),
  fits=lapply(c("minimum","one_se"),function(rule) {
    j<-match(full_tune$choices[rule],LAMBDAS)
    list(rule=rule,lambda=LAMBDAS[j],coefficients=as.list(full_path$beta[,j]))
  }),prediction="Construct the same factor dummies, eta=intercept+sum(x*coefficient), probability=plogis(eta). Coefficients are already back-transformed; do not standardize again.")
jsonlite::write_json(model_export,file.path(out,"LASSO_model.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
saveRDS(model_export,file.path(out,"LASSO_model_parameters.rds"))
show_table(selection,6)
show_table(coefficients[c("method","term","factor","coefficient","selected")],6)

# ---- selection-stability ----
outer_coefs<-do.call(rbind,fold_coefficients)
term_owner<-function(term) {
  ifelse(grepl("^female_",term),"female",ifelse(grepl("^race_",term),"race",
    ifelse(grepl("^education_",term),"education",term)))
}
outer_coefs$factor<-term_owner(outer_coefs$term)
stability_rows<-list()
for(rule in c("minimum","one_se"))for(factor in FACTORS) {
  z<-outer_coefs[outer_coefs$rule==rule & outer_coefs$factor==factor,]
  retained<-vapply(1:5,function(k)any(abs(z$coefficient[z$fold==k])>1e-6),logical(1))
  stability_rows[[length(stability_rows)+1]]<-data.frame(rule=rule,factor=factor,
    retained_folds=sum(retained),outer_folds=5,fraction=mean(retained))
}
selection_stability<-do.call(rbind,stability_rows)
contrast_stability<-do.call(rbind,lapply(split(outer_coefs[outer_coefs$term!="Intercept",],
  interaction(outer_coefs$rule[outer_coefs$term!="Intercept"],outer_coefs$term[outer_coefs$term!="Intercept"],drop=TRUE)),function(z) {
    data.frame(rule=z$rule[1],term=z$term[1],retained_folds=sum(abs(z$coefficient)>1e-6),
      outer_folds=length(unique(z$fold)),minimum_coefficient=min(z$coefficient),maximum_coefficient=max(z$coefficient))
}))
save_csv(selection_stability,"LASSO_outer_factor_retention.csv")
save_csv(contrast_stability,"LASSO_outer_contrast_retention.csv")
show_table(selection_stability,2)

# ---- repeated-nested-functions ----
random_psu_folds<-function(frame,seed) {
  set.seed(seed)
  units<-sort(unique(as.character(frame$PSU_ID)))
  counts<-table(factor(as.character(frame$PSU_ID),levels=units))
  units<-units[sample.int(length(units))]
  load<-numeric(5);allocation<-setNames(integer(length(units)),units)
  for(unit in units) {
    available<-which(load==min(load))
    j<-available[sample.int(length(available),1)]
    allocation[unit]<-j;load[j]<-load[j]+counts[unit]
  }
  fold<-unname(allocation[as.character(frame$PSU_ID)])
  assert_folds(frame,fold)
  fold
}
partition_signature<-function(PSU_ID,fold) {
  z<-unique(data.frame(PSU_ID=as.character(PSU_ID),fold=fold))
  z<-z[order(z$PSU_ID),]
  stopifnot(!anyDuplicated(z$PSU_ID))
  paste(match(z$fold,unique(z$fold)),collapse=",")
}
run_lasso_repeat<-function(repeat_id,frame,base_seed) {
  outer_seed<-as.integer(base_seed+repeat_id)
  outer_fold<-random_psu_folds(frame,outer_seed)
  rows<-list();coefs<-list();flags<-list();tuning_rows<-list();audit_rows<-list()
  inner_ledger<-list();outer_preprocessing<-list();inner_preprocessing<-list()
  outer_ledger<-unique(data.frame(repeat_id=repeat_id,outer_seed=outer_seed,
    PSU_ID=as.character(frame$PSU_ID),outer_fold=outer_fold))
  for(k in 1:5) {
    tr<-frame[outer_fold!=k,];te<-frame[outer_fold==k,]
    inner_seed<-as.integer(base_seed+1000L*repeat_id+k)
    inner<-random_psu_folds(tr,inner_seed)
    stopifnot(!length(intersect(tr$PSU_ID,te$PSU_ID)))
    inner_ledger[[k]]<-unique(data.frame(repeat_id=repeat_id,outer_fold=k,inner_seed=inner_seed,
      PSU_ID=as.character(tr$PSU_ID),inner_fold=inner))
    tune<-tune_lasso(tr,inner,audit=TRUE)
    path<-lasso_path(tr,te)
    tuning_rows[[k]]<-cbind(repeat_id=repeat_id,outer_fold=k,tune$table)
    inner_preprocessing[[k]]<-cbind(repeat_id=repeat_id,outer_fold=k,tune$preprocessing)
    outer_preprocessing[[k]]<-data.frame(repeat_id=repeat_id,outer_fold=k,
      term=canonical(names(path$training_mean)),training_mean=unname(path$training_mean),
      training_SD=unname(path$training_SD),training_n=nrow(tr),test_n=nrow(te))
    for(rule in c("minimum","one_se")) {
      j<-match(tune$choices[rule],LAMBDAS)
      beta<-path$beta[,j]
      reconstructed<-plogis(drop(cbind(1,x_matrix(te))%*%beta))
      error<-max(abs(reconstructed-path$p[,j]))
      stopifnot(error<1e-9)
      rows[[length(rows)+1]]<-data.frame(repeat_id=repeat_id,outer_fold=k,rule=rule,
        SEQN=te$SEQN,PSU_ID=as.character(te$PSU_ID),y=te$y,
        weight=te$WTMECPRP,probability=path$p[,j])
      coefs[[length(coefs)+1]]<-data.frame(repeat_id=repeat_id,outer_fold=k,rule=rule,
        lambda=LAMBDAS[j],term=canonical(names(beta)),factor=c("Intercept",owners),
        coefficient=unname(beta),selected=c(NA,abs(beta[-1])>1e-6))
      for(factor in FACTORS)flags[[length(flags)+1]]<-data.frame(repeat_id=repeat_id,
        outer_fold=k,rule=rule,factor=factor,selected=any(abs(beta[-1][owners==factor])>1e-6))
      audit_rows[[length(audit_rows)+1]]<-data.frame(repeat_id=repeat_id,outer_fold=k,rule=rule,
        outer_seed=outer_seed,inner_seed=inner_seed,lambda=LAMBDAS[j],
        train_n=nrow(tr),test_n=nrow(te),train_PSUs=length(unique(tr$PSU_ID)),
        test_PSUs=length(unique(te$PSU_ID)),PSU_overlap=length(intersect(tr$PSU_ID,te$PSU_ID)),
        KKT_max=path$KKT[j],raw_coefficient_prediction_error=error,
        nonzero_columns=sum(abs(beta[-1])>1e-6),
        selected_factors=paste(FACTORS[FACTORS %in% owners[abs(beta[-1])>1e-6]],collapse=" | "))
    }
  }
  pred<-do.call(rbind,rows)
  performance<-do.call(rbind,lapply(c("minimum","one_se"),function(rule) {
    z<-pred[pred$rule==rule,]
    stopifnot(nrow(z)==nrow(frame),!anyDuplicated(z$SEQN),setequal(z$SEQN,frame$SEQN))
    cbind(data.frame(repeat_id=repeat_id,rule=rule,n=nrow(z)),
      as.data.frame(as.list(scores(z$y,z$probability,z$weight))))
  }))
  list(predictions=pred,performance=performance,coefficients=do.call(rbind,coefs),
    factor_flags=do.call(rbind,flags),fit_audit=do.call(rbind,audit_rows),
    tuning=do.call(rbind,tuning_rows),outer_ledger=outer_ledger,
    inner_ledger=do.call(rbind,inner_ledger),
    outer_preprocessing=do.call(rbind,outer_preprocessing),
    inner_preprocessing=do.call(rbind,inner_preprocessing))
}
validate_repeat<-function(result,repeat_id,frame,base_seed) {
  stopifnot(setequal(result$predictions$repeat_id,repeat_id),
    nrow(result$predictions)==2*nrow(frame),nrow(result$coefficients)==2*5*13,
    nrow(result$factor_flags)==2*5*length(FACTORS),
    all(result$fit_audit$PSU_overlap==0),all(result$fit_audit$KKT_max<2e-5),
    all(result$fit_audit$raw_coefficient_prediction_error<1e-9))
  expected_outer<-random_psu_folds(frame,as.integer(base_seed+repeat_id))
  m<-result$outer_ledger
  stopifnot(!anyDuplicated(m$PSU_ID),setequal(m$PSU_ID,frame$PSU_ID),
    identical(as.integer(m$outer_fold[match(frame$PSU_ID,m$PSU_ID)]),as.integer(expected_outer)))
  for(k in 1:5) {
    tr<-frame[expected_outer!=k,];test_units<-unique(frame$PSU_ID[expected_outer==k])
    inner<-result$inner_ledger[result$inner_ledger$outer_fold==k,]
    stopifnot(!anyDuplicated(inner$PSU_ID),setequal(inner$PSU_ID,tr$PSU_ID),
      !length(intersect(inner$PSU_ID,test_units)))
    expected_inner<-random_psu_folds(tr,as.integer(base_seed+1000L*repeat_id+k))
    stopifnot(identical(as.integer(inner$inner_fold[match(tr$PSU_ID,inner$PSU_ID)]),as.integer(expected_inner)))
  }
  for(rule in c("minimum","one_se")) {
    z<-result$predictions[result$predictions$rule==rule,]
    stopifnot(!anyDuplicated(z$SEQN),setequal(z$SEQN,frame$SEQN))
    at<-match(z$SEQN,frame$SEQN)
    stopifnot(all(z$y==frame$y[at]),all(z$weight==frame$WTMECPRP[at]),
      all(z$outer_fold==expected_outer[at]),all(z$PSU_ID==frame$PSU_ID[at]))
  }
  invisible(TRUE)
}
text_md5<-function(x) {
  path<-tempfile(fileext=".txt")
  on.exit(unlink(path))
  writeLines(x,path,useBytes=TRUE)
  unname(tools::md5sum(path))
}

# ---- repeated-nested-CV ----
REPEATS<-as.integer(analysis_params$repeated_cv_repeats)
BASE_SEED<-as.integer(analysis_params$repeated_cv_seed)
stopifnot(REPEATS>=2,REPEATS<=100,BASE_SEED>0,
  BASE_SEED+1000*REPEATS+5<.Machine$integer.max)
repeat_out<-file.path(out,"repeated_nested_CV")
dir.create(repeat_out,recursive=TRUE,showWarnings=FALSE)
function_hash<-text_md5(c(deparse(body(lasso_path)),deparse(body(tune_lasso)),
  deparse(body(scores)),deparse(body(roc_points)),deparse(body(assert_folds)),
  deparse(body(x_matrix)),deparse(body(canonical)),
  deparse(body(random_psu_folds)),deparse(body(run_lasso_repeat))))
repeat_manifest<-list(method="Repeated nested PSU-grouped weighted LASSO logistic regression",
  repeats=REPEATS,outer_folds=5,inner_folds=5,base_seed=BASE_SEED,
  data_md5=unname(tools::md5sum(file.path(out,"cleaned_diabetes_complete_cases.csv"))),
  function_md5=function_hash,lambda_grid=as.list(LAMBDAS),formula=paste(deparse(FORM),collapse=" "),
  predictors=as.list(FACTORS),selection_tolerance=1e-6,
  factor_levels=list(female=levels(g$female),race=levels(g$race),education=levels(g$education)),
  matrix_columns=as.list(canonical(colnames(x_matrix(g)))),matrix_owners=as.list(owners),
  glmnet_control=list(thresh=1e-14,maxit=1000000),
  glmnet_version=as.character(packageVersion("glmnet")),R_version=R.version.string,
  splitting="Random PSU order; greedy balance by participant count; random ties; outcomes unused",
  preprocessing="Weighted means/SDs fitted separately on every training sample",
  tuning="Mean of five inner weighted log losses; minimum and one-SE rules, as in original analysis",
  outer_seed_rule="base_seed + repeat_id",inner_seed_rule="base_seed + 1000*repeat_id + outer_fold",
  interpretation="Repeated partitions of the same cohort, not independent cohorts or significance tests")
encoded_manifest<-jsonlite::toJSON(repeat_manifest,auto_unbox=TRUE,pretty=TRUE,digits=16)
cache_id<-text_md5(encoded_manifest)
cache_dir<-file.path(repeat_out,"checkpoints",cache_id)
dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
writeLines(encoded_manifest,file.path(repeat_out,"run_manifest.json"))
repeated_results<-vector("list",REPEATS)
for(repeat_id in seq_len(REPEATS)) {
  checkpoint<-file.path(cache_dir,sprintf("repeat_%03d.rds",repeat_id))
  if(file.exists(checkpoint)) {
    result<-readRDS(checkpoint)
  } else {
    result<-run_lasso_repeat(repeat_id,g,BASE_SEED)
    validate_repeat(result,repeat_id,g,BASE_SEED)
    saveRDS(result,checkpoint)
  }
  validate_repeat(result,repeat_id,g,BASE_SEED)
  repeated_results[[repeat_id]]<-result
  cat("Completed repeat",repeat_id,"of",REPEATS,"\n")
  flush.console()
}
combined<-lapply(names(repeated_results[[1]]),function(name)
  do.call(rbind,lapply(repeated_results,function(z)z[[name]])))
names(combined)<-names(repeated_results[[1]])
original_signature<-partition_signature(g$PSU_ID,outer)
signatures<-vapply(repeated_results,function(z)partition_signature(z$outer_ledger$PSU_ID,z$outer_ledger$outer_fold),character(1))
stopifnot(!anyDuplicated(signatures),all(signatures!=original_signature))
for(name in names(combined))write.csv(combined[[name]],file.path(repeat_out,paste0(name,".csv")),row.names=FALSE,na="")
repeat_performance<-combined$performance
performance_summary<-do.call(rbind,lapply(c("minimum","one_se"),function(rule) {
  z<-repeat_performance[repeat_performance$rule==rule,]
  do.call(rbind,lapply(c("AUC","log_loss","Brier"),function(metric) {
    value<-z[[metric]]
    data.frame(rule=rule,metric=metric,repeats=length(value),mean=mean(value),SD=sd(value),
      minimum=min(value),maximum=max(value),
      original=CV[[metric]][CV$method==paste0("LASSO_",rule)])
  }))
}))
factor_stability<-do.call(rbind,lapply(c("minimum","one_se"),function(rule) {
  do.call(rbind,lapply(FACTORS,function(factor) {
    z<-combined$factor_flags[combined$factor_flags$rule==rule & combined$factor_flags$factor==factor,]
    per_repeat<-aggregate(selected~repeat_id,z,mean)$selected
    data.frame(rule=rule,factor=factor,retained_fits=sum(z$selected),total_fits=nrow(z),
      retention_percent=100*mean(z$selected),minimum_repeat_percent=100*min(per_repeat),
      maximum_repeat_percent=100*max(per_repeat))
  }))
}))
repeat_contrasts<-combined$coefficients[combined$coefficients$term!="Intercept",]
contrast_stability_repeated<-do.call(rbind,lapply(split(repeat_contrasts,
  interaction(repeat_contrasts$rule,repeat_contrasts$term,drop=TRUE)),function(z) {
    data.frame(rule=z$rule[1],term=z$term[1],factor=z$factor[1],retained_fits=sum(z$selected),
      total_fits=nrow(z),retention_percent=100*mean(z$selected),
      negative_fits=sum(z$coefficient< -1e-6),positive_fits=sum(z$coefficient>1e-6),
      minimum_coefficient=min(z$coefficient),maximum_coefficient=max(z$coefficient))
}))
factor_sets<-aggregate(repeat_id~rule+selected_factors,combined$fit_audit,length)
names(factor_sets)[3]<-"fits"
factor_sets$percent<-100*factor_sets$fits/(REPEATS*5)
lambda_frequency<-aggregate(repeat_id~rule+lambda,combined$fit_audit,length)
names(lambda_frequency)[3]<-"fits"
lambda_frequency$percent<-100*lambda_frequency$fits/(REPEATS*5)
mean_size<-do.call(rbind,lapply(c("minimum","one_se"),function(rule) {
  z<-combined$fit_audit[combined$fit_audit$rule==rule,]
  data.frame(rule=rule,mean_nonzero_columns=mean(z$nonzero_columns),
    minimum_nonzero_columns=min(z$nonzero_columns),maximum_nonzero_columns=max(z$nonzero_columns))
}))
for(pair in list(list(performance_summary,"performance_summary.csv"),
  list(factor_stability,"factor_retention_summary.csv"),
  list(contrast_stability_repeated,"contrast_retention_summary.csv"),
  list(factor_sets,"selected_factor_sets.csv"),list(lambda_frequency,"lambda_frequency.csv"),
  list(mean_size,"model_size_summary.csv")))
  write.csv(pair[[1]],file.path(repeat_out,pair[[2]]),row.names=FALSE,na="")
repeat_verification<-list(status="passed",planned_repeats=REPEATS,completed_repeats=length(repeated_results),
  new_unique_outer_partitions=length(unique(signatures)),original_partition_reused=FALSE,
  outer_training_fits=REPEATS*5,inner_training_paths=REPEATS*5*5,scored_rules=2,
  cohort_n=nrow(g),PSUs=length(unique(g$PSU_ID)),cases=sum(g$y),
  heldout_predictions_per_rule_per_repeat=nrow(g),total_saved_predictions=nrow(combined$predictions),
  outer_or_inner_PSU_overlap=0,all_observed_outcomes_weights_preserved=TRUE,
  preprocessing_rows_audited=nrow(combined$outer_preprocessing)+nrow(combined$inner_preprocessing),
  maximum_outer_fit_KKT=max(combined$fit_audit$KKT_max),
  maximum_raw_coefficient_prediction_error=max(combined$fit_audit$raw_coefficient_prediction_error),
  manifest_id=cache_id,original_fit_and_lambda_rules_unchanged=TRUE,
  uncertainty="SD/range across dependent data partitions, not population confidence intervals",
  selection="Descriptive coefficient/group retention; not significance or formal stability selection")
jsonlite::write_json(repeat_verification,file.path(repeat_out,"verification.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
show_table(performance_summary,5)
show_table(factor_stability[c("rule","factor","retained_fits","total_fits","retention_percent")],1)
show_table(contrast_stability_repeated[c("rule","term","retained_fits","total_fits","retention_percent")],1)

# ---- repeated-nested-figures ----
draw_repeated_results<-function() {
  op<-par(mfrow=c(1,2),mar=c(4.5,6,3,1));on.exit(par(op))
  colors<-c(minimum="#416b91",one_se="#af7948")
  plot(c(1,2),range(repeat_performance$AUC),xlim=c(.7,2.3),type="n",xaxt="n",xlab="LASSO rule",ylab="Held-out AUC",
    main=paste(REPEATS,"complete nested-CV repeats"))
  axis(1,at=1:2,labels=c("Minimum loss","One SE"))
  for(j in 1:2) {
    rule<-names(colors)[j];z<-repeat_performance[repeat_performance$rule==rule,]
    points(j+.08*sin(z$repeat_id*1.618),z$AUC,pch=16,cex=.7,col=adjustcolor(colors[rule],alpha.f=.7))
    segments(j-.18,mean(z$AUC),j+.18,mean(z$AUC),lwd=2,col=colors[rule])
  }
  pretty_names<-c(BMI_5="BMI",age_10="Age",female="Sex",race="Race",education="Education")
  plot(c(0,100),c(.7,5.3),type="n",yaxt="n",xlab="Retained outer-training fits (%)",ylab="",
    main=paste(REPEATS*5,"fits per rule"))
  axis(2,at=1:5,labels=pretty_names[FACTORS],las=1)
  abline(v=c(0,25,50,75,100),col="grey90",lty=3)
  for(j in 1:2) {
    rule<-names(colors)[j];z<-factor_stability[factor_stability$rule==rule,]
    at<-match(FACTORS,z$factor)
    points(z$retention_percent[at],seq_along(FACTORS)+if(j==1).08 else -.08,
      pch=if(j==1)16 else 17,cex=.9,col=colors[rule])
  }
  legend("topleft",legend=c("Minimum loss","One SE"),pch=c(16,17),col=colors,bty="n",cex=.8)
}
png(file.path(repeat_out,"repeated_nested_CV.png"),width=1650,height=750,res=150)
draw_repeated_results();dev.off()
if(interactive())draw_repeated_results()

# ---- heldout-classification ----
classification<-do.call(rbind,lapply(unique(predictions$method),function(method) {
  z<-predictions[predictions$method==method,]
  positive<-z$probability>=.5;w<-z$weight;y<-z$y
  TP<-sum(w[positive&y==1]);FN<-sum(w[!positive&y==1])
  TN<-sum(w[!positive&y==0]);FP<-sum(w[positive&y==0])
  result<-data.frame(method=method,threshold=.5,n=nrow(z),
    sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),
    precision=if(TP+FP>0)TP/(TP+FP)else NA_real_,accuracy=(TP+TN)/sum(w),
    mean_binomial_deviance=2*unname(scores(y,z$probability,w)["log_loss"]),
    always_predict_no_accuracy=weighted.mean(y==0,w),
    raw_TP=sum(positive&y==1),raw_FN=sum(!positive&y==1),
    raw_TN=sum(!positive&y==0),raw_FP=sum(positive&y==0))
  stopifnot(result$raw_TP+result$raw_FN+result$raw_TN+result$raw_FP==nrow(z))
  result
}))
save_csv(classification,"LASSO_heldout_classification_at_half.csv")
show_table(classification[c("method","threshold","sensitivity","specificity","accuracy","mean_binomial_deviance")],4)

# ---- LASSO-residual-checks ----
j<-match(full_tune$choices["minimum"],LAMBDAS)
pfit<-full_path$p[,j]
bounded_p<-pmin(pmax(pfit,1e-12),1-1e-12)
diagnostics<-data.frame(SEQN=g$SEQN,y=g$y,fitted_probability=pfit,
  raw_residual=g$y-pfit,
  Pearson_residual=(g$y-bounded_p)/sqrt(bounded_p*(1-bounded_p)),
  deviance_residual=sign(g$y-bounded_p)*sqrt(pmax(-2*(g$y*log(bounded_p)+(1-g$y)*log1p(-bounded_p)),0)))
stopifnot(nrow(diagnostics)==8038,!anyDuplicated(diagnostics$SEQN),
  all(is.finite(diagnostics$raw_residual)),all(is.finite(diagnostics$Pearson_residual)),
  all(is.finite(diagnostics$deviance_residual)),
  abs(weighted.mean(diagnostics$deviance_residual^2,g$WTMECPRP)-
    2*unname(scores(g$y,pfit,g$WTMECPRP)["log_loss"]))<1e-10)
save_csv(diagnostics,"LASSO_residual_diagnostics.csv")
draw_residuals<-function() {
  op<-par(mfrow=c(1,3),mar=c(4.5,4.5,3,1));on.exit(par(op))
  labels<-c(raw_residual="Raw residual: y - p",Pearson_residual="Pearson residual",deviance_residual="Signed deviance residual")
  titles<-c(raw_residual="Raw logistic residuals",Pearson_residual="Pearson logistic residuals",deviance_residual="Binomial deviance residuals")
  for(v in names(labels)) {
    plot(diagnostics$fitted_probability,diagnostics[[v]],pch=16,cex=.3,
      col=adjustcolor("#416b91",alpha.f=.2),xlab="Fitted probability",
      ylab=labels[v],main=titles[v])
    abline(h=0,lty=2,col="grey50")
  }
}
png(file.path(out,"LASSO_residual_diagnostics.png"),width=1800,height=675,res=150)
draw_residuals();dev.off()
if(interactive())draw_residuals()
verification$residual_checks<-list(status="passed",n=nrow(diagnostics),
  scope="Descriptive raw/Pearson/deviance residuals for the selected penalized binomial fit; no Gaussian normality test or post-selection significance claim.")

# ---- figures ----
draw_figures<-function() {
  op<-par(mfrow=c(1,3),mar=c(4.5,4.5,3,1));on.exit(par(op))
  plot(g$y+.12*sin(g$SEQN*1.618),g$BMXBMI,pch=16,cex=.3,
    col=adjustcolor("#416b91",alpha.f=.2),xaxt="n",xlab="Observed diagnosed diabetes",
    ylab=expression(BMI~(kg/m^2)),main="All 8,038 complete cases")
  axis(1,at=0:1,labels=c("No","Yes"))
  tab<-full_tune$table;lx<-log10(tab$lambda)
  plot(lx,tab$mean_log_loss,type="n",ylim=range(c(tab$mean_log_loss+tab$SE,tab$mean_log_loss-tab$SE,full_tune$fold_loss)),
    xlab=expression(log[10](lambda)),ylab="Weighted CV log loss",main="Full-data penalty tuning")
  for(k in 1:5)points(lx,full_tune$fold_loss[k,],pch=16,cex=.35,col=adjustcolor("grey40",alpha.f=.35))
  arrows(lx,tab$mean_log_loss-tab$SE,lx,tab$mean_log_loss+tab$SE,angle=90,code=3,length=.025,col="grey55")
  lines(lx,tab$mean_log_loss,type="o",pch=16,cex=.5,col="#416b91")
  abline(v=log10(full_tune$choices),lty=c(2,3),col=c("#416b91","#af7948"))
  legend("topleft",legend=c("Minimum","One SE"),lty=c(2,3),col=c("#416b91","#af7948"),bty="n",cex=.7)
  plot(0:1,0:1,type="n",xlab="False positive rate",ylab="True positive rate",main="Held-out LASSO ROC")
  abline(0,1,lty=2,col="grey65")
  colors<-c(LASSO_minimum="#416b91",LASSO_one_se="#af7948")
  for(method in names(colors)) {
    z<-predictions[predictions$method==method,];r<-roc_points(z$y,z$probability,z$weight)
    lines(r$FPR,r$TPR,col=colors[method],lwd=1.5)
    points(r$FPR,r$TPR,pch=16,cex=.15,col=adjustcolor(colors[method],alpha.f=.2))
  }
  legend("bottomright",legend=c("Minimum","One SE"),lty=1,col=colors,bty="n",cex=.7)
}
png(file.path(out,"Cleaning_and_LASSO.png"),width=1800,height=675,res=150)
draw_figures();dev.off()
if(interactive())draw_figures()

# ---- verification ----
stopifnot(identical(unname(tools::md5sum(paths)),hashes$md5))
verification$raw_files_unchanged<-TRUE
verification$lecture_review<-list(date="2026-11-12",status="passed",
  sources=c("Lecture18_model_comparison.html","Lecture19_CV.html","Lecture20_modelselection.pdf"),
  cleaning_rules_changed=FALSE,primary_lambda_grid_changed=FALSE,
  additions=c("Included/excluded summaries","Existing outer-fold factor/contrast retention",
    "Held-out classification and mean binomial deviance","Signed deviance residuals"),
  additional_random_partitions_executed=TRUE)
verification$repeated_nested_CV<-repeat_verification
verification$scope<-"Data cleaning and weighted LASSO logistic regression only."
verification$R_version<-R.version.string
verification$package_versions<-as.list(vapply(required,function(x)as.character(packageVersion(x)),character(1)))
verification$script_md5<-if(length(script_file))unname(tools::md5sum(script_file))else NA_character_
verification$status<-"passed"
jsonlite::write_json(verification,file.path(out,"verification.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
show_table(data.frame(check=c("Cleaned inputs versus previous","Nested LASSO validation",
  "Final lambdas and selected factors",paste(REPEATS,"new nested PSU partitions"),
  "Original raw files"),status=c("passed","passed","passed","passed","unchanged")))
