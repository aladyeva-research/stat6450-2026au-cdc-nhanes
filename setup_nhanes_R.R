# Install the three packages needed by NHANES_Cleaning_LASSO.R.
# Run from the existing project folder containing data/raw/ and inputs/.
# Packages are installed in that project's .r-library, not system-wide.
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
dir.create(task_library,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(normalizePath(task_library,winslash="/"),.libPaths()))
packages<-c("foreign","glmnet","jsonlite")
missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,lib=task_library,repos="https://cloud.r-project.org",
  type=if(.Platform$OS.type=="windows")"binary"else"source")
stopifnot(all(vapply(packages,requireNamespace,logical(1),quietly=TRUE)))
print(vapply(packages,function(x)as.character(packageVersion(x)),character(1)))
