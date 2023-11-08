#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
path <- as.character(args[1])
source(
  file.path(
    path,
    ".setup",
    "scripts",
    "project.R"
  )
)
dot_library_folder <- file.path(
  path,
  ".library"
)
.libPaths(
  dot_library_folder
)
Sys.setenv(
  R_LIBS_USER = dot_library_folder
)
rProject::PkgDevel(
  path = path
)
rProject::PkgDevelGitHub(
  path = path
)
rProject::PkgProject(
  path = path,
  pkg = pkg_cran
)
if (length(pkg_github_ref) == 0) {
  rProject::PkgProjectGitHub(
    path = path,
    pkg = pkg_github
  )
} else {
  rProject::PkgProjectGitHub(
    path = path,
    pkg = pkg_github,
    ref = pkg_github_ref
  )
}
rProject::PkgProjectVersion(
  path = path,
  pkg = pkg_ver,
  ver = ver
)
warnings()
