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
rProject::ConfigFiles(
  git_user = git_user
)
warnings()
