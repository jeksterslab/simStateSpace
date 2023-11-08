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
secrets_path <- file.path(
  path,
  ".setup",
  "scripts",
  "secrets.R"
)
if (file.exists(secrets_path)) {
  source(secrets_path)
  keys <- names(tokens)
  env <- rep(x = "", times = length(keys))
  for (i in seq_along(tokens)) {
    env[i] <- paste0(
      "export",
      " ",
      keys[i],
      "=",
      tokens[i]
    )
  }
  writeLines(env, "~/.bash-secrets")
  Sys.setenv(
    GITHUB_PAT = tokens["GITHUB_PAT"]
  )
} else {
  message(
    paste0(
      "\n\n",
      secrets_path,
      " not found.\n\n"
    )
  )
  tokens <- NULL
}
dot_library_folder <- file.path(
  path,
  ".library"
)
dir.create(
  dot_library_folder,
  showWarnings = FALSE,
  recursive = TRUE
)
dot_gitignore_file <- file.path(
  dot_library_folder,
  ".gitignore"
)
if (!file.exists(dot_gitignore_file)) {
  dot_gitignore <- c(
    "*",
    "*/"
  )
  con <- file(dot_gitignore_file)
  writeLines(
    text = dot_gitignore,
    con = con,
    sep = "\n"
  )
  close(con)
}
.libPaths(
  dot_library_folder
)
Sys.setenv(
  R_LIBS_USER = dot_library_folder
)
installed <- installed.packages()
pkg_installed <- installed[, "Package"]
if (!("remotes" %in% pkg_installed)) {
  install.packages(
    "remotes",
    repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/latest"),
    lib = dot_library_folder,
    quiet = TRUE
  )
}
if (project == "rProject") {
  if (!("rProject" %in% pkg_installed)) {
    remotes::install_github(
      "ijapesigan/rProject",
      quiet = TRUE,
      lib = dot_library_folder
    )
  }
} else {
  remotes::install_github(
    "ijapesigan/rProject",
    quiet = TRUE,
    lib = dot_library_folder
  )
}
rProject::Project(
  path = path
)
rProject::Profile(
  path = path,
  project = project
)
rProject::Environment(
  path = path,
  project = project,
  tokens = tokens
)
rProject::BuildIgnore(
  path = path,
  project = project,
  add = ignore
)
rProject::Binary(
  path = path
)
rProject::License(
  path = path,
  type = license
)
warnings()
