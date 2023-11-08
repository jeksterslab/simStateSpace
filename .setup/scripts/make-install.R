#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
path <- as.character(args[1])
if (
  file.exists(
    file.path(
      path,
      "DESCRIPTION"
    )
  )
) {
  deps <- list.files(
    path = file.path(
      path,
      ".setup",
      "r-dependencies"
    ),
    pattern = ".*\\.R",
    full.names = TRUE, recursive = TRUE
  )
  if (length(deps) > 0) {
    file.copy(
      from = deps,
      to = file.path(path, "R")
    )
  }
  cpp <- list.files(
    path = file.path(path, ".setup", "cpp"),
    pattern = "^.*\\.cpp",
    full.names = TRUE
  )
  if (length(cpp) > 0) {
    dir.create(
      path = file.path(path, "src"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    unlink(
      file.path(
        path,
        "R",
        "RcppExports.R"
      )
    )
    unlink(
      file.path(
        path,
        "src",
        "RcppExports.cpp"
      )
    )
    namespace <- file.path(
      path,
      "NAMESPACE"
    )
    if (!file.exists(namespace)) {
      file.create(namespace)
    }
    cpp_file_fn <- file.path(
      path,
      "src",
      "source.cpp"
    )
    cpp_file <- file(
      description = cpp_file_fn,
      open = "w"
    )
    for (i in seq_along(cpp)) {
      writeLines(
        text = readLines(con = cpp[i]),
        con = cpp_file
      )
    }
    file.copy(
      from = file.path(path, ".setup", "cpp", "Makevars"),
      to = file.path(path, "src")
    )
    file.copy(
      from = file.path(path, ".setup", "cpp", "Makevars.win"),
      to = file.path(path, "src")
    )
    close(cpp_file)
    Rcpp::compileAttributes(pkgdir = path)
    roxygen2::roxygenize(
      package.dir = path,
      roclets = "rd"
    )
    unlink(namespace)
  }
  devtools::document(path)
  devtools::install(path, dependencies = FALSE)
}
warnings()
