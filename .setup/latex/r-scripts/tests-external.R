# find root directory
rproj <- rprojroot::is_rstudio_project
# load functions
source(
  rproj$find_file(
    ".setup",
    "load",
    "load.R"
  )
)
Load(rproj)
# load data
x <- list.files(
  path = rproj$find_file(
    "data"
  ),
  pattern = "\\.rda$",
  full.names = TRUE,
  all.files = TRUE,
  recursive = TRUE
)
x <- c(
  x,
  list.files(
    path = rproj$find_file(
      ".setup",
      "data-dependencies"
    ),
    pattern = "\\.rda$",
    full.names = TRUE,
    all.files = TRUE,
    recursive = TRUE
  )
)
if (length(x) > 0) {
  for (i in seq_along(x)) {
    load(x[i])
  }
  rm(i)
}
rm(x)
# run tests_external
lapply(
  X = rproj$find_file(
    ".setup",
    "tests-external"
  ),
  FUN = function(x) {
    x <- list.files(
      path = x,
      pattern = "\\.R$",
      full.names = TRUE,
      all.files = TRUE,
      recursive = TRUE
    )
    if (length(x) > 0) {
      lapply(
        X = x,
        FUN = source
      )
    }
  }
)
rm(rproj)
