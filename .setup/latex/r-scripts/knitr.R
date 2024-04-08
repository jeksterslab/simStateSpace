# find root directory
rproj <- rprojroot::is_rstudio_project
source(
  rproj$find_file(
    ".setup",
    "latex",
    "r-scripts",
    "knitr-options-default.R"
  )
)
source(
  rproj$find_file(
    ".setup",
    "latex",
    "r-scripts",
    "knitr-options-custom.R"
  )
)
knitr::opts_chunk$set(
  modifyList(
    KnitrOptionsDefault(rproj),
    KnitrOptionsCustom()
  )
)
knitr::knit_hooks$set(
  document = function(x) {
    sub(
      "\\usepackage[]{color}",
      "\\usepackage{xcolor}",
      x,
      fixed = TRUE
    )
  }
)
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
x <- c(
  x,
  list.files(
    path = rproj$find_file(
      "R"
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
rm(
  KnitrOptionsDefault,
  KnitrOptionsCustom,
  rproj
)
