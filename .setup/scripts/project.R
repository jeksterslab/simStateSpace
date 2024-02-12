#' Define project related objects here.
#'
#' | Object Name      | Description                                         |
#' |:-----------------|:----------------------------------------------------|
#' | `project`        | Project name.                                       |
#' | `pkg_cran`       | CRAN packages to install on top of                  |
#' |                  | package dependencies in `DESCRIPTION`.              |
#' | `pkg_github`     | GitHub packages to install.                         |
#' | `pkg_github_ref` | GitHub branch corresponding to packages             |
#' |                  | in `pkg_github`. `if (length(pkg_github_ref) == 0)` |
#' |                  | use the `HEAD` branch.                              |
#' | `pkg_ver`        | Packages with specific version.                     |
#' | `ver`            | Version corresponding to packages in `pkg_ver`.     |
#' | `ignore`         | Items to add to `.Rbuildignore`.                    |
#' | `license`        | Options are `"mit"`, `"gpl3"`, or `NULL`            |
#'

project <- "simStateSpace"

pkg_cran <- c(
  "Rcpp",
  "RcppArmadillo",
  "expm",
  "lavaan",
  "betaDelta"
)

# dynr has to be installed which requires more work to install

pkg_github <- c()

pkg_github_ref <- c()

pkg_ver <- c()

ver <- c()

ignore <- "^vignettes$"

license <- "gpl3"

git_user <- "jeksterslab"
git_email <- "learn.jeksterslab@gmail.com"
r_email <- "r.jeksterslab@gmail.com"
