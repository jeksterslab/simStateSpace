#' Load R Scripts in `R`, `.setup/r-dependencies`, `.setup/r-miscellaneous`
#'
#' @param rproj Object.
#'   Output of `rprojroot::is_rstudio_project`
Load <- function(rproj) {
  # sources the R scripts in the project's R, .setup/r-dependencies, and
  # .setup/r-miscellaneous folders
  x <- list.files(
    path = rproj$find_file(
      "R"
    ),
    pattern = "\\.R$",
    full.names = TRUE,
    all.files = TRUE,
    recursive = TRUE
  )
  x <- c(
    x,
    list.files(
      path = rproj$find_file(
        ".setup",
        "r-dependencies"
      ),
      pattern = "\\.R$",
      full.names = TRUE,
      all.files = TRUE,
      recursive = TRUE
    )
  )
  x <- c(
    x,
    list.files(
      path = rproj$find_file(
        ".setup",
        "r-miscellaneous"
      ),
      pattern = "\\.R$",
      full.names = TRUE,
      all.files = TRUE,
      recursive = TRUE
    )
  )
  if (length(x) > 0) {
    invisible(
      lapply(
        X = x,
        FUN = source
      )
    )
  }
}
