#' @param rproj Object.
#'   Output of `rprojroot::is_rstudio_project`
KnitrOptionsDefault <- function(rproj) {
  path <- dirname(
    rproj$find_file(
      "project.Rproj"
    )
  )
  source(
    rproj$find_file(
      ".setup",
      "latex",
      "r-scripts",
      "detritus.R"
    )
  )
  return(
    list(
      collapse = TRUE,
      comment = "#>",
      fig.width = 6,
      fig.height = 6,
      fig.path = Detritus(path = path)$tex_figures_folder,
      out.width = "100%",
      dev = "pdf",
      dev.args = list(
        pdf = list(
          family = "Times"
        )
      )
    )
  )
}
