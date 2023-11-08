#' validate email
Validate <- function() {
  root <- rprojroot::is_rstudio_project
  path <- dirname(
    root$find_file(
      "project.Rproj"
    )
  )
  source(
    file.path(
      path,
      ".setup",
      "scripts",
      "project.R"
    )
  )
  rhub::validate_email(email = r_email)
}
