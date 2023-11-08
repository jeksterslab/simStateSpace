#' Create `detritus/tex/tex-fig-*`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams LibPaths
#' @export
Detritus <- function(path) {
  detritus_folder <- file.path(
    path,
    "detritus"
  )
  tex_folder <- file.path(
    detritus_folder,
    "tex"
  )
  tex_figures_folder <- file.path(
    tex_folder,
    paste0(
      "tex-fig-",
      paste0(
        sample(
          letters,
          size = 8,
          replace = TRUE
        ),
        collapse = ""
      )
    )
  )
  dir.create(
    tex_figures_folder,
    showWarnings = FALSE,
    recursive = TRUE
  )
  dot_gitignore_file <- file.path(
    detritus_folder,
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
  return(
    list(
      detritus_folder = detritus_folder,
      tex_folder = tex_folder,
      tex_figures_folder = tex_figures_folder
    )
  )
}
