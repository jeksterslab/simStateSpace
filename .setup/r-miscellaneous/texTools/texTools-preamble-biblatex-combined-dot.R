#' Create Biblatex Preamble
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @returns Returns a character string.
#'
#' @param path Character string.
#'   Path to `*.bib` files.
#' @param output_path Character string.
#'   Output path.
#' @param style Character string.
#'   `"apa"` or `"reading"`.
#' @param sortcites Logical.
#'   Sort citations.
#' @param sorting Character string.
#'   `"nty"` sort by name, title, year;
#'   `"nyt"` sort by name, year, title;
#'   `"nyvt"` sort by name, year, volume, title;
#'   `"anyt "` sort by alphabetic label, name, year, title;
#'   `"anyvt"` sort by alphabetic label, name, year, volume, title;
#'   `"ydnt"` sort by year (descending), name, title; and
#'   `"none"` entries are processed in citation order.
#' @param map Logical.
#'   Null mapping for `addendum`, `note`, and `annotation`.
#' @param fn_bib Character string.
#'   Biblatex filename.
#'
#' @family Bibliography Functions
#' @keywords texTools biblatex internal
#' @noRd
.PreambleBiblatexCombined <- function(path,
                                      output_path,
                                      style = "apa",
                                      sortcites = TRUE,
                                      sorting = "nyt",
                                      map = TRUE,
                                      fn_bib = "bib.bib") {
  if (dir.exists(path)) {
    bibs <- list.files(
      path = path,
      pattern = "\\.bib$",
      full.names = TRUE,
      all.files = TRUE
    )
    if (length(bibs) > 0) {
      if (sortcites) {
        sortcites <- "true"
      } else {
        sortcites <- "false"
      }
      bibs <- unlist(
        lapply(
          X = bibs,
          FUN = readLines
        )
      )
      dir.create(
        path = output_path,
        showWarnings = FALSE
      )
      output_file <- file.path(
        output_path,
        fn_bib
      )
      con <- file(output_file)
      writeLines(
        text = bibs,
        con = con
      )
      close(con)
      biblatex <- paste0(
        "\n",
        "\\usepackage[",
        "style=",
        style,
        ",",
        "sortcites=",
        sortcites,
        ",",
        "sorting=",
        sorting,
        ",",
        "backend=biber",
        ",",
        "labeldate=year",
        "]{biblatex}",
        "\n",
        "\\DeclareLanguageMapping{american}{american-apa}",
        "\n",
        collapse = ""
      )
      if (map) {
        map <- paste0(
          "\n",
          "\\DeclareSourcemap{",
          "\n",
          "\\maps[datatype = bibtex]{",
          "\n",
          "\\map{",
          "\n",
          "\\step[fieldset = addendum, null]",
          "\n",
          "\\step[fieldset = note, null]",
          "\n",
          "\\step[fieldset = annotation, null]",
          "\n",
          "}",
          "\n",
          "}",
          "\n",
          "}",
          collapse = ""
        )
      } else {
        map <- ""
      }
      return(
        paste0(
          biblatex,
          paste0(
            "\n",
            "\\addbibresource{",
            output_file,
            "}",
            "\n",
            collapse = ""
          ),
          map,
          collapse = ""
        )
      )
    } else {
      return("")
    }
  } else {
    return("")
  }
}
