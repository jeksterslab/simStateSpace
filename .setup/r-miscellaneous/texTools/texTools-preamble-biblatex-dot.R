#' Create Biblatex Preamble
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @returns Returns a character string.
#'
#' @param path Character string.
#'   Path to `*.bib` files.
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
#'
#' @family Bibliography Functions
#' @keywords texTools biblatex internal
#' @noRd
.PreambleBiblatex <- function(path,
                              style = "apa",
                              sortcites = TRUE,
                              sorting = "nyt",
                              map = TRUE) {
  if (dir.exists(path)) {
    bibs <- file.path(path, "bib.bib")
    if (file.exists(bibs)) {
      if (sortcites) {
        sortcites <- "true"
      } else {
        sortcites <- "false"
      }
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
            bibs,
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
