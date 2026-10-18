#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
path <- as.character(args[1])
source(
  file.path(
    path,
    ".setup",
    "scripts",
    "project.R"
  )
)
# Secrets are read from the local machine, never from the project.
# The first file that defines `tokens` wins:
#   1. Sys.getenv("RPROJECT_SECRETS")
#   2. ~/.config/rproject/<project>.R
#   3. ~/.config/rproject/secrets.R
#   4. .setup/scripts/secrets.R (legacy; prints a warning)
secrets_legacy <- file.path(
  path,
  ".setup",
  "scripts",
  "secrets.R"
)
secrets_candidates <- c(
  Sys.getenv("RPROJECT_SECRETS"),
  file.path(
    "~",
    ".config",
    "rproject",
    paste0(project, ".R")
  ),
  file.path(
    "~",
    ".config",
    "rproject",
    "secrets.R"
  ),
  secrets_legacy
)
secrets_candidates <- secrets_candidates[
  nzchar(secrets_candidates) & file.exists(secrets_candidates)
]
tokens <- NULL
secrets_path <- NULL
for (secrets_candidate in secrets_candidates) {
  secrets_env <- new.env(parent = globalenv())
  sys.source(
    file = secrets_candidate,
    envir = secrets_env
  )
  tokens <- get0(
    x = "tokens",
    envir = secrets_env,
    inherits = FALSE
  )
  if (length(tokens) > 0) {
    secrets_path <- secrets_candidate
    break
  }
}
rm(list = ls(pattern = "^secrets_(env|candidate|candidates)$"))
if (!is.null(secrets_path)) {
  if (secrets_path == secrets_legacy) {
    warning(
      paste0(
        "\n\nSecrets were read from ",
        secrets_legacy,
        ".\nMove them to ~/.config/rproject/secrets.R (chmod 600)",
        " and remove the tokens from Overleaf.\n\n"
      ),
      call. = FALSE
    )
  } else {
    secrets_mode <- file.info(secrets_path)$mode
    if (bitwAnd(as.integer(secrets_mode), strtoi("077", base = 8L)) != 0L) {
      Sys.chmod(
        paths = secrets_path,
        mode = "600"
      )
      message(
        paste0(
          "\n\nPermissions of ",
          secrets_path,
          " changed from ",
          format(secrets_mode),
          " to 600.\n\n"
        )
      )
    }
  }
  message(
    paste0(
      "\n\nSecrets read from ",
      secrets_path,
      ": ",
      paste(names(tokens), collapse = ", "),
      "\n\n"
    )
  )
  keys <- names(tokens)
  env <- rep(x = "", times = length(keys))
  for (i in seq_along(tokens)) {
    env[i] <- paste0(
      "export",
      " ",
      keys[i],
      "=",
      tokens[i]
    )
  }
  writeLines(env, "~/.bash-secrets")
  Sys.setenv(
    GITHUB_PAT = tokens["GITHUB_PAT"]
  )
} else {
  message(
    paste0(
      "\n\n",
      "No secrets found.",
      " Put tokens in ~/.config/rproject/secrets.R",
      " (not needed in GitHub Actions).",
      "\n\n"
    )
  )
  tokens <- NULL
}
dot_library_folder <- file.path(
  path,
  ".library"
)
dir.create(
  dot_library_folder,
  showWarnings = FALSE,
  recursive = TRUE
)
dot_gitignore_file <- file.path(
  dot_library_folder,
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
.libPaths(
  dot_library_folder
)
Sys.setenv(
  R_LIBS_USER = dot_library_folder
)
installed <- installed.packages()
pkg_installed <- installed[, "Package"]
if (!("pak" %in% pkg_installed)) {
  install.packages(
    "pak",
    repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/latest"),
    lib = dot_library_folder,
    quiet = TRUE
  )
}
# install rProject alone into .library: pak looks only in `lib` for what is
# installed, so with dependencies it would rebuild from source every package
# the image already has (spec 2026-10-13-make-project-pak); the dependencies
# come from the image, which .libPaths() keeps after .library. In GitHub
# Actions, make cleanpkg has removed the image's packages (the container runs
# as root), so there pak installs the dependencies too
rproject_deps <- identical(Sys.getenv("GITHUB_ACTIONS"), "true")
rproject_install <- function(ref) {
  message(
    paste0(
      "\n\nInstalling rProject ",
      ref,
      " into .library",
      ifelse(
        rproject_deps,
        " with its dependencies (GitHub Actions).\n\n",
        " (dependencies from the image).\n\n"
      )
    )
  )
  rproject_missing <- paste0(
    "\n\nIf a dependency is missing, add it to the image",
    " (docker-docs or the project's .sif);",
    " make project does not install dependencies into .library.\n\n"
  )
  tryCatch(
    expr = pak::pkg_install(
      pkg = ref,
      lib = dot_library_folder,
      dependencies = if (rproject_deps) NA else FALSE
    ),
    error = function(e) {
      stop(
        paste0(
          "\n\nrProject could not be installed in .library:\n",
          conditionMessage(e),
          rproject_missing
        ),
        call. = FALSE
      )
    }
  )
  # load it in a fresh R process; a missing dependency belongs in the image
  rproject_load <- suppressWarnings(
    system2(
      command = file.path(R.home("bin"), "Rscript"),
      args = c(
        "--vanilla",
        "-e",
        shQuote(
          paste0(
            ".libPaths(c(\"",
            dot_library_folder,
            "\", .libPaths())); loadNamespace(\"rProject\")"
          )
        )
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  if (!is.null(attr(rproject_load, "status"))) {
    stop(
      paste0(
        "\n\nrProject was installed in .library but does not load:\n",
        paste(rproject_load, collapse = "\n"),
        rproject_missing
      ),
      call. = FALSE
    )
  }
}
# every project, rProject included: install rProject when it is missing,
# and upgrade it when the installed copy is older than the GitHub HEAD
# (or than rproject_ver)
if (is.null(rproject_ver)) {
  # pak does not upgrade an installed GitHub package,
  # so install the current commit when the project's copy is older
  rproject_head <- suppressWarnings(
    system2(
      command = "git",
      args = c(
        "ls-remote",
        "https://github.com/ijapesigan/rProject.git",
        "HEAD"
      ),
      stdout = TRUE,
      stderr = FALSE,
      env = c(
        "GIT_CONFIG_GLOBAL=/dev/null",
        "GIT_TERMINAL_PROMPT=0"
      )
    )
  )
  rproject_head <- substr(
    x = c(rproject_head, NA_character_)[1],
    start = 1,
    stop = 40
  )
  rproject_sha <- suppressWarnings(
    utils::packageDescription(
      pkg = "rProject",
      lib.loc = dot_library_folder,
      fields = "RemoteSha"
    )
  )
  if (is.na(rproject_head)) {
    if (is.na(rproject_sha)) {
      rproject_install(ref = "ijapesigan/rProject")
    } else {
      message(
        "\n\nCould not reach GitHub; keeping the installed rProject.\n\n"
      )
    }
  } else if (!identical(rproject_sha, rproject_head)) {
    rproject_install(
      ref = paste0(
        "ijapesigan/rProject",
        "@",
        rproject_head
      )
    )
  }
} else {
  rproject_install(
    ref = paste0(
      "ijapesigan/rProject",
      "@",
      rproject_ver
    )
  )
}
rProject::Project(
  path = path
)
rProject::Profile(
  path = path,
  project = project
)
rProject::Environment(
  path = path,
  project = project,
  tokens = tokens
)
rProject::BuildIgnore(
  path = path,
  project = project,
  add = ignore
)
rProject::Binary(
  path = path
)
rProject::License(
  path = path,
  type = license
)
# per-project AI instructions (rProject 0.0.34 or later)
if (exists("AI", envir = asNamespace("rProject"))) {
  rProject::AI(
    path = path
  )
}
warnings()
