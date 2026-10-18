# simStateSpace

This project follows the house rules in your global instructions
(`~/.claude/CLAUDE.md`, `~/.codex/AGENTS.md`,
`~/.copilot/copilot-instructions.md`). This file adds only what is
specific to this project.

**Read `.setup/ai/project.md` first.** It holds the project facts that
the user maintains (aims, data, notation, co-authors, open issues).

## Project

- Name: `simStateSpace`
- Type: package
- Private repository: `jeksterslab-psu/simStateSpace` (Overleaf pushes
  here)
- Overleaf project: not recorded (`.setup/repo/private-overleaf.md`)
- Public repository: `jeksterslab/simStateSpace` (`make auto` pushes
  here)
- These come from `.setup/repo/` (private). O2G checks the private
  repository; PUB checks the public one.
- Image: **the user chooses it** and runs
  `apptainer exec <image> make <target>`; R is not installed on the
  host. Do not pick an image. When a command needs one, ask, or give the
  command with `<image>` for the user to fill in. `.setup/ai/project.md`
  may note the image the project usually uses.

## Folder map

| Folder | Contents |
|----|----|
| `R/` | package functions, one per file |
| `tests/testthat/` | tests |
| `.setup/scripts/` | `project.R` (project settings) and the make scripts |
| `.setup/data-raw/` | shareable raw data |
| `.setup/data-private/` | restricted data: never public |
| `.setup/ai/` | specs, the plan, `project.md` (private) |
| `.setup/logs/` | one log per Overleaf push, and `index.md` (private) |
| `.setup/repo/` | where the project lives (private) |

## Make targets

| Target | Does |
|----|----|
| `make all` | build, docs, LaTeX, pkgdown, quarto, coverage |
| `make lint` | style and lint the R code, and report Overleaf limits |
| `make sizecheck` | check Overleaf limits; run before uploading |
| `make spec TOPIC=<topic>` | start a spec for a change set |
| `make log TOPIC=<topic> SUMMARY="..."` | write the change-set log and its file list (`UPDATE=TRUE` refreshes the list) |
| `make verify RECORD=TRUE` | after the user pushes from Overleaf (O2G): verify and record the push. On Roar, run `bash .setup/scripts/verify.sh` on the host first to pull. |
| `make auto` | publish to the public repository (checks the export first) |

## Generated files: never edit by hand

`NAMESPACE`, `man/`, `README.md`, `CITATION.cff`, `LICENSE*`,
`.Rbuildignore`, `project.Rproj`, `R/RcppExports.R`,
`src/RcppExports.cpp`, `*/bib.bib`, `vignettes/*.Rmd` built from
`.setup/vignettes/`, `.setup/latex/pdf/`, `.setup/latex/figures/`,
`data/*.rda`. Change their sources and run `make`.

## Size

Overleaf’s GitHub sync stops above 7 MB of editable text and caps text
files at 2 MB. Keep large data compressed (binary, e.g. `.Rds` or `.xz`)
and never store the same file twice. `make sizecheck` reports problems.

## Package

- Every exported function has roxygen documentation and tests.
- `make all` runs the checks; fix every error, warning, and note.
- Keep the R style rules from your global instructions (BigCamelCase,
  one function per file, internal `.Foo` functions in
  `R/<pkg>-<name>-dot.R`).
