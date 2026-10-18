# simStateSpace

See `AGENTS.md` and `.setup/ai/project.md` for this project (package).

- R: BigCamelCase functions; snake_case objects; `pkg::fn()` with named
  arguments; no explicit
  [`return()`](https://rdrr.io/r/base/function.html); one function per
  file.
- Never edit generated files (`NAMESPACE`, `man/`, `README.md`,
  `.setup/latex/pdf/`).
- Run code only inside the image the user chooses:
  `apptainer exec <image> make <target>`.
