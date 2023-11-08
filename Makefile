.PHONY: all build rproject clean cleanall cleanpkg cleanproj cleanpush cleantinytex bib bibrproject clone coverage data dependencies readme docs dotfiles install latex lint local localforce pdf pkg pkgdown project push quarto style tinytex tinytexforce vignettes

push: build docs latex coverage cleanpush

project:
	@echo "\n\nBuilding project...\n\n"
	@Rscript .setup/scripts/make-project.R ${PWD}

pkg:
	@echo "\n\nInstalling packages...\n\n"
	@Rscript .setup/scripts/make-packages.R ${PWD}

dotfiles:
	@echo "\n\nBuilding dotfiles...\n\n"
	@Rscript .setup/scripts/make-config.R ${PWD}

clean:
	@echo "\n\nCleaning...\n\n"
	@Rscript -e "rProject::Clean(\"${PWD}\")"
	@rm -rf "${PWD}/TEMPDIR.*"

cleanpush:
	@echo "\n\nCleaning in prep for push...\n\n"
	@Rscript -e "rProject::Clean(\"${PWD}\", push = TRUE)"

cleanpkg:
	@echo "\n\nCleaning packages...\n\n"
	@Rscript -e "rProject::CleanPkg(\"${PWD}\")"

cleantinytex:
	@echo "\n\nCleaning TinyTex...\n\n"
	@Rscript -e "rProject::CleanTinyTex(\"${PWD}\")"

cleanproj:
	@echo "\n\nCleaning project...\n\n"
	@Rscript -e "rProject::CleanProj(\"${PWD}\")"

bib:
	@echo "\n\nCleaning project...\n\n"
	@Rscript -e "rProject::Bib(\"${PWD}\", bib_lib = FALSE)"

bibrproject:
	@echo "\n\nCleaning project...\n\n"
	@Rscript -e "rProject::Bib(\"${PWD}\", bib_lib = TRUE)"

tinytex: 
	@echo "\n\nTinyTex...\n\n"
	@Rscript -e "rProject::TinyTex(\"${PWD}\", force = FALSE)"

tinytexforce: 
	@echo "\n\nTinyTex...\n\n"
	@Rscript -e "rProject::TinyTex(\"${PWD}\", force = TRUE)"

style:
	@echo "\n\nStyling...\n\n"
	@Rscript -e "rProject::Style(\"${PWD}\")"

lint: style
	@echo "\n\nLinting...\n\n"
	@Rscript -e "rProject::Lint(\"${PWD}\")"

data:
	@echo "\n\nBuilding data...\n\n"
	@Rscript -e "rProject::DataProcess(\"${PWD}\")"

dependencies:
	@echo "\n\nBuilding dependencies...\n\n"
	@Rscript -e "rProject::Dependencies(\"${PWD}\")"

vignettes:
	@echo "\n\nInitial build...\n\n"
	@Rscript -e "rProject::Build(\"${PWD}\")"
	@Rscript -e "rProject::DataAnalysis(\"${PWD}\")"
	@echo "\n\nPrecompiling vignettes...\n\n"
	@Rscript -e "rProject::VignettesPrecompile(\"${PWD}\")"

build: project pkg dotfiles tinytex lint data dependencies bib vignettes
	@echo "\n\nBuilding package...\n\n"
	@Rscript -e "rProject::Build(\"${PWD}\")"

coverage:
	@echo "\n\nCode coverage...\n\n"
	@Rscript -e "rProject::Coverage(\"${PWD}\")"

latex:
	@echo "\n\nCompiling latex...\n\n"
	@Rscript -e "rProject::LatexMake(\"${PWD}\")"

pdf:
	@echo "\n\nCompiling latex...\n\n"
	@Rscript -e "rProject::LatexMake(\"${PWD}\", clean = TRUE)"

pkgdown:
	@echo "\n\nBuilding pkgdown website...\n\n"
	@Rscript -e "rProject::Site(\"${PWD}\")"

quarto:
	@echo "\n\nRendering quarto...\n\n"
	@Rscript -e "rProject::Quarto(\"${PWD}\")"

readme:
	@echo "\n\nBuilding README.md...\n\n"
	@Rscript -e "rProject::ReadMe(\"${PWD}\")"

docs: readme
	@echo "\n\nBuilding manual...\n\n"
	@Rscript -e "rProject::Manual(\"${PWD}\", project = Sys.getenv(\"PROJECT\"))"
	@echo "\n\nBuilding CITATION.cff...\n\n"
	@Rscript -e "rProject::CFF(\"${PWD}\")"

all: build docs latex pkgdown quarto coverage

cleanall: clean cleanpush cleanpkg cleantinytex

install: project pkg data clean
	@echo "\n\nInstalling...\n\n"
	@Rscript .setup/scripts/make-install.R ${PWD}

local: project dotfiles
	@echo "\n\nInstalling local applications...\n\n"
	@Rscript -e "rProject::InstallLocal(all = TRUE)"
	@Rscript .setup/scripts/make-config.R ${PWD}

localforce: project dotfiles
	@echo "\n\nInstalling local applications...\n\n"
	@Rscript -e "rProject::InstallLocal(all = TRUE, force = TRUE)"
	@Rscript .setup/scripts/make-config.R ${PWD}

clone:
	@bash .setup/scripts/clone.sh

auto: clone
	@git add --all
	@git commit -m "Automated build."
	@git push origin main
