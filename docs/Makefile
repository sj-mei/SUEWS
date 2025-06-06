# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = SUEWS
SOURCEDIR     = source
BUILDDIR      = build
SUPYDIR      = ./source/related-softwares/supy

# Variables for test builds
TESTSRC ?= testsource
TESTMASTER ?= input_files/yaml_input

LIVEHTML_OPTS  = \
	-j auto \
	--re-ignore 'supy\.util\..*' \
	--ignore "source/SuPy.log" \
	--ignore "source/yaml_input/*.rst"

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile pip testyaml config-ui schema

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


pip:
	pip install -e ../src/supy

proc-csv:
	python $(SUPYDIR)/proc_var_info/gen_rst.py


livehtml:
	sphinx-autobuild $(LIVEHTML_OPTS) \
		-b html \
		-a "$(SOURCEDIR)" "$(BUILDDIR)"\
		$(SPHINXOPTS)

schema:
	@echo "Generating JSON schema from SUEWSConfig Pydantic model..."
	python gen_schema.py

config-ui:
	@echo "Starting SUEWS Configuration UI server..."
	@echo "This will serve the config UI at http://localhost:8080"
	@echo "Press Ctrl+C to stop the server"
	@echo ""
	cd source/_static && python run_server.py