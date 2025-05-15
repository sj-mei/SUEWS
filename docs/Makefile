# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = SUEWS
SOURCEDIR     = source
BUILDDIR      = build
SUPYDIR      = ./source/related-softwares/supy



LIVEHTML_OPTS  = \
  --re-ignore 'supy\.util\..*' \
  --ignore "source/SuPy.log" \
  --ignore "source/yaml_input/*.rst" \
  --fresh-env

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile pip

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