.PHONY: install docs
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
else
SEDI=sed -i
endif

venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --python=python3
	${IN_VENV} && pip install pip --upgrade

minimap2/libminimap2.a:
	${SEDI} 's/int\ mm_verbose\ =\ 3;/int\ mm_verbose\ =\ 2;/' minimap2/misc.c
	${SEDI} 's/CFLAGS=.*/CFLAGS=-g\ -Wc++-compat -Wall\ -Wno-unused-function\ -O2\ -fPIC/' minimap2/Makefile
	cd minimap2 && make libminimap2.a 

install: venv minimap2/libminimap2.a
	${IN_VENV} && pip install -r requirements.txt && python setup.py install


# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

DOCSRC = docs

docs: install # we install the package to ensure imports
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll
