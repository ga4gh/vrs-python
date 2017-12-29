# Makefile for Python project

.DELETE_ON_ERROR:
.PHONY: FORCE
.PRECIOUS :
.SUFFIXES:

SHELL:=/bin/bash -e -o pipefail
SELF:=$(firstword $(MAKEFILE_LIST))

PKG=vmc
PKGD=$(subst .,/,${PKG})
VEDIR=venv/3.6

export SEQREPO_ROOT_DIR=tests/_data/seqrepo


############################################################################
#= BASIC USAGE
default: help

#=> help -- display this help message
help:
	@sbin/makefile-extract-documentation "${SELF}"


############################################################################
#= SETUP, INSTALLATION, PACKAGING

#=> venv: make a Python virtual environment
.PHONY: venv/2.7
venv/2.7: venv/%:
	virtualenv -p $$(type -p python$*) $@; \
	source $@/bin/activate; \
	pip install --upgrade pip setuptools

#=> venv: make a Python 3 virtual environment
.PHONY: venv/3.5 venv/3.6
venv/3.5 venv/3.6: venv/%:
	pyvenv-$* $@; \
	source $@/bin/activate; \
	python -m ensurepip --upgrade; \
	pip install --upgrade pip setuptools

#=> setup: setup/upgrade packages *in current environment*
.PHONY: setup
setup: etc/develop.reqs etc/install.reqs
	if [ -s $(word 1,$^) ]; then pip install --upgrade -r $(word 1,$^); fi
	if [ -s $(word 2,$^) ]; then pip install --upgrade -r $(word 2,$^); fi

#=> devready: create venv, install prerequisites, install pkg in develop mode
.PHONY: devready
devready:
	make ${VEDIR} && source ${VEDIR}/bin/activate && make setup develop
	@echo '#################################################################################'
	@echo '###  Do not forget to `source ${VEDIR}/bin/activate` to use this environment  ###'
	@echo '#################################################################################'

#=> develop: install package in develop mode
.PHONY: develop
develop:
	pip install -e .

#=> install: install package
#=> bdist bdist_egg bdist_wheel build sdist: distribution options
.PHONY: bdist bdist_egg bdist_wheel build build_sphinx sdist install
bdist bdist_egg bdist_wheel build sdist install: %:
	python setup.py $@

#=> upload: upload to pypi
#=> upload_*: upload to named pypi service (requires config in ~/.pypirc)
.PHONY: upload upload_%
upload: upload_pypi
upload_%:
	python setup.py bdist_egg bdist_wheel sdist upload -r $*


############################################################################
#= TESTING
# see test configuration in setup.cfg

#=> test: execute tests
.PHONY: test
test:
	type -p python
	python setup.py pytest --addopts="--cov=${PKG} ${PKGD} tests"

#=> tox: execute tests via tox
.PHONY: tox
tox:
	tox


############################################################################
#= UTILITY TARGETS

# N.B. Although code is stored in github, I use hg and hg-git on the command line
#=> reformat: reformat code with yapf and commit
.PHONY: reformat
reformat:
	@if ! git diff --cached --exit-code; then echo "Repository not clean" 1>&2; exit 1; fi
	yapf -i -r "${PKGD}" tests
	git commit -a -m "reformatted with yapf"

#=> docs -- make sphinx docs
.PHONY: doc docs
doc docs: develop
	# RTD makes json. Build here to ensure that it works.
	make -C doc html json

############################################################################
#= CLEANUP

#=> clean: remove temporary and backup files
.PHONY: clean
clean:
	find . \( -name \*~ -o -name \*.bak \) -print0 | xargs -0r rm

#=> cleaner: remove files and directories that are easily rebuilt
.PHONY: cleaner
cleaner: clean
	rm -f devready.log
	rm -fr .cache *.egg-info build dist doc/_build htmlcov
	find . \( -name \*.pyc -o -name \*.orig -o -name \*.rej \) -print0 | xargs -0r rm
	find . -name __pycache__ -print0 | xargs -0r rm -fr

#=> cleanest: remove files and directories that require more time/network fetches to rebuild
.PHONY: cleanest distclean
cleanest distclean: cleaner
	rm -fr .eggs .tox venv


## <LICENSE>
## Copyright 2016 Source Code Committers
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
