# Makefile for Python project

.DELETE_ON_ERROR:
.PHONY: FORCE
.PRECIOUS:
.SUFFIXES:

SHELL:=/bin/bash -e -o pipefail
SELF:=$(firstword $(MAKEFILE_LIST))

UNAME = $(shell uname)
ifeq (${UNAME},Darwin)
    _XRM_R:=
else
    _XRM_R:=r
endif
XRM=xargs -0${_XRM_R} rm

PKG=ga4gh.vrs
PKGD=$(subst .,/,${PKG})
PYV:=3.12
VEDIR=venv/${PYV}

export SEQREPO_ROOT_DIR=tests/data/seqrepo/latest


############################################################################
#= BASIC USAGE
default: help

#=> help -- display this help message
help:
	@sbin/makefile-extract-documentation "${SELF}"


############################################################################
#= SETUP, INSTALLATION, PACKAGING

#=> venv: make a Python 3 virtual environment
.PHONY: venv/%
venv/%:
	python$* -m venv $@; \
	source $@/bin/activate; \
	python -m ensurepip --upgrade; \
	pip install --upgrade pip setuptools

#=> develop: install package in develop mode
.PHONY: develop setup
develop setup:
	pip install -e .[dev,extras,notebooks]

#=> devready: create venv, install prerequisites, install pkg in develop mode
.PHONY: devready
devready:
	make ${VEDIR} && source ${VEDIR}/bin/activate && make develop
	@echo '#################################################################################'
	@echo '###  Do not forget to `source ${VEDIR}/bin/activate` to use this environment  ###'
	@echo '#################################################################################'

#=> nbready: create venv, install prerequisites, install pkg in notebook mode
.PHONY: nbready
nbready:
	make ${VEDIR} && source ${VEDIR}/bin/activate && pip install -e '.[extras,notebooks]'
	@echo '#################################################################################'
	@echo '###  Do not forget to `source ${VEDIR}/bin/activate` to use this environment  ###'
	@echo '#################################################################################'

#=> install: install package
.PHONY: install-extras
install-extras:
	pip install -e .[extras]


############################################################################
#= TESTING
# see test configuration in pyproject.toml

#=> test: execute tests
.PHONY: test
test:
	pytest

#=> doctest: execute documentation tests (requires extra data)
.PHONY: doctest
doctest:
	pytest --doctest-modules

############################################################################
#= UTILITY TARGETS

# N.B. Although code is stored in github, I use hg and hg-git on the command line
#=> reformat: reformat code with yapf and commit
.PHONY: reformat
reformat:
	@if ! git diff --cached --exit-code; then echo "Repository not clean" 1>&2; exit 1; fi
	yapf -i -r "${PKGD}" tests
	git commit -a -m "reformatted with yapf"

#=> lint -- static analysis check
.PHONY: lint
lint:
	pylint src/ga4gh/{core,vrs} | tee $@

#=> docs -- make sphinx docs
.PHONY: docs
docs: develop
	# RTD makes json. Build here to ensure that it works.
	make -C doc html json

############################################################################
#= CLEANUP

#=> clean: remove temporary and backup files
.PHONY: clean
clean:
	find . \( -name \*~ -o -name \*.bak \) -print0 | ${XRM}

#=> cleaner: remove files and directories that are easily rebuilt
.PHONY: cleaner
cleaner: clean
	rm -fr .cache *.egg-info .pytest_cache build dist doc/_build htmlcov
	find . \( -name \*.pyc -o -name \*.orig -o -name \*.rej \) -print0 | ${XRM}
	find . -name __pycache__ -print0 | ${XRM} -fr

#=> cleanest: remove files and directories that require more time/network fetches to rebuild
.PHONY: cleanest
cleanest: cleaner
	rm -fr .eggs venv


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
