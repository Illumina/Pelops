.PHONY: clean clean-build clean-pyc clean-test coverage dist docs test tests help devinstall install lint lint/flake8 lint/black
.DEFAULT_GOAL := help

PYTHON=python3

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := $(PYTHON) -c "$$BROWSER_PYSCRIPT"

VERSION := $(shell cat VERSION)
PACKAGE_NAME := ilmn-pelops
PACKAGE_FILE := dist/$(PACKAGE_NAME)-$(VERSION).tar.gz


help:
	@$(PYTHON) -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	rm -fr Versionfile
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache
	rm -fr .mypy_cache

lint/black: ## check style with black
	black --check ilmn tests

lint: lint/black ## check style

tests: typecheck

fulltest: versionfile ## run full tox tests
	tox run-parallel

typecheck: ## run pytest with additional typechecking feature
	mypy ilmn \
		--strict \
		--config-file mypy.ini \
		--explicit-package-bases
	pytest --verbose --mypy-config-file=mypy.ini tests

unittest: ## run tests quickly with the default Python
	pytest --verbose tests

coverage: ## check code coverage quickly with the default Python
	pytest --verbose \
		--cov=ilmn/pelops \
		--cov=tests \
		--cov-report term-missing \
		tests

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/dux4r_caller.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ dux4r_caller
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

devinstall: versionfile dev-prerequisites ## install the package in developer mode with test requirements
	$(PYTHON) -m pip install -e .

dev-prerequisites: requirements-dev.txt
	$(PYTHON) -m pip install -r requirements-dev.txt

install: $(PACKAGE_FILE) ## install the package to the active Python's site-packages
	$(PYTHON) -m pip install $<

dist: $(PACKAGE_FILE) ## build package

$(PACKAGE_FILE): versionfile ## create distribution package
	$(PYTHON) -m pip install --upgrade build
	$(PYTHON) -m build

versionfile:
	@echo $(VERSION) >> Versionfile

version:
	@echo $(VERSION)
