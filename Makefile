PYTHON ?= python
PYTEST ?= pytest
CTAGS ?= ctags

all: clean inplace test

inplace: 
	$(PYTHON) setup.py build_ext -i

test: test-code

test-code: inplace
	$(PYTEST) --showlocals -v akasthesia --durations=20

test-coverage:
	rm -rf coverage .coverage
	$(PYTEST) akasthesia --doctest-modules --showlocals -v --cov=akasthesia

clean-ctags:
	rm -f tags

clean: clean-ctags
	$(PYTHON) setup.py clean
	rm -rf dist
	rm -rf build

trailing-spaces:
	find akasthesia -name "*.py" -exec perl -pi -e 's/[ \t]*$$//' {} \;
