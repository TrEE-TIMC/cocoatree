#!/usr/bin/env bash
# Fail on non-zero exit and echo the commands
set -ev

python -m pip list

(cd .. && pytest --doctest-modules --cov=akasthesia --pyargs akasthesia)
(cd ..\\akasthesia && nose2)
flake8 --exit-zero akasthesia examples

set +ev
