#!/bin/bash

PYTHON_CMD=${PYTHON_VERSION:-python3}
ENVPATH="../.venv"

if ! command -v $PYTHON_CMD &> /dev/null; then
    echo "Error: $PYTHON_CMD not found. Please install Python or set PYTHON_VERSION environment variable."
    exit 1
fi

echo "Using Python: $(which $PYTHON_CMD)"
echo "Python version: $($PYTHON_CMD --version)"

$PYTHON_CMD -m venv $ENVPATH

source $ENVPATH/bin/activate
which python

pip install --upgrade pip
pip install numpy galois
