#!/usr/bin/env bash
set -e

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
PY310="/usr/local/bin/python3.10"

echo "Project: $PROJECT_DIR"
echo "Using Python: $PY310"

if [ ! -x "$PY310" ]; then
  echo "ERROR: Python 3.10 not found at $PY310"
  exit 1
fi

cd "$PROJECT_DIR"

if [ -d "venv" ]; then
  echo "Removing old virtual environment"
  rm -rf venv
fi

echo "Creating virtual environment"
"$PY310" -m venv venv

echo "Activating virtual environment"
source venv/bin/activate

echo "Python in venv:"
python --version

echo "Upgrading pip"
pip install -U pip

if [ -f "requirements.txt" ]; then
  echo "Installing requirements.txt"
  pip install -r requirements.txt
else
  echo "Installing otree"
  pip install otree
fi

echo "Starting oTree devserver"
otree devserver
