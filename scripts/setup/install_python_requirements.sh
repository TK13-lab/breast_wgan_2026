#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
python_bin="${PYTHON_BIN:-python3}"
venv_dir="${project_root}/.venv"

"${python_bin}" -m venv "${venv_dir}"
"${venv_dir}/bin/pip" install --upgrade pip
"${venv_dir}/bin/pip" install numpy pandas scikit-learn matplotlib torch

echo "Python environment ready at ${venv_dir}"
echo "Interpreter: ${venv_dir}/bin/python"
