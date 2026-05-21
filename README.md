# breast_wgan_2026

Public analysis scripts for the paper:

Translating lineage-resolved single-cell programs to bulk clinical prognosis: adversarial generative learning reveals lineage-informed immune - stromal - epithelial risk signatures

This repository contains a script-focused release for model reruns, benchmarks, and figure panel generation.

## Setup

```bash
conda env create -f environment.yml
conda activate paper3-breast-wgan

# Use native/system R 4.5.3 on PATH for the manuscript R stack
Rscript scripts/setup/install_packages.R
bash scripts/setup/install_python_requirements.sh
```

`environment.yml` provisions the Conda-side Python/helper tooling only. The manuscript R stack is intentionally installed with native `Rscript` through `scripts/setup/install_packages.R`, which pins the R and Bioconductor package versions used for the workflow.
