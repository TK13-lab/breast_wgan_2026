# breast_wgan_2026

Public analysis scripts for the paper:

Translating lineage-resolved single-cell programs to bulk clinical prognosis: adversarial generative learning reveals lineage-informed immune - stromal - epithelial risk signatures

This repository contains a script-focused release for model reruns, benchmarks, and figure panel generation.



`environment.yml` provisions the Conda-side Python/helper tooling only. The manuscript R stack is intentionally installed with native `Rscript` through `scripts/setup/install_packages.R`, which pins the R and Bioconductor package versions used for the workflow.
