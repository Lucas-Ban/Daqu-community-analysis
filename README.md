# iCAMP analysis pipeline (16S / ITS) + robustness simulations

This repository provides a reproducible R workflow to run iCAMP-based community assembly analysis
and (optionally) robustness/removal simulations for microbial ecological networks.

## Directory structure

- `data/` : input data (not tracked by git by default; see `.gitignore`)
- `R/` : core R functions
- `scripts/` : entry scripts
- `config/` : YAML configuration files
- `output/` : generated results (not tracked by git)

## Requirements

- R (>= 4.1 recommended)
- R packages:
  - `iCAMP`, `ape`, `yaml`
  - (optional, for robustness) `foreach`, `doParallel`

Install packages in R:

```r
install.packages(c("yaml", "ape", "foreach", "doParallel"))
# iCAMP is on CRAN for most environments; if not, install from the official source you used.
install.packages("iCAMP")
