# iCAMP analysis pipeline (16S / ITS) + robustness simulations

This repository provides a reproducible R workflow to run iCAMP-based community assembly analysis
and (optionally) robustness/removal simulations for microbial ecological networks.  
It also includes additional downstream analyses commonly used in the revision workflow:
Cohesion, MEN/network construction, Procrustes, Mantel-link correlation plots, NMDS, and
Aitchison distance + alpha diversity indices.

## Directory structure

- `data/` : input data (not tracked by git by default; see `.gitignore`)
- `R/` : core R functions (reusable helpers)
- `scripts/` : entry scripts (one-click runners)
- `config/` : YAML configuration files
- `output/` or `out/` : generated results (not tracked by git)

> Note: Some scripts write to `out/` by default. You can change the output folder in each script.

## Requirements

- R (>= 4.1 recommended)

### Core packages (iCAMP)
- `iCAMP`, `ape`, `yaml`

Install packages in R:

```r
install.packages(c("yaml", "ape"))
install.packages("iCAMP")
