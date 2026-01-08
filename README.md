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
install.packages(c(
  "vegan", "ggplot2", "patchwork", "tidyverse",
  "compositions", "linkET", "foreach", "doParallel"
))
install.packages(c(
  "vegan", "ggplot2", "patchwork", "tidyverse",
  "compositions", "linkET", "foreach", "doParallel"
))
Rscript scripts/run_iCAMP.R
2) Cohesion

Computes positive/negative connectedness and cohesion based on (observed - expected) correlations.

Input:

otu_all.txt (or configured path inside the script)
Rscript scripts/run_cohesion.R
3) Procrustes: BGC vs SM species

Performs Procrustes (ordination alignment) and exports a cross-matrix correlation heatmap.

Inputs:

BGC.txt (sample x BGC features)

SM_species.txt (sample x secondary-metabolite species)
Rscript scripts/procrustes_BGC_vs_SM.R --bgc BGC.txt --sm SM_species.txt --outdir out
Outputs:

out/procrustes_plot.png

out/crosscorr_heatmap.png (Top 40 SM species prioritized; color range [-0.6, 0.6])

4) Mantel test correlation plot (Bacteria / Fungi)

Produces heatmap (process-process correlations) + Mantel links from env factors to processes.

Inputs:

meta.txt (Sample, Temperature, Fermentation_time)

process_bacteria.txt (Sample + 5 process proportions)

process_fungi.txt (Sample + 5 process proportions)
Rscript scripts/run_mantel_linkET.R
Outputs:

out/mantel_bacteria_fungi.png
5) NMDS (one species table, two groupings)

Runs NMDS twice using the same species.csv but different grouping files.

Inputs:

species.csv

group_stime.csv

group_temperature.csv
Rscript scripts/NMDS_time_vs_temperature.R
Outputs:

out/NMDS_time.png

out/NMDS_temperature.png
Outputs:

out/NMDS_time.png

out/NMDS_temperature.png
Rscript scripts/aitchison_simpson_from_species_csv.R
Outputs:

out/Aitchison_distance_matrix.tsv

out/Aitchison_adjacent_distance.tsv

out/Simpson_indices.tsv

Reproducibility tips

Prefer running entry scripts from the project root.

If you run from inside R/ or scripts/, add an auto ¡°setwd to root¡± snippet:

detect /R or /scripts and setwd(dirname(getwd())).

Keep raw input data under data/ and avoid committing large generated outputs.