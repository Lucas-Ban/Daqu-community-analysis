#!/usr/bin/env Rscript

# Usage (in R): source("scripts/run_robustness.R")

source("R/00_utils.R")
source("R/run_robustness.R")

cfg <- read_config("config/config.yml")

# robustness output directory
ensure_dir(cfg$robustness$out_dir)

message("[run_robustness] Starting ...")
run_robustness_pipeline(cfg)
message("[run_robustness] Finished. Outputs are in: ", cfg$robustness$out_dir)
