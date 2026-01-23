#!/usr/bin/env Rscript

# Entry point: run the iCAMP pipeline
# Usage (in R): source("scripts/run_iCAMP.R")

# If executed from scripts/, go back to project root
wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (grepl("/scripts$", wd)) setwd(dirname(wd))

source("R/00_utils.R")
source("R/02_run_icamp.R")

cfg <- read_config("config/config.yml")
ensure_dir(cfg$output$out_dir)

message("[run_iCAMP] Starting iCAMP pipeline ...")
run_icamp_pipeline(cfg)
message("[run_iCAMP] Finished. Outputs are in: ", cfg$output$out_dir)
