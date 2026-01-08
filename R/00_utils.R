# =========================================================
# One-click runner for ALL figure scripts
# Put this file under: R/00_run_all_figures.R
#
# Usage:
#   source("R/00_run_all_figures.R")
# or
#   Rscript R/00_run_all_figures.R
# =========================================================

# ---- 1) always setwd to project root ----
wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
# if running inside .../R or .../scripts, go back to root
if (grepl("/R$", wd) || grepl("/scripts$", wd)) setwd(dirname(wd))

cat("Project root:", normalizePath(getwd(), winslash = "/"), "\n")

# ---- 2) (optional) global options ----
options(stringsAsFactors = FALSE)

# ---- 3) define figure scripts to run (in order) ----
# IMPORTANT:
# - If your figure scripts are stored in the repo root (same level as R/), use "Figure_*.R"
# - If stored inside R/, use "R/Figure_*.R"
# Please adjust paths to match your repo.

FIG_SCRIPTS <- c(
  # Main figures
  "Figure_1A-B.R",
  "Figure_1C.R",
  "Figure_1D-E.R",
  "Figure_3E-F.R",
  "Figure_4C.R",
  "Figure_5A-B.R",
  
  # Supplementary
  "Figure_S4.R"
)

# If your figure scripts are actually under R/ folder, uncomment this:
# FIG_SCRIPTS <- file.path("R", FIG_SCRIPTS)

# ---- 4) run with logging ----
log_dir <- file.path("out", "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

run_one <- function(path){
  if(!file.exists(path)){
    stop("File not found: ", path, "\n(Current wd: ", getwd(), ")")
  }
  
  cat("\n=============================\n")
  cat("Running:", path, "\n")
  cat("=============================\n")
  
  # capture console output
  log_file <- file.path(log_dir, paste0(gsub("[/\\\\: ]", "_", path), ".log"))
  zz <- file(log_file, open = "wt")
  sink(zz, type = "output")
  sink(zz, type = "message")
  
  ok <- TRUE
  err_msg <- NULL
  t0 <- Sys.time()
  
  tryCatch({
    source(path, local = new.env(parent = globalenv()))
  }, error = function(e){
    ok <<- FALSE
    err_msg <<- conditionMessage(e)
  })
  
  t1 <- Sys.time()
  
  sink(type = "message")
  sink(type = "output")
  close(zz)
  
  if(ok){
    cat("OK:", path, " | time:", round(as.numeric(difftime(t1, t0, units = "secs")), 1), "s\n")
  } else {
    cat("FAILED:", path, "\n  Error:", err_msg, "\n  Log:", log_file, "\n")
    stop("Stop at failed script: ", path)
  }
}

for(f in FIG_SCRIPTS){
  run_one(f)
}

cat("\nAll figure scripts finished.\nLogs in:", normalizePath(log_dir, winslash="/"), "\n")
