# =========================================================
# Mantel test + correlation heatmap + links (Bacteria & Fungi)
# Style: heatmap + Mantel links; P color; R width
#
# Inputs:
#   meta.txt              : Sample, Temperature, Fermentation_time
#   process_bacteria.txt  : Sample, (5 process proportion columns)
#   process_fungi.txt     : Sample, (5 process proportion columns)
#
# Outputs:
#   out/mantel_bacteria_fungi.png
#   out/mantel_bacteria.png
#   out/mantel_fungi.png
# =========================================================

# ---------- packages ----------
pkgs <- c("tidyverse", "vegan", "patchwork")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if(length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)

# linkET
if(!requireNamespace("linkET", quietly = TRUE)){
  install.packages("linkET", repos = "https://cloud.r-project.org")
}
library(linkET)

# ---------- user params ----------
META_FILE <- "meta.txt"
BAC_FILE  <- "process_bacteria.txt"
FUN_FILE  <- "process_fungi.txt"

OUTDIR <- "out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Mantel settings
N_PERM <- 9999
SPEC_DIST <- "bray"       # distance for process composition (proportions)
ENV_DIST  <- "euclidean"  # distance for env variables (univariate)
COR_METHOD <- "pearson"   # heatmap correlations among processes

# Plot settings
PANEL_W <- 12
PANEL_H <- 5
DPI <- 300

# ---------- helpers ----------
read_table_auto <- function(path){
  # Try tab first; if fails, try csv
  x <- tryCatch({
    read.table(path, header = TRUE, sep = "\t", quote = "", check.names = FALSE, comment.char = "")
  }, error = function(e) NULL)
  
  if(is.null(x) || ncol(x) <= 1){
    x <- read.csv(path, header = TRUE, check.names = FALSE)
  }
  x
}

as_sample_row_matrix <- function(df, sample_col = "Sample"){
  if(sample_col %in% colnames(df)){
    rn <- df[[sample_col]]
    df[[sample_col]] <- NULL
    mat <- as.matrix(df)
    rownames(mat) <- rn
  } else {
    mat <- as.matrix(df)
  }
  storage.mode(mat) <- "numeric"
  mat
}

prep_inputs <- function(meta_df, proc_df,
                        sample_col = "Sample",
                        temp_col = "Temperature",
                        time_col = "Fermentation_time",
                        process_cols = NULL){
  meta_mat <- as_sample_row_matrix(meta_df, sample_col = sample_col)
  
  stopifnot(temp_col %in% colnames(meta_mat))
  stopifnot(time_col %in% colnames(meta_mat))
  
  proc_mat0 <- as_sample_row_matrix(proc_df, sample_col = sample_col)
  
  if(is.null(process_cols)){
    process_cols <- colnames(proc_mat0)
  } else {
    stopifnot(all(process_cols %in% colnames(proc_mat0)))
  }
  proc_mat <- proc_mat0[, process_cols, drop = FALSE]
  
  # align samples
  common <- intersect(rownames(meta_mat), rownames(proc_mat))
  if(length(common) < 4) stop("Common samples < 4, please check sample IDs.")
  meta_mat <- meta_mat[common, c(temp_col, time_col), drop = FALSE]
  proc_mat <- proc_mat[common, , drop = FALSE]
  
  # drop all-zero process columns
  proc_mat <- proc_mat[, colSums(proc_mat, na.rm = TRUE) > 0, drop = FALSE]
  
  list(env = meta_mat, proc = proc_mat)
}

make_panel <- function(env_mat, proc_mat, title_text){
  # Heatmap: correlations among processes
  proc_df <- as.data.frame(proc_mat)
  
  p_cor <- linkET::quickcor(proc_df, method = COR_METHOD, type = "upper") +
    geom_square() +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson's R") +
    theme(
      legend.position = "left",
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # Mantel: env (Temperature, time) vs process composition (all processes together)
  env_df <- as.data.frame(env_mat)
  proc_df2 <- as.data.frame(proc_mat)
  
  mres <- linkET::mantel_test(
    spec = proc_df2,
    env  = env_df,
    spec.dist.method = SPEC_DIST,
    env.dist.method  = ENV_DIST,
    permutations = N_PERM
  )
  
  # Bin P & R like your legend style
  mres <- mres %>%
    mutate(
      `Mantel's P` = case_when(
        p.value < 0.01 ~ "< 0.01",
        p.value < 0.05 ~ "0.01 – 0.05",
        TRUE ~ ">= 0.05"
      ),
      `Mantel's R` = case_when(
        statistic < 0.2 ~ "< 0.2",
        statistic < 0.4 ~ "0.2 – 0.4",
        TRUE ~ ">= 0.4"
      )
    )
  
  # Add links to heatmap
  p <- p_cor +
    anno_link(
      data = mres,
      aes(colour = `Mantel's P`, size = `Mantel's R`),
      curvature = 0.2
    ) +
    scale_colour_manual(
      values = c("< 0.01" = "#e74c3c", "0.01 – 0.05" = "#f39c12", ">= 0.05" = "#7fbf7b"),
      name = "Mantel's P"
    ) +
    scale_size_manual(
      values = c("< 0.2" = 0.8, "0.2 – 0.4" = 1.8, ">= 0.4" = 3.0),
      name = "Mantel's R"
    ) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.box = "vertical"
    )
  
  p
}

# ---------- main ----------
meta_df <- read_table_auto(META_FILE)
bac_df  <- read_table_auto(BAC_FILE)
fun_df  <- read_table_auto(FUN_FILE)

# If your sample column is not "Sample", change here:
SAMPLE_COL <- "Sample"
TEMP_COL   <- "Temperature"
TIME_COL   <- "Fermentation_time"

# If you want to lock the five process columns, set them explicitly:
# PROCESS_COLS <- c("Selection","Dispersal_limitation","Homogenizing_dispersal","Drift","Undominated")
PROCESS_COLS <- NULL

bac_in <- prep_inputs(meta_df, bac_df,
                      sample_col = SAMPLE_COL, temp_col = TEMP_COL, time_col = TIME_COL,
                      process_cols = PROCESS_COLS)

fun_in <- prep_inputs(meta_df, fun_df,
                      sample_col = SAMPLE_COL, temp_col = TEMP_COL, time_col = TIME_COL,
                      process_cols = PROCESS_COLS)

p_bac <- make_panel(bac_in$env, bac_in$proc, "Bacteria")
p_fun <- make_panel(fun_in$env, fun_in$proc, "Fungi")

p_all <- p_bac + p_fun + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "left")

ggsave(file.path(OUTDIR, "mantel_bacteria_fungi.png"),
       p_all, width = PANEL_W, height = PANEL_H, dpi = DPI)

ggsave(file.path(OUTDIR, "mantel_bacteria.png"),
       p_bac, width = PANEL_W/2, height = PANEL_H, dpi = DPI)

ggsave(file.path(OUTDIR, "mantel_fungi.png"),
       p_fun, width = PANEL_W/2, height = PANEL_H, dpi = DPI)

cat("Done.\nOutputs:\n",
    "- out/mantel_bacteria_fungi.png\n",
    "- out/mantel_bacteria.png\n",
    "- out/mantel_fungi.png\n", sep = "")
