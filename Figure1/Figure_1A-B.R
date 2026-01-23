# =========================================================
# NMDS + ANOSIM using ONE species table, TWO grouping files
#   1) group_stime.csv        (time grouping)
#   2) group_temperature.csv  (temperature grouping)
#
# Input (CSV):
#   - species.csv: rows = samples, cols = taxa/features
#       first column is sample ID (ID/Sample/first col)
#   - group_*.csv: two columns: sample ID + group
#
# Output:
#   out/NMDS_time.png, out/NMDS_temperature.png
#   out/NMDS_time_points.csv, out/NMDS_temperature_points.csv
#   out/ANOSIM_time.txt, out/ANOSIM_temperature.txt
# =========================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# files
# -----------------------------
species_file <- "species.csv"
group_time_file <- "group_stime.csv"
group_temp_file <- "group_temperature.csv"

OUTDIR <- "out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# params
# -----------------------------
DIST_METHOD <- "bray"
PERM <- 999
K <- 2
TRYMAX <- 200
ELLIPSE_LEVEL <- 0.50  # like your example

# -----------------------------
# helpers
# -----------------------------
read_species_csv <- function(path){
  df <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  df <- na.omit(df)
  
  # find ID column
  id_col <- NULL
  for(cn in c("ID", "Sample", "sample", "id")){
    if(cn %in% colnames(df)) { id_col <- cn; break }
  }
  if(is.null(id_col)) id_col <- colnames(df)[1]
  
  rownames(df) <- df[[id_col]]
  df[[id_col]] <- NULL
  
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  
  # drop all-zero rows/cols
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  mat
}

read_group_csv <- function(path){
  df <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  df <- na.omit(df)
  
  # id column
  id_col <- NULL
  for(cn in c("ID", "Sample", "sample", "id")){
    if(cn %in% colnames(df)) { id_col <- cn; break }
  }
  if(is.null(id_col)) id_col <- colnames(df)[1]
  
  # group column
  grp_col <- NULL
  for(cn in c("group", "Group", "grp", "Grp")){
    if(cn %in% colnames(df)) { grp_col <- cn; break }
  }
  if(is.null(grp_col)) grp_col <- setdiff(colnames(df), id_col)[1]
  
  out <- df[, c(id_col, grp_col)]
  colnames(out) <- c("sample", "group")
  out
}

align_species_group <- function(species_mat, group_df){
  common <- intersect(rownames(species_mat), group_df$sample)
  if(length(common) < 4) stop("Common samples < 4. Check sample IDs in species.csv and group file.")
  
  species2 <- species_mat[common, , drop = FALSE]
  group2 <- group_df %>% filter(sample %in% common)
  
  # reorder rows to match group order
  group2 <- group2 %>% arrange(match(sample, rownames(species2)))
  species2 <- species2[group2$sample, , drop = FALSE]
  
  grp <- factor(group2$group)
  list(species = species2, group = grp, meta = group2)
}

run_nmds <- function(species_mat, seed = 1){
  set.seed(seed)
  d <- vegdist(species_mat, method = DIST_METHOD)
  fit <- metaMDS(d, k = K, trymax = TRYMAX, autotransform = FALSE, trace = 0)
  pts <- as.data.frame(scores(fit, display = "sites"))
  pts$sample <- rownames(pts)
  list(dist = d, fit = fit, points = pts)
}

run_anosim <- function(dist_obj, grp){
  anosim(dist_obj, grp, permutations = PERM)
}

save_anosim <- function(an, file){
  sink(file)
  cat("ANOSIM\n======\n")
  print(summary(an))
  cat("\nR statistic: ", an$statistic, "\n", sep = "")
  cat("P value:     ", an$signif, "\n", sep = "")
  sink()
}

plot_nmds <- function(points_df, title_text){
  ggplot(points_df, aes(NMDS1, NMDS2, colour = group, shape = group)) +
    stat_ellipse(level = ELLIPSE_LEVEL, linewidth = 0.8) +
    geom_point(size = 3.5) +
    theme_bw(base_size = 12) +
    labs(title = title_text, colour = NULL, shape = NULL)
}

do_one <- function(species_mat, group_file, tag){
  gdf <- read_group_csv(group_file)
  al <- align_species_group(species_mat, gdf)
  
  nm <- run_nmds(al$species, seed = 1)
  an <- run_anosim(nm$dist, al$group)
  
  pts <- nm$points
  pts$group <- al$group
  
  # add stress to title
  stress <- nm$fit$stress
  p <- plot_nmds(pts, sprintf("NMDS - %s (Bray), stress=%.3f", tag, stress))
  
  # outputs
  write.csv(pts, file.path(OUTDIR, paste0("NMDS_", tag, "_points.csv")),
            row.names = FALSE, quote = FALSE)
  
  ggsave(file.path(OUTDIR, paste0("NMDS_", tag, ".png")),
         p, width = 6.5, height = 5.5, dpi = 300)
  
  save_anosim(an, file.path(OUTDIR, paste0("ANOSIM_", tag, ".txt")))
  
  cat("[OK] ", tag, "\n", sep = "")
}

# -----------------------------
# main
# -----------------------------
if(!file.exists(species_file)) stop("species.csv not found: ", species_file)
species_mat <- read_species_csv(species_file)

do_one(species_mat, group_time_file, "time")
do_one(species_mat, group_temp_file, "temperature")

cat("Done. Check out/ folder.\n")
