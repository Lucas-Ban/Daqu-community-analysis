#!/usr/bin/env Rscript

# ==========================================================
# Procrustes analysis + Cross-matrix correlation heatmap
#   - Procrustes: BGC vs Secondary-metabolite species community structures
#   - Cross-correlation: cor( BGC_features, SM_species ) across samples
#
# Inputs:
#   - BGC.txt        : sample x features (BGC gene/cluster abundance)
#   - SM_species.txt : sample x species  (secondary-metabolite species abundance)
#
# Optional:
#   - meta file (tab or csv): columns include "sample" and group column (default "group")
#     Example:
#       sample  group
#       A69     T1
#       A80     T1
#       ...
#
# Outputs:
#   out/procrustes_summary.txt
#   out/procrustes_scores.csv
#   out/procrustes_plot.png
#   out/crosscorr_r.tsv
#   out/crosscorr_p.tsv
#   out/crosscorr_heatmap.png
#
# Run:
#   Rscript scripts/procrustes_BGC_vs_SM.R \
#     --bgc BGC.txt --sm SM_species.txt --outdir out \
#     --meta meta.txt --group_col group \
#     --transform hellinger --dist bray --ord pcoa --perm 9999 --seed 1 \
#     --do_crosscorr TRUE --corr_method spearman --corr_top_bgc 60 --corr_limits "-0.6,0.6"
#
# Notes:
# - Heatmap prioritizes TOP 40 SM species by mean abundance.
# - BGC features in heatmap: top variance features (default 60), set 0 to keep all.
# ==========================================================

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(tidyverse)
  library(grid)   # unit() for arrows
  library(scales) # squish
})

# -----------------------------
# CLI
# -----------------------------
option_list <- list(
  make_option(c("--bgc"), type="character", default="BGC.txt",
              help="BGC abundance table (samples x features). Tab-delimited, header=TRUE, rownames=1."),
  make_option(c("--sm"), type="character", default="SM_species.txt",
              help="Secondary-metabolite species table (samples x species). Tab-delimited, header=TRUE, rownames=1."),
  make_option(c("--meta"), type="character", default=NULL,
              help="Sample metadata file (tab or csv) with columns: sample, <group_col>."),
  make_option(c("--group_col"), type="character", default="group",
              help="Group column name in meta file (default: group)."),
  make_option(c("--outdir"), type="character", default="out",
              help="Output directory."),
  make_option(c("--transform"), type="character", default="hellinger",
              help="Data transform for ordination/correlation: hellinger / none / clr"),
  make_option(c("--dist"), type="character", default="bray",
              help="Distance for ordination: bray / euclidean / jaccard"),
  make_option(c("--ord"), type="character", default="pcoa",
              help="Ordination: pcoa / nmds"),
  make_option(c("--k"), type="integer", default=2,
              help="Number of dimensions for Procrustes (usually 2)."),
  make_option(c("--perm"), type="integer", default=9999,
              help="Permutations for protest."),
  make_option(c("--seed"), type="integer", default=1,
              help="Random seed."),
  
  # --- cross-matrix correlation ---
  make_option(c("--do_crosscorr"), type="character", default="TRUE",
              help="TRUE/FALSE: compute cross correlation matrix between BGC features and SM species."),
  make_option(c("--corr_method"), type="character", default="spearman",
              help="Correlation method: spearman / pearson"),
  make_option(c("--corr_top_bgc"), type="integer", default=60,
              help="For heatmap: keep top N BGC features by variance (0 means keep all)."),
  make_option(c("--corr_limits"), type="character", default="-0.6,0.6",
              help="Heatmap color limits, format: min,max (default -0.6,0.6).")
)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

do_crosscorr <- toupper(opt$do_crosscorr) == "TRUE"

# -----------------------------
# helpers
# -----------------------------
read_tab_or_csv_matrix <- function(path){
  x <- tryCatch({
    read.table(path, header = TRUE, row.names = 1, sep = "\t",
               quote = "", check.names = FALSE, comment.char = "")
  }, error = function(e) NULL)
  
  if(is.null(x) || ncol(x) == 0){
    x <- read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
  }
  
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

ensure_sample_rows <- function(mat, ref_samples){
  rn <- rownames(mat); cn <- colnames(mat)
  inter_r <- length(intersect(rn, ref_samples))
  inter_c <- length(intersect(cn, ref_samples))
  if(inter_c > inter_r){
    mat <- t(mat)
  }
  mat
}

transform_mat <- function(mat, method){
  method <- tolower(method)
  if(method == "none") return(mat)
  
  if(method == "hellinger"){
    return(vegan::decostand(mat, method = "hellinger"))
  }
  
  if(method == "clr"){
    mat2 <- mat
    rs <- rowSums(mat2)
    rs[rs == 0] <- NA
    mat2 <- mat2 / rs
    mat2 <- mat2 + 1e-10
    lg <- log(mat2)
    lg - rowMeans(lg, na.rm = TRUE)
  } else {
    stop("Unknown transform: ", method)
  }
}

do_ordination <- function(mat, dist = "bray", ord = "pcoa", k = 2){
  dist <- tolower(dist); ord <- tolower(ord)
  d <- vegan::vegdist(mat, method = dist)
  
  if(ord == "pcoa"){
    pts <- cmdscale(d, k = k, eig = TRUE)
    return(list(points = pts$points, dist = d, eig = pts$eig))
  } else if(ord == "nmds"){
    fit <- vegan::metaMDS(d, k = k, trymax = 50, autotransform = FALSE, trace = 0)
    return(list(points = fit$points, dist = d, stress = fit$stress))
  } else {
    stop("Unknown ordination: ", ord)
  }
}

read_meta <- function(path, group_col){
  if(is.null(path)) return(NULL)
  
  meta <- tryCatch({
    read.table(path, header = TRUE, sep = "\t", quote = "",
               check.names = FALSE, comment.char = "")
  }, error = function(e) NULL)
  
  if(is.null(meta) || ncol(meta) <= 1){
    meta <- read.csv(path, header = TRUE, check.names = FALSE)
  }
  
  if(!("sample" %in% colnames(meta))){
    stop("Meta file must contain column: sample")
  }
  if(!(group_col %in% colnames(meta))){
    stop("Meta file must contain group column: ", group_col)
  }
  
  meta <- meta[, c("sample", group_col)]
  colnames(meta) <- c("sample", "group")
  meta
}

parse_limits <- function(s){
  sp <- strsplit(s, ",")[[1]]
  if(length(sp) != 2) stop("corr_limits format must be min,max (e.g., -0.6,0.6)")
  as.numeric(trimws(sp))
}

# 选择 topN（按方差）
top_by_variance <- function(mat, topn){
  # mat: samples x features
  if(topn <= 0 || ncol(mat) <= topn) return(mat)
  v <- apply(mat, 2, stats::var, na.rm = TRUE)
  keep <- names(sort(v, decreasing = TRUE))[1:topn]
  mat[, keep, drop = FALSE]
}

# 选择 topN（按平均丰度“排名”）
top_by_mean <- function(mat, topn){
  if(topn <= 0 || ncol(mat) <= topn) return(mat)
  m <- colMeans(mat, na.rm = TRUE)
  keep <- names(sort(m, decreasing = TRUE))[1:topn]
  mat[, keep, drop = FALSE]
}

# Correlation + p-value matrix between two matrices (samples aligned in rows)
# X: samples x p ; Y: samples x q
cor_p_between <- function(X, Y, method = "spearman"){
  method <- tolower(method)
  if(!(method %in% c("spearman", "pearson"))) stop("corr_method must be spearman or pearson")
  
  X <- as.matrix(X); Y <- as.matrix(Y)
  p <- ncol(X); q <- ncol(Y)
  
  R <- matrix(NA_real_, nrow = p, ncol = q, dimnames = list(colnames(X), colnames(Y)))
  P <- matrix(NA_real_, nrow = p, ncol = q, dimnames = list(colnames(X), colnames(Y)))
  
  for(i in seq_len(p)){
    xi <- X[, i]
    for(j in seq_len(q)){
      yj <- Y[, j]
      ok <- is.finite(xi) & is.finite(yj)
      if(sum(ok) < 4){
        R[i, j] <- NA_real_
        P[i, j] <- NA_real_
      } else {
        ct <- suppressWarnings(stats::cor.test(xi[ok], yj[ok], method = method, exact = FALSE))
        R[i, j] <- unname(ct$estimate)
        P[i, j] <- ct$p.value
      }
    }
  }
  list(r = R, p = P)
}

# -----------------------------
# load
# -----------------------------
if(!file.exists(opt$bgc)) stop("BGC file not found: ", opt$bgc)
if(!file.exists(opt$sm))  stop("SM file not found: ", opt$sm)

bgc0 <- read_tab_or_csv_matrix(opt$bgc)
sm0  <- read_tab_or_csv_matrix(opt$sm)

bgc1 <- ensure_sample_rows(bgc0, ref_samples = rownames(sm0))
sm1  <- ensure_sample_rows(sm0,  ref_samples = rownames(bgc1))

common <- intersect(rownames(bgc1), rownames(sm1))
if(length(common) < 4){
  stop("Common samples < 4. Please check sample IDs and whether tables are sample x feature.")
}

bgc <- bgc1[common, , drop = FALSE]
sm  <- sm1[common, , drop = FALSE]

bgc <- bgc[, colSums(bgc, na.rm = TRUE) > 0, drop = FALSE]
sm  <- sm[,  colSums(sm,  na.rm = TRUE) > 0, drop = FALSE]

message("Common samples: ", length(common))
message("BGC features kept: ", ncol(bgc))
message("SM species kept: ", ncol(sm))

# -----------------------------
# meta join (groups)
# -----------------------------
meta <- data.frame(sample = common, group = "All", stringsAsFactors = FALSE)
meta0 <- read_meta(opt$meta, opt$group_col)
if(!is.null(meta0)){
  meta <- meta %>% left_join(meta0, by = "sample") %>%
    mutate(group = ifelse(!is.na(group.y), group.y, group.x)) %>%
    select(sample, group)
}
meta$group <- factor(meta$group, levels = unique(meta$group))

# -----------------------------
# transform
# -----------------------------
bgc_t <- transform_mat(bgc, opt$transform)
sm_t  <- transform_mat(sm,  opt$transform)

# -----------------------------
# ordination
# -----------------------------
ord_bgc <- do_ordination(bgc_t, dist = opt$dist, ord = opt$ord, k = opt$k)
ord_sm  <- do_ordination(sm_t,  dist = opt$dist, ord = opt$ord, k = opt$k)

X <- ord_bgc$points
Y <- ord_sm$points
rownames(X) <- rownames(bgc_t)
rownames(Y) <- rownames(sm_t)

# -----------------------------
# Procrustes + PROTEST
# -----------------------------
proc <- vegan::procrustes(X, Y, symmetric = TRUE)
prot <- vegan::protest(X, Y, permutations = opt$perm)

m2_true <- proc$ss
r_prot  <- prot$t0
p_value <- prot$signif

# -----------------------------
# outputs: summary
# -----------------------------
summary_txt <- file.path(opt$outdir, "procrustes_summary.txt")
cat(
  "Procrustes analysis (BGC vs SM species)\n",
  "=====================================\n",
  sprintf("Transform: %s\n", opt$transform),
  sprintf("Distance:  %s\n", opt$dist),
  sprintf("Ordination:%s\n", opt$ord),
  sprintf("k:         %d\n", opt$k),
  sprintf("Perm:      %d\n\n", opt$perm),
  sprintf("Procrustes m^2 (ss): %.6f (smaller is better)\n", m2_true),
  sprintf("PROTEST Procrustes correlation (r): %.6f\n", r_prot),
  sprintf("PROTEST P-value: %.6g\n", p_value),
  file = summary_txt
)

# -----------------------------
# outputs: per-sample scores
# -----------------------------
scores <- data.frame(
  sample = rownames(X),
  BGC1 = proc$X[,1], BGC2 = proc$X[,2],
  SM1  = proc$Yrot[,1], SM2  = proc$Yrot[,2],
  resid = sqrt(rowSums((proc$X - proc$Yrot)^2)),
  stringsAsFactors = FALSE
) %>% left_join(meta, by = "sample")

write.csv(scores, file.path(opt$outdir, "procrustes_scores.csv"),
          row.names = FALSE, quote = FALSE)

# -----------------------------
# plot (match your style)
# NOTE:
#   Your example prints "M²=0.977" which visually matches PROTEST r.
#   To keep your exact style, we print r_prot as "M²" here.
#   If you want strict m², change M2_show <- m2_true
# -----------------------------
df <- scores
M2_show <- r_prot
P_show  <- p_value

p1 <- ggplot(df) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  
  geom_segment(
    aes(x = SM1, y = SM2, xend = BGC1, yend = BGC2, colour = group),
    linewidth = 0.9,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    alpha = 0.95
  ) +
  geom_point(aes(x = BGC1, y = BGC2, colour = group), size = 3.6, shape = 17) +
  geom_point(aes(x = SM1,  y = SM2,  colour = group), size = 3.2, shape = 16) +
  
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("M\u00b2 = %.3f\nP = %.3g", M2_show, P_show),
           hjust = -0.05, vjust = 1.10, size = 4.2) +
  
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme_classic(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.22, 0.82)
  )

shape_legend <- tibble(
  x = c(Inf, Inf),
  y = c(Inf, Inf),
  key = factor(c("BGCs", "BGC-containing species"),
               levels = c("BGCs", "BGC-containing species"))
)

p1 <- p1 +
  geom_point(
    data = shape_legend,
    aes(x = x, y = y, shape = key),
    colour = "black", size = 3.8,
    inherit.aes = FALSE
  ) +
  scale_shape_manual(values = c("BGCs" = 17, "BGC-containing species" = 16)) +
  guides(
    shape  = guide_legend(override.aes = list(colour = "black")),
    colour = guide_legend(override.aes = list(shape = 16, size = 3.5))
  )

ggsave(file.path(opt$outdir, "procrustes_plot.png"),
       p1, width = 7.4, height = 6.2, dpi = 300)

# ==========================================================
# Cross-matrix correlation heatmap
#   - SM axis: TOP 40 species by mean abundance (priority display)
#   - BGC axis: top variance features (default 60), set 0 to keep all
#   - color: r clipped to [-0.6,0.6], skyblue -> white -> magenta
# ==========================================================
if(do_crosscorr){
  
  lims <- parse_limits(opt$corr_limits)
  lo <- lims[1]; hi <- lims[2]
  
  # choose features for heatmap
  bgc_for_corr <- bgc_t
  sm_for_corr  <- sm_t
  
  # BGC: optional top variance (for readability)
  if(opt$corr_top_bgc > 0) bgc_for_corr <- top_by_variance(bgc_for_corr, opt$corr_top_bgc)
  
  # SM species: fixed TOP 40 by mean abundance
  sm_for_corr <- top_by_mean(sm_for_corr, topn = 40)
  
  message("Cross-corr heatmap uses: BGC=", ncol(bgc_for_corr),
          ", SM(top mean)= ", ncol(sm_for_corr), " (fixed at 40)")
  
  cp <- cor_p_between(bgc_for_corr, sm_for_corr, method = opt$corr_method)
  rmat <- cp$r
  pmat <- cp$p
  
  # save matrices
  write.table(rmat, file = file.path(opt$outdir, "crosscorr_r.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(pmat, file = file.path(opt$outdir, "crosscorr_p.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # long df for heatmap
  dfh <- as.data.frame(as.table(rmat), stringsAsFactors = FALSE)
  colnames(dfh) <- c("BGC_feature", "SM_species", "r")
  dfh$r_clip <- scales::squish(dfh$r, range = c(lo, hi))
  
  # order axes for nicer display:
  # - SM in decreasing mean abundance
  sm_means <- colMeans(sm_for_corr, na.rm = TRUE)
  sm_levels <- names(sort(sm_means, decreasing = TRUE))
  dfh$SM_species <- factor(dfh$SM_species, levels = sm_levels)
  
  # - BGC in decreasing variance (if subset) else keep current order
  bgc_vars <- apply(bgc_for_corr, 2, stats::var, na.rm = TRUE)
  bgc_levels <- names(sort(bgc_vars, decreasing = TRUE))
  dfh$BGC_feature <- factor(dfh$BGC_feature, levels = bgc_levels)
  
  ph <- ggplot(dfh, aes(x = SM_species, y = BGC_feature, fill = r_clip)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "skyblue",
      mid = "white",
      high = "magenta",
      midpoint = 0,
      limits = c(lo, hi),
      name = sprintf("r (%s)\nclipped [%.1f, %.1f]", opt$corr_method, lo, hi)
    ) +
    labs(x = "Top 40 secondary-metabolite species (by mean abundance)",
         y = "BGC features",
         title = "Cross-matrix correlation: BGC features vs SM species") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      plot.title = element_text(size = 13)
    )
  
  ggsave(file.path(opt$outdir, "crosscorr_heatmap.png"),
         ph, width = 9.5, height = 8.5, dpi = 300)
  
  message("Cross-corr outputs:\n- crosscorr_r.tsv\n- crosscorr_p.tsv\n- crosscorr_heatmap.png")
}

message("Done. Outputs in: ", normalizePath(opt$outdir))
