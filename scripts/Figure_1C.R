# =========================================================
# Aitchison distance + Simpson index
# Input : species.csv (rows=samples, cols=taxa/features)
# Output:
#   out/Aitchison_distance_matrix.tsv
#   out/Aitchison_adjacent_distance.tsv
#   out/Simpson_indices.tsv
# =========================================================

suppressPackageStartupMessages({
  library(compositions)
  library(vegan)
})

# -----------------------------
# user settings
# -----------------------------
IN_FILE <- "species.csv"   # CSV
OUTDIR <- "out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

PSEUDOCOUNT <- 1e-6        # avoid log(0) before CLR
DO_CLOSURE  <- TRUE        # recommended for compositional data
ORDER_FILE  <- NULL        # optional: file listing sample order (one sample per line)

# -----------------------------
# read data (CSV)
# -----------------------------
df <- read.csv(IN_FILE, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
df <- na.omit(df)

# detect ID column (prefer ID/Sample; otherwise first column)
id_col <- NULL
for(cn in c("ID", "Sample", "sample", "id")){
  if(cn %in% colnames(df)) { id_col <- cn; break }
}
if(is.null(id_col)) id_col <- colnames(df)[1]

rownames(df) <- df[[id_col]]
df[[id_col]] <- NULL

x <- as.matrix(df)
storage.mode(x) <- "numeric"

# remove all-zero samples / taxa
x <- x[rowSums(x, na.rm = TRUE) > 0, , drop = FALSE]
x <- x[, colSums(x, na.rm = TRUE) > 0, drop = FALSE]

# optional reorder samples
if(!is.null(ORDER_FILE) && file.exists(ORDER_FILE)){
  ord <- readLines(ORDER_FILE)
  ord <- ord[ord %in% rownames(x)]
  if(length(ord) >= 2){
    x <- x[ord, , drop = FALSE]
  }
}

# -----------------------------
# Aitchison distance
# Steps: (optional closure) -> add pseudocount -> CLR -> Euclidean distance
# -----------------------------
x_comp <- x

if(DO_CLOSURE){
  rs <- rowSums(x_comp)
  rs[rs == 0] <- NA
  x_comp <- x_comp / rs
}

x_comp <- x_comp + PSEUDOCOUNT

clr_x <- compositions::clr(x_comp)

d <- dist(clr_x, method = "euclidean")
dmat <- as.matrix(d)
rownames(dmat) <- rownames(x)
colnames(dmat) <- rownames(x)

write.table(dmat,
            file = file.path(OUTDIR, "Aitchison_distance_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# adjacent (i vs i-1) according to current row order
n <- nrow(dmat)
adj <- data.frame(
  sample = rownames(dmat)[2:n],
  prev_sample = rownames(dmat)[1:(n-1)],
  aitchison_distance = dmat[cbind(2:n, 1:(n-1))],
  stringsAsFactors = FALSE
)

write.table(adj,
            file = file.path(OUTDIR, "Aitchison_adjacent_distance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Simpson index
# vegan::diversity(x, "simpson") = 1 - sum(p^2)  (Gini-Simpson)
# also output classic Simpson D = sum(p^2) and InvSimpson
# -----------------------------
p <- x / rowSums(x)

simpson_gini <- vegan::diversity(x, index = "simpson")      # 1 - sum(p^2)
simpson_D    <- rowSums(p^2, na.rm = TRUE)                  # sum(p^2)
inv_simpson  <- vegan::diversity(x, index = "invsimpson")   # 1 / sum(p^2)

simpson_out <- data.frame(
  sample = rownames(x),
  Simpson_Gini = simpson_gini,
  Simpson_D = simpson_D,
  InvSimpson = inv_simpson,
  stringsAsFactors = FALSE
)

write.table(simpson_out,
            file = file.path(OUTDIR, "Simpson_indices.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\nOutputs:\n",
    "- out/Aitchison_distance_matrix.tsv\n",
    "- out/Aitchison_adjacent_distance.tsv\n",
    "- out/Simpson_indices.tsv\n", sep = "")
