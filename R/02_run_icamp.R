# Main iCAMP pipeline
# - Reads comm/tree/clas from cfg
# - treat_file/env_file can be null (skipped)
# - Matches sample IDs (only if treat provided)
# - Matches taxa IDs across comm/clas/tree
# - Computes pdist.big (if pd.desc not exist)
# - Runs icamp.big and saves outputs

run_icamp_pipeline <- function(cfg) {
  # ---- required packages ----
  pkgs <- c("iCAMP", "ape")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "),
         "\nInstall them first (e.g., install.packages('iCAMP')).")
  }
  
  # ---- paths ----
  com_file   <- cfg$paths$com_file
  tree_file  <- cfg$paths$tree_file
  clas_file  <- cfg$paths$clas_file
  treat_file <- cfg$paths$treat_file
  env_file   <- cfg$paths$env_file
  
  out_dir <- cfg$output$out_dir
  ensure_dir(out_dir)
  
  # required inputs
  stop_if_missing(c(com_file, tree_file, clas_file), label = "input")
  
  # optional inputs
  if (!is.null(treat_file) && !identical(treat_file, "null")) stop_if_missing(treat_file, label = "treat")
  if (!is.null(env_file) && !identical(env_file, "null")) stop_if_missing(env_file, label = "env")
  
  # ---- parameters ----
  prefix         <- cfg$project$prefix
  rand.time      <- cfg$project$rand_time
  nworker        <- cfg$project$nworker
  memory.G       <- cfg$project$memory_G
  ds             <- cfg$project$ds
  bin.size.limit <- cfg$project$bin_size_limit
  sig.index      <- cfg$project$sig_index
  
  t0 <- Sys.time()
  
  # ---- read data ----
  comm <- t(read_tsv_table(com_file, row_names = TRUE)) # samples x taxa
  tree <- ape::read.tree(file = tree_file)
  clas <- read_tsv_table(clas_file, row_names = TRUE)
  
  treat <- NULL
  if (!is.null(treat_file) && !identical(treat_file, "null")) {
    treat <- read_tsv_table(treat_file, row_names = TRUE)
  }
  
  env <- NULL
  if (!is.null(env_file) && !identical(env_file, "null")) {
    env <- read_tsv_table(env_file, row_names = TRUE)
  }
  
  # ---- match sample IDs (only if treat exists) ----
  if (!is.null(treat)) {
    sampid.check <- iCAMP::match.name(rn.list = list(comm = comm, treat = treat))
    comm  <- sampid.check$comm
    treat <- sampid.check$treat
  }
  
  # remove taxa that become all-zero after any filtering
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  
  # ---- match taxa IDs across comm / clas / tree ----
  spid.check <- iCAMP::match.name(
    cn.list   = list(comm = comm),
    rn.list   = list(clas = clas),
    tree.list = list(tree = tree)
  )
  comm <- spid.check$comm
  clas <- spid.check$clas
  tree <- spid.check$tree
  
  # ---- phylogenetic distance big matrix ----
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)
  
  if (!file.exists("pd.desc")) {
    message("[iCAMP] Computing phylogenetic distance matrix (pdist.big) ...")
    pd.big <- iCAMP::pdist.big(tree = tree, wd = out_dir, nworker = nworker, memory.G = memory.G)
  } else {
    message("[iCAMP] Found pd.desc; reusing existing pd files.")
    pd.big <- list()
    pd.big$tip.label <- read.csv(file.path(out_dir, "pd.taxon.name.csv"),
                                 row.names = 1, stringsAsFactors = FALSE)[, 1]
    pd.big$pd.wd <- out_dir
    pd.big$pd.file <- "pd.desc"
    pd.big$pd.name.file <- "pd.taxon.name.csv"
  }
  
  # ---- run icamp.big ----
  message("[iCAMP] Running icamp.big ...")
  icres <- iCAMP::icamp.big(
    comm = comm,
    pd.desc = pd.big$pd.file,
    pd.spname = pd.big$tip.label,
    pd.wd = pd.big$pd.wd,
    rand = rand.time,
    tree = tree,
    prefix = prefix,
    ds = ds,
    pd.cut = NA,
    sp.check = TRUE,
    phylo.rand.scale = "within.bin",
    taxa.rand.scale = "across.all",
    phylo.metric = "bMPD",
    sig.index = sig.index,
    bin.size.limit = bin.size.limit,
    nworker = nworker,
    memory.G = memory.G,
    rtree.save = FALSE,
    detail.save = TRUE,
    qp.save = FALSE,
    detail.null = FALSE,
    ignore.zero = TRUE,
    output.wd = out_dir,
    correct.special = TRUE,
    unit.sum = rowSums(comm),
    special.method = "depend",
    ses.cut = 1.96,
    rc.cut = 0.95,
    conf.cut = 0.975,
    omit.option = "no",
    meta.ab = NULL
  )
  
  # ---- save outputs ----
  rds_path <- file.path(out_dir, paste0(prefix, "_icamp_result.rds"))
  save_rds(icres, rds_path)
  
  # optional export
  if (!is.null(icres$CbMPDiCBraya)) {
    write.csv(icres$CbMPDiCBraya,
              file = file.path(out_dir, paste0(prefix, "_CbMPDiCBraya.csv")),
              row.names = FALSE)
  }
  
  message("[iCAMP] Saved: ", rds_path)
  message("[iCAMP] Time elapsed: ", format(Sys.time() - t0))
  
  invisible(list(icres = icres, comm = comm, tree = tree, clas = clas, treat = treat, env = env))
}
