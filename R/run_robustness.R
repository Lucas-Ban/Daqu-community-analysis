# Robustness / node-removal simulation pipeline
# - Reproduces your original logic with config-driven paths
# - Parallelizes across MEN windows (1..N)
# - Saves one combined CSV (and optionally per-window CSVs)

`%||%` <- function(x, y) if (is.null(x)) y else x

rand.remov.once <- function(netRaw, rm.count, sp.ra, abundance.weighted = TRUE) {
  id.rm <- sample(seq_len(nrow(netRaw)), rm.count)
  
  net.Raw <- netRaw
  net.Raw[id.rm, ] <- 0
  net.Raw[, id.rm] <- 0
  
  if (abundance.weighted) {
    # multiply each column by species relative abundance (sp.ra aligned to columns)
    net.strength <- sweep(net.Raw, 2, sp.ra, `*`)
  } else {
    net.strength <- net.Raw
  }
  
  sp.meanInteraction <- colMeans(net.strength)
  id.rm2 <- which(sp.meanInteraction <= 0)
  remain.percent <- (nrow(netRaw) - length(id.rm2)) / nrow(netRaw)
  remain.percent
}

rmsimu <- function(netRaw, sp.ra, abundance.weighted = TRUE, nperm = 100) {
  results <- numeric(0)
  rm.count <- 1
  
  while (TRUE) {
    remains <- vapply(seq_len(nperm), function(m) {
      rand.remov.once(netRaw = netRaw, rm.count = rm.count, sp.ra = sp.ra,
                      abundance.weighted = abundance.weighted)
    }, numeric(1))
    
    remain.mean <- mean(remains)
    results[rm.count] <- remain.mean
    
    if (remain.mean == 0) break
    rm.count <- rm.count + 1
  }
  results
}

# Your original correlation logic (kept as-is)
# For each pair:
# - if both are 0 in a sample => NA (ignored)
# - else replace zeros with pseudocount (default 0.01), then cor(log(x), log(y))
calc_cor_matrix <- function(comm, pseudocount = 0.01) {
  p <- ncol(comm)
  cormatrix <- matrix(0, p, p)
  colnames(cormatrix) <- colnames(comm)
  rownames(cormatrix) <- colnames(comm)
  
  for (j in seq_len(p)) {
    for (k in j:p) {
      x <- comm[, j]
      y <- comm[, k]
      
      # samples where at least one is >0
      keep <- (x > 0) | (y > 0)
      if (!any(keep)) {
        corij <- NA_real_
      } else {
        x2 <- x[keep]
        y2 <- y[keep]
        x2[x2 <= 0] <- pseudocount
        y2[y2 <= 0] <- pseudocount
        corij <- suppressWarnings(cor(log(x2), log(y2)))
      }
      
      cormatrix[j, k] <- corij
      cormatrix[k, j] <- corij
    }
  }
  
  cormatrix[is.na(cormatrix)] <- 0
  diag(cormatrix) <- 0
  cormatrix
}

build_network <- function(comm, corr_cutoff = 0.80, pseudocount = 0.01) {
  # comm: samples x species
  corm <- calc_cor_matrix(comm, pseudocount = pseudocount)
  
  adj <- corm * (abs(corm) >= corr_cutoff)
  adj[is.na(adj)] <- 0
  diag(adj) <- 0
  
  # drop isolates
  keep <- colSums(abs(adj)) > 0
  adj2 <- adj[keep, keep, drop = FALSE]
  
  list(adj = adj2, keep_species = colnames(comm)[keep])
}

run_robustness_pipeline <- function(cfg) {
  # ---- packages ----
  pkgs <- c("foreach", "doParallel", "parallel")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop("Missing package(s): ", paste(miss, collapse = ", "),
         "\nInstall first: install.packages(c('foreach','doParallel'))")
  }
  
  # ---- config ----
  rb <- cfg$robustness
  if (is.null(rb)) stop("Missing 'robustness' section in config/config.yml")
  
  input_dir <- rb$input_dir
  selected_rows_file <- rb$selected_rows_file
  node_attribute_file <- rb$node_attribute_file
  out_dir <- rb$out_dir
  
  corr_cutoff <- rb$corr_cutoff %||% 0.80
  pseudocount <- rb$pseudocount %||% 0.01
  nperm <- rb$nperm %||% 100
  abundance_weighted <- rb$abundance_weighted %||% TRUE
  n_windows <- rb$n_windows %||% 76
  file_pattern <- rb$file_pattern %||% "{i}.txt"   # e.g., "1.txt"
  sep <- rb$sep %||% "\t"
  
  ensure_dir(out_dir)
  
  stop_if_missing(c(input_dir, selected_rows_file, node_attribute_file), label = "robustness input")
  
  selected_rows <- readLines(selected_rows_file, warn = FALSE)
  species_to_remove <- readLines(node_attribute_file, warn = FALSE)
  
  # ---- parallel ----
  numCores <- rb$nworker %||% max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  start_time <- Sys.time()
  message("[run_robustness] windows=", n_windows, " cores=", numCores,
          " cutoff=", corr_cutoff, " nperm=", nperm)
  
  # ---- loop across MEN windows ----
  # IMPORTANT: return a data.frame per window; combine with rbind
  combined_results <- foreach::foreach(i = seq_len(n_windows),
                                       .combine = rbind,
                                       .packages = character(0),
                                       .errorhandling = "stop") %dopar% {
                                         
                                         # build file path
                                         file_name <- gsub("\\{i\\}", as.character(i), file_pattern)
                                         file_path <- file.path(input_dir, file_name)
                                         if (!file.exists(file_path)) {
                                           stop("Missing OTU table: ", file_path)
                                         }
                                         
                                         otutab <- read.table(file_path, header = TRUE, row.names = 1, sep = sep,
                                                              check.names = FALSE, comment.char = "", quote = "")
                                         otutab[is.na(otutab)] <- 0
                                         
                                         # keep selected rows (species)
                                         otutab <- otutab[rownames(otutab) %in% selected_rows, , drop = FALSE]
                                         
                                         # remove target species (all at once, consistent with your current code)
                                         otutab2 <- otutab[!rownames(otutab) %in% species_to_remove, , drop = FALSE]
                                         
                                         # transpose to samples x species
                                         comm <- t(otutab2)
                                         
                                         # species relative abundance (aligned to columns)
                                         total_counts <- colSums(comm)
                                         sp.ra <- total_counts / sum(total_counts)
                                         
                                         # build network adjacency
                                         net <- build_network(comm = comm,
                                                              corr_cutoff = corr_cutoff,
                                                              pseudocount = pseudocount)
                                         
                                         # align sp.ra to kept species
                                         if (nrow(net$adj) == 0) {
                                           # no network left after filtering
                                           return(data.frame(
                                             file = i,
                                             weighted = c("weighted", "unweighted"),
                                             Proportion.removed = 1,
                                             remain.mean = 0,
                                             removed_species = paste(species_to_remove, collapse = ", "),
                                             stringsAsFactors = FALSE
                                           ))
                                         }
                                         
                                         keep_sp <- rownames(net$adj)
                                         sp.ra2 <- sp.ra[keep_sp]
                                         sp.ra2[is.na(sp.ra2)] <- 0
                                         
                                         # simulate
                                         Weighted.simu <- rmsimu(netRaw = net$adj, sp.ra = sp.ra2, abundance.weighted = TRUE, nperm = nperm)
                                         Unweighted.simu <- rmsimu(netRaw = net$adj, sp.ra = sp.ra2, abundance.weighted = FALSE, nperm = nperm)
                                         
                                         dat <- data.frame(
                                           file = rep(i, length(Weighted.simu) * 2),
                                           weighted = rep(c("weighted", "unweighted"), each = length(Weighted.simu)),
                                           Proportion.removed = rep(seq_along(Weighted.simu), 2),
                                           remain.mean = c(Weighted.simu, Unweighted.simu),
                                           removed_species = paste(species_to_remove, collapse = ", "),
                                           stringsAsFactors = FALSE
                                         )
                                         
                                         dat
                                       }
  
  # ---- save outputs ----
  out_csv <- file.path(out_dir, "robustness_removed_species_results.csv")
  write.csv(combined_results, out_csv, row.names = FALSE)
  
  message("[run_robustness] Saved: ", out_csv)
  message("[run_robustness] Time elapsed: ", format(Sys.time() - start_time))
  
  invisible(combined_results)
}
