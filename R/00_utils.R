# Utilities: config, IO checks, small helpers

read_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install it with: install.packages('yaml')")
  }
  if (!file.exists(path)) {
    stop("Config file not found: ", path)
  }
  yaml::read_yaml(path)
}

ensure_dir <- function(path) {
  if (is.null(path) || identical(path, "") || identical(path, "null")) return(invisible(FALSE))
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

stop_if_missing <- function(paths, label = "file") {
  paths <- paths[!is.null(paths)]
  paths <- paths[!identical(paths, "null")]
  paths <- paths[!is.na(paths)]
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing ", label, "(s):\n- ", paste(missing, collapse = "\n- "))
  }
}

read_tsv_table <- function(path, row_names = TRUE) {
  read.table(
    file = path,
    header = TRUE,
    sep = "\t",
    row.names = if (row_names) 1 else NULL,
    as.is = TRUE,
    stringsAsFactors = FALSE,
    comment.char = "",
    check.names = FALSE,
    quote = ""
  )
}

save_rds <- function(obj, path) {
  ensure_dir(dirname(path))
  saveRDS(obj, file = path)
}
