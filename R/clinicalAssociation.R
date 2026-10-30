# Standalone R implementation of clinicalAssociation.cc.py
#
# For each gene signature:
#   1. Per-gene z-score computed across samples that have a given clinical variable
#   2. Activity score = sum(z-scores of up-regulated genes) - sum(z-scores of down-regulated genes)
#   3. Numerical variables  -> Pearson correlation + clustered heatmap
#   4. Categorical variables -> Kruskal-Wallis + BH FDR + clustered heatmap
#
# Required packages: pheatmap, ggplot2
# Install if needed: install.packages(c("pheatmap", "ggplot2"))

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
})

# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

NA_STRINGS   <- c("", "NA", "[Not Available]", "[Not Applicable]", "[Unknown]", "Unknown")
PVALUE_FLOOR <- 1e-100   # clamp before -log10 to avoid Inf

# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS: DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

# Load all signature files matching the glob pattern.
# Signature file format (tab-separated):
#   Line 1  : signature name
#   Lines 2+: GENE_SYMBOL <TAB> log2FC  (positive = up, negative = down)
# Returns a named list; each element has $name, $up_genes, $down_genes.
load_signatures <- function(sig_folder, sig_pattern) {
  sig_files <- Sys.glob(file.path(sig_folder, sig_pattern))
  if (length(sig_files) == 0)
    stop("No signature files found matching: ", file.path(sig_folder, sig_pattern))

  message(sprintf("[Signatures] Found %d file(s)", length(sig_files)))
  signatures <- list()

  for (f in sig_files) {
    lines      <- readLines(f, warn = FALSE)
    lines      <- lines[nchar(trimws(lines)) > 0]
    sig_name   <- trimws(lines[1])
    up_genes   <- character(0)
    down_genes <- character(0)

    for (gl in lines[-1]) {
      parts  <- strsplit(gl, "\t")[[1]]
      if (length(parts) != 2) { warning("Malformed line in ", basename(f), ": '", gl, "' — skipping"); next }
      gene   <- toupper(trimws(parts[1]))
      log2fc <- suppressWarnings(as.numeric(parts[2]))
      if (is.na(log2fc)) next
      if (log2fc > 0) up_genes   <- c(up_genes,   gene)
      else            down_genes <- c(down_genes,  gene)
    }

    signatures[[sig_name]] <- list(name = sig_name, up_genes = up_genes, down_genes = down_genes)
    message(sprintf("  '%s': %d up-genes, %d down-genes", sig_name, length(up_genes), length(down_genes)))
  }
  return(signatures)
}

# Load expression CSV.
# First column = gene symbols (becomes rownames); remaining columns = samples.
# Returns a numeric matrix (genes x samples).
load_expression <- function(expr_file) {
  df  <- data.table::fread(expr_file, header = TRUE, data.table = FALSE,
                           check.names = FALSE)
  genes <- as.character(df[[1]])
  df[[1]] <- NULL
  mat <- as.matrix(df)
  rm(df)
  mode(mat) <- "numeric"
  rownames(mat) <- toupper(genes)
  message(sprintf("[Expression] %d genes x %d samples", nrow(mat), ncol(mat)))
  return(mat)
}

# Detect whether a character vector is numeric (after removing NA-like strings).
is_numeric_col <- function(values) {
  clean <- as.character(values)
  clean <- clean[!is.na(clean) & !(clean %in% NA_STRINGS)]
  if (length(clean) == 0) return(FALSE)
  !any(is.na(suppressWarnings(as.numeric(clean))))
}

# Load clinical CSV and auto-split into numerical and categorical data frames.
# First column must be sample IDs.
# Returns list($sample_col, $numerical, $categorical) — either df or NULL.
load_and_split_clinical <- function(clinical_file) {
  df         <- read.csv(clinical_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  sample_col <- colnames(df)[1]
  var_cols   <- colnames(df)[-1]

  is_num   <- sapply(var_cols, function(v) is_numeric_col(as.character(df[[v]])))
  num_vars <- var_cols[ is_num]
  cat_vars <- var_cols[!is_num]

  message(sprintf("[Clinical] %d samples | %d numerical, %d categorical variable(s)",
                  nrow(df), length(num_vars), length(cat_vars)))
  if (length(num_vars) > 0) message("  Numerical:   ", paste(num_vars, collapse = ", "))
  if (length(cat_vars) > 0) message("  Categorical: ", paste(cat_vars, collapse = ", "))

  rownames(df) <- df[[sample_col]]
  list(
    sample_col  = sample_col,
    numerical   = if (length(num_vars) > 0) df[, c(sample_col, num_vars), drop = FALSE] else NULL,
    categorical = if (length(cat_vars) > 0) df[, c(sample_col, cat_vars), drop = FALSE] else NULL
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS: Z-SCORE AND ACTIVITY SCORE
# ══════════════════════════════════════════════════════════════════════════════

# Per-gene z-score across samples.
# Genes with sd == 0 are set to 0 (matches Python behaviour).
compute_zscore_per_gene <- function(expr_mat) {
  zscore_mat <- t(apply(expr_mat, 1, function(row) {
    s <- sd(row, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(row)))
    (row - mean(row, na.rm = TRUE)) / s
  }))
  dimnames(zscore_mat) <- dimnames(expr_mat)
  return(zscore_mat)
}

# Compute signature activity scores across all samples.
# activity(sample) = Σ z(up-gene) − Σ z(down-gene)
# Returns matrix: signatures x samples.
compute_activity_scores <- function(zscore_mat, signatures) {
  genes_in_mat <- rownames(zscore_mat)
  sig_names    <- names(signatures)
  activity_mat <- matrix(0,
                         nrow     = length(sig_names),
                         ncol     = ncol(zscore_mat),
                         dimnames = list(sig_names, colnames(zscore_mat)))
  for (sig_name in sig_names) {
    sig        <- signatures[[sig_name]]
    up_found   <- intersect(sig$up_genes,   genes_in_mat)
    down_found <- intersect(sig$down_genes, genes_in_mat)
    if (length(up_found)   > 0)
      activity_mat[sig_name, ] <- activity_mat[sig_name, ] + colSums(zscore_mat[up_found,   , drop = FALSE])
    if (length(down_found) > 0)
      activity_mat[sig_name, ] <- activity_mat[sig_name, ] - colSums(zscore_mat[down_found, , drop = FALSE])
  }
  return(activity_mat)
}

# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS: STATISTICS
# ══════════════════════════════════════════════════════════════════════════════

# Pearson correlation between each signature's activity scores and a numeric vector.
# Returns list($all_pcc, $sig_pcc, $neg_log10_pval).
correlate_with_numerical <- function(activity_mat, var_values) {
  sig_names      <- rownames(activity_mat)
  all_pcc        <- setNames(numeric(length(sig_names)), sig_names)
  sig_pcc        <- setNames(numeric(length(sig_names)), sig_names)
  neg_log10_pval <- setNames(numeric(length(sig_names)), sig_names)

  for (sig in sig_names) {
    scores <- activity_mat[sig, ]
    if (length(unique(scores)) == 1 || length(unique(var_values)) == 1) next
    ct                  <- cor.test(var_values, scores, method = "pearson")
    pcc                 <- as.numeric(ct$estimate)
    all_pcc[sig]        <- pcc
    neg_log10_pval[sig] <- -log10(max(ct$p.value, PVALUE_FLOOR))
    sig_pcc[sig]        <- if (ct$p.value < 0.05) pcc else 0
  }
  list(all_pcc = all_pcc, sig_pcc = sig_pcc, neg_log10_pval = neg_log10_pval)
}

# Kruskal-Wallis test per signature + BH FDR across all signatures.
# Returns list($pvalues, $qvalues).
associate_with_categorical <- function(activity_mat, group_labels) {
  sig_names <- rownames(activity_mat)
  pvalues   <- setNames(rep(1, length(sig_names)), sig_names)
  for (sig in sig_names) {
    tryCatch({
      kw           <- kruskal.test(activity_mat[sig, ] ~ factor(group_labels))
      pvalues[sig] <- kw$p.value
    }, error = function(e) { pvalues[sig] <<- 1 })
  }
  list(pvalues = pvalues, qvalues = p.adjust(pvalues, method = "BH"))
}

# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS: PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

# Save a clustered heatmap PDF.
# Clustering is disabled when fewer than 2 rows or columns (mirrors Python).
save_heatmap <- function(mat, pdf_file, fig_size, color_pal, breaks, title = "") {
  mat        <- as.matrix(mat)
  do_cluster <- nrow(mat) >= 2 && ncol(mat) >= 2
  pdf(pdf_file, width = fig_size, height = fig_size)
  # silent=TRUE suppresses auto-drawing; capture and explicitly print to device
  ph <- pheatmap(mat,
                 color                    = color_pal,
                 breaks                   = breaks,
                 cluster_rows             = do_cluster,
                 cluster_cols             = do_cluster,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method        = "complete",
                 main                     = title,
                 fontsize                 = 8,
                 border_color             = "black",
                 silent                   = TRUE)
  print(ph)
  dev.off()
  message("  Saved: ", pdf_file)
}

# One-page-per-signature scatterplots: activity score vs numerical variable.
plot_numerical_support <- function(activity_mat, var_values, var_name, pdf_file,
                                   pcc_vec, neg_log10_pval_vec) {
  pdf(pdf_file, width = 8, height = 8)
  for (sig in rownames(activity_mat)) {
    df_plot <- data.frame(x = var_values, y = activity_mat[sig, ])
    p <- ggplot(df_plot, aes(x = x, y = y)) +
      geom_point(color = "red", alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7, formula = y ~ x) +
      labs(x     = paste("Clinical variable:", var_name),
           y     = paste(sig, "activity score"),
           title = sprintf("r = %.3f   -log10(p) = %.2f", pcc_vec[sig], neg_log10_pval_vec[sig])) +
      theme_classic(base_size = 11)
    print(p)
  }
  dev.off()
  message("  Saved: ", pdf_file)
}

# One-page-per-signature boxplots: activity score by categorical group.
plot_categorical_support <- function(activity_mat, group_labels, var_name, pdf_file, kw_pvalues) {
  pdf(pdf_file, width = 8, height = 8)
  for (sig in rownames(activity_mat)) {
    df_plot <- data.frame(group = factor(group_labels), score = activity_mat[sig, ])
    kw_p    <- signif(kw_pvalues[sig], 3)
    anova_p <- tryCatch(
      signif(summary(aov(score ~ group, data = df_plot))[[1]][["Pr(>F)"]][1], 3),
      error = function(e) NA
    )
    p <- ggplot(df_plot, aes(x = group, y = score)) +
      geom_boxplot(outlier.shape = NA, fill = "steelblue", alpha = 0.6) +
      labs(x     = paste("Clinical variable:", var_name),
           y     = paste(sig, "activity score"),
           title = sprintf("ANOVA p = %s   K-W p = %s", anova_p, kw_p)) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    print(p)
  }
  dev.off()
  message("  Saved: ", pdf_file)
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN FUNCTION
# ══════════════════════════════════════════════════════════════════════════════
#
# Arguments:
#   sig_folder    Path to folder containing signature files
#   sig_pattern   Glob pattern within sig_folder (default "*.txt")
#   expression    Path to expression CSV (genes x samples, first col = gene symbols)
#   clinical      Path to clinical CSV (samples x variables, first col = sample IDs)
#   output_root   File path prefix for all outputs  (e.g. "results/my_project")
#   sample_list   Optional path to a text file of sample IDs to include (one per line)
#   min_samples   Minimum valid samples per clinical variable (default 10)
#   support_plots Generate per-signature scatterplots / boxplots (default FALSE)
#   figure_size   Width and height of summary heatmaps in inches (default 20)
#
run_clinical_association <- function(
    sig_folder,
    sig_pattern   = "*.txt",
    expression,
    clinical,
    output_root,
    sample_list   = NULL,
    min_samples   = 6,
    support_plots = FALSE,
    figure_size   = 20,
    save_temp     = TRUE
) {
  # ── Output directory setup (mirrors Python script structure) ─────────────────
  # output_root is a file-path prefix, e.g. "results/my_project"
  #   Main results     → results/
  #   pdf-numerical.my_project/   per-variable scatterplots
  #   pdf-categorical.my_project/ per-variable boxplots
  #   support-clinical.my_project/ per-variable activity + clinical joined CSVs
  #   temp-files.my_project/      intermediate summaries (all PCCs, q-values, filtered expr)
  out_dir   <- dirname(normalizePath(output_root, mustWork = FALSE))
  base_name <- basename(output_root)

  dir_pdf_num <- file.path(out_dir, paste0("pdf-numerical.",    base_name))
  dir_pdf_cat <- file.path(out_dir, paste0("pdf-categorical.",  base_name))
  dir_support <- file.path(out_dir, paste0("support-clinical.", base_name))
  dir_temp    <- file.path(out_dir, paste0("temp-files.",       base_name))

  dirs_to_create <- c(out_dir, dir_pdf_num, dir_pdf_cat, dir_support)
  if (isTRUE(save_temp)) dirs_to_create <- c(dirs_to_create, dir_temp)
  for (d in dirs_to_create)
    dir.create(d, showWarnings = FALSE, recursive = TRUE)

  # ── Load all inputs ──────────────────────────────────────────────────────────
  signatures <- load_signatures(sig_folder, sig_pattern)
  expr_mat   <- load_expression(expression)
  clinical   <- load_and_split_clinical(clinical)

  all_sig_genes <- unique(toupper(c(
    unlist(lapply(signatures, `[[`, "up_genes")),
    unlist(lapply(signatures, `[[`, "down_genes"))
  )))

  # ── Filter expression to signature genes only ────────────────────────────────
  genes_found <- intersect(all_sig_genes, rownames(expr_mat))
  message(sprintf("[Filter] %d/%d signature genes found in expression matrix",
                  length(genes_found), length(all_sig_genes)))
  if (length(genes_found) == 0) stop("No signature genes found in expression matrix.")
  expr_mat <- expr_mat[genes_found, , drop = FALSE]

  # Save filtered expression matrix to temp-files
  if (isTRUE(save_temp))
    write.csv(as.data.frame(expr_mat),
              file.path(dir_temp, paste0(base_name, ".expr.consistent.csv")),
              quote = FALSE)

  # ── Optional sample whitelist ────────────────────────────────────────────────
  if (!is.null(sample_list)) {
    keep     <- trimws(readLines(sample_list, warn = FALSE))
    keep     <- intersect(keep, colnames(expr_mat))
    message(sprintf("[Samples] Using %d samples from sample list", length(keep)))
    expr_mat <- expr_mat[, keep, drop = FALSE]
  }

  # Color palettes
  bwr_pal      <- colorRampPalette(c("blue", "white", "red"))(101)
  bwr_breaks   <- seq(-1, 1, length.out = 102)
  # First color slot = grey for 0; remaining 100 slots = white (minimal) → red (max)
  rainbow_pal  <- c("grey80", colorRampPalette(c("white", "red"))(100))
  rainbow_brks <- seq(0, 4, length.out = 102)

  # ══════════════════════════════════════════════════════════════════════════════
  # NUMERICAL VARIABLES
  # ══════════════════════════════════════════════════════════════════════════════
  if (!is.null(clinical$numerical)) {
    message("\n── Numerical variables ──────────────────────────────────────────────────")
    num_df       <- clinical$numerical
    sample_col   <- clinical$sample_col
    num_var_cols <- setdiff(colnames(num_df), sample_col)

    all_pcc_list  <- list()
    sig_pcc_list  <- list()
    pval_list     <- list()
    raw_pval_list <- list()   # raw p-values for temp-files

    for (var_name in num_var_cols) {
      var_data             <- num_df[, c(sample_col, var_name), drop = FALSE]
      var_data[[var_name]] <- as.character(var_data[[var_name]])
      var_data             <- var_data[!(var_data[[var_name]] %in% NA_STRINGS) &
                                         !is.na(var_data[[var_name]]), , drop = FALSE]
      rownames(var_data)   <- var_data[[sample_col]]
      valid_samples        <- intersect(rownames(var_data), colnames(expr_mat))

      if (length(valid_samples) < min_samples) {
        message(sprintf("  SKIP '%s': %d valid samples (min = %d)", var_name, length(valid_samples), min_samples))
        next
      }
      message(sprintf("  '%s' — %d samples", var_name, length(valid_samples)))

      # Z-scores recomputed on this variable's sample subset
      expr_sub   <- expr_mat[, valid_samples, drop = FALSE]
      zscore_sub <- compute_zscore_per_gene(expr_sub)
      activity   <- compute_activity_scores(zscore_sub, signatures)

      var_values <- as.numeric(var_data[valid_samples, var_name])
      result     <- correlate_with_numerical(activity, var_values)

      all_pcc_list[[var_name]]  <- result$all_pcc
      sig_pcc_list[[var_name]]  <- result$sig_pcc
      pval_list[[var_name]]     <- result$neg_log10_pval
      raw_pval_list[[var_name]] <- 10^(-result$neg_log10_pval)   # back-transform for temp CSV

      # temp-files: per-variable activity score matrix (signatures x samples)
      if (isTRUE(save_temp))
        write.csv(as.data.frame(activity),
                  file.path(dir_temp, paste0(base_name, ".Numerical_variable.", make.names(var_name),
                                             ".activity_score_all.csv")),
                  quote = FALSE)

      # support-clinical: activity scores joined with clinical variable (samples x cols)
      ac_df        <- as.data.frame(t(activity))
      ac_df[[var_name]] <- var_values
      write.csv(ac_df,
                file.path(dir_support, paste0(base_name, ".Numerical_variable.", make.names(var_name),
                                              ".activity_and_clinvar.csv")),
                quote = FALSE)

      if (isTRUE(support_plots))
        plot_numerical_support(
          activity, var_values, var_name,
          file.path(dir_pdf_num, paste0(base_name, ".numerical.", make.names(var_name), ".scatterplots.pdf")),
          result$all_pcc, result$neg_log10_pval
        )
    }

    if (length(all_pcc_list) > 0) {
      all_pcc_mat  <- do.call(rbind, all_pcc_list)   # variables x signatures
      sig_pcc_mat  <- do.call(rbind, sig_pcc_list)
      raw_pval_mat <- do.call(rbind, raw_pval_list)

      # temp-files: all PCC, raw p-values, all-correlations heatmap
      if (isTRUE(save_temp)) {
        write.csv(all_pcc_mat,
                  file.path(dir_temp, paste0(base_name, ".numerical_var_correlations_all_pcc.csv")),
                  quote = FALSE)
        write.csv(raw_pval_mat,
                  file.path(dir_temp, paste0(base_name, ".numerical_var_correlations_pcc_pvalue.csv")),
                  quote = FALSE)
        save_heatmap(all_pcc_mat,
                     file.path(dir_temp, paste0(base_name, ".numerical_var_all_correlations.euclidean.pdf")),
                     figure_size, bwr_pal, bwr_breaks,
                     title = "All Pearson correlations: signature activity vs numerical variables")
      }

      # Primary results (at output_root level)
      write.csv(sig_pcc_mat,
                file.path(out_dir, paste0(base_name, ".numerical_var_correlations_significant_pcc.csv")),
                quote = FALSE)
      save_heatmap(sig_pcc_mat,
                   file.path(out_dir, paste0(base_name, ".numerical_var_significant_correlations.euclidean.pdf")),
                   figure_size, bwr_pal, bwr_breaks,
                   title = "Significant Pearson correlations (p < 0.05; others = 0)")
    } else {
      message("  No numerical variables passed the minimum-sample filter.")
    }
  }

  # ══════════════════════════════════════════════════════════════════════════════
  # CATEGORICAL VARIABLES
  # ══════════════════════════════════════════════════════════════════════════════
  if (!is.null(clinical$categorical)) {
    message("\n── Categorical variables ────────────────────────────────────────────────")
    cat_df       <- clinical$categorical
    sample_col   <- clinical$sample_col
    cat_var_cols <- setdiff(colnames(cat_df), sample_col)

    pvalue_list <- list()
    qvalue_list <- list()

    for (var_name in cat_var_cols) {
      var_data             <- cat_df[, c(sample_col, var_name), drop = FALSE]
      var_data[[var_name]] <- as.character(var_data[[var_name]])
      var_data             <- var_data[!(var_data[[var_name]] %in% NA_STRINGS) &
                                         !is.na(var_data[[var_name]]), , drop = FALSE]
      rownames(var_data)   <- var_data[[sample_col]]
      valid_samples        <- intersect(rownames(var_data), colnames(expr_mat))

      if (length(valid_samples) < min_samples) {
        message(sprintf("  SKIP '%s': %d valid samples (min = %d)", var_name, length(valid_samples), min_samples))
        next
      }

      group_labels <- var_data[valid_samples, var_name]
      n_groups     <- length(unique(group_labels))
      if (n_groups < 2) {
        message(sprintf("  SKIP '%s': only %d group(s)", var_name, n_groups))
        next
      }
      message(sprintf("  '%s' — %d samples, %d groups", var_name, length(valid_samples), n_groups))

      expr_sub   <- expr_mat[, valid_samples, drop = FALSE]
      zscore_sub <- compute_zscore_per_gene(expr_sub)
      activity   <- compute_activity_scores(zscore_sub, signatures)
      result     <- associate_with_categorical(activity, group_labels)

      pvalue_list[[var_name]] <- -log10(pmax(result$pvalues, PVALUE_FLOOR))
      qvalue_list[[var_name]] <- -log10(pmax(result$qvalues, PVALUE_FLOOR))

      # temp-files: per-variable activity score matrix (signatures x samples)
      if (isTRUE(save_temp))
        write.csv(as.data.frame(activity),
                  file.path(dir_temp, paste0(base_name, ".Categorical_variable.", make.names(var_name),
                                             ".activity_score_all.csv")),
                  quote = FALSE)

      # support-clinical: activity scores joined with group label (samples x cols)
      ac_df             <- as.data.frame(t(activity))
      ac_df[[var_name]] <- group_labels
      write.csv(ac_df,
                file.path(dir_support, paste0(base_name, ".Categorical_variable.", make.names(var_name),
                                              ".activity_and_clinvar.csv")),
                quote = FALSE)

      if (isTRUE(support_plots))
        plot_categorical_support(
          activity, group_labels, var_name,
          file.path(dir_pdf_cat, paste0(base_name, ".categorical.", make.names(var_name), ".boxplots.pdf")),
          result$pvalues
        )
    }

    if (length(pvalue_list) > 0) {
      pval_mat <- do.call(rbind, pvalue_list)
      qval_mat <- do.call(rbind, qvalue_list)

      # Primary results (at output_root level)
      write.csv(pval_mat,
                file.path(out_dir, paste0(base_name, ".categorical_variables.minusLog10_pvalues.csv")),
                quote = FALSE)
      save_heatmap(pval_mat,
                   file.path(out_dir, paste0(base_name, ".categorical_variables_pvalues.euclidean.pdf")),
                   figure_size, rainbow_pal, rainbow_brks,
                   title = "Kruskal-Wallis -log10(p-value): signature activity vs categorical variables")

      # temp-files: q-value CSV and heatmap
      if (isTRUE(save_temp)) {
        write.csv(qval_mat,
                  file.path(dir_temp, paste0(base_name, ".categorical_variables.minusLog10_qvalues.csv")),
                  quote = FALSE)
        save_heatmap(qval_mat,
                     file.path(dir_temp, paste0(base_name, ".categorical_variables_qvalues.euclidean.pdf")),
                     figure_size, rainbow_pal, rainbow_brks,
                     title = "Kruskal-Wallis -log10(BH q-value): signature activity vs categorical variables")
      }
    } else {
      message("  No categorical variables passed the minimum-sample and group filters.")
    }
  }

  message(sprintf("\nDone. Outputs written to:"))
  message(sprintf("  Main results : %s/", out_dir))
  message(sprintf("  PDFs (num)   : %s/", dir_pdf_num))
  message(sprintf("  PDFs (cat)   : %s/", dir_pdf_cat))
  message(sprintf("  Support files: %s/", dir_support))
  if (isTRUE(save_temp))
    message(sprintf("  Temp files   : %s/", dir_temp))
}

# ══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLE — edit paths and run in RStudio
# ══════════════════════════════════════════════════════════════════════════════

# source("clinicalAssociation.R")
#
# run_clinical_association(
#   sig_folder    = "path/to/signatures",
#   sig_pattern   = "*.txt",
#   expression    = "path/to/expression.csv",
#   clinical      = "path/to/clinical.csv",
#   output_root   = "results/my_project",
#   sample_list   = NULL,          # or "path/to/samples.txt"
#   min_samples   = 10,
#   support_plots = TRUE,
#   figure_size   = 20
# )
