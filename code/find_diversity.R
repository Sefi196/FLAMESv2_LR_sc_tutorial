find_diversity <- function(obj,
                           genes,
                           assay = "iso",
                           slot.counts = "counts",
                           slot.data = "data",
                           ident.name = NULL,
                           alpha = 0,
                           min_counts_per_cell = 10,
                           isoform_min_pct_cells = 0.05,
                           isoform_cumulative_pct = 0.95,
                           min_cell_fraction = 0.20,
                           p = NULL) {
  if (!is.null(ident.name)) {
    Idents(obj) <- ident.name
  }
  
  if (is.null(p)) {
    p <- function(x) { message(x) }  # simple fallback that just prints messages
  }
  
  cells <- colnames(obj)
  mat_counts <- GetAssayData(obj, assay = assay, slot = slot.counts)[, cells, drop = FALSE]
  mat_data   <- GetAssayData(obj, assay = assay, slot = slot.data)
  
  norm_entropy_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(cells))
  raw_entropy_mat  <- matrix(NA_real_, nrow = length(genes), ncol = length(cells))
  rownames(norm_entropy_mat) <- genes
  colnames(norm_entropy_mat) <- cells
  rownames(raw_entropy_mat) <- genes
  colnames(raw_entropy_mat) <- cells
  
  isoform_summary_list <- list()
  
  for (gene in genes) {
    p(sprintf("Processing %s", gene))
    
    # all isoforms matching gene pattern
    isoforms <- grep(paste0("(^|-|\\b)", gene, "($|\\b)"), rownames(mat_data), value = TRUE)
    total_isoforms <- length(isoforms)
    if (total_isoforms < 2) {
      message("Skipping ", gene, ": fewer than 2 isoforms (found ", total_isoforms, ")")
      next
    }
    
    submat <- mat_counts[isoforms, , drop = FALSE]
    gene_total_counts <- colSums(submat)
    keep_cells <- names(gene_total_counts)[gene_total_counts >= min_counts_per_cell]
    if (length(keep_cells) == 0) {
      message("Skipping ", gene, ": no cells with at least ", min_counts_per_cell, " gene counts")
      next
    }
    
    submat_f <- submat[, keep_cells, drop = FALSE]
    
    # filter by percent of cells expressing each isoform
    pct_cells_expressed <- rowMeans(submat_f > 0)
    iso_keep_pct <- names(pct_cells_expressed)[pct_cells_expressed >= isoform_min_pct_cells]
    isoforms_passing_pct <- iso_keep_pct
    submat_f <- submat_f[iso_keep_pct, , drop = FALSE]
    if (nrow(submat_f) == 0) {
      message("Skipping ", gene, ": no isoforms pass pct filter")
      next
    }
    
    # cumulative proportion filtering per cell, then intersect
    filter_top_isoforms <- function(counts, threshold = 0.9) { 
      proportions <- counts / sum(counts)
      sorted <- sort(proportions, decreasing = TRUE)
      keep <- which(cumsum(sorted) <= threshold)
      if (length(keep) == 0) keep <- 1
      names(sorted)[c(keep, length(keep) + 1)]
    }
    
    iso_keep_cum <- unique(unlist(apply(submat_f, 2, filter_top_isoforms, threshold = isoform_cumulative_pct)))
    iso_keep_final <- base::intersect(rownames(submat_f), iso_keep_cum)
    n_isoforms_used <- length(iso_keep_final)
    if (n_isoforms_used < 2) {
      message("Skipping ", gene, ": fewer than 2 isoforms remain after cumulative filtering (kept ", n_isoforms_used, ")")
      next
    }
    
    submat_final <- submat_f[iso_keep_final, , drop = FALSE]
    mat_smoothed <- submat_final + alpha
    p_mat <- sweep(mat_smoothed, 2, colSums(mat_smoothed), FUN = "/")  # proportions
    
    # raw entropy: -sum(p log2 p)
    entropy_raw <- vapply(seq_len(ncol(p_mat)), function(i) {
      p_vec <- p_mat[, i]
      valid <- p_vec > 0
      -sum(p_vec[valid] * log2(p_vec[valid]))
    }, numeric(1))
    
    # normalized: divide by log2(number of nonzero isoforms) per cell
    denom <- vapply(seq_len(ncol(p_mat)), function(i) {
      p_vec <- p_mat[, i]
      sum(p_vec > 0)
    }, numeric(1))
    norm_factor <- log2(denom)
    entropy_norm <- entropy_raw / norm_factor
    
    norm_entropy_mat[gene, colnames(p_mat)] <- entropy_norm
    raw_entropy_mat[gene, colnames(p_mat)]  <- entropy_raw
    
    # record cells used (where entropy was computed)
    cells_used <- colnames(p_mat)
    n_cells_used <- length(cells_used)
    
    isoform_summary_list[[gene]] <- list(
      gene = gene,
      isoforms_all = paste(isoforms, collapse = ";"),
      total_isoforms = total_isoforms,
      isoforms_passing_pct = paste(isoforms_passing_pct, collapse = ";"),
      n_isoforms_passing_pct = length(isoforms_passing_pct),
      isoforms_used = paste(iso_keep_final, collapse = ";"),
      n_isoforms_used = n_isoforms_used,
      cells_used = paste(cells_used, collapse = ";"),
      n_cells_used = n_cells_used
    )
  }
  
  # assemble summary df
  isoform_summary_df <- do.call(rbind, lapply(isoform_summary_list, function(x) {
    data.frame(
      gene = x$gene,
      isoforms_all = x$isoforms_all,
      total_isoforms = x$total_isoforms,
      isoforms_passing_pct = x$isoforms_passing_pct,
      n_isoforms_passing_pct = x$n_isoforms_passing_pct,
      isoforms_used = x$isoforms_used,
      n_isoforms_used = x$n_isoforms_used,
      cells_used = x$cells_used,
      n_cells_used = x$n_cells_used,
      stringsAsFactors = FALSE
    )
  }))
  rownames(isoform_summary_df) <- isoform_summary_df$gene
  
  # added this to deal with Nan values thata rise from 0 entropy values 
  nan_idx <- is.nan(norm_entropy_mat)
  norm_entropy_mat[nan_idx & raw_entropy_mat == 0] <- 0
  
  # filter genes by min_cell_fraction (based on normalized entropy presence)
  keep_genes <- rownames(norm_entropy_mat)[rowMeans(!is.na(norm_entropy_mat)) >= min_cell_fraction]
  norm_entropy_mat <- norm_entropy_mat[keep_genes, , drop = FALSE]
  raw_entropy_mat  <- raw_entropy_mat[keep_genes, , drop = FALSE]
  isoform_summary_df <- isoform_summary_df[keep_genes, , drop = FALSE]
  
  return(list(
    normalized_entropy = norm_entropy_mat,
    raw_entropy = raw_entropy_mat,
    isoform_summary = isoform_summary_df
  ))
}
