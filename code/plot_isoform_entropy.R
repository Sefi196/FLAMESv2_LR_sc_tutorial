plot_isoform_entropy <- function(
    obj,
    genes,
    ident.name             = "BroadType",
    cell.type              = NULL,
    assay                  = "iso",
    slot.counts            = "counts",
    slot.data              = "data",
    reduction              = "umap",
    min_gene_counts        = 10,
    isoform_min_pct_cells  = 0.05,
    isoform_cumulative_pct = 0.95,
    alpha                  = 0.1
) {
  library(Seurat); library(dplyr); library(tidyr)
  library(ggplot2); library(ggbeeswarm); library(ggpubr)

  # 1) Set identities & optionally subset cells
  Idents(obj) <- ident.name
  cells <- if (is.null(cell.type)) {
    colnames(obj)
  } else {
    WhichCells(obj, ident = cell.type)
  }

  all_plots <- list()

  for (gene in genes) {
    # 2) Grab raw counts & isoform names
    mat_iso <- GetAssayData(obj, assay = assay, slot = slot.counts)[, cells, drop = FALSE]
    isoforms <- grep(paste0("(^|-|\\b)", gene, "($|\\b)"),
                     rownames(GetAssayData(obj, assay = assay, slot = slot.data)),
                     value = TRUE)
    if (length(isoforms) < 2) {
      message(gene, ": fewer than 2 isoforms found—skipping.")
      next
    }

    # 3) Subset to that gene’s isoforms × selected cells
    submat <- mat_iso[isoforms, , drop = FALSE]

    # 4) Filter out cells with low total gene UMI
    keep_cells <- colnames(submat)[colSums(submat) >= min_gene_counts]
    submat_f   <- submat[, keep_cells, drop = FALSE]
    if (ncol(submat_f) == 0) {
      message(gene, ": no cells pass min_gene_counts=", min_gene_counts, "—skipping.")
      next
    }

    # 5) Filter isoforms by proportion of cells expressing (same as isoform_min_pct_cells logic)
    pct_cells_expressed <- rowMeans(submat_f > 0)
    iso_keep_pct <- names(pct_cells_expressed)[pct_cells_expressed >= isoform_min_pct_cells]
    submat_f <- submat_f[iso_keep_pct, , drop = FALSE]
    if (nrow(submat_f) < 2) {
      message(gene, ": fewer than 2 isoforms remain after isoform_min_pct_cells; skipping.")
      next
    }

    # 6) Cumulative proportion filtering per cell and intersection
    filter_top_isoforms <- function(counts, threshold = 0.9) {
      proportions <- counts / sum(counts)
      sorted <- sort(proportions, decreasing = TRUE)
      keep <- which(cumsum(sorted) <= threshold)
      if (length(keep) == 0) keep <- 1
      names(sorted)[c(keep, length(keep) + 1)]
    }
    iso_keep_cum <- unique(unlist(apply(submat_f, 2, filter_top_isoforms, threshold = isoform_cumulative_pct)))
    iso_keep_final <- intersect(rownames(submat_f), iso_keep_cum)
    if (length(iso_keep_final) < 2) {
      message(gene, ": fewer than 2 isoforms remain after cumulative filtering; skipping.")
      next
    }
    submat_final <- submat_f[iso_keep_final, , drop = FALSE]

    # 7) Additive smoothing
    mat_sm <- submat_final + alpha
    p_mat  <- sweep(mat_sm, 2, colSums(mat_sm), FUN = "/")

    # 8) Entropy: raw and per-cell normalized
    entropy_raw <- apply(p_mat, 2, function(p) {
      valid <- p > 0
      -sum(p[valid] * log2(p[valid]))
    })
    n_active <- apply(p_mat, 2, function(p) sum(p > 0))
    entropy_norm <- entropy_raw / log2(n_active)

    # Safe fix: only replace NaN if raw entropy == 0
    nan_idx <- is.nan(entropy_norm)
    entropy_norm[nan_idx & entropy_raw == 0] <- 0


    # Summary stats
    median_raw <- median(entropy_raw, na.rm = TRUE)
    median_norm <- median(entropy_norm, na.rm = TRUE)  # on 0-1 scale
    max_H <- 1  # normalized entropy
    norm_med <- median_norm  # already scaled

    cat(glue::glue("
------ {gene} ------
Cells used: {length(entropy_norm)}
Isoforms used: {nrow(p_mat)}
Median raw H = {round(median_raw,3)}
Normalized median = {round(norm_med,2)} → {if(norm_med>0.5) 'HIGH' else 'LOW'} diversity
"))

    # 9a) Beeswarm plot with normalized entropy
    df <- tibble(
      gene    = gene,
      cell    = names(entropy_norm),
      entropy = unname(entropy_norm),
      cluster = Idents(obj)[names(entropy_norm)]
    )
    p_bees <- ggplot(df, aes(x = cluster, y = entropy)) +
      geom_quasirandom(alpha = 0.4, varwidth = TRUE) +
      geom_jitter(width = 0.2, height = 0, alpha = 0.4) +
      #geom_jitter(width = 0.25, height = 0.03, alpha = 0.3) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.6, fatten = 2, color = "red") +
      labs(
        title    = paste0("Entropy: ", gene,
                          if (!is.null(cell.type)) paste0(" (", cell.type, ")")),
        subtitle = paste0("Norm. median H = ", round(norm_med,2)),
        x        = "Cluster",
        y        = "Normalized Shannon entropy"
      ) +
      theme_minimal() +
      #ylim(-0.01,1)
      theme(
        axis.text.x   = element_text(angle = 45, hjust = 1),
        plot.subtitle = element_text(size = 10, color = "gray40")
      )

    # 9b) UMAP overlay with normalized entropy
    obj$.__tmp_entropy <- NA_real_
    obj$.__tmp_entropy[names(entropy_norm)] <- entropy_norm
    p_umap <- FeaturePlot(
      obj,
      features  = ".__tmp_entropy",
      reduction = reduction,
      cols      = c("lightgray","navy")
    ) +
      labs(title = paste0("UMAP: ", gene, " normalized entropy")) +
      scale_color_viridis_c(option = "C") +
      theme_minimal() +
      theme(legend.position = "right")
    obj$.__tmp_entropy <- NULL  # clean up

    all_plots[[gene]] <- list(
      beeswarm = p_bees,
      umap     = p_umap
    )
  }

  return(all_plots)
}
