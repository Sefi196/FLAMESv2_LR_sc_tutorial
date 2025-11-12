plotIsoformPieCells <- function(
    seurat_obj,
    gene,
    assay = "iso",
    slot = "counts",
    top_n = 3,
    n_cells = 12,
    plot_col=5
) {
  # Load required packages
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("scales")
  requireNamespace("stringr")

  # 1. Extract count matrix and identify matching isoforms
  mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  isoform_names <- grep(paste0("(^|-|\\b)", gene, "($|\\b)"), rownames(mat), value = TRUE)

  if (length(isoform_names) < 2) {
    stop("Fewer than 2 isoforms found for this gene.")
  }

  # 2. Subset matrix and filter to expressing cells
  mat_gene <- mat[isoform_names, , drop = FALSE]
  gene_counts <- Matrix::colSums(mat_gene)
  expressing_cells <- names(gene_counts[gene_counts > 0])
  mat_gene <- mat_gene[, expressing_cells, drop = FALSE]

  if (ncol(mat_gene) == 0) {
    stop("No cells express this gene.")
  }

  # 3. Convert to long format and normalize
  df_long <- as.data.frame(as.table(as.matrix(mat_gene))) %>%
    dplyr::rename(
      transcript_id = Var1,
      cell          = Var2,
      expression    = Freq
    ) %>%
    dplyr::filter(expression > 0)

  # 4. Compute proportions and clean isoform names
  df_long <- df_long %>%
    dplyr::group_by(cell) %>%
    dplyr::mutate(
      proportion = expression / sum(expression),
      stripped_id = sub("\\..*$", "", transcript_id)
    ) %>%
    dplyr::ungroup()

  # 5. Choose top N cells based on total counts
  top_cells <- df_long %>%
    dplyr::group_by(cell) %>%
    dplyr::summarise(total_expression = sum(expression)) %>%
    dplyr::arrange(dplyr::desc(total_expression)) %>%
    dplyr::slice_head(n = n_cells)

  df_top <- df_long %>% dplyr::filter(cell %in% top_cells$cell)

  # 6. Identify top N isoforms by max proportion across cells
  top_isoforms <- df_top %>%
    dplyr::group_by(stripped_id) %>%
    dplyr::summarise(peak_prop = max(proportion), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(peak_prop)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(stripped_id)

  # 7. Label isoforms as top or "Other isoforms"
  df_top <- df_top %>%
    dplyr::mutate(
      isoform_group = ifelse(stripped_id %in% top_isoforms, stripped_id, "Other isoforms")
    )

  # 8. Aggregate again by cell + isoform_group
  df_summed <- df_top %>%
    dplyr::group_by(cell, isoform_group) %>%
    dplyr::summarise(
      proportion = sum(proportion),
      .groups    = "drop"
    )

  # Add total counts to each cell's label
  cell_labels <- top_cells %>%
    dplyr::mutate(
      cell_label = paste0(cell, "\n(counts = ", scales::comma(total_expression), ")")
    )

  df_summed <- df_summed %>%
    dplyr::left_join(cell_labels, by = "cell")

  # 9. Set color palette
  hues <- scales::hue_pal()(length(top_isoforms))
  color_map <- c(stats::setNames(hues, sort(top_isoforms)), "Other isoforms" = "grey50")

  # — NOW SWAP the blue & green entries —
  color_map[c("ENST00000422006", "ENST00000417292")] <-
  color_map[c("ENST00000417292", "ENST00000422006")]


  color_map["ENST00000417292"] <- "#C77CFF"

  # 10. Plot
  p <- ggplot2::ggplot(df_summed, ggplot2::aes(
    x = "", y = proportion, fill = isoform_group
  )) +
    ggplot2::geom_col(color = "black", size = 0.2, width = 1) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~cell_label, ncol = plot_col, strip.position = "bottom") +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(proportion > 0.05, paste0(round(proportion * 100, 1), "%"), "")
      ),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white",
      size = 3
    ) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::labs(
      title = paste0("Isoform Proportions per Cell for ", gene),
      fill = "Isoform"
    ) +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      panel.spacing = ggplot2::unit(1, "lines")
    )

  print(p)

  p2 <-   # Now build *just* the “no strip labels” version:
    p <- ggplot2::ggplot(df_summed, ggplot2::aes(x = "", y = proportion, fill = isoform_group)) +
    ggplot2::geom_col(color = "black", size = 0.2, width = 1) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(proportion > 0.05,
                                  paste0(round(proportion * 100, 1), "%"),
                                  "")
      ),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white", size = 3
    ) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::labs(
      title = paste0("Isoform Proportions per Cell for ", gene),
      fill  = "Isoform"
    ) +
    ggplot2::facet_wrap(~cell, ncol = plot_col) +  # still one pie per cell
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      strip.text       = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),

      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(face = "bold", size = 10),
      legend.text      = ggplot2::element_text(size = 8),

      panel.spacing    = ggplot2::unit(0.5, "lines"),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    )

  #print(p2)
}
