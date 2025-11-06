plot_isoforms_exp_umap <- function(obj, gene, reduction = "umap.harm", slot = "data", n_features = NULL, assay="iso", ncol=3) {
  DefaultAssay(obj) <- assay

  # Ensure layers are joined
  obj <- JoinLayers(obj, overwrite = TRUE)

  # Get expression matrix
  expression_matrix <- GetAssayData(obj, assay = assay, slot = slot)

  # Find features matching the gene name
  matching_features <- grep(paste0("(^|-|\\b)", gene, "($|\\b)"),
                            rownames(expression_matrix), value = TRUE)

  if (length(matching_features) == 0) {
    warning("No matching features found for gene: ", gene)
    return(NULL)
  }

  # Subset and rank by total expression
  subset_expression <- expression_matrix[matching_features, , drop = FALSE]
  total_expression <- Matrix::rowSums(subset_expression)
  top_features <- names(sort(total_expression, decreasing = TRUE))

  # Limit number of features if specified
  if (!is.null(n_features)) {
    top_features <- head(top_features, n_features)
  }

  # Print feature ranking
  print(data.frame(Feature = top_features, Expression = total_expression[top_features]))

  # Plot features
  FeaturePlot(obj, features = top_features, reduction = reduction, ncol = ncol)
}



library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

plot_isoforms_exp_umap_publication <- function(
    obj,
    gene,
    assay = "iso",
    slot = "data",
    reduction = "umap",
    n_features = NULL,
    title = NULL,
    ncol = 2,
    pt.size = 1.2,
    show_legend = TRUE,
    color_mode = c("seurat", "viridis", "custom"),
    viridis_option = "inferno",
    clip_quant = 0.99,
    gene_assay = "RNA",
    gene_slot = "data",
    include_gene_plot = TRUE
) {
  color_mode <- match.arg(color_mode)

  require(Seurat)
  require(Matrix)
  require(ggplot2)
  require(patchwork)
  if (color_mode == "viridis") require(viridisLite)

  DefaultAssay(obj) <- assay
  obj <- JoinLayers(obj, overwrite = TRUE)

  expr_mat <- GetAssayData(obj, assay = assay, slot = slot)
  matching_features <- grep(paste0("(^|-)", gene, "($|\\b)"), rownames(expr_mat), value = TRUE)
  if (length(matching_features) == 0) {
    warning("No matching features found for gene: ", gene)
    return(NULL)
  }

  total_expression <- Matrix::rowSums(expr_mat[matching_features, , drop = FALSE])
  top_features <- names(sort(total_expression, decreasing = TRUE))
  if (!is.null(n_features)) top_features <- head(top_features, n_features)

  message("Top features:\n", paste(
    sprintf("%s  (sum=%.2f)", top_features, total_expression[top_features]),
    collapse = "\n"
  ))

  feature_max <- vapply(top_features, function(f) {
    x <- expr_mat[f, ]
    if (all(x == 0)) return(0)
    as.numeric(quantile(x, probs = clip_quant, na.rm = TRUE))
  }, numeric(1))
  feature_max[feature_max == 0] <- 1

  plots <- list()

  # 1. Add gene-level plot if requested
  if (include_gene_plot && gene %in% rownames(GetAssayData(obj, assay = gene_assay))) {
    gene_expr <- GetAssayData(obj, assay = gene_assay, slot = gene_slot)[gene, ]
    gene_clip <- quantile(gene_expr[gene_expr > 0], clip_quant, na.rm = TRUE)
    if (is.na(gene_clip) || gene_clip == 0) gene_clip <- 1

    p_gene <- FeaturePlot(
      obj,
      features = gene,
      reduction = reduction,
      pt.size = pt.size
    )

    if (color_mode == "viridis") {
      p_gene <- p_gene + scale_color_viridis_c(
        option   = viridis_option,
        limits   = c(0, gene_clip),
        oob      = scales::squish,
        name     = "Expression",
        na.value = "grey90"
      )
    } else if (color_mode == "custom") {
      p_gene <- p_gene + scale_color_gradientn(
        colours = c("grey95", "#cbc9e2", "#9e9ac8", "#756bb1", "#54278f"),
        limits  = c(0, gene_clip),
        oob     = scales::squish,
        name    = "Expression",
        na.value = "grey90"
      )
    }

    p_gene <- p_gene +
      ggtitle(paste0(gene, " (gene-level)")) +
      theme_void(base_size = 14) +
      coord_fixed() +
      theme(
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = if (show_legend) "right" else "none",
        legend.title    = element_text(face = "bold", size = 9),
        legend.text     = element_text(size = 8),
        plot.margin     = margin(2, 2, 2, 2, "pt")
      )

    plots <- c(plots, list(p_gene))
  }

  # 2. Add isoform plots
  isoform_plots <- lapply(top_features, function(f) {
    clip_val <- feature_max[f]
    p <- FeaturePlot(
      obj,
      features = f,
      reduction = reduction,
      pt.size = pt.size
    )

    if (color_mode == "viridis") {
      p <- p + scale_color_viridis_c(
        option   = viridis_option,
        limits   = c(0, clip_val),
        oob      = scales::squish,
        name     = "Expression",
        na.value = "grey90"
      )
    } else if (color_mode == "custom") {
      p <- p + scale_color_gradientn(
        colours = c("grey95", "#cbc9e2", "#9e9ac8", "#756bb1", "#54278f"),
        limits  = c(0, clip_val),
        oob     = scales::squish,
        name    = "Expression",
        na.value = "grey90"
      )
    }

    p +
      theme_void(base_size = 14) +
      coord_fixed() +
      theme(
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = if (show_legend) "right" else "none",
        legend.title    = element_text(face = "bold", size = 9),
        legend.text     = element_text(size = 8),
        plot.margin     = margin(2, 2, 2, 2, "pt")
      ) +
      ggtitle(f)
  })

  plots <- c(plots, isoform_plots)

  combined <- patchwork::wrap_plots(plots, ncol = ncol) &
    theme(plot.margin = margin(5, 5, 5, 5, "pt"))

  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      )
    )
  }

  return(combined)
}



