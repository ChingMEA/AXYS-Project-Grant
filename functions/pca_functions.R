# ---- plot_pca_scores ----
plot_pca_scores <- function(data, group_var = NULL, title = "", pc1_var, pc2_var, nudge_y = 0.5) {
  
  p <- ggplot(data, aes(x = PC1, y = PC2)) +
    geom_text_repel(
      aes(label = DONOR),
      size = 3, color = "black", max.overlaps = 10, nudge_y = nudge_y,
      segment.color = "grey70", segment.alpha = 0.5
    ) +
    labs(
      title = title,
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14),
      legend.position = "right"
    )
  
  if (!is.null(group_var)) {
    group_var <- rlang::ensym(group_var)
    
    p <- p +
      geom_point(aes(colour = !!group_var), size = 3) +
      stat_ellipse(aes(group = !!group_var, colour = !!group_var), level = 0.95, linetype = "dashed")
  } else {
    p <- p +
      geom_point(size = 3) +
      stat_ellipse(aes(group = 1), level = 0.95, linetype = "dashed")
  }
  
  return(p)
}

# ---- make_pca_df ----
make_pca_df <- function(scores, meta, donor_col = "Donor", sex_col = "Sex", bmi_col = "BMI_category", cit_col = "coldischemiatime", dt_col = "digesttime", purity_col = "puritypercentage", ipi_col = "isletparticleindex") {
  data.frame(
    PC1 = scores[, 1],
    PC2 = scores[, 2],
    DONOR = meta[[donor_col]],
    Sex = meta[[sex_col]],
    BMI = meta[[bmi_col]],
    ColdIschemiaTime = meta[[cit_col]],
    DigestTime = meta[[dt_col]],
    Purity = meta[[purity_col]],
    IsletParticleIndex = meta[[ipi_col]]
  )
}

# ---- as_prcomp_pm ----
as_prcomp_pm <- function(pca_obj, original_data) {
  if (!inherits(pca_obj, "pcaRes")) stop("Not a 'pcaRes' object")
  
  prcomp_like <- list(
    sdev     = sqrt(colSums(pca_obj@scores^2) / (nrow(pca_obj@scores) - 1)),
    rotation = pca_obj@loadings,
    x        = pca_obj@scores,
    center   = colMeans(original_data),
    scale    = apply(original_data, 2, sd)
  )
  
  class(prcomp_like) <- "prcomp"
  prcomp_like
}

# ---- plot_biplot_pm  ----
plot_biplot_pm <- function(pca_obj, metadata,
                           donor_col = "DONOR",
                           group_col = NULL,
                           pc_x = 1,
                           pc_y = 2,
                           explained_var = NULL,
                           ellipse = TRUE,
                           alpha = 0.6,
                           point_size = 3,
                           dot_label_size = 3,
                           arrow_label_size = 3,
                           point_colors = NULL,
                           dot_label_color = "black",
                           arrow_label_color = "black",
                           arrow_color = "gray30",
                           arrow_label_nudge_x = 0.2,
                           arrow_label_nudge_y = 0.2,
                           title = "PCA Biplot") {
  # Extract scores and loadings
  scores <- as.data.frame(pca_obj@scores)
  loadings <- as.data.frame(pca_obj@loadings)
  
  # Rename columns
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  colnames(loadings) <- paste0("PC", seq_len(ncol(loadings)))
  
  # Add metadata and donor labels if donor_col not NULL
  if (!is.null(donor_col)) {
    scores$DONOR <- metadata[[donor_col]]
  }
  if (!is.null(group_col)) {
    scores$Group <- metadata[[group_col]]
  }
  
  # Select PCs
  pcx <- paste0("PC", pc_x)
  pcy <- paste0("PC", pc_y)
  
  # Axis labels with % variance if provided
  xlab <- if (!is.null(explained_var)) {
    paste0(pcx, " (", round(100 * explained_var[pc_x], 1), "%)")
  } else pcx
  
  ylab <- if (!is.null(explained_var)) {
    paste0(pcy, " (", round(100 * explained_var[pc_y], 1), "%)")
  } else pcy
  
  # Automatic scaling of loadings
  x_range <- range(scores[[pcx]], na.rm = TRUE)
  y_range <- range(scores[[pcy]], na.rm = TRUE)
  loadings_range <- range(c(loadings[[pcx]], loadings[[pcy]]), na.rm = TRUE)
  scaling_factor <- 0.8 * min(diff(x_range), diff(y_range)) / max(abs(loadings_range))
  
  # Scale loadings
  loadings_scaled <- loadings[, c(pcx, pcy)] * scaling_factor
  colnames(loadings_scaled) <- c("PCx", "PCy")
  loadings_scaled$Variable <- rownames(loadings)
  
  # Position arrow labels with nudges
  loadings_scaled$LabelX <- loadings_scaled$PCx + arrow_label_nudge_x
  loadings_scaled$LabelY <- loadings_scaled$PCy + arrow_label_nudge_y
  
  # Begin plot
  p <- ggplot()
  
  # Sample points with or without group colors
  if (!is.null(group_col)) {
    if (!is.null(point_colors)) {
      p <- p + geom_point(data = scores, aes_string(x = pcx, y = pcy, color = "Group"), size = point_size, alpha = alpha) +
        scale_color_manual(values = point_colors)
    } else {
      p <- p + geom_point(data = scores, aes_string(x = pcx, y = pcy, color = "Group"), size = point_size, alpha = alpha)
    }
  } else {
    if (!is.null(point_colors)) {
      p <- p + geom_point(data = scores, aes_string(x = pcx, y = pcy), size = point_size, alpha = alpha, color = point_colors)
    } else {
      p <- p + geom_point(data = scores, aes_string(x = pcx, y = pcy), size = point_size, alpha = alpha)
    }
  }
  
  # Donor labels (only if donor_col not NULL)
  if (!is.null(donor_col)) {
    p <- p + geom_text_repel(
      data = scores,
      aes_string(x = pcx, y = pcy, label = "DONOR"),
      size = dot_label_size,
      color = dot_label_color,
      segment.color = "grey70",
      segment.alpha = 0.5,
      max.overlaps = 10
    )
  }
  
  # Confidence ellipses
  if (ellipse && !is.null(group_col)) {
    p <- p + stat_ellipse(data = scores,
                          aes_string(x = pcx, y = pcy, group = "Group", color = "Group"),
                          level = 0.95, linetype = "dashed", alpha = 0.2)
  } else if (ellipse && is.null(group_col)) {
    # Ellipse around all points (no group)
    p <- p + stat_ellipse(data = scores,
                          aes_string(x = pcx, y = pcy, group = 1),
                          level = 0.95, linetype = "dashed", alpha = 0.2, color = "black")
  }
  
  # Loadings (arrows)
  p <- p + geom_segment(
    data = loadings_scaled,
    aes(x = 0, y = 0, xend = PCx, yend = PCy),
    arrow = arrow(length = unit(0.2, "cm")),
    color = arrow_color
  )
  
  # Variable labels (arrow labels)
  p <- p + geom_text_repel(
    data = loadings_scaled,
    aes(x = LabelX, y = LabelY, label = Variable),
    size = arrow_label_size,
    color = arrow_label_color,
    max.overlaps = 20
  )
  
  # Final layout
  p <- p +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.position = ifelse(is.null(group_col), "none", "right")
    )
  
  return(p)
}