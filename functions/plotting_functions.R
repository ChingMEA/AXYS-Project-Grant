# plotting_functions.r

# ----- Outline -----
# col, custom_palette: colour palettes for correlation matrix
# plot_correlation
# plot_all_correlations
# save_corrplot
# export_all_corrplots
# save_correlation_matrices

# ----- Colour palette for correlation matrix -----
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
custom_palette <- c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")

# ----- plot_correlation -----
plot_correlation <- function(corrplot_data, index) {
  corrplot(corrplot_data[[1]], type = "upper", order = "original", method = "circle", cl.align.text = "l", 
           diag = FALSE, addCoef.col = 'grey20', sig.level = 1, insig = 'blank', number.cex = 0.5, 
           p.mat = corrplot_data[[2]], tl.col = "black", tl.srt = 45, tl.cex = 0.75)
}

# ----- plot_all_correlations -----
plot_all_correlations <- function(df_list, output_dir = NULL, 
                                  width = 8, height = 6, res = 300) {
  # Create output directory if saving is enabled
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each dataframe in the list
  for (i in seq_along(df_list)) {
    df_name <- names(df_list)[i]
    corrplot_data <- df_list[[i]]
    
    if (!is.null(output_dir)) {
      png(file.path(output_dir, paste0(df_name, "_corrplot.png")), 
          width = width, height = height, units = "in", res = res)
    }
    
    corrplot::corrplot(corrplot_data[[1]], 
                       type = "upper", order = "original", method = "circle", cl.align.text = "l", 
                       diag = FALSE, addCoef.col = 'grey20', sig.level = 1, insig = 'blank', number.cex = 0.5, 
                       p.mat = corrplot_data[[2]], tl.col = "black", tl.srt = 45, tl.cex = 0.75, 
                       title = df_name, mar = c(0, 0, 2, 0))
    
    if (!is.null(output_dir)) dev.off()
  }
}

# ----- save_corrplot -----
save_corrplot <- function(filename, corrplot_data, width, height, show_legend = TRUE) {
  png(filename = filename, width = width, height = height, units = "in", res = 600)
  
  # Conditional cl.pos based on show_legend argument
  legend_position <- ifelse(show_legend, "l", "n")
  
  corrplot(corrplot_data[[1]], type = "upper", order = "original", method = "circle", cl.align.text = legend_position, 
           diag = FALSE, addCoef.col = 'gray20', sig.level = 1, insig = 'blank', col = col(10), number.cex = 0.5, 
           p.mat = corrplot_data[[2]], tl.col = "black", tl.srt = 45, tl.cex = 0.75, cl.pos = legend_position)
  
  dev.off()
}

# ----- export_all_corrplots -----
export_all_corrplots <- function(corr_results, 
                                 output_dir = "figures/correlation_plots",
                                 width = 6, 
                                 height = 6,
                                 show_legend = TRUE,
                                 res = 600) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Loop over all nested correlation results
  purrr::walk2(names(corr_results), corr_results, function(group_name, result_list) {
    purrr::imap(result_list, function(corr_data, sub_name) {
      filename <- file.path(output_dir, paste0(group_name, "_", sub_name, "_corrplot.png"))
      save_corrplot(
        filename = filename,
        corrplot_data = corr_data,
        width = width,
        height = height,
        show_legend = show_legend
      )
    })
  })
}

# ----- save_correlation_matrices -----
save_correlation_matrices <- function(
    correlation_list,
    date_today = Sys.Date(),
    parent_folder = ".",
    suffix = ".csv",
    method = "Kendall's tau",
    notes = NULL,
    write_metadata = TRUE
) {
  purrr::imap(correlation_list, function(sublist, level1_name) {
    # Create folder for this level1 group
    group_folder <- file.path(parent_folder, paste0(date_today, "_", level1_name))
    if (!dir.exists(group_folder)) dir.create(group_folder, recursive = TRUE)
    
    # Save CSVs using level2 name in filename
    filenames <- purrr::imap_chr(sublist, function(df, level2_name) {
      file_name <- paste0(date_today, "_", level2_name, "_", suffix)
      file_path <- file.path(group_folder, file_name)
      write.csv(df, file = file_path, row.names = TRUE)
      return(file_name)
    })
    
    # Optionally write metadata
    if (write_metadata) {
      note_text <- if (!is.null(notes) && !is.null(notes[[level1_name]])) notes[[level1_name]] else "None"
      
      metadata_path <- file.path(group_folder, "metadata.txt")
      metadata_contents <- paste0(
        "Export Date: ", date_today, "\n",
        "Group Folder: ", level1_name, "\n",
        "Correlation Method: ", method, "\n",
        "Exported Files:\n  - ", paste(filenames, collapse = "\n  - "), "\n",
        "Notes: ", note_text, "\n"
      )
      writeLines(metadata_contents, con = metadata_path)
    }
  })
}
