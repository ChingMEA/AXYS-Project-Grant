# correlation_functions.r

# ----- Outline -----
# cor.mtest
# cor.flat
# cor.flat.ci
# get.upper.r
# apply_r
# cor.mh, cor.jk, cor.adi, cor.meta

# ----- cor.mtest -----
cor.mtest <- function(df, method = "kendall", ...) {
  n <- ncol(df)
  mat <- as.matrix(df[, 2:n])
  cor <- psych::corTest(mat, method = method, adjust = "fdr", use = "pairwise")
  cormat <- as.matrix(cor$r)
  pmat <- as.matrix(cor$p)
  return(list(cormat, pmat, cor))
}

# ----- cor.flat -----
cor.flat <- function(df, method = "kendall", ...) {
  n <- ncol(df)
  mat <- as.matrix(df[, 2:n])
  cor <- psych::corTest(mat, method = method, adjust = "fdr", use = "pairwise")
  
  # Convert correlation (r) matrix to long format
  r_df <- as.data.frame(cor$r)
  r_df$row <- rownames(r_df)
  melted_r <- reshape2::melt(r_df, id.vars = "row", variable.name = "column", value.name = "r")
  
  # Convert p-value matrix to long format
  p_df <- as.data.frame(cor$p)
  p_df$row <- rownames(p_df)
  melted_p <- reshape2::melt(p_df, id.vars = "row", variable.name = "column", value.name = "p")
  
  # Merge and adjust p-values
  final_data <- dplyr::left_join(melted_r, melted_p, by = c("row", "column")) %>%
    dplyr::distinct(row, column, .keep_all = TRUE) %>%
    dplyr::mutate(p_adj = p.adjust(p, method = "fdr"))
  
  return(final_data)
}

# ----- cor.flat.ci -----
cor.flat.ci <- function(df, method = "kendall", ...) {
  n <- ncol(df)
  mat <- as.matrix(df[, 2:n])

  cor <- tryCatch({
    psych::corTest(mat, method = method, adjust = "fdr", use = "pairwise")
  }, error = function(e) return(NULL))

  if (is.null(cor)) return(NULL)

  # Melt r and p matrices
  r_df <- as.data.frame(cor$r)
  r_df$row <- rownames(r_df)
  melted_r <- reshape2::melt(r_df, id.vars = "row", variable.name = "column", value.name = "r")

  p_df <- as.data.frame(cor$p)
  p_df$row <- rownames(p_df)
  melted_p <- reshape2::melt(p_df, id.vars = "row", variable.name = "column", value.name = "p")

  final_data <- dplyr::left_join(melted_r, melted_p, by = c("row", "column")) %>%
    dplyr::filter(row != column) %>%
    dplyr::mutate(p_adj = p.adjust(p, method = "fdr"))

  # Add n, standard error, and confidence intervals
  final_data <- final_data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n = sum(complete.cases(mat[, row], mat[, column])),
      se = sqrt((2 * (2 * n + 5)) / (9 * n * (n - 1))),
      r_lower = max(-1, r - 1.96 * se),
      r_upper = min(1, r + 1.96 * se)
    ) %>%
    dplyr::ungroup()

  return(final_data)
}

# ----- get.upper.r -----
get.upper.r <- function(df, method = "kendall", ...) {
  n <- ncol(df)
  mat <- as.matrix(df[, 2:n])
  cor_result <- psych::corTest(mat, method = method, use = "pairwise")
  
  # Extract Kendall's r matrix
  r_values <- cor_result$r
  
  # Keep only the upper triangle and diagonal elements
  upper_triangle <- upper.tri(r_values, diag = TRUE)
  r_values[!upper_triangle] <- NA
  
  # Convert matrix to data frame
  r_df <- as.data.frame(r_values)
  
  return(r_df)
}

# ----- apply_r -----
apply_r <- function(list_of_dfs, output = c("ci", "flat", "upper")) {
  output <- match.arg(output)
  method_fun <- switch(output,
                       ci    = cor.flat.ci,
                       flat  = cor.flat,
                       upper = get.upper.r)
  purrr::map(list_of_dfs, method_fun)
}

# ----- Correlation filtering functions: filter rows and columns by lab or data type -----
cor.mh <- function(df) {
  df %>%
    cor.flat.ci() %>%
    dplyr::filter(!str_detect(row, "MH"), str_detect(column, "MH"))
}

cor.jk <- function(df) {
  df %>%
    cor.flat.ci() %>%
    dplyr::filter(!str_detect(row, "JK"), str_detect(column, "JK"))
}

cor.adi <- function(df) {
  df %>%
    cor.flat.ci() %>%
    dplyr::filter(!str_detect(row, "IsletCore"), str_detect(column, "IsletCore"))
}

cor.meta <- function(df) {
  df %>%
    cor.flat.ci() %>%
    dplyr::filter((row == "Age" | row == "BMI" | row == "Hba1c") &
                    !(column == "Age" | column == "BMI" | column == "Hba1c"))
}