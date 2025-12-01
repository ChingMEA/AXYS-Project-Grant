# partial_correlation_functions.r

# ----- Outline -----
# partial_ppcor_pops
# cor.mtest.partial
# partial_ppcor_analyte_gsis

# ----- Packages -----
library(dplyr)
library(tidyr)
library(ppcor)
library(purrr)

# ---- partial_ppcor_pops ----
partial_ppcor_pops <- function(pops_list,
                               analyte_groups,
                               tissue_label = "Pancreas",
                               ...) {
  purrr::imap_dfr(pops_list, function(df_joined, nm) {
    # guess pollutant class from element name
    cls <- names(analyte_groups)[
      vapply(names(analyte_groups),
             function(k) grepl(k, nm, ignore.case = TRUE),
             logical(1))
    ]
    if (length(cls) == 0) return(tibble())
    cls <- cls[[1]]
    
    # select analyte columns for this class
    analytes <- intersect(names(df_joined), analyte_groups[[cls]])
    if (length(analytes) < 2) return(tibble())  # need at least 2 for correlations
    
    res <- partial_ppcor_analyte_gsis(
      df           = df_joined,
      analyte_cols = analytes,
      gsis_cols    = analytes,  # same set, so you get analyte × analyte
      ...
    )
    
    dplyr::mutate(res$long,
                  tissue = tissue_label,
                  pollutant_class = cls)
  })
}

# ----- cor.mtest.partial -----
cor.mtest.partial <- function(
    df,
    analyte_cols,
    gsis_cols,
    covariates = NULL,
    method = "kendall",
    adjust = "fdr",
    min_n = 5,
    conf  = 0.95,
    donor_col = NULL,
    drop_all_na = TRUE,
    symmetric = FALSE          
) {
  if (!is.null(donor_col) && donor_col %in% names(df)) {
    df <- dplyr::distinct(df, .data[[donor_col]], .keep_all = TRUE)
  }
  
  analyte_cols <- intersect(analyte_cols, names(df))
  gsis_cols    <- intersect(gsis_cols,    names(df))
  if (length(analyte_cols) == 0 || length(gsis_cols) == 0) {
    return(list(
      cormat = matrix(NA_real_, 0, 0),
      pmat   = matrix(NA_real_, 0, 0),
      raw_results = tibble(
        var1=character(), var2=character(), n=integer(),
        r=double(), p=double(), se=double(), r_lower=double(), r_upper=double(),
        p_adj=double()
      )
    ))
  }
  
  to_num <- function(x) suppressWarnings(as.numeric(x))
  df[analyte_cols] <- lapply(df[analyte_cols], to_num)
  df[gsis_cols]    <- lapply(df[gsis_cols],    to_num)
  
  if (drop_all_na) {
    analyte_cols <- analyte_cols[vapply(df[analyte_cols], function(v) !all(is.na(v)), logical(1))]
    gsis_cols    <- gsis_cols[   vapply(df[gsis_cols],    function(v) !all(is.na(v)), logical(1))]
    if (length(analyte_cols) == 0 || length(gsis_cols) == 0) {
      return(list(
        cormat = matrix(NA_real_, 0, 0),
        pmat   = matrix(NA_real_, 0, 0),
        raw_results = tibble(
          var1=character(), var2=character(), n=integer(),
          r=double(), p=double(), se=double(), r_lower=double(), r_upper=double(),
          p_adj=double()
        )
      ))
    }
  }
  
  if (!is.null(covariates) && length(covariates) > 0) {
    if (!all(covariates %in% names(df))) stop("Some covariates are not present in the dataframe.")
    df[covariates] <- lapply(df[covariates], to_num)
    has_var <- vapply(df[covariates], function(v) stats::sd(v, na.rm = TRUE) > 0, logical(1))
    covariates <- covariates[has_var]
  }
  
  # ----- pairs grid (toggle)
  if (symmetric) {
    # full square minus diagonal across analytes ∪ gsis
    all_vars <- c(analyte_cols, gsis_cols)
    pairs <- tidyr::expand_grid(var1 = all_vars, var2 = all_vars) %>%
      dplyr::filter(var1 != var2)
  } else {
    # cross-product: analytes (rows) × gsis (cols) — no duplicates
    pairs <- tidyr::expand_grid(var1 = analyte_cols, var2 = gsis_cols)
  }
  
  zcrit <- stats::qnorm(1 - (1 - conf)/2)
  
  compute_one <- function(v1, v2) {
    d <- dplyr::select(df,
                       x = dplyr::all_of(v1),
                       y = dplyr::all_of(v2),
                       dplyr::any_of(covariates)) %>%
      tidyr::drop_na()
    n <- nrow(d)
    if (n < min_n) {
      return(tibble(var1=v1, var2=v2, n=n, r=NA_real_, p=NA_real_,
                    se=NA_real_, r_lower=NA_real_, r_upper=NA_real_))
    }
    
    if (is.null(covariates) || length(covariates) == 0) {
      ct  <- suppressWarnings(stats::cor.test(d$x, d$y, method = method))
      est <- unname(ct$estimate); pval <- ct$p.value
    } else {
      if (length(covariates) == 0) {
        ct  <- suppressWarnings(stats::cor.test(d$x, d$y, method = method))
        est <- unname(ct$estimate); pval <- ct$p.value
      } else {
        res <- suppressWarnings(ppcor::pcor.test(d$x, d$y,
                                                 d[, covariates, drop = FALSE],
                                                 method = method))
        est <- unname(res$estimate); pval <- res$p.value
      }
    }
    
    if (method == "pearson") {
      se <- NA_real_
      z  <- atanh(est); se_z <- 1 / sqrt(n - 3)
      lo <- tanh(z - zcrit * se_z); hi <- tanh(z + zcrit * se_z)
    } else {
      se <- sqrt((2 * (2 * n + 5)) / (9 * n * (n - 1)))
      lo <- max(-1, est - zcrit * se); hi <- min(1, est + zcrit * se)
    }
    
    tibble(var1=v1, var2=v2, n=n, r=est, p=pval, se=se, r_lower=lo, r_upper=hi)
  }
  
  results <- purrr::map2_dfr(pairs$var1, pairs$var2, compute_one) %>%
    dplyr::mutate(p_adj = ifelse(is.na(.data$p), NA_real_, p.adjust(.data$p, method = adjust)))
  
  # Only build matrices in the cross-product case (symmetric isn't rectangular)
  if (!symmetric) {
    rmat <- matrix(NA_real_, nrow = length(analyte_cols), ncol = length(gsis_cols),
                   dimnames = list(analyte_cols, gsis_cols))
    pmat <- matrix(NA_real_, nrow = length(analyte_cols), ncol = length(gsis_cols),
                   dimnames = list(analyte_cols, gsis_cols))
    if (nrow(results) > 0) {
      for (i in seq_len(nrow(results))) {
        a <- results$var1[i]; g <- results$var2[i]
        if (a %in% rownames(rmat) && g %in% colnames(rmat)) {
          rmat[a, g] <- results$r[i]
          pmat[a, g] <- results$p_adj[i]
        }
      }
    }
  } else {
    rmat <- matrix(NA_real_, 0, 0)
    pmat <- matrix(NA_real_, 0, 0)
  }
  
  list(cormat = rmat, pmat = pmat, raw_results = results)
}

# ----- partial_ppcor_analyte_gsis -----
partial_ppcor_analyte_gsis <- function(
    df,
    analyte_cols,
    gsis_cols,
    covar_sets = list(
      Unadjusted = NULL,
      Age        = "Age",
      BMI        = "BMI",
      Age_BMI    = c("Age","BMI")
    ),
    method = "kendall",
    adjust = "fdr",
    min_n = 5,
    conf  = 0.95,
    fdr_scope = c("within","overall"),
    donor_col = NULL,
    symmetric = FALSE      
) {
  fdr_scope <- match.arg(fdr_scope)
  
  pieces <- purrr::imap(covar_sets, function(covs, lab) {
    out <- cor.mtest.partial(
      df           = df,
      analyte_cols = analyte_cols,
      gsis_cols    = gsis_cols,
      covariates   = covs,
      method       = method,
      adjust       = adjust,
      min_n        = min_n,
      conf         = conf,
      donor_col    = donor_col,
      symmetric    = FALSE
    )
    dplyr::mutate(out$raw_results, adjustment = lab)
  })
  
  long <- dplyr::bind_rows(pieces)
  if (nrow(long) == 0) return(list(long = long, wide = tibble::tibble()))
  
  # recompute p_adj at desired scope across all pairs
  if (fdr_scope == "within") {
    long <- long %>% dplyr::group_by(adjustment) %>% dplyr::mutate(p_adj = p.adjust(p, method = adjust)) %>% dplyr::ungroup()
  } else {
    long <- long %>% dplyr::mutate(p_adj = p.adjust(p, method = adjust))
  }
  
  # Wide table only meaningful for cross-product (non-symmetric)
  if (symmetric) {
    wide <- tibble::tibble()
  } else {
    wide <- long %>%
      dplyr::select(var1, var2, adjustment, r, p, p_adj, r_lower, r_upper, n) %>%
      tidyr::pivot_wider(
        names_from  = adjustment,
        values_from = c(r, p, p_adj, r_lower, r_upper, n),
        names_glue  = "{.value}_{adjustment}"
      )
  }
  
  list(long = long, wide = wide)
}

# ----- partial_ppcor_over_list -----
# pops_list: nested list (class -> gsis_type -> df)
# analyte_groups: list(dioxin = dioxin_list, pcb = pcb_list, ocp = ocp_list)
# gsis_groups:    list(core = gsis_adi_list, jk_glu = glu_peri_jk_list, ...)
# gsis_suffix: TRUE if your joined GSIS columns are suffixed with _{gsis_type}; FALSE if not
partial_ppcor_over_list <- function(pops_list, analyte_groups, gsis_groups,
                                    tissue_label = "Pancreas",
                                    gsis_suffix = FALSE,
                                    ...) {
  
  # sanity: require names on both levels
  if (is.null(names(pops_list))) stop("pops_list must be a *named* list by class.")
  if (any(vapply(pops_list, function(x) is.null(names(x)), logical(1)))) {
    stop("Each sublist in pops_list must be *named* by GSIS type.")
  }
  
  purrr::imap_dfr(pops_list, function(by_type_list, cls) {
    # look up analyte vector for this class
    analyte_vec <- analyte_groups[[cls]]
    if (is.null(analyte_vec)) {
      message("No analyte list for class '", cls, "'. Skipping.")
      return(tibble())
    }
    
    purrr::imap_dfr(by_type_list, function(df_joined, gsis_type) {
      gsis_base <- gsis_groups[[gsis_type]]
      if (is.null(gsis_base)) {
        message("No GSIS list for type '", gsis_type, "'. Skipping.")
        return(tibble())
      }
      
      # pick columns present in the dataframe
      analytes  <- intersect(names(df_joined), analyte_vec)
      gsis_cols <- if (gsis_suffix) {
        intersect(names(df_joined), paste0(gsis_base, "_", gsis_type))
      } else {
        intersect(names(df_joined), gsis_base)
      }
      
      if (length(analytes) == 0 || length(gsis_cols) == 0) {
        # quiet skip; return empty tibble for binding
        return(tibble())
      }
      
      res <- partial_ppcor_analyte_gsis(
        df           = df_joined,
        analyte_cols = analytes,
        gsis_cols    = gsis_cols,
        ...          # pass method = "kendall", fdr_scope = "within", donor_col = "DONOR", etc.
      )
      
      dplyr::mutate(res$long,
                    tissue = tissue_label,
                    pollutant_class = cls,
                    gsis_type = gsis_type)
    })
  })
}