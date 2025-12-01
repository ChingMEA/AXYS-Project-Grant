library(dplyr)
library(purrr)
library(rlang)

# ----- kendall_tau_by_compound -----

kendall_tau_by_compound <- function(
    df,
    compound_col = "COMPOUND",
    adipose_col   = "adipose",
    pancreas_col  = "pancreas",
    min_pairs = 3,
    ci = c("asymptotic","bootstrap"),
    conf = 0.95,
    nboot = 2000,
    seed = NULL
) {
  ci <- match.arg(ci)
  
  # helper: bootstrap CI for Kendall's tau on a 2-col dataframe
  boot_tau_ci <- function(dat, conf = 0.95, B = 2000, seed = NULL) {
    if (!requireNamespace("boot", quietly = TRUE))
      stop("Please install the 'boot' package for bootstrap CIs (install.packages('boot')).")
    if (!is.null(seed)) set.seed(seed)
    
    stat_fun <- function(d, idx) {
      x <- d[idx, 1, drop = TRUE]
      y <- d[idx, 2, drop = TRUE]
      suppressWarnings(cor(x, y, method = "kendall", use = "pairwise"))
    }
    
    b <- boot::boot(dat, statistic = stat_fun, R = B)
    q <- boot::boot.ci(b, type = "perc", conf = conf)
    if (is.null(q$percent)) return(c(NA_real_, NA_real_))
    c(q$percent[4], q$percent[5])  # percentile CI
  }
  
  # symbols
  comp_sym <- sym(compound_col); a_sym <- sym(adipose_col); p_sym <- sym(pancreas_col)
  
  out <- df %>%
    group_by(!!comp_sym) %>%
    group_modify(~{
      # use fully-qualified dplyr::select to avoid masking issues
      d <- .x %>%
        dplyr::select(adipose = !!a_sym, pancreas = !!p_sym) %>%
        mutate(
          adipose  = suppressWarnings(as.numeric(adipose)),
          pancreas = suppressWarnings(as.numeric(pancreas))
        ) %>%
        filter(complete.cases(adipose, pancreas))
      
      n <- nrow(d)
      if (n < min_pairs) {
        return(tibble(n = n, r = NA_real_, p = NA_real_,
                      se = NA_real_, r_lower = NA_real_, r_upper = NA_real_,
                      ci_method = ci))
      }
      
      ct <- suppressWarnings(cor.test(d$adipose, d$pancreas, method = "kendall"))
      r  <- unname(ct$estimate); p <- ct$p.value
      if (is.na(r)) {
        return(tibble(n = n, r = NA_real_, p = NA_real_,
                      se = NA_real_, r_lower = NA_real_, r_upper = NA_real_,
                      ci_method = ci))
      }
      
      if (ci == "asymptotic") {
        # Your SE + normal-approx CI (no-ties approximation)
        z <- qnorm(1 - (1 - conf)/2)
        se <- sqrt((2 * (2 * n + 5)) / (9 * n * (n - 1)))
        lo <- max(-1, r - z * se)
        hi <- min( 1, r + z * se)
      } else {
        se <- NA_real_
        ci_vec <- boot_tau_ci(as.data.frame(d[, c("adipose","pancreas")]), conf = conf, B = nboot, seed = seed)
        lo <- ci_vec[1]; hi <- ci_vec[2]
      }
      
      tibble(n = n, r = r, p = p, se = se,
             r_lower = lo, r_upper = hi,
             ci_method = ci)
    }) %>%
    ungroup() %>%
    rename(COMPOUND = !!comp_sym) %>%
    mutate(p_adj = p.adjust(p, method = "fdr")) %>%
    arrange(COMPOUND)
  
  out
}

# ----- partial_kendall_ppcor_multi -----
library(dplyr)
library(tidyr)
library(ppcor)
library(rlang)
library(purrr)

partial_kendall_ppcor_multi <- function(df,
                                        donor_col    = "DONOR",
                                        compound_col = "COMPOUND",
                                        adipose_col  = "adipose",
                                        pancreas_col = "pancreas",
                                        covar_sets = list(
                                          Unadjusted = character(0),
                                          Age        = "Age",
                                          BMI        = "BMI",
                                          Age_BMI    = c("Age","BMI")
                                        ),
                                        method = "kendall",          
                                        min_pairs = 3,
                                        conf = 0.95,
                                        ci = c("asymptotic","bootstrap"),
                                        nboot = 2000,
                                        seed = NULL) {
  ci <- match.arg(ci)
  donor_sym <- sym(donor_col); comp_sym <- sym(compound_col)
  a_sym <- sym(adipose_col);   p_sym <- sym(pancreas_col)
  zcrit <- qnorm(1 - (1 - conf)/2)
  
  # helper: estimate + p
  est_p <- function(d, covs) {
    if (length(covs) == 0) {
      ct  <- suppressWarnings(cor.test(d$adipose, d$pancreas, method = method))
      list(est = unname(ct$estimate), p = ct$p.value)
    } else {
      res <- suppressWarnings(ppcor::pcor.test(d$adipose, d$pancreas,
                                               d[, covs, drop = FALSE], method = method))
      list(est = unname(res$estimate), p = res$p.value)
    }
  }
  
  # bootstrap CI (resample donors)
  boot_ci <- function(d, covs) {
    if (!is.null(seed)) set.seed(seed)
    B <- nboot; n <- nrow(d); reps <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n, replace = TRUE)
      db <- d[idx, , drop = FALSE]
      reps[b] <- est_p(db, covs)$est
    }
    stats::quantile(reps, probs = c((1-conf)/2, 1-(1-conf)/2), na.rm = TRUE, names = FALSE)
  }
  
  results <- imap_dfr(covar_sets, function(covs, label) {
    df %>%
      group_by(!!comp_sym) %>%
      group_modify(~{
        # ensure unique donors within each compound slice
        d <- .x %>%
          distinct(!!donor_sym, .keep_all = TRUE) %>%
          dplyr::select(adipose = !!a_sym, pancreas = !!p_sym, all_of(covs)) %>%
          tidyr::drop_na()
        
        n <- nrow(d)
        if (n < min_pairs) {
          return(tibble(n = n, r = NA_real_, p = NA_real_,
                        se = NA_real_, r_lower = NA_real_, r_upper = NA_real_))
        }
        
        ep <- est_p(d, covs)
        est <- ep$est; pval <- ep$p
        
        if (ci == "bootstrap") {
          se <- NA_real_
          qi <- boot_ci(d, covs)
          lo <- max(-1, qi[1]); hi <- min(1, qi[2])
        } else {
          if (method == "pearson") {
            se <- NA_real_
            z  <- atanh(est); se_z <- 1 / sqrt(n - 3)
            lo <- tanh(z - zcrit * se_z); hi <- tanh(z + zcrit * se_z)
          } else {
            # your large-sample SE for Kendall/Spearman
            se <- sqrt((2 * (2 * n + 5)) / (9 * n * (n - 1)))
            lo <- max(-1, est - zcrit * se); hi <- min(1, est + zcrit * se)
          }
        }
        
        tibble(n, r = est, p = pval, se, r_lower = lo, r_upper = hi)
      }) %>%
      ungroup() %>%
      mutate(adjustment = label)
  }) %>%
    group_by(adjustment) %>%
    mutate(p_adj = p.adjust(p, "fdr")) %>%
    ungroup() %>%
    rename(COMPOUND = !!comp_sym)
  
  # side-by-side comparison
  results_wide <- results %>%
    dplyr::select(COMPOUND, adjustment, r, p, p_adj, r_lower, r_upper) %>%
    pivot_wider(
      names_from = adjustment,
      values_from = c(r, p, p_adj, r_lower, r_upper),
      names_glue = "{.value}_{adjustment}"
    )
  
  list(long = results, wide = results_wide)
}
