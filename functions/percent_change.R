# percent_change.r

# ----- Outline -----
# combine_correlations
# generate_rchange_lists
# run_all_gsis_correlations

# ----- combine_correlations -----
combine_correlations <- function(initial_corr, new_corr) {
  # Left join the two correlation tables based on 'row' and 'column'
  combined_corr <- dplyr::left_join(initial_corr, new_corr, by = c("row", "column"), suffix = c("_initial", "_new"))
  
  # Calculate the absolute percent change in r and add it as a new column
  combined_corr <- combined_corr %>%
    mutate(r_difference = round(abs(r_new - r_initial), 2),
           percent_change_r = round(abs((r_new - r_initial) / r_initial) * 100, 2),
           percent_change_r_dir = round((r_new - r_initial) / abs(r_initial) * 100, 2))
  
  return(combined_corr)
}

# ---- Generate and compare POPs vs POPs correlations: generate_rchange_lists ----
generate_rchange_lists <- function(
    nested_corr_list,
    tissue = c("panc", "adip"),
    reference = "orig",
    comparisons = NULL,
    save_csv = FALSE,
    out_dir = paste0(Sys.Date(), "_pops_vs_pops_rchange")
) {
  tissue <- match.arg(tissue)
  
  # Filter by tissue prefix
  tissue_keys <- names(nested_corr_list)[grepl(paste0("^", tissue, "_"), names(nested_corr_list))]
  tissue_list <- nested_corr_list[tissue_keys]
  
  # Extract reference list (e.g., panc_orig or adip_orig)
  ref_key <- paste0(tissue, "_", reference)
  if (!ref_key %in% names(tissue_list)) {
    stop("Reference key not found: ", ref_key)
  }
  ref_list <- tissue_list[[ref_key]]
  
  # Determine comparison keys
  if (is.null(comparisons)) {
    comparisons <- setdiff(names(tissue_list), ref_key)
  }
  
  # Run comparisons
  rchange_results <- purrr::map(comparisons, function(comp_key) {
    comp_list <- tissue_list[[comp_key]]
    purrr::map2(ref_list, comp_list, combine_correlations)
  })
  names(rchange_results) <- comparisons
  
  # Optionally write to CSV
  if (save_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    purrr::iwalk(rchange_results, function(comp_list, comp_name) {
      purrr::iwalk(comp_list, function(tbl, pol_name) {
        out_file <- file.path(out_dir, paste0(comp_name, "_", pol_name, "_rchange.csv"))
        write.csv(tbl, out_file, row.names = FALSE)
      })
    })
  }
  
  return(rchange_results)
}

# ---- Generate and compare GSIS correlations: run_all_gsis_correlations ----
run_all_gsis_correlations <- function(
    tissues = c("panc", "adip"),
    gsis_types = c("core", "jk_glu", "jk_leu", "jk_fa", "mh_glu"),
    variants = c("orig", "females", "males", "imputed", "rm_ind", "rm_389", "rm_419", "rm_421", "rm_448"),
    compare_to = c("imputed", "rm_ind", "rm_389", "rm_419", "rm_421", "rm_448"),
    suffix = "_list",
    cor_fun_map = list(core = cor.adi, jk_glu = cor.jk, jk_leu = cor.jk, jk_fa = cor.jk, mh_glu = cor.mh),
    comb_fun = combine_correlations,
    save_csv = FALSE,
    out_dir = "combined_gsis_correlations",
    env = .GlobalEnv
) {
  # Step 1: Correlation for each variant
  for (tissue in tissues) {
    for (gsis in gsis_types) {
      cor_fun <- cor_fun_map[[gsis]]
      for (variant in variants) {
        
        # Skip invalid combos
        if ((tissue == "adip" && variant %in% c("rm_389", "rm_448")) ||
            (tissue == "panc" && variant == "rm_421")) next
        
        list_name <- glue::glue("axys_pops_with_sums_{tissue}_{gsis}_gsis_{variant}{suffix}")
        tbl_name <- glue::glue("corr_axys_{tissue}_{gsis}_gsis_{variant}_tbl")
        
        if (exists(list_name, envir = env)) {
          datalist <- get(list_name, envir = env)
          cor_result <- purrr::map(datalist, cor_fun)
          assign(tbl_name, cor_result, envir = env)
        }
      }
    }
  }
  
  # Step 2: Combine correlation results
  for (tissue in tissues) {
    for (gsis in gsis_types) {
      for (variant in compare_to) {
        if ((tissue == "adip" && variant %in% c("rm_389", "rm_448")) ||
            (tissue == "panc" && variant == "rm_421")) next
        
        orig_name <- glue::glue("corr_axys_{tissue}_{gsis}_gsis_orig_tbl")
        alt_name <- glue::glue("corr_axys_{tissue}_{gsis}_gsis_{variant}_tbl")
        out_name <- glue::glue("corr_axys_{tissue}_{gsis}_gsis_orig_vs_{variant}")
        
        if (exists(orig_name, envir = env) && exists(alt_name, envir = env)) {
          orig_tbl <- get(orig_name, envir = env)
          alt_tbl <- get(alt_name, envir = env)
          combined_tbl <- purrr::map2(orig_tbl, alt_tbl, comb_fun)
          assign(out_name, combined_tbl, envir = env)
          
          if (save_csv) {
            path <- file.path(out_dir, tissue, gsis, paste0("orig_vs_", variant))
            if (!dir.exists(path)) dir.create(path, recursive = TRUE)
            
            purrr::iwalk(combined_tbl, function(df, name) {
              file <- file.path(path, paste0(name, "_combined_corr.csv"))
              readr::write_csv(df, file)
            })
          }
        }
      }
    }
  }
  
  message("✓ Correlation + combination finished. CSVs saved: ", save_csv)
}