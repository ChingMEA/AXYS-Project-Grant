# correlations_gsis.R

# ----- Outline -----
# run_corr_analysis
# generate_all_corr_lists_and_tests

# ---- helper function: run_corr_analysis ----
run_corr_analysis <- function(df, method) {
  switch(method,
         "cor.mtest" = {
           res <- cor.mtest(df)
           list(r = res[[1]], p = res[[2]])  # Convert to named list
         },
         "cor.flat.ci" = cor.flat.ci(df),
         "get.upper.r" = get.upper.r(df),
         stop(paste("Unknown method:", method)))
}

# ---- generate_gsis_pollutant_lists ----
generate_gsis_pollutant_lists <- function(
    tissue = c("panc", "adip", "adip-to-panc"),
    pollutants = c("dioxin", "pcb", "ocp"),
    gsis_types = c("core", "jk_glu", "jk_leu", "jk_fa", "mh_glu"),
    variants = c("orig", "females", "males", "imputed", "rm_ind", "rm_389", "rm_419", "rm_421", "rm_448"),
    suffix = "_list",
    env = .GlobalEnv,
    log_env = .GlobalEnv
) {
  tissue <- match.arg(tissue)
  build_skipped <- character()
  
  for (v in variants) {
    for (g in gsis_types) {
      list_name <- glue::glue("axys_pops_with_sums_{tissue}_{g}_gsis_{v}{suffix}")
      list_obj <- list()
      
      for (pol in pollutants) {
        obj_name <- glue::glue("axys_pops_{tissue}_{pol}_{v}_{g}")
        if (exists(obj_name, envir = env)) {
          list_obj[[obj_name]] <- get(obj_name, envir = env)
        } else {
          build_skipped <- c(build_skipped, obj_name)
        }
      }
      
      if (length(list_obj) > 0) {
        assign(list_name, list_obj, envir = env)
      }
    }
  }
  
  assign(glue::glue("skipped_data_vars_{tissue}"), build_skipped, envir = log_env)
  message("✓ List generation completed for: ", tissue, " | Skipped data vars: ", length(build_skipped))
  invisible(NULL)
}

# ---- generate_gsis_pollutant_lists ----
run_all_corr_tests_on_lists <- function(
    tissue = c("panc", "adip", "adip-to-panc"),
    gsis_types = c("core", "jk_glu", "jk_leu", "jk_fa", "mh_glu"),
    variants = c("orig", "females", "males", "imputed", "rm_ind", "rm_389", "rm_419", "rm_421", "rm_448"),
    methods = c("cor.mtest", "get.upper.r", "cor.flat.ci"),
    suffix = "_list",
    save_csv = FALSE,
    out_dir = "correlation_results",
    env = .GlobalEnv,
    export_r_matrix = TRUE,
    export_p_matrix = TRUE,
    log_env = .GlobalEnv
) {
  tissue <- match.arg(tissue)
  cor_skipped <- character()
  corr_results <- list()
  
  for (v in variants) {
    for (g in gsis_types) {
      list_name <- glue::glue("axys_pops_with_sums_{tissue}_{g}_gsis_{v}{suffix}")
      if (!exists(list_name, envir = env)) {
        cor_skipped <- c(cor_skipped, list_name)
        next
      }
      
      list_obj <- get(list_name, envir = env)
      
      for (method in methods) {
        corr_name <- glue::glue("corr_{tissue}_{g}_gsis_{v}_{method}")
        corr_out <- purrr::map(list_obj, ~ run_corr_analysis(.x, method))
        assign(corr_name, corr_out, envir = env)
        corr_results[[corr_name]] <- corr_out
        
        if (save_csv) {
          out_path <- file.path(out_dir, paste0(tissue, "_", g, "_", v, "_", method))
          if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
          
          purrr::iwalk(corr_out, function(df, name) {
            if (is.data.frame(df)) {
              write.csv(df, file = file.path(out_path, paste0(name, ".csv")), row.names = FALSE)
            } else if (is.list(df) && all(c("r", "p") %in% names(df))) {
              if (export_r_matrix) {
                write.csv(df$r, file = file.path(out_path, paste0(name, "_r_matrix.csv")))
              }
              if (export_p_matrix) {
                write.csv(df$p, file = file.path(out_path, paste0(name, "_p_matrix.csv")))
              }
            }
          })
        }
      }
    }
  }
  
  assign(glue::glue("skipped_corr_lists_{tissue}"), cor_skipped, envir = log_env)
  message("✓ Correlation tests completed for: ", tissue, " | Skipped correlation lists: ", length(cor_skipped))
  invisible(corr_results)
}

# ---- main function: generate_all_corr_lists_and_tests ----
generate_all_corr_lists_and_tests <- function(
    tissue = c("panc", "adip", "adip-to-panc"),
    pollutants = c("dioxin", "pcb", "ocp"),
    gsis_types = c("core", "jk_glu", "jk_leu", "jk_fa", "mh_glu"),
    variants = c("orig", "females", "males", "imputed", "rm_ind", "rm_389", "rm_419", "rm_421", "rm_448"),
    methods = c("cor.mtest", "get.upper.r", "cor.flat.ci"),
    save_csv = FALSE,
    out_dir = "correlation_results",
    suffix = "_list",
    env = .GlobalEnv,
    log_env = .GlobalEnv,
    export_r_matrix = TRUE,
    export_p_matrix = TRUE
) {
  tissue <- match.arg(tissue)
  build_skipped <- character()
  cor_skipped <- character()
  corr_results <- list()
  
  for (v in variants) {
    for (g in gsis_types) {
      list_name <- glue::glue("axys_pops_with_sums_{tissue}_{g}_gsis_{v}{suffix}")
      list_obj <- list()
      
      for (pol in pollutants) {
        obj_name <- glue::glue("axys_pops_{tissue}_{pol}_{v}_{g}")
        if (exists(obj_name, envir = env)) {
          list_obj[[obj_name]] <- get(obj_name, envir = env)
        } else {
          build_skipped <- c(build_skipped, obj_name)
        }
      }
      
      if (length(list_obj) > 0) {
        assign(list_name, list_obj, envir = env)
        
        for (method in methods) {
          corr_name <- glue::glue("corr_{tissue}_{g}_gsis_{v}_{method}")
          corr_out <- purrr::map(list_obj, ~ run_corr_analysis(.x, method))
          assign(corr_name, corr_out, envir = env)
          corr_results[[corr_name]] <- corr_out
          
          if (save_csv) {
            out_path <- file.path(out_dir, paste0(tissue, "_", g, "_", v, "_", method))
            if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
            
            purrr::iwalk(corr_out, function(df, name) {
              if (is.data.frame(df)) {
                write.csv(df, file = file.path(out_path, paste0(name, ".csv")), row.names = FALSE)
              } else if (is.list(df) && all(c("r", "p") %in% names(df))) {
                if (export_r_matrix) {
                  write.csv(df$r, file = file.path(out_path, paste0(name, "_r_matrix.csv")))
                }
                if (export_p_matrix) {
                  write.csv(df$p, file = file.path(out_path, paste0(name, "_p_matrix.csv")))
                }
              }
            })
          }
        }
        
      } else {
        cor_skipped <- c(cor_skipped, list_name)
      }
    }
  }
  
  assign(glue::glue("skipped_data_vars_{tissue}"), build_skipped, envir = log_env)
  assign(glue::glue("skipped_corr_lists_{tissue}"), cor_skipped, envir = log_env)
  message("✓ Completed for: ", tissue,
          " | Skipped data vars: ", length(build_skipped),
          " | Skipped correlation lists: ", length(cor_skipped))
  
  invisible(corr_results)
}