#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(terra)
  library(geodata)
  library(Biostrings)
  library(glmmTMB)
  library(ggVennDiagram)
  library(Cairo)
})

# ============================================================
# TE residuals, climate associations, TFBS enrichment,
# and spatial gradients in Malus sylvestris
#
# 
#
# Main analyses:
#   1. Load TE counts, residual TE list, metadata, and TE FASTA
#   2. Restrict to M. sylvestris individuals with coordinates
#   3. Extract WorldClim bioclim variables
#   4. Run TE ~ climate GLMs within groups
#   5. Build Venn diagrams and common TE sets
#   6. Run spatial negative-binomial models on common TEs
#   7. Optionally summarize TFBS/FIMO results
#
# Example:
#   Rscript scripts/11_run_te_climate_spatial_tfbs.R \
#     --sample-info data/metadata/sample_info.csv \
#     --te-counts data/processed/te_counts.csv \
#     --residuals data/processed/te_residuals.tsv \
#     --te-fasta data/reference/te_database.fasta \
#     --outdir results/s11_te_climate_spatial \
#     --tf-family data/reference/Mdo_TF_list.txt \
#     --fimo-copia results/fimo/fimo_copia.tsv \
#     --fimo-gypsy results/fimo/fimo_gypsy.tsv \
#     --fimo-ltr results/fimo/fimo_LTR.tsv \
#     --fimo-nocat results/fimo/fimo_nocat.tsv \
#     --fimo-tir results/fimo/fimo_TIR.tsv
# ============================================================

# ----------------------------
# CLI arguments
# ----------------------------
option_list <- list(
  make_option("--sample-info", type = "character", help = "Sample metadata CSV"),
  make_option("--te-counts", type = "character", help = "TE count matrix CSV"),
  make_option("--residuals", type = "character", help = "Residual TE list TSV"),
  make_option("--te-fasta", type = "character", help = "TE FASTA file"),
  make_option("--outdir", type = "character", default = "results/s11_te_climate_spatial",
              help = "Output directory [default: %default]"),

  make_option("--tf-family", type = "character", default = NULL,
              help = "Optional TF family table TSV"),
  make_option("--fimo-copia", type = "character", default = NULL,
              help = "Optional FIMO TSV for Copia"),
  make_option("--fimo-gypsy", type = "character", default = NULL,
              help = "Optional FIMO TSV for Gypsy"),
  make_option("--fimo-ltr", type = "character", default = NULL,
              help = "Optional FIMO TSV for LTR"),
  make_option("--fimo-nocat", type = "character", default = NULL,
              help = "Optional FIMO TSV for NoCat"),
  make_option("--fimo-tir", type = "character", default = NULL,
              help = "Optional FIMO TSV for TIR"),

  make_option("--worldclim-dir", type = "character", default = tempdir(),
              help = "Directory used to cache WorldClim downloads [default: tempdir()]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------
# Logging helpers
# ----------------------------
log_info <- function(...) message("[INFO] ", paste0(..., collapse = ""))
log_ok   <- function(...) message("[OK] ", paste0(..., collapse = ""))
log_warn <- function(...) warning(paste0("[WARN] ", paste0(..., collapse = "")), call. = FALSE)

stop_if_missing <- function(path, label) {
  if (is.null(path) || !nzchar(path)) {
    stop("Missing required argument: ", label, call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("File not found for ", label, ": ", path, call. = FALSE)
  }
}

safe_write_lines <- function(x, path) {
  if (length(x) == 0) {
    file.create(path)
  } else {
    writeLines(x, path)
  }
}

# ----------------------------
# Validate required inputs
# ----------------------------
stop_if_missing(opt$`sample-info`, "--sample-info")
stop_if_missing(opt$`te-counts`,   "--te-counts")
stop_if_missing(opt$residuals,     "--residuals")
stop_if_missing(opt$`te-fasta`,    "--te-fasta")

# ----------------------------
# Output directories
# ----------------------------
dir_main         <- opt$outdir
dir_te_results   <- file.path(dir_main, "te_results")
dir_ggvenn       <- file.path(dir_te_results, "ggvenn")
dir_maps         <- file.path(dir_main, "maps_true_vs_sylvestris")
dir_vs_lists     <- file.path(dir_main, "vs_lists_by_threshold")
dir_tfbs         <- file.path(dir_main, "tfbs")

dir.create(dir_main, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_te_results, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_ggvenn, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_maps, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_vs_lists, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tfbs, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helpers
# ----------------------------
clean_sample_ids <- function(x) {
  x <- gsub("^A_|^B_", "", x)
  x <- gsub("_?list_sorting\\.txt$", "", x)
  x
}

extract_te_family <- function(te_id) {
  dplyr::case_when(
    str_detect(te_id, "_LTR(_|$)") ~ "LTR",
    str_detect(te_id, "_LINE(_|$)") ~ "LINE",
    str_detect(te_id, "_TIR(_|$)") ~ "TIR",
    str_detect(te_id, "_SINE(_|$)") ~ "SINE",
    str_detect(te_id, "_SSR(_|$)") ~ "SSR",
    str_detect(te_id, "_Helitron(_|$)") ~ "Helitron",
    str_detect(te_id, "_mixture(_|$)") ~ "mixture",
    str_detect(te_id, "_pararetrovirus(_|$)") ~ "pararetrovirus",
    str_detect(te_id, "_noCat(_|$)") ~ "noCat",
    str_detect(te_id, "_MITE(_|$)") ~ "MITE",
    TRUE ~ "Unknown"
  )
}

extract_te_superfamily <- function(te_id) {
  parts <- stringr::str_split(te_id, "_", simplify = TRUE)
  out <- rep(NA_character_, nrow(parts))
  if (ncol(parts) >= 4) out <- parts[, 4]
  out
}

fit_quasipoisson_glm <- function(df) {
  if (nrow(df) < 3) return(list(valid = FALSE, estimate = NA_real_, p_value = NA_real_))
  if (length(unique(df$te_count)) <= 1) return(list(valid = FALSE, estimate = NA_real_, p_value = NA_real_))
  if (length(unique(df$climate)) <= 1) return(list(valid = FALSE, estimate = NA_real_, p_value = NA_real_))

  model <- try(glm(te_count ~ climate, data = df, family = quasipoisson()), silent = TRUE)
  if (inherits(model, "try-error")) {
    return(list(valid = FALSE, estimate = NA_real_, p_value = NA_real_))
  }

  coef_table <- summary(model)$coefficients
  if (!"climate" %in% rownames(coef_table)) {
    return(list(valid = FALSE, estimate = NA_real_, p_value = NA_real_))
  }

  list(
    valid = TRUE,
    estimate = coef_table["climate", 1],
    p_value = coef_table["climate", 4]
  )
}

direction_from_betas <- function(beta_lat, beta_long) {
  dplyr::case_when(
    beta_lat > 0 & beta_long > 0 ~ "NE",
    beta_lat > 0 & beta_long < 0 ~ "NW",
    beta_lat < 0 & beta_long > 0 ~ "SE",
    beta_lat < 0 & beta_long < 0 ~ "SW",
    TRUE ~ "Undetermined"
  )
}

# ----------------------------
# Load data
# ----------------------------
log_info("Loading input files...")

sample_info <- readr::read_csv(opt$`sample-info`, show_col_types = FALSE)
count_te_raw <- readr::read_csv(opt$`te-counts`, show_col_types = FALSE)
residuals_tbl <- readr::read_tsv(opt$residuals, col_names = FALSE, show_col_types = FALSE)
colnames(residuals_tbl)[1] <- "TE_ID"
fasta_te <- Biostrings::readDNAStringSet(opt$`te-fasta`)

log_ok("Loaded files successfully.")
log_info("Metadata rows: ", nrow(sample_info))
log_info("TE count matrix dimensions: ", nrow(count_te_raw), " x ", ncol(count_te_raw))
log_info("Residual TE count: ", nrow(residuals_tbl))
log_info("FASTA sequence count: ", length(fasta_te))

# ----------------------------
# Clean TE count matrix
# ----------------------------
count_te_clean <- count_te_raw %>%
  tibble::as_tibble()

colnames(count_te_clean)[1] <- "Sample_ID_raw"

count_te_clean <- count_te_clean %>%
  dplyr::mutate(Sample_ID_clean = clean_sample_ids(Sample_ID_raw)) %>%
  dplyr::select(-Sample_ID_raw) %>%
  dplyr::select(Sample_ID_clean, dplyr::everything())

sequence_names <- names(fasta_te)
if (length(sequence_names) == (ncol(count_te_clean) - 1)) {
  colnames(count_te_clean)[-1] <- sequence_names
  log_ok("TE column names replaced using FASTA headers.")
} else {
  log_warn("TE matrix and FASTA do not contain the same number of entries; original TE column names were kept.")
}

te_residual_names <- residuals_tbl$TE_ID
te_residuals_present <- intersect(te_residual_names, colnames(count_te_clean))

if (length(te_residuals_present) == 0) {
  stop("No residual TE IDs from the residual list were found in the TE count matrix.", call. = FALSE)
}

log_info("Residual TEs found in matrix: ", length(te_residuals_present))

# ----------------------------
# Keep only M. sylvestris with coordinates
# ----------------------------
required_meta_cols <- c("Sample_ID", "Species", "latitude", "longitude", "group_structure_CORE_admix_K6")
missing_meta_cols <- setdiff(required_meta_cols, names(sample_info))
if (length(missing_meta_cols) > 0) {
  stop("Missing required metadata column(s): ", paste(missing_meta_cols, collapse = ", "), call. = FALSE)
}

meta_with_coords <- sample_info %>%
  dplyr::filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
  dplyr::select(Sample_ID, latitude, longitude, group_structure_CORE_admix_K6)

if (nrow(meta_with_coords) == 0) {
  stop("No M. sylvestris samples with non-missing latitude and longitude were found.", call. = FALSE)
}

all_te_data_filtered <- count_te_clean %>%
  dplyr::rename(Sample_ID = Sample_ID_clean) %>%
  dplyr::filter(Sample_ID %in% meta_with_coords$Sample_ID)

if (nrow(all_te_data_filtered) == 0) {
  stop("No TE count rows matched the filtered metadata. Check Sample_ID formatting.", call. = FALSE)
}

log_info("Individuals retained after coordinate filtering: ", nrow(all_te_data_filtered))

all_te_final <- all_te_data_filtered %>%
  dplyr::select(Sample_ID, dplyr::all_of(te_residuals_present))

readr::write_csv(all_te_final, file.path(dir_main, "coord_syl_te.csv"))

# ----------------------------
# Climate data and grouping
# ----------------------------
metadata <- meta_with_coords %>%
  dplyr::mutate(
    group_label = dplyr::case_when(
      group_structure_CORE_admix_K6 %in% c("Q5_Admix", "Q1_Admix") ~ "SYL_Admix",
      group_structure_CORE_admix_K6 == "Q5" ~ "SYLE",
      group_structure_CORE_admix_K6 == "Q1" ~ "SYLW",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group_label))

te_meta <- all_te_final %>%
  dplyr::left_join(metadata, by = "Sample_ID") %>%
  dplyr::filter(!is.na(latitude) & !is.na(longitude))

log_info("Group counts used for climate models:")
print(table(te_meta$group_label, useNA = "ifany"))

log_info("Downloading/loading WorldClim bioclim variables...")
bio <- geodata::worldclim_global(var = "bio", res = 10, path = opt$`worldclim-dir`)

coords <- te_meta %>%
  dplyr::select(Sample_ID, longitude, latitude) %>%
  tidyr::drop_na()

points <- terra::vect(coords, geom = c("longitude", "latitude"), crs = terra::crs(bio[[1]]))

clim_data <- terra::extract(bio, points) %>%
  tibble::as_tibble()

clim_data$Sample_ID <- coords$Sample_ID
clim_data <- clim_data %>% dplyr::select(-ID)

full_data <- te_meta %>% dplyr::left_join(clim_data, by = "Sample_ID")

bio_vars <- names(clim_data)[grepl("^bio", names(clim_data))]
te_vars <- names(all_te_final)[-1]

if (length(bio_vars) == 0) {
  stop("No bioclimatic variables were extracted from WorldClim.", call. = FALSE)
}

log_info("Bioclimatic variables extracted: ", length(bio_vars))

# ----------------------------
# Climate GLMs
# ----------------------------
groups <- c("SYLW", "SYLE", "SYL_Admix")
sig_results <- list()
total_tests <- 0L
total_skipped <- 0L
total_models_run <- 0L
all_glm_results <- list()

for (grp in groups) {
  log_info("Running climate GLMs for group: ", grp)
  group_data <- dplyr::filter(full_data, group_label == grp)

  sig_results[[grp]] <- list(
    "p<0.05" = character(0),
    "p<0.01" = character(0),
    "p<0.001" = character(0)
  )

  for (bio_var in bio_vars) {
    for (te in te_vars) {
      total_tests <- total_tests + 1L

      df <- group_data %>%
        dplyr::select(
          te_count = dplyr::all_of(te),
          climate = dplyr::all_of(bio_var)
        ) %>%
        tidyr::drop_na()

      fit <- fit_quasipoisson_glm(df)

      if (fit$valid) {
        total_models_run <- total_models_run + 1L
      } else {
        total_skipped <- total_skipped + 1L
      }

      all_glm_results[[length(all_glm_results) + 1]] <- tibble::tibble(
        group = grp,
        bio_var = bio_var,
        TE_ID = te,
        n = nrow(df),
        estimate = fit$estimate,
        p_value = fit$p_value
      )

      if (!is.na(fit$p_value)) {
        if (fit$p_value < 0.05)  sig_results[[grp]][["p<0.05"]]  <- c(sig_results[[grp]][["p<0.05"]], te)
        if (fit$p_value < 0.01)  sig_results[[grp]][["p<0.01"]]  <- c(sig_results[[grp]][["p<0.01"]], te)
        if (fit$p_value < 0.001) sig_results[[grp]][["p<0.001"]] <- c(sig_results[[grp]][["p<0.001"]], te)
      }
    }
  }
}

glm_results <- dplyr::bind_rows(all_glm_results)
readr::write_csv(glm_results, file.path(dir_main, "glm_te_vs_climate_all_results_raw_p.csv"))

log_info("====== Climate GLM summary ======")
log_info("Total TE x climate combinations tested: ", total_tests)
log_info("Valid models run: ", total_models_run)
log_info("Skipped models (constant or insufficient variation): ", total_skipped)

for (thresh in c("p<0.05", "p<0.01", "p<0.001")) {
  union_te <- unique(unlist(lapply(sig_results, function(x) x[[thresh]])))
  safe_write_lines(
    sort(union_te),
    file.path(dir_te_results, paste0("significant_te_", gsub("<", "p", thresh), ".txt"))
  )
}

for (thresh in c("p<0.05", "p<0.01", "p<0.001")) {
  sets <- list(
    SYLW = unique(sig_results[["SYLW"]][[thresh]]),
    SYLE = unique(sig_results[["SYLE"]][[thresh]]),
    SYL_Admix = unique(sig_results[["SYL_Admix"]][[thresh]])
  )

  sets_nonempty <- sets[lengths(sets) > 0]

  if (length(sets_nonempty) >= 2) {
    p_venn <- ggVennDiagram(sets_nonempty, label_alpha = 0) +
      ggtitle(paste("Venn diagram of significant TEs for", thresh)) +
      theme(legend.position = "none")

    ggsave(
      filename = file.path(dir_ggvenn, paste0("ggvenn_", gsub("<", "p", thresh), ".pdf")),
      plot = p_venn,
      width = 6,
      height = 5
    )
  }
}

common_te <- Reduce(intersect, list(
  unique(sig_results[["SYLW"]][["p<0.01"]]),
  unique(sig_results[["SYLE"]][["p<0.01"]]),
  unique(sig_results[["SYL_Admix"]][["p<0.01"]])
))

safe_write_lines(common_te, file.path(dir_te_results, "common_te_p0.01.txt"))
log_info("Common TEs across SYLW, SYLE, and SYL_Admix at p < 0.01: ", length(common_te))

te_matrix <- tibble::tibble(TE_ID = common_te)

for (grp in groups) {
  group_data <- dplyr::filter(full_data, group_label == grp)

  for (bio_var in bio_vars) {
    col_name <- paste0(grp, "_", bio_var)
    te_matrix[[col_name]] <- NA

    for (te in common_te) {
      df <- group_data %>%
        dplyr::select(
          te_count = dplyr::all_of(te),
          climate = dplyr::all_of(bio_var)
        ) %>%
        tidyr::drop_na()

      fit <- fit_quasipoisson_glm(df)
      if (!is.na(fit$p_value) && fit$p_value < 0.01) {
        te_matrix[te_matrix$TE_ID == te, col_name] <- TRUE
      }
    }
  }
}

readr::write_csv(te_matrix, file.path(dir_te_results, "common_te_p0.01_wide_format.csv"))

fasta_filtered <- fasta_te[names(fasta_te) %in% common_te]
Biostrings::writeXStringSet(fasta_filtered, filepath = file.path(dir_te_results, "common_te_p0.01.fasta"))

# ----------------------------
# Spatial model
# ----------------------------
if (length(common_te) > 0) {
  meta_spatial <- sample_info %>%
    dplyr::filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
    dplyr::filter(group_structure_CORE_admix_K6 != "Q2_Admix") %>%
    dplyr::select(Sample_ID, latitude, longitude)

  all_te_spatial <- count_te_clean %>%
    dplyr::rename(Sample_ID = Sample_ID_clean) %>%
    dplyr::filter(Sample_ID %in% meta_spatial$Sample_ID) %>%
    dplyr::select(Sample_ID, dplyr::all_of(common_te)) %>%
    dplyr::left_join(meta_spatial, by = "Sample_ID")

  df_geo <- all_te_spatial %>%
    tidyr::pivot_longer(cols = dplyr::all_of(common_te), names_to = "TE_ID", values_to = "Count") %>%
    dplyr::mutate(
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      sample_id = as.factor(Sample_ID)
    )

  df_nested <- df_geo %>%
    dplyr::group_by(TE_ID) %>%
    tidyr::nest()

  log_info("Fitting spatial negative-binomial models for common TEs...")

  vs_results <- df_nested %>%
    dplyr::mutate(
      model = purrr::map(data, ~ try(
        glmmTMB::glmmTMB(
          Count ~ latitude + longitude,
          data = .x,
          family = glmmTMB::nbinom2()
        ),
        silent = TRUE
      )),
      valid = purrr::map_lgl(model, ~ !inherits(.x, "try-error")),
      beta_lat = purrr::map_dbl(model, ~ {
        if (inherits(.x, "try-error")) return(NA_real_)
        fx <- try(glmmTMB::fixef(.x)$cond, silent = TRUE)
        if (!inherits(fx, "try-error") && "latitude" %in% names(fx)) fx[["latitude"]] else NA_real_
      }),
      beta_long = purrr::map_dbl(model, ~ {
        if (inherits(.x, "try-error")) return(NA_real_)
        fx <- try(glmmTMB::fixef(.x)$cond, silent = TRUE)
        if (!inherits(fx, "try-error") && "longitude" %in% names(fx)) fx[["longitude"]] else NA_real_
      }),
      angle_rad = atan2(beta_lat, beta_long),
      projection = purrr::map2(data, angle_rad, ~ .x$longitude * cos(.y) + .x$latitude * sin(.y)),
      true_vs = purrr::map2_dbl(data, projection, ~ suppressWarnings(cor(.x$Count, .y, use = "complete.obs")))
    ) %>%
    dplyr::filter(valid, !is.na(beta_lat), !is.na(beta_long), !is.na(true_vs)) %>%
    dplyr::mutate(
      angle_rad_corrected = ifelse(true_vs < 0, angle_rad + pi, angle_rad),
      true_vs_corrected = abs(true_vs)
    ) %>%
    dplyr::arrange(dplyr::desc(abs(true_vs_corrected)))

  spatial_summary <- vs_results %>%
    dplyr::select(TE_ID, beta_lat, beta_long, angle_rad, angle_rad_corrected, true_vs, true_vs_corrected) %>%
    dplyr::mutate(direction = direction_from_betas(beta_lat, beta_long))

  readr::write_csv(spatial_summary, file.path(dir_main, "supplementary_table_glmm_vs_tes.csv"))

  vs_results_export <- vs_results %>%
    dplyr::transmute(
      TE_ID = TE_ID,
      n_samples = purrr::map_int(data, nrow),
      beta_lat = beta_lat,
      beta_long = beta_long,
      angle_rad = angle_rad,
      angle_rad_corrected = angle_rad_corrected,
      true_vs = true_vs,
      true_vs_corrected = true_vs_corrected
    )

  readr::write_csv(vs_results_export, file.path(dir_main, "vs_results_summary.csv"))

  thresholds <- seq(0.05, 1.0, by = 0.05)
  for (thresh in thresholds) {
    file_name <- file.path(dir_vs_lists, paste0("te_list_vs_gt_", sprintf("%.2f", thresh), ".txt"))
    te_ids <- vs_results %>%
      dplyr::filter(abs(true_vs_corrected) > thresh) %>%
      dplyr::pull(TE_ID) %>%
      unique() %>%
      sort()
    safe_write_lines(te_ids, file_name)
  }

  te_vs_fasta_ids <- vs_results$TE_ID
  fasta_subset <- fasta_te[names(fasta_te) %in% te_vs_fasta_ids]
  Biostrings::writeXStringSet(fasta_subset, filepath = file.path(dir_vs_lists, "te_sequences_used_in_vs.fasta"))

  p_vs <- spatial_summary %>%
    dplyr::mutate(
      family = stringr::str_replace(TE_ID, "^TE_[0-9]+_", ""),
      TE_ID = factor(TE_ID, levels = rev(TE_ID[order(true_vs_corrected)]))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = TE_ID, y = true_vs_corrected)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::labs(
      x = NULL,
      y = "Corrected vector strength",
      title = "Spatial vector strength per TE"
    )

  ggplot2::ggsave(file.path(dir_main, "vs_per_te.pdf"), p_vs, width = 9, height = 10)
  ggplot2::ggsave(file.path(dir_main, "vs_per_te.svg"), p_vs, width = 9, height = 10)

  world_map <- ggplot2::map_data("world")

  for (i in seq_len(nrow(vs_results))) {
    te_id <- vs_results$TE_ID[i]
    vs_value <- round(vs_results$true_vs_corrected[i], 3)
    model <- vs_results$model[[i]]
    df_te <- vs_results$data[[i]] %>%
      dplyr::mutate(predicted = stats::predict(model, type = "response"))

    coefs <- glmmTMB::fixef(model)$cond
    beta_lat <- coefs[["latitude"]]
    beta_long <- coefs[["longitude"]]
    vs_norm <- sqrt(beta_lat^2 + beta_long^2)

    arrow_data <- df_te %>%
      dplyr::summarise(
        mean_lat = mean(latitude),
        mean_long = mean(longitude),
        arrow_lat = mean_lat + (beta_lat / ifelse(vs_norm == 0, 1, vs_norm)) * 5,
        arrow_long = mean_long + (beta_long / ifelse(vs_norm == 0, 1, vs_norm)) * 5
      ) %>%
      dplyr::mutate(label = paste0(te_id, "\nVS = ", vs_value))

    p_map <- ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = world_map,
        ggplot2::aes(x = long, y = lat, group = group),
        fill = "gray95", color = "gray70"
      ) +
      ggplot2::geom_point(
        data = df_te,
        ggplot2::aes(x = longitude, y = latitude, color = predicted),
        size = 2.5, alpha = 0.85
      ) +
      ggplot2::geom_segment(
        data = arrow_data,
        ggplot2::aes(x = mean_long, y = mean_lat, xend = arrow_long, yend = arrow_lat),
        arrow = grid::arrow(length = grid::unit(0.25, "cm")),
        color = "black", linewidth = 1.2
      ) +
      ggplot2::geom_text(
        data = arrow_data,
        ggplot2::aes(x = mean_long, y = mean_lat, label = label),
        vjust = -1, hjust = 0.5, size = 4, fontface = "italic"
      ) +
      ggplot2::scale_color_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = mean(df_te$predicted, na.rm = TRUE),
        name = "Predicted"
      ) +
      ggplot2::coord_cartesian(xlim = c(-10, 40), ylim = c(35, 70)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0(vs_value, " | TE: ", te_id),
        subtitle = "Arrow = geographic gradient; color = spatial model prediction",
        x = "Longitude",
        y = "Latitude"
      )

    ggplot2::ggsave(
      filename = file.path(dir_maps, paste0(sprintf("%04d_", i), vs_value, "_", te_id, ".svg")),
      plot = p_map,
      width = 10,
      height = 7,
      device = Cairo::CairoSVG
    )
  }
} else {
  log_info("Skipping spatial model block: no common TE at p < 0.01 across the three groups.")
}

# ----------------------------
# Optional TFBS block
# ----------------------------
fimo_files <- c(
  Copia = opt$`fimo-copia`,
  Gypsy = opt$`fimo-gypsy`,
  LTR   = opt$`fimo-ltr`,
  NoCat = opt$`fimo-nocat`,
  TIR   = opt$`fimo-tir`
)

tfbs_inputs_available <- !is.null(opt$`tf-family`) &&
  file.exists(opt$`tf-family`) &&
  all(!is.na(fimo_files)) &&
  all(file.exists(fimo_files))

if (tfbs_inputs_available) {
  log_info("Running optional TFBS block...")

  tf_fam <- readr::read_tsv(opt$`tf-family`, show_col_types = FALSE)

  fimo_all <- purrr::map_dfr(names(fimo_files), function(type) {
    readr::read_tsv(fimo_files[[type]], show_col_types = FALSE) %>%
      dplyr::mutate(TE_type = type)
  })

  fimo_sig <- fimo_all %>%
    dplyr::filter(`q-value` < 0.05) %>%
    dplyr::mutate(neglog10_q = -log10(`q-value`)) %>%
    dplyr::left_join(tf_fam, by = c("motif_id" = "TF_ID")) %>%
    dplyr::filter(!is.na(Family))

  heat_data <- fimo_sig %>%
    dplyr::group_by(Family, sequence_name, TE_type) %>%
    dplyr::summarise(median_neglog10_q = median(neglog10_q), .groups = "drop")

  readr::write_csv(heat_data, file.path(dir_tfbs, "tfbs_vs_te_summary.csv"))

  p_tfbs <- ggplot2::ggplot(
    heat_data,
    ggplot2::aes(
      x = sequence_name,
      y = forcats::fct_reorder(Family, median_neglog10_q, max, na.rm = TRUE)
    )
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill = median_neglog10_q), color = "black", linewidth = 0.1) +
    ggplot2::scale_fill_viridis_c(name = "Median -log10(q)") +
    ggplot2::facet_grid(. ~ TE_type, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = "TE", y = "TF family") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.text.y = ggplot2::element_text(size = 7)
    )

  ggplot2::ggsave(file.path(dir_tfbs, "tfbs_vs_te_heatmap.pdf"), p_tfbs, width = 14, height = 7)
  ggplot2::ggsave(file.path(dir_tfbs, "tfbs_vs_te_heatmap.svg"), p_tfbs, width = 14, height = 7)
} else {
  log_info("Skipping TFBS block: optional TFBS/FIMO inputs not fully provided.")
}

# ----------------------------
# Session info
# ----------------------------
writeLines(capture.output(sessionInfo()), file.path(dir_main, "session_info.txt"))
log_ok("Pipeline completed successfully.")