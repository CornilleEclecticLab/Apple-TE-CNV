# =============================================================================
# Script:  03_spatial_vs.R
# Purpose: Re-run the manuscript's Model B (the SPATIAL model) on the final
#          38-TE candidate set, so that Figure 3A and 3B can be regenerated.
#
#          NAMING: this script implements "Model B" in the manuscript's own
#          sense - TE copy number as a function of geographical coordinates,
#          fitted with a negative-binomial GLM (glmmTMB, nbinom2):
#
#              Count ~ latitude + longitude
#
#          from which the spatial vector strength (VS) of Fig. 3A-B is
#          derived. It is NOT the PC-covariate climate model; that one is
#          called "Model A-PC" everywhere in this folder.
#
#          Why it needs re-running: the candidate list changed from the
#          published 34 (Model A alone) to the final non-redundant union of
#          38 produced by 02_candidate_TEs.R. The four
#          TEs recovered only by the structure-controlled climate model have
#          no VS yet, so Figure 3A-B cannot include them.
#
#          The script first refits the published 34 and checks the results
#          against the archived VS table. If they do not reproduce, it stops
#          rather than emitting VS values computed by a different route.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   data/coord_syl_TE.csv              residual TE counts, M. sylvestris
#   data/250923_sample_list_TE.csv     coordinates and K6 clusters
#   data/published_VS_table.csv        archived VS values for the 34
#   outputs/02_candidate_TEs/tables/
#       FINAL_union_candidates.csv     the 38 final candidates
#
# Outputs (all under outputs/03_spatial_vs/):
#   tables/validation_VS_reproduction.csv  reproduction gate on the 34
#   tables/VS_table_38_candidates.csv      VS for the final 38 (Fig. 3A-B)
#   tables/VS_new_candidates_only.csv      the TEs added since the 34
#   figures/*.pdf + *.svg
#   ../../logs/03_spatial_vs_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(glmmTMB)   # negative-binomial GLM, as in the published pipeline
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
# The published spatial model excluded the single Q2_Admix individual.
EXCLUDE_CLUSTER <- "Q2_Admix"
VS_STRONG       <- 0.5   # threshold called "strong" in the Results

dir_out  <- file.path(dir_outputs, "03_spatial_vs")
dir_tab  <- file.path(dir_out, "tables")
dir_fig  <- file.path(dir_out, "figures")
for (d in c(dir_out, dir_tab, dir_fig)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p,
                  width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p,
                  width = width, height = height, device = svglite::svglite)
}

message("=== 03_spatial_vs.R ===")

# =============================================================================
# 1. Data, exactly as assembled by the published spatial model
# =============================================================================
counts <- readr::read_csv(file.path(dir_data, "coord_syl_TE.csv"),
                          show_col_types = FALSE) %>%
  dplyr::mutate(Sample_ID = toupper(Sample_ID))

meta_spatial <- readr::read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                                show_col_types = FALSE) %>%
  dplyr::filter(Species == "M. sylvestris",
                !is.na(latitude), !is.na(longitude),
                group_structure_CORE_admix_K6 != EXCLUDE_CLUSTER) %>%
  dplyr::transmute(Sample_ID = toupper(Sample_ID),
                   latitude  = as.numeric(latitude),
                   longitude = as.numeric(longitude))

message(sprintf("Individuals available for the spatial model: %d",
                nrow(dplyr::semi_join(meta_spatial, counts, by = "Sample_ID"))))

#' Fit the spatial model for one TE and derive its vector strength.
#'
#' Faithful to the published pipeline: negative-binomial GLM of count on
#' latitude and longitude; the gradient angle is atan2(beta_lat, beta_long);
#' VS is the correlation between the counts and the individual coordinates
#' projected onto that angle; a negative correlation flips the angle by pi and
#' the reported VS is its absolute value.
fit_vs <- function(te_id) {
  d <- counts %>%
    dplyr::select(Sample_ID, Count = dplyr::all_of(te_id)) %>%
    dplyr::inner_join(meta_spatial, by = "Sample_ID") %>%
    tidyr::drop_na(Count, latitude, longitude)

  fit <- try(glmmTMB::glmmTMB(Count ~ latitude + longitude, data = d,
                              family = glmmTMB::nbinom2()), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble::tibble(TE_ID = te_id, beta_lat = NA_real_, beta_long = NA_real_,
                          angle_rad = NA_real_, true_vs = NA_real_))
  }
  fx <- try(glmmTMB::fixef(fit)$cond, silent = TRUE)
  if (inherits(fx, "try-error") ||
      !all(c("latitude", "longitude") %in% names(fx))) {
    return(tibble::tibble(TE_ID = te_id, beta_lat = NA_real_, beta_long = NA_real_,
                          angle_rad = NA_real_, true_vs = NA_real_))
  }
  b_lat  <- fx[["latitude"]]; b_long <- fx[["longitude"]]
  ang    <- atan2(b_lat, b_long)
  proj   <- d$longitude * cos(ang) + d$latitude * sin(ang)
  vs     <- suppressWarnings(stats::cor(d$Count, proj, use = "complete.obs"))

  tibble::tibble(TE_ID = te_id, beta_lat = b_lat, beta_long = b_long,
                 angle_rad = ang, true_vs = vs)
}

#' Add the derived columns the published table carries.
finalise_vs <- function(df) {
  df %>%
    dplyr::filter(!is.na(beta_lat), !is.na(beta_long), !is.na(true_vs)) %>%
    dplyr::mutate(
      angle_rad_corrected = ifelse(true_vs < 0, angle_rad + pi, angle_rad),
      true_vs_corrected   = abs(true_vs),
      direction = dplyr::case_when(
        beta_lat > 0 & beta_long > 0 ~ "NE",
        beta_lat > 0 & beta_long < 0 ~ "NW",
        beta_lat < 0 & beta_long > 0 ~ "SE",
        beta_lat < 0 & beta_long < 0 ~ "SW",
        TRUE ~ "Undetermined"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(true_vs_corrected))
}

# =============================================================================
# 2. Validation gate - reproduce the published VS values for the 34
# =============================================================================
published_vs <- readr::read_csv(file.path(dir_data, "published_VS_table.csv"),
                                show_col_types = FALSE)
te_34 <- sort(unique(published_vs$TE_ID))
message(sprintf("Refitting the %d published candidates for validation ...",
                length(te_34)))

vs_34 <- finalise_vs(purrr::map_dfr(te_34, fit_vs))

check <- vs_34 %>%
  dplyr::select(TE_ID, vs_refit = true_vs_corrected, lat_refit = beta_lat,
                long_refit = beta_long) %>%
  dplyr::inner_join(
    published_vs %>%
      dplyr::select(TE_ID, vs_pub = true_vs_corrected, lat_pub = beta_lat,
                    long_pub = beta_long),
    by = "TE_ID"
  )

validation <- tibble::tibble(
  n_refit          = nrow(vs_34),
  n_published      = nrow(published_vs),
  n_matched        = nrow(check),
  vs_pearson_r     = stats::cor(check$vs_refit, check$vs_pub, use = "complete.obs"),
  vs_max_abs_diff  = max(abs(check$vs_refit - check$vs_pub), na.rm = TRUE),
  beta_lat_max_diff  = max(abs(check$lat_refit - check$lat_pub), na.rm = TRUE),
  beta_long_max_diff = max(abs(check$long_refit - check$long_pub), na.rm = TRUE)
)
readr::write_csv(validation, file.path(dir_tab, "validation_VS_reproduction.csv"))
message("Validation of the spatial-model reproduction:")
print(as.data.frame(validation))

if (is.na(validation$vs_max_abs_diff) || validation$vs_max_abs_diff > 1e-4) {
  stop("The refitted VS values do not reproduce the published table ",
       "(max abs difference = ", signif(validation$vs_max_abs_diff, 3),
       "); aborting rather than mixing two spatial models in one figure.")
}
message("  -> published VS values reproduced; the spatial model is faithful.")

# =============================================================================
# 3. Model B on the final 38 candidates
# =============================================================================
final_union <- readr::read_csv(
  file.path(dir_outputs, "02_candidate_TEs", "tables",
            "FINAL_union_candidates.csv"),
  show_col_types = FALSE
)
te_38 <- sort(unique(final_union$TE_ID))
message(sprintf("Fitting the spatial model on the final %d candidates ...",
                length(te_38)))

vs_38 <- finalise_vs(purrr::map_dfr(te_38, fit_vs)) %>%
  dplyr::left_join(
    final_union %>% dplyr::select(TE_ID, superfamily, support,
                                  in_model_A, in_model_A_PC),
    by = "TE_ID"
  ) %>%
  dplyr::mutate(
    in_published34 = TE_ID %in% te_34,
    strong_VS      = true_vs_corrected > VS_STRONG
  )

readr::write_csv(vs_38, file.path(dir_tab, "VS_table_38_candidates.csv"))

new_only <- vs_38 %>% dplyr::filter(!in_published34)
readr::write_csv(new_only, file.path(dir_tab, "VS_new_candidates_only.csv"))

message(sprintf("VS computed for %d of the %d candidates", nrow(vs_38), length(te_38)))
message(sprintf("  TEs with strong VS (> %.1f): %d of %d",
                VS_STRONG, sum(vs_38$strong_VS), nrow(vs_38)))
message("  newly added candidates:")
print(as.data.frame(new_only %>%
  dplyr::select(TE_ID, superfamily, true_vs_corrected, direction, support)))
message("  direction breakdown:")
print(as.data.frame(vs_38 %>% dplyr::count(direction, name = "n")))

# =============================================================================
# 4. Figures - the Fig. 3B panel, regenerated on the final set
# =============================================================================
p_vs <- vs_38 %>%
  dplyr::mutate(TE_ID = stats::reorder(TE_ID, true_vs_corrected)) %>%
  ggplot2::ggplot(ggplot2::aes(x = TE_ID, y = true_vs_corrected,
                               fill = superfamily)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::geom_hline(yintercept = VS_STRONG, linetype = 2, colour = "grey40") +
  ggplot2::coord_flip() +
  # Palette taken verbatim from the Figure 1 heatmaps
  # (1-PCA/250925_PCA_TE.rmd, ann_colors), so a superfamily keeps one colour
  # across the whole paper. Same specification as the map panel.
  ggplot2::scale_fill_manual(values = c(LTR_Copia = "lightgreen", LTR_Gypsy = "lightblue",
                                        LTR = "orange", noCat = "grey",
                                        TIR = "purple"),
                             labels = c(LTR_Copia = "Copia", LTR_Gypsy = "Gypsy",
                                        LTR = "LTR", noCat = "noCat", TIR = "TIR")) +
  ggplot2::labs(x = NULL, y = "Spatial vector strength (VS)", fill = "TE superfamily",
                title = sprintf("Spatial vector strength of the %d final candidates",
                                nrow(vs_38)),
                subtitle = "Dashed line = VS 0.5, the threshold described as strong in the Results") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(legend.position = "bottom",
                 axis.text.y = ggplot2::element_text(size = 6))
save_plot(p_vs, "VS_barplot_38_candidates", width = 7, height = 8)
# Saved for the Fig. 3 assembly (17_figure3_final.R), so the published panel
# and this script can never drift apart.
saveRDS(p_vs, file.path(dir_out, "panel_B_vs.rds"))

# Highlight which candidates are new relative to the published 34.
p_new <- vs_38 %>%
  dplyr::mutate(status = ifelse(in_published34, "In the published 34",
                                "Added by the union")) %>%
  ggplot2::ggplot(ggplot2::aes(x = status, y = true_vs_corrected,
                               colour = status)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.8,
                       show.legend = FALSE) +
  ggplot2::geom_hline(yintercept = VS_STRONG, linetype = 2, colour = "grey40") +
  ggplot2::scale_colour_manual(values = c("In the published 34" = "#4C72B0",
                                          "Added by the union" = "#C44E52")) +
  ggplot2::labs(x = NULL, y = "Spatial vector strength (VS)",
                title = "Do the newly added candidates show comparable spatial structure?") +
  ggplot2::theme_classic(base_size = 10)
save_plot(p_new, "VS_new_vs_published", width = 5.5, height = 5)

# =============================================================================
# 5. Session log
# =============================================================================
log_session("03_spatial_vs")
message("=== 03_spatial_vs.R finished ===")
