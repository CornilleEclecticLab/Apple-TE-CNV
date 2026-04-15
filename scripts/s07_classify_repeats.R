#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(data.table)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(plotly)
  library(FactoMineR)
  library(randomForest)
  library(pls)
  library(RColorBrewer)
  library(tibble)
  library(VennDiagram)
  library(grid)
})

# ============================================================
# Repeat classification
#
#
# Main analyses:
#   1. Compute ANOVA p-value, entropy, and FST per TE
#   2. Compute random forest importance
#   3. Build TE lists across threshold stringencies
#   4. Compare TE subsets to SNP structure using PLS
#   5. Run analyses for:
#        - population structure groups
#        - geographic groups
#   6. Export selected TE sets, summary tables, PCA, and Venn diagram
#
# Example:
#   Rscript scripts/08_classify_repeats.R \
#     --sample-info data/metadata/sample_info.csv \
#     --te-counts data/processed/te_counts_cleaned.csv \
#     --dosage-rda data/processed/dosage_data.rda \
#     --outdir results/s08_repeat_classification \
#     --species-threshold-index 8 \
#     --geo-threshold-index 6
# ============================================================

# ----------------------------
# CLI arguments
# ----------------------------
option_list <- list(
  make_option("--sample-info", type = "character", help = "Sample metadata CSV/TSV"),
  make_option("--te-counts", type = "character", help = "TE count matrix with sample IDs in first column"),
  make_option("--dosage-rda", type = "character", help = "RData file containing object 'dosage_data'"),
  make_option("--outdir", type = "character", default = "results/s08_repeat_classification",
              help = "Output directory [default: %default]"),
  make_option("--species-threshold-index", type = "integer", default = 8,
              help = "Selected threshold index for structure/species TE set [default: %default]"),
  make_option("--geo-threshold-index", type = "integer", default = 6,
              help = "Selected threshold index for geography-associated TE set [default: %default]")
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

# ----------------------------
# Validate inputs
# ----------------------------
stop_if_missing(opt$`sample-info`, "--sample-info")
stop_if_missing(opt$`te-counts`, "--te-counts")
stop_if_missing(opt$`dosage-rda`, "--dosage-rda")

# ----------------------------
# Output directories
# ----------------------------
dir_main <- opt$outdir
dir_species <- file.path(dir_main, "species_structure")
dir_geo <- file.path(dir_main, "geography")
dir_pca <- file.path(dir_main, "pca")
dir_lists <- file.path(dir_main, "lists")
dir_plots <- file.path(dir_main, "plots")

dir.create(dir_main, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_species, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_geo, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_pca, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_lists, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Parameters
# ----------------------------
fraction_threshold <- tibble::tibble(
  index = seq_len(13),
  values = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 7.5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 9.5e-1)
)

# ----------------------------
# Helper functions
# ----------------------------
compute_anova_entropy_fst <- function(column, groups) {
  df <- tibble(value = as.numeric(column), group = as.factor(groups)) %>%
    tidyr::drop_na()

  if (nrow(df) < 3 || length(unique(df$group)) < 2 || length(unique(df$value)) < 2) {
    return(c(p_value = NA_real_, entropy = NA_real_, fst = NA_real_))
  }

  anova_p <- tryCatch({
    summary(aov(value ~ group, data = df))[[1]]$"Pr(>F)"[1]
  }, error = function(e) NA_real_)

  freq <- table(df$value) / nrow(df)
  entropy <- -sum(freq * log2(freq), na.rm = TRUE)

  group_means <- tapply(df$value, df$group, mean, na.rm = TRUE)
  group_vars  <- tapply(df$value, df$group, var, na.rm = TRUE)
  overall_mean <- mean(df$value, na.rm = TRUE)

  if (length(group_means) < 2) {
    fst <- NA_real_
  } else {
    va <- sum((group_means - overall_mean)^2, na.rm = TRUE) / (length(group_means) - 1)
    vw <- mean(group_vars, na.rm = TRUE)
    fst <- ifelse(is.finite(va + vw) && (va + vw) > 0, va / (va + vw), 0)
  }

  c(p_value = anova_p, entropy = entropy, fst = fst)
}

run_random_forest <- function(data_tbl, groups, seed = 53820) {
  set.seed(seed)

  x <- data_tbl[, -1, drop = FALSE]
  y <- as.factor(groups)

  tuned <- tuneRF(
    x = x,
    y = y,
    ntreeTry = 500,
    stepFactor = 1.5,
    improve = 0.01,
    trace = TRUE,
    plot = FALSE
  )

  best_mtry <- tuned[tuned[, 2] == min(tuned[, 2]), 1][1]

  rf <- randomForest(
    x = x,
    y = y,
    importance = TRUE,
    ntree = 500,
    mtry = best_mtry
  )

  list(model = rf, best_mtry = best_mtry, tuning = tuned)
}

build_te_lists <- function(fraction_threshold, results_stats, type_label) {
  out <- list()

  for (i in seq_len(nrow(fraction_threshold))) {
    th <- fraction_threshold$values[i]
    idx <- fraction_threshold$index[i]
    n <- nrow(results_stats)

    ranked_stats <- results_stats %>%
      mutate(
        rank_p_value = ntile(p_value, n),
        rank_fst = ntile(FST, n),
        rank_gini = ntile(MeanDecreaseGini, n)
      ) %>%
      mutate(
        category = case_when(
          rank_p_value <= n * th &
            rank_fst >= (n - n * th) &
            rank_gini >= (n - n * th) ~ "Structured",

          rank_p_value > n * th &
            rank_fst < (n - n * th) &
            rank_gini < (n - n * th) ~ "Not_structured",

          TRUE ~ "Intermediate"
        )
      )

    te_selected <- ranked_stats %>%
      filter(category == "Structured") %>%
      pull(TE)

    out[[paste0("te_", type_label, "_", idx)]] <- te_selected
  }

  out
}

remove_zero_variance_cols <- function(df) {
  if (ncol(df) == 0) return(df)
  keep <- vapply(df, function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    stats::var(x_num, na.rm = TRUE) > 0
  }, logical(1))
  df[, keep, drop = FALSE]
}

pls_cca_function <- function(list_te_lists, te_data, snp_data) {
  canonical_corr <- rep(NA_real_, length(list_te_lists))
  names(canonical_corr) <- names(list_te_lists)

  if (!"sample_id" %in% names(te_data)) {
    stop("Missing 'sample_id' column in te_data.", call. = FALSE)
  }
  if (!"sample_id" %in% names(snp_data)) {
    stop("Missing 'sample_id' column in snp_data.", call. = FALSE)
  }

  for (i in seq_along(list_te_lists)) {
    dataset_name <- names(list_te_lists)[i]
    te_names <- list_te_lists[[i]]

    if (is.null(te_names) || length(te_names) == 0 || all(is.na(te_names))) {
      log_warn("No usable TE set for ", dataset_name)
      next
    }

    te_names <- intersect(te_names, names(te_data))
    if (length(te_names) == 0) {
      log_warn("No valid TE columns found for ", dataset_name)
      next
    }

    te_df <- te_data %>% select(sample_id, all_of(te_names))
    merged_df <- inner_join(te_df, snp_data, by = "sample_id")

    if (nrow(merged_df) == 0) {
      log_warn("No shared samples for ", dataset_name)
      next
    }

    te_cols <- setdiff(names(te_df), "sample_id")
    snp_cols <- setdiff(names(snp_data), "sample_id")

    x <- as.matrix(merged_df[, te_cols, drop = FALSE])
    y <- as.matrix(merged_df[, snp_cols, drop = FALSE])

    if (ncol(x) > 0) {
      var_x <- apply(x, 2, var, na.rm = TRUE)
      x <- x[, var_x > 0, drop = FALSE]
    }
    if (ncol(y) > 0) {
      var_y <- apply(y, 2, var, na.rm = TRUE)
      y <- y[, var_y > 0, drop = FALSE]
    }

    if (ncol(x) == 0 || ncol(y) == 0) {
      log_warn("No variable TE or SNP columns for ", dataset_name)
      next
    }

    x <- log1p(x)
    x_scaled <- scale(x)
    y_scaled <- scale(y)

    if (any(!is.finite(x_scaled)) || any(!is.finite(y_scaled))) {
      log_warn("Non-finite scaled values for ", dataset_name)
      next
    }

    tryCatch({
      pls_model <- plsr(y_scaled ~ x_scaled, ncomp = 1, scale = FALSE, validation = "none")
      pls1_te <- scores(pls_model)[, 1]
      pls1_snp <- pls_model$Yscores[, 1]
      canonical_corr[i] <- cor(pls1_te, pls1_snp, use = "complete.obs")
    }, error = function(e) {
      log_warn("PLS failed for ", dataset_name, ": ", e$message)
      canonical_corr[i] <- NA_real_
    })
  }

  canonical_corr
}

plot_rf_error <- function(rf_model) {
  rf_model$err.rate %>%
    as.data.frame() %>%
    rownames_to_column("tree_number") %>%
    mutate(tree_number = as.numeric(tree_number)) %>%
    pivot_longer(-tree_number, names_to = "category", values_to = "rate") %>%
    filter(category == "OOB") %>%
    ggplot(aes(x = tree_number, y = rate)) +
    geom_line() +
    labs(x = "Number of trees", y = "OOB error rate") +
    theme_classic()
}

plot_rf_confusion <- function(rf_model, n_groups) {
  rf_model$confusion %>%
    as.data.frame() %>%
    rownames_to_column("population") %>%
    pivot_longer(cols = 2:(n_groups + 1), names_to = "attribution", values_to = "n") %>%
    ggplot(aes(x = population, y = attribution)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = n)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

plot_summary_histograms <- function(results_tbl) {
  p1 <- results_tbl %>%
    ggplot(aes(x = p_value)) +
    geom_histogram() +
    scale_x_log10() +
    labs(x = "ANOVA p-value") +
    theme_minimal()

  p2 <- results_tbl %>%
    ggplot(aes(x = FST)) +
    geom_histogram() +
    labs(x = "FST") +
    theme_minimal()

  p3 <- results_tbl %>%
    ggplot(aes(x = MeanDecreaseGini)) +
    geom_histogram() +
    scale_x_log10() +
    labs(x = "Mean decrease in Gini") +
    theme_minimal()

  ggpubr::ggarrange(p1, p2, p3, nrow = 1)
}

plot_threshold_correlation <- function(canonical_corr, fraction_threshold, selected_threshold, title_text) {
  tibble(
    threshold_name = names(canonical_corr),
    corr = as.numeric(canonical_corr)
  ) %>%
    mutate(
      index = str_remove(threshold_name, "^te_[A-Za-z_]+_"),
      index = as.integer(index)
    ) %>%
    left_join(fraction_threshold, by = "index") %>%
    ggplot(aes(x = values, y = corr)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    geom_vline(xintercept = fraction_threshold$values[fraction_threshold$index == selected_threshold], linetype = 2) +
    geom_text(aes(x = values, y = corr + 0.01, label = index), size = 3, na.rm = TRUE) +
    scale_x_log10() +
    labs(
      x = "Threshold value",
      y = "PLS correlation",
      title = title_text
    ) +
    theme_classic()
}

run_feature_screening <- function(data_tbl, groups_vec) {
  results_a <- t(apply(data_tbl[, -1, drop = FALSE], 2, compute_anova_entropy_fst, groups = groups_vec))
  results_a <- as.data.frame(results_a)
  colnames(results_a) <- c("p_value", "Entropy", "FST")
  results_a$TE <- colnames(data_tbl)[-1]

  rf_out <- run_random_forest(data_tbl, groups_vec)

  te_importance <- as.data.frame(rf_out$model$importance)
  te_importance$TE <- rownames(te_importance)

  results_b <- full_join(results_a, te_importance, by = "TE")

  list(
    stats = results_b,
    rf = rf_out$model,
    best_mtry = rf_out$best_mtry,
    tuning = rf_out$tuning
  )
}

make_pca_plotly <- function(pca_obj, groups_df, color_col, symbol_col, title_text) {
  pca_res <- data.frame(pca_obj$ind$coord[, 1:3], groups_df)

  num_colors <- length(unique(groups_df[[color_col]]))
  color_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(max(3, num_colors))

  symbol_list <- c("circle", "square", "diamond", "circle-open", "diamond-open", "x", "square-open")

  plot_ly(pca_res) %>%
    add_trace(
      x = ~Dim.1,
      y = ~Dim.2,
      z = ~Dim.3,
      symbol = as.formula(paste0("~`", symbol_col, "`")),
      symbols = symbol_list,
      color = as.formula(paste0("~`", color_col, "`")),
      colors = color_palette,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 6, opacity = 0.8),
      showlegend = TRUE
    ) %>%
    layout(
      title = title_text,
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      ),
      legend = list(title = list(text = ""))
    )
}

save_te_list <- function(te_vec, file_path) {
  if (length(te_vec) == 0) {
    file.create(file_path)
  } else {
    writeLines(te_vec, file_path)
  }
}

# ----------------------------
# Load inputs
# ----------------------------
log_info("Loading metadata and TE count matrix...")

sample_info <- suppressMessages(fread(opt$`sample-info`)) %>% as_tibble()
te_counts <- suppressMessages(fread(opt$`te-counts`)) %>% as.data.frame()

if (!"V1" %in% names(te_counts)) {
  names(te_counts)[1] <- "V1"
}

load(opt$`dosage-rda`)
if (!exists("dosage_data")) {
  stop("Object 'dosage_data' not found in ", opt$`dosage-rda`, call. = FALSE)
}

dosage_data <- dosage_data %>%
  as.data.frame() %>%
  dplyr::select(where(~ !any(is.na(.)))) %>%
  tibble::rownames_to_column("sample_id")

log_ok("Inputs loaded.")

# ----------------------------
# Required metadata columns
# ----------------------------
required_cols <- c("Sample_ID", "Code_group_geo_residual")
missing_cols <- setdiff(required_cols, names(sample_info))
if (length(missing_cols) > 0) {
  stop("Missing required metadata column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

structure_col <- names(sample_info)[str_detect(names(sample_info), "group_structure_CORE_admix")]
if (length(structure_col) == 0) {
  stop("Could not find a metadata column matching 'group_structure_CORE_admix'.", call. = FALSE)
}
structure_col <- structure_col[1]

# ============================================================
# 1. Structure/species-associated TE classification
# ============================================================
log_info("Running structure/species-associated TE classification...")

groups_species <- sample_info %>%
  select(Sample_ID, cluster = all_of(structure_col)) %>%
  filter(!grepl("Admix", cluster))

data_species <- te_counts %>%
  filter(V1 %in% groups_species$Sample_ID)

groups_species <- groups_species %>% arrange(Sample_ID)
data_species <- data_species %>% arrange(V1)

common_species <- intersect(groups_species$Sample_ID, data_species$V1)
groups_species <- groups_species %>% filter(Sample_ID %in% common_species) %>% arrange(Sample_ID)
data_species <- data_species %>% filter(V1 %in% common_species) %>% arrange(V1)

stopifnot(all(groups_species$Sample_ID == data_species$V1))

screen_species <- run_feature_screening(data_species, groups_species$cluster)

readr::write_csv(screen_species$stats, file.path(dir_species, "te_statistics_species.csv"))

p_rf_error_species <- plot_rf_error(screen_species$rf)
ggsave(file.path(dir_species, "rf_error_species.pdf"), p_rf_error_species, width = 6, height = 4)

n_species_groups <- length(unique(groups_species$cluster))
p_rf_conf_species <- plot_rf_confusion(screen_species$rf, n_species_groups)
ggsave(file.path(dir_species, "rf_confusion_species.pdf"), p_rf_conf_species, width = 7, height = 5)

p_hist_species <- plot_summary_histograms(screen_species$stats)
ggsave(file.path(dir_species, "summary_histograms_species.pdf"), p_hist_species, width = 12, height = 4)

dosage_species <- dosage_data %>%
  filter(sample_id %in% groups_species$Sample_ID) %>%
  arrange(sample_id)

te_lists_species <- build_te_lists(
  fraction_threshold = fraction_threshold,
  results_stats = screen_species$stats,
  type_label = "species"
)

corr_species <- pls_cca_function(
  list_te_lists = te_lists_species,
  te_data = data_species %>% rename(sample_id = V1),
  snp_data = dosage_species
)

corr_species_tbl <- tibble(
  threshold_name = names(corr_species),
  correlation = as.numeric(corr_species)
)
readr::write_csv(corr_species_tbl, file.path(dir_species, "pls_correlations_species.csv"))

p_corr_species <- plot_threshold_correlation(
  canonical_corr = corr_species,
  fraction_threshold = fraction_threshold,
  selected_threshold = opt$`species-threshold-index`,
  title_text = "PLS correlation across threshold stringencies (species/structure)"
)

ggsave(file.path(dir_species, "pls_correlation_species.pdf"), p_corr_species, width = 7, height = 5)

te_structure <- te_lists_species[[paste0("te_species_", opt$`species-threshold-index`)]]
save_te_list(te_structure, file.path(dir_lists, "te_structure.txt"))

# ============================================================
# 2. Geography-associated TE classification
# ============================================================
log_info("Running geography-associated TE classification...")

groups_geo <- sample_info %>%
  group_by(Code_group_geo_residual) %>%
  filter(n() >= 3) %>%
  ungroup()

n_geo_groups <- length(unique(groups_geo$Code_group_geo_residual))

data_geo <- te_counts %>%
  filter(V1 %in% groups_geo$Sample_ID) %>%
  select(V1, -all_of(intersect(te_structure, names(.))))

groups_geo <- groups_geo %>% arrange(Sample_ID)
data_geo <- data_geo %>% arrange(V1)

common_geo <- intersect(groups_geo$Sample_ID, data_geo$V1)
groups_geo <- groups_geo %>% filter(Sample_ID %in% common_geo) %>% arrange(Sample_ID)
data_geo <- data_geo %>% filter(V1 %in% common_geo) %>% arrange(V1)

stopifnot(all(groups_geo$Sample_ID == data_geo$V1))

screen_geo <- run_feature_screening(data_geo, groups_geo$Code_group_geo_residual)

readr::write_csv(screen_geo$stats, file.path(dir_geo, "te_statistics_geography.csv"))

p_rf_error_geo <- plot_rf_error(screen_geo$rf)
ggsave(file.path(dir_geo, "rf_error_geography.pdf"), p_rf_error_geo, width = 6, height = 4)

p_rf_conf_geo <- plot_rf_confusion(screen_geo$rf, n_geo_groups)
ggsave(file.path(dir_geo, "rf_confusion_geography.pdf"), p_rf_conf_geo, width = 8, height = 5)

p_hist_geo <- plot_summary_histograms(screen_geo$stats)
ggsave(file.path(dir_geo, "summary_histograms_geography.pdf"), p_hist_geo, width = 12, height = 4)

dosage_geo <- dosage_data %>%
  filter(sample_id %in% groups_geo$Sample_ID) %>%
  arrange(sample_id)

te_lists_geo <- build_te_lists(
  fraction_threshold = fraction_threshold,
  results_stats = screen_geo$stats,
  type_label = "geo"
)

corr_geo <- pls_cca_function(
  list_te_lists = te_lists_geo,
  te_data = data_geo %>% rename(sample_id = V1),
  snp_data = dosage_geo
)

corr_geo_tbl <- tibble(
  threshold_name = names(corr_geo),
  correlation = as.numeric(corr_geo)
)
readr::write_csv(corr_geo_tbl, file.path(dir_geo, "pls_correlations_geography.csv"))

p_corr_geo <- plot_threshold_correlation(
  canonical_corr = corr_geo,
  fraction_threshold = fraction_threshold,
  selected_threshold = opt$`geo-threshold-index`,
  title_text = "PLS correlation across threshold stringencies (geography)"
)

ggsave(file.path(dir_geo, "pls_correlation_geography.pdf"), p_corr_geo, width = 7, height = 5)

te_geo <- te_lists_geo[[paste0("te_geo_", opt$`geo-threshold-index`)]]
save_te_list(te_geo, file.path(dir_lists, "te_geography.txt"))

# ============================================================
# 3. PCA on selected TE subsets
# ============================================================
log_info("Running PCA on selected TE subsets...")

te_structure <- readLines(file.path(dir_lists, "te_structure.txt"), warn = FALSE)
te_geo <- readLines(file.path(dir_lists, "te_geography.txt"), warn = FALSE)

pca_groups <- sample_info %>%
  group_by(Code_group_geo_residual) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  arrange(Sample_ID)

pca_data <- te_counts %>%
  filter(V1 %in% pca_groups$Sample_ID) %>%
  arrange(V1)

stopifnot(all(pca_groups$Sample_ID == pca_data$V1))

te_structure_only <- setdiff(te_structure, te_geo)
te_geo_only <- setdiff(te_geo, te_structure)

te_top_struct <- pca_data %>%
  select(V1, any_of(te_structure_only))

te_top_geo <- pca_data %>%
  select(V1, any_of(te_geo_only))

if (ncol(te_top_struct) > 1) {
  mat_struct <- te_top_struct %>% select(-V1) %>% remove_zero_variance_cols()
  if (ncol(mat_struct) >= 2) {
    pca_struct <- FactoMineR::PCA(mat_struct, graph = FALSE)
    htmlwidgets::saveWidget(
      as_widget(make_pca_plotly(
        pca_obj = pca_struct,
        groups_df = pca_groups,
        color_col = "Code_group_geo_residual",
        symbol_col = structure_col,
        title_text = "PCA based on species/structure-associated TE"
      )),
      file = file.path(dir_pca, "pca_species_structure_te.html"),
      selfcontained = TRUE
    )
  }
}

if (ncol(te_top_geo) > 1) {
  mat_geo <- te_top_geo %>% select(-V1) %>% remove_zero_variance_cols()
  if (ncol(mat_geo) >= 2) {
    pca_geo <- FactoMineR::PCA(mat_geo, graph = FALSE)
    htmlwidgets::saveWidget(
      as_widget(make_pca_plotly(
        pca_obj = pca_geo,
        groups_df = pca_groups,
        color_col = "Code_group_geo_residual",
        symbol_col = structure_col,
        title_text = "PCA based on geography-associated TE"
      )),
      file = file.path(dir_pca, "pca_geography_te.html"),
      selfcontained = TRUE
    )
  }
}

# ============================================================
# 4. Venn diagram
# ============================================================
log_info("Building Venn diagram...")

venn.plot <- venn.diagram(
  x = list(
    Species = te_structure,
    Species_geo = te_geo
  ),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  cat.cex = 1,
  cat.pos = c(-10, 10)
)

pdf(file.path(dir_plots, "venn_species_vs_geography.pdf"), width = 6, height = 6)
grid::grid.draw(venn.plot)
dev.off()

# ============================================================
# 5. Export overlap summary
# ============================================================
overlap_tbl <- tibble(
  category = c("species_only", "geography_only", "shared"),
  n = c(
    length(setdiff(te_structure, te_geo)),
    length(setdiff(te_geo, te_structure)),
    length(intersect(te_structure, te_geo))
  )
)

readr::write_csv(overlap_tbl, file.path(dir_main, "overlap_summary.csv"))

save_te_list(setdiff(te_structure, te_geo), file.path(dir_lists, "te_species_only.txt"))
save_te_list(setdiff(te_geo, te_structure), file.path(dir_lists, "te_geography_only.txt"))
save_te_list(intersect(te_structure, te_geo), file.path(dir_lists, "te_shared.txt"))

# ============================================================
# 6. Session info
# ============================================================
writeLines(capture.output(sessionInfo()), file.path(dir_main, "session_info.txt"))

log_ok("S08 repeat classification completed successfully.")