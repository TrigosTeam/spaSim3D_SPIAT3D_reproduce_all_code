# SPIAT3D ----
calculate_all_gradient_cc_metrics3D <- function(spe,
                                                reference_cell_type,
                                                target_cell_types,
                                                radii,
                                                feature_colname,
                                                plot_image = T) {
  
  # Define constants
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  mixing_score_df_colnames <- c("ref_cell_type",
                                "tar_cell_type",
                                "n_ref_cells",
                                "n_tar_cells",
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions",
                                "mixing_score",
                                "normalised_mixing_score")
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "neighbourhood_counts" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "neighbourhood_entropy" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["neighbourhood_counts"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["neighbourhood_entropy"]]) <- target_cell_types
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define individual data frames for mixing_score and cross_G
  for (target_cell_type in target_cell_types) {
    if (reference_cell_type != target_cell_type) {
      result[["mixing_score"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(mixing_score_df_colnames)))
      colnames(result[["mixing_score"]][[target_cell_type]]) <- mixing_score_df_colnames
    }
    result[["cross_G"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(cross_G_df_colnames)))
    colnames(result[["cross_G"]][[target_cell_type]]) <- cross_G_df_colnames
  }
  
  # Get gradient results for each metric
  for (i in seq(length(radii))) {
    df <- calculate_all_single_radius_cc_metrics3D(spe,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["neighbourhood_counts"]]$ref_cell_id <- NULL
    
    result[["neighbourhood_counts"]][i, ] <- apply(df[["neighbourhood_counts"]], 2, mean)
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["neighbourhood_entropy"]][i, ] <- apply(df[["neighbourhood_entropy"]][ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
    result[["cross_K"]][i, ] <- df[["cross_K"]]
    result[["cross_L"]][i, ] <- df[["cross_L"]]
    result[["co_occurrence"]][i, ] <- df[["co_occurrence"]]
    
    for (target_cell_type in names(df[["mixing_score"]])) {
      result[["mixing_score"]][[target_cell_type]][i, ] <- df[["mixing_score"]][[target_cell_type]]
    }
    for (target_cell_type in names(df[["cross_G"]])) {
      result[["cross_G"]][[target_cell_type]][i, ] <- df[["cross_G"]][[target_cell_type]]
    }
  }
  
  # Add radius column to each data frame
  result[["neighbourhood_counts"]]$radius <- radii
  result[["cells_in_neighbourhood"]]$radius <- radii
  result[["neighbourhood_entropy"]]$radius <- radii
  result[["cross_K"]]$radius <- radii
  result[["cross_L"]]$radius <- radii
  result[["co_occurrence"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  for (target_cell_type in names(df[["cross_G"]])) {
    result[["cross_G"]][[target_cell_type]]$radius <- radii
  }
  
  
  ## Plot
  if (plot_image) {
    fig_ANC <- plot_neighbourhood_counts_gradient3D(result[["neighbourhood_counts"]], reference_cell_type)
    methods::show(fig_ANC)
    
    fig_ACIN <- plot_cells_in_neighbourhood_gradient3D(result[["cells_in_neighbourhood"]], reference_cell_type)
    methods::show(fig_ACIN)
    
    fig_ANE <- plot_neighbourhood_entropy_gradient3D(result[["neighbourhood_entropy"]], reference_cell_type)
    methods::show(fig_ANE)
    
    for (target_cell_type in names(result[["mixing_score"]])) {
      fig_NMS <- plot_mixing_scores_gradient3D(result[["mixing_score"]][[target_cell_type]], "NMS")
      fig_MS <- plot_mixing_scores_gradient3D(result[["mixing_score"]][[target_cell_type]], "MS")
      fig_NMS_MS <- cowplot::plot_grid(fig_NMS, fig_MS, nrow = 2)
      methods::show(fig_NMS_MS)
    }
    fig_CK <- plot_cross_K_gradient3D(result[["cross_K"]])
    methods::show(fig_CK)
    
    fig_CL <- plot_cross_L_gradient3D(result[["cross_L"]])
    methods::show(fig_CL)
    
    for (target_cell_type in names(result[["cross_G"]])) {
      fig_CG <- plot_cross_G_gradient3D(result[["cross_G"]][[target_cell_type]], reference_cell_type, target_cell_type)
      methods::show(fig_CG)
    }
    
    fig_COO <- plot_co_occurrence_gradient3D(result[["co_occurrence"]])
    methods::show(fig_COO)
  }
  
  return(result)
}

calculate_all_single_radius_cc_metrics3D <- function(spe,
                                                     reference_cell_type,
                                                     target_cell_types,
                                                     radius,
                                                     feature_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop(paste(radius, " is not a positive numeric."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Define result
  result <- list("neighbourhood_counts" = list(),
                 "cells_in_neighbourhood" = list(),
                 "neighbourhood_entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list(),
                 "cross_L" = list(),
                 "cross_G" = list(),
                 "co_occurrence" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type",
                                "tar_cell_type",
                                "n_ref_cells",
                                "n_tar_cells",
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions",
                                "mixing_score",
                                "normalised_mixing_score")
  
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  # Get rough dimensions of window for cross_K
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  length <- (max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- (max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- (max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # All single radius cc metrics stem from calculate_neighbourhood_entropy3D function
  neighbourhood_entropy_df <- calculate_neighbourhood_entropy3D(spe,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radius,
                                                                feature_colname)
  
  ## Neighbourhood counts ----------
  result[["neighbourhood_counts"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_prop", sep = ""))]
  
  ## Neighbourhood_entropy --------------
  result[["neighbourhood_entropy"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_entropy", sep = ""))]
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(neighbourhood_entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(neighbourhood_entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
  }
  
  ## Cross_K ---------------------
  cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
  colnames(cross_K_df) <- cross_K_df_colnames
  cross_K_df$reference <- reference_cell_type
  cross_K_df$expected <- (4 / 3) * pi * radius^3
  
  for (target_cell_type in target_cell_types) {
    cross_K_df[[target_cell_type]] <- (((volume * sum(neighbourhood_entropy_df[[target_cell_type]])) / sum(spe[[feature_colname]] == reference_cell_type)) / sum(spe[[feature_colname]] == target_cell_type))
  }
  result[["cross_K"]] <- cross_K_df
  
  ## Cross_L ---------------------
  cross_L_df <- cross_K_df
  cross_L_df[ , c("expected", target_cell_types)] <- (cross_L_df[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  result[["cross_L"]] <- cross_L_df
  
  ## Cross_G ---------------------
  for (target_cell_type in target_cell_types) {
    cross_G_df <- data.frame(matrix(nrow = 1, ncol = length(cross_G_df_colnames)))
    colnames(cross_G_df) <- cross_G_df_colnames
    
    reference_target_interactions <- neighbourhood_entropy_df[[target_cell_type]]
    n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / volume
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  }
  
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spe[[feature_colname]])
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  n_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe),
                                                  eps = radius,
                                                  query = SpatialExperiment::spatialCoords(spe),
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    
    # Account for case where target cell type is not present in data
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      co_occurrence_df[[target_cell_type]] <- NA
      next
    }
    
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe)[spe[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = SpatialExperiment::spatialCoords(spe),
                                                           sort = FALSE)
    n_target_cells_in_all_cell_type_radius <- sum(rapply(n_target_cells_in_all_cell_type_radius$id, base::length))
    
    target_cell_type_proportion_in_all_cell_type_radius <- n_target_cells_in_all_cell_type_radius / n_cells_in_all_cell_type_radius
    
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_all_cell_type_radius
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}

calculate_cell_proportion_grid_metrics3D <- function(spe,
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname,
                                                     plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics3D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^3
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  # Fill in the result data frame
  if (length(reference_cell_types) == 1) {
    result$reference <- grid_prism_cell_matrix[[reference_cell_types]]
  }
  else {
    result$reference <- rowSums(grid_prism_cell_matrix[ , reference_cell_types])
  }
  if (length(target_cell_types) == 1) {
    result$target <- grid_prism_cell_matrix[[target_cell_types]]
  }
  else {
    result$target <- rowSums(grid_prism_cell_matrix[ , target_cell_types])
  }
  result$total <- result$reference + result$target
  result$proportion <- result$target / result$total
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}

calculate_cell_proportions3D <- function(spe,
                                         cell_types_of_interest = NULL,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (ncol(spe) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spe[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Subset for cell types chosen by user
    cell_proportions <- cell_proportions[(cell_proportions$cell_type %in% cell_types_of_interest), ]
    
  }
  
  # Get frequency total for all cells
  cell_type_frequency_total <- sum(cell_proportions$frequency)
  
  # Get proportions and percentages
  cell_proportions$proportion <- cell_proportions$frequency / cell_type_frequency_total
  cell_proportions$percentage <- cell_proportions$proportion * 100
  
  # Order the cell types by proportion (highest cell proportion is first)
  cell_proportions <- cell_proportions[rev(order(cell_proportions$proportion)), ]
  rownames(cell_proportions) <- seq(nrow(cell_proportions))
  
  # Plot
  if (plot_image) {
    
    labels <- paste(round(cell_proportions$percentage, 1), "%", sep = "")
    
    fig <- ggplot(cell_proportions, aes(x = factor(cell_type, cell_type), y = percentage, fill = cell_type)) +
      geom_bar(stat='identity') +
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(fig)
  }
  
  return(cell_proportions)
}

calculate_cells_in_neighbourhood_gradient3D <- function(spe,
                                                        reference_cell_type,
                                                        target_cell_types,
                                                        radii,
                                                        feature_colname,
                                                        plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                                        reference_cell_type,
                                                                                        target_cell_types,
                                                                                        radii[i],
                                                                                        feature_colname)
    
    if (is.null(cell_proportions_neighbourhood_proportions_df)) return(NULL)
    
    result[i, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cells_in_neighbourhood3D <- function(spe,
                                               reference_cell_type,
                                               target_cell_types,
                                               radius,
                                               feature_colname) {
  
  ## Get cells in neighbourhood df
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              c(reference_cell_type, target_cell_types),
                                                              radius,
                                                              feature_colname)
  
  if (is.null(neighbourhood_counts_df)) return(NULL)
  
  neighbourhood_counts_df[ , paste(target_cell_types, "_prop", sep = "")] <-
    neighbourhood_counts_df[ , target_cell_types] / (neighbourhood_counts_df[ , target_cell_types] + neighbourhood_counts_df[ , reference_cell_type])
  
  # If reference cell type is in target cell types, proportion should be 1
  if (reference_cell_type %in% target_cell_types) {
    neighbourhood_counts_df[neighbourhood_counts_df[[reference_cell_type]] != 0, paste(reference_cell_type, "_prop", sep = "")] <- 1
  }
  
  return(neighbourhood_counts_df)
}

calculate_co_occurrence_gradient3D <- function(spe,
                                               reference_cell_type,
                                               target_cell_types,
                                               radii,
                                               feature_colname,
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types) + 1))
  colnames(result) <- c("reference", target_cell_types)
  
  for (i in seq(length(radii))) {
    co_occurrence_df <- calculate_co_occurrence3D(spe,
                                                  reference_cell_type,
                                                  target_cell_types,
                                                  radii[i],
                                                  feature_colname)
    
    result[i, ] <- co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_co_occurrence_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_co_occurrence3D <- function(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radius,
                                      feature_colname) {
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get all cell types in spe
  all_cell_types <- unique(spe[[feature_colname]])
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  # Get total number of cells in radius around all cell types
  n_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe),
                                                  eps = radius,
                                                  query = SpatialExperiment::spatialCoords(spe),
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    
    # Account for case where target cell type is not present in data
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      co_occurrence_df[[target_cell_type]] <- NA
      next
    }
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get total number of target cells in radius around all cell types
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe)[spe[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = SpatialExperiment::spatialCoords(spe),
                                                           sort = FALSE)
    n_target_cells_in_all_cell_type_radius <- sum(rapply(n_target_cells_in_all_cell_type_radius$id, base::length))
    
    # Get proportion of target cells in radius around all cell types
    target_cell_type_proportion_in_all_cell_type_radius <- n_target_cells_in_all_cell_type_radius / n_cells_in_all_cell_type_radius
    
    # Get co-occurence value for target cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_all_cell_type_radius
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}

calculate_cross_G_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_type,
                                         radii,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check if radii input is valid
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 2))
  colnames(result) <- c("observed_cross_G",
                        "expected_cross_G")
  
  # Iterate through each radii
  for (i in seq(length(radii))) {
    cross_G_df <- calculate_cross_G3D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_G_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cross_G_gradient3D(result, reference_cell_type, target_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cross_G3D <- function(spe,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname) {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              target_cell_type,
                                                              radius,
                                                              feature_colname)
  
  reference_target_interactions <- neighbourhood_counts_df[[target_cell_type]]
  
  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  
  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  length <- (max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- (max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- (max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # Get the number of target cells
  n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}

calculate_cross_K_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check if radii input is valid
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  # Iterate through each radii
  for (i in seq(length(radii))) {
    cross_K_df <- calculate_cross_K3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cross_K_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cross_K3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Check if inputs are valid
  if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  
  ## Get expected cross K-function
  expected_cross_K <- (4/3) * pi * radius^3
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  length <- (max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- (max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- (max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              target_cell_types,
                                                              radius,
                                                              feature_colname)
  
  # Calculate cross K-fucnction for each target cell type
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
    
    ## Get observed cross K-function
    if (n_tar_cells == 0) {
      observed_cross_K <- NA
    }
    else {
      observed_cross_K <- (volume * n_ref_tar_interactions) / (n_ref_cells * n_tar_cells)
    }
    result[[target_cell_type]] <- observed_cross_K
  }
  
  return(result)
}

calculate_cross_L_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_L_df <- calculate_cross_L3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_L_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cross_L_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cross_L3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Calculate cross K first
  result <- calculate_cross_K3D(spe = spe,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  # Linearise cross K to get cross L
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  
  return(result)
}

calculate_entropy_background3D <- function(spe,
                                           cell_types_of_interest,
                                           feature_colname) {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions3D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy)
}

calculate_entropy_grid_metrics3D <- function(spe,
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname,
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spe[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics3D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^3
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  for (cell_type in cell_types_of_interest) {
    result[[cell_type]] <- grid_prism_cell_matrix[[cell_type]]
  }
  result$total <- rowSums(result)
  
  ## Get data frame containing proportions for cell_types_of_interest
  df_props <- result[ , cell_types_of_interest] / result$total
  
  ## Use proportion data frame to get entropy
  calculate_entropy <- function(x) {
    entropy <- -1 * sum(x * ifelse(is.infinite(log(x, 2)), 0, log(x, 2))) / log(length(x), 2)
    return(entropy)
  }
  result$entropy <- apply(df_props, 1, calculate_entropy)
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}

calculate_minimum_distances_between_cell_types3D <- function(spe,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname,
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Get spe coords
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spe[[feature_colname]] == cell_type1) == 0 || sum(spe[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spe_coords[spe[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spe_coords[spe[[feature_colname]] == cell_type2, ]
    
    # Find all of closest points
    # For each cell of cell_type1, find the closest cell of cell_type2
    if (cell_type1 != cell_type2) {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords,
                                      query = cell_type1_coords,
                                      k = 1)
    }
    # If we are comparing the same cell_type, and there is only one of this cell type, move on
    else if (nrow(cell_type1_coords) == 1) {
      warning("There is only 1 '", cell_type1, "' cell in your data. It has no nearest neighbour of the same cell type.", sep = "")
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    # If we are comparing the same cell_type, use the second closest neighbour
    else {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords,
                                      query = cell_type1_coords,
                                      k = 2)
      
      # Get cell IDs and distances
      nearest_neighbours[['nn.idx']] <- nearest_neighbours[['nn.idx']][ , 2]
      nearest_neighbours[['nn.dists']] <- nearest_neighbours[['nn.dists']][ , 2]
    }
    
    # Create the data frame containing the chosen cells and their ids, as well as the nearest cell to them and their ids, and the distance between
    df <- data.frame(
      ref_cell_id = cell_type_ids[[cell_type1]],
      ref_cell_type = cell_type1,
      nearest_cell_id = cell_type_ids[[cell_type2]][c(nearest_neighbours$nn.idx)],
      nearest_cell_type = cell_type2,
      distance = nearest_neighbours$nn.dists
    )
    result <- rbind(result, df)
  }
  
  # Add cell pair column
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types3D(result))
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin3D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_mixing_scores_gradient3D <- function(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii,
                                               feature_colname,
                                               plot_image = TRUE) {
  
  # Check radii input
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 8))
  colnames(result) <- c("ref_cell_type",
                        "tar_cell_type",
                        "n_ref_cells",
                        "n_tar_cells",
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions",
                        "mixing_score",
                        "normalised_mixing_score")
  
  # Iterate through each radius
  for (i in seq(length(radii))) {
    mixing_scores <- calculate_mixing_scores3D(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii[i],
                                               feature_colname)
    
    result[i, ] <- mixing_scores
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_mixing_scores_gradient3D(result, "NMS")
    fig2 <- plot_mixing_scores_gradient3D(result, "MS")
    combined_fig <- cowplot::plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}

calculate_mixing_scores3D <- function(spe,
                                      reference_cell_types,
                                      target_cell_types,
                                      radius,
                                      feature_colname) {
  
  # Define result
  result <- data.frame()
  
  # Iterate through each reference cell type
  for (reference_cell_type in reference_cell_types) {
    
    # Iterate through each target cell type
    for (target_cell_type in target_cell_types) {
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Get number of reference cells and target cells
      n_ref <- sum(spe[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spe[[feature_colname]] == target_cell_type)
      
      
      # Can't get mixing scores if there are 0 or 1 reference cells
      if (n_ref == 0 || n_ref == 1) {
        result <-  rbind(result,
                         c(reference_cell_type,
                           target_cell_type,
                           n_ref,
                           n_tar,
                           0,
                           0,
                           NA,
                           NA))
      }
      
      
      ## Get cells in neighbourhood df
      neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                                  reference_cell_type,
                                                                  c(reference_cell_type, target_cell_type),
                                                                  radius,
                                                                  feature_colname)
      
      # Get number of ref-ref interactions
      # Halve it to avoid counting each ref-ref interaction twice
      n_ref_ref_interactions <- 0.5 * sum(neighbourhood_counts_df[[reference_cell_type]])
      
      # Get number of ref-tar interactions
      n_ref_tar_interactions <- sum(neighbourhood_counts_df[[target_cell_type]])
      
      
      # Can't get mixing scores if there are no target cells
      if (n_tar == 0) {
        
        result <-  rbind(result,
                         c(reference_cell_type,
                           target_cell_type,
                           n_ref,
                           0,
                           0,
                           n_ref_ref_interactions,
                           NA,
                           NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        
        if (n_ref_ref_interactions != 0) {
          mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
          normalised_mixing_score <- 0.5 * mixing_score * n_ref / n_tar
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        result <-  rbind(result,
                         c(reference_cell_type,
                           target_cell_type,
                           n_ref,
                           n_tar,
                           n_ref_tar_interactions,
                           n_ref_ref_interactions,
                           mixing_score,
                           normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  colnames(result) <- c("ref_cell_type",
                        "tar_cell_type",
                        "n_ref_cells",
                        "n_tar_cells",
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions",
                        "mixing_score",
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}

calculate_neighbourhood_counts_gradient3D <- function(spe,
                                                      reference_cell_type,
                                                      target_cell_types,
                                                      radii,
                                                      feature_colname,
                                                      plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radii[i],
                                                                feature_colname)
    
    if (is.null(neighbourhood_counts_df)) return(NULL)
    
    neighbourhood_counts_df$ref_cell_id <- NULL
    result[i, ] <- apply(neighbourhood_counts_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_neighbourhood_counts_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_neighbourhood_counts3D <- function(spe,
                                             reference_cell_type,
                                             target_cell_types,
                                             radius,
                                             feature_colname) {
  
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # Get spe coords
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spe_coords[spe[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(ref_cell_id = spe$Cell.ID[spe[[feature_colname]] == reference_cell_type])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spe_coords[spe[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells specified distance for each reference cell
    ref_tar_result <- dbscan::frNN(target_cell_type_coords,
                                   eps = radius,
                                   query = reference_cell_type_coords,
                                   sort = FALSE)
    
    n_targets <- rapply(ref_tar_result$id, length)
    
    
    # Don't want to include the reference cell as one of the target cells
    if (reference_cell_type == target_cell_type) n_targets <- n_targets - 1
    
    ## Add to data frame
    result[[target_cell_type]] <- n_targets
  }
  
  return(result)
}

calculate_neighbourhood_entropy_gradient3D <- function(spe,
                                                       reference_cell_type,
                                                       target_cell_types,
                                                       radii,
                                                       feature_colname,
                                                       plot_image = TRUE) {
  
  # Check radii input
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  # Iterate through each radii
  for (i in seq(length(radii))) {
    neighbourhood_entropy_df <- calculate_neighbourhood_entropy3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radii[i],
                                                                  feature_colname)
    
    if (is.null(neighbourhood_entropy_df)) return(NULL)
    
    result[i, ] <- apply(neighbourhood_entropy_df[ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_neighbourhood_entropy_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_neighbourhood_entropy3D <- function(spe,
                                              reference_cell_type,
                                              target_cell_types,
                                              radius,
                                              feature_colname) {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  ## Get neighbourhood_entropy for target_cell_type
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_entropy", sep = "")] <-
    -1 *
    (cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] * log(cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")], 2) +
       (1 - cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")]) * log(1 - (cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")]), 2))
  
  return(cells_in_neighbourhood_df)
}

calculate_pairwise_distances_between_cell_types3D <- function(spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname,
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(SpatialExperiment::spatialCoords(spe))
  rownames(distance_matrix) <- spe$Cell.ID
  colnames(distance_matrix) <- spe$Cell.ID
  
  result <- data.frame()
  
  for (i in seq(length(cell_types_of_interest))) {
    
    for (j in i:length(cell_types_of_interest)) {
      
      # Get current cell types and cell ids
      cell_type1 <- names(cell_type_ids)[i]
      cell_type2 <- names(cell_type_ids)[j]
      
      cell_type1_ids <- cell_type_ids[[cell_type1]]
      cell_type2_ids <- cell_type_ids[[cell_type2]]
      
      ## Don't have a cell type, or the same cell type with only one cell
      if (length(cell_type1_ids) == 0 || length(cell_type2_ids) == 0) {
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      ## Same cell type only one cell
      if (cell_type1 == cell_type2 && length(cell_type1_ids) == 1) {
        warning("There is only 1 '", cell_type1, "' cell in your data. It has no pair of the same cell type.", sep = "")
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      # Subset distance_matrix for current cell types
      distance_matrix_subset <- distance_matrix[rownames(distance_matrix) %in% cell_type1_ids,
                                                colnames(distance_matrix) %in% cell_type2_ids]
      
      ## Different cell types, each only has one cell
      if (length(cell_type1_ids) == 1 && length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        rownames(distance_matrix_subset) <- cell_type1_ids
        colnames(distance_matrix_subset) <- cell_type2_ids
      }
      ## Different cell types, only one cell of cell_type1
      else if (length(cell_type1_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type1_ids
      }
      ## Different cell types, only one cell of cell_type2
      else if (length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type2_ids
      }
      ## Same cell type, only need part of the matrix (make irrelevant part of matrix equal to NA)
      if (cell_type1 == cell_type2) distance_matrix_subset[upper.tri(distance_matrix_subset, diag = TRUE)] <- NA
      
      # Convert distance_matrix_subset to a data frame
      df <- reshape2::melt(distance_matrix_subset, na.rm = TRUE)
      df$cell_type1 <- cell_type1
      df$cell_type2 <- cell_type2
      df$pair <- paste(cell_type1, cell_type2, sep="/")
      
      result <- rbind(result, df)
    }
  }
  
  # Rearrange columns
  colnames(result)[c(1, 2, 3)] <- c("cell_type1_id", "cell_type2_id", "distance")
  result <- result[ , c("cell_type1_id", "cell_type1", "cell_type2_id", "cell_type2", "distance", "pair")]
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types3D(result))
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin3D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_prevalence_gradient3D <- function(grid_metrics,
                                            metric_colname,
                                            show_AUC = TRUE,
                                            plot_image = TRUE) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!is.logical(show_AUC)) {
    stop("`show_AUC` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  # Get prevalences for each threshold
  result$prevalence <- sapply(thresholds, function(threshold) {
    calculate_prevalence3D(grid_metrics, metric_colname, threshold)
  })
  
  # Show AUC of prevalence gradient graph
  if (show_AUC) {
    AUC_value <- sum(diff(result$threshold) * (head(result$prevalence, -1) + tail(result$prevalence, -1)) / 2)
    print(paste("AUC:", round(AUC_value, 2)))
  }
  
  # Plot
  if (plot_image) {
    fig <- ggplot(result, aes(threshold, prevalence)) +
      geom_line() +
      theme_bw() +
      labs(x = "Threshold",
           y = "Prevalence",
           title = paste("Prevalence vs Threshold (", metric_colname, ")", sep = "")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 100)
    methods::show(fig)
  }
  
  return(result)
}

calculate_prevalence3D <- function(grid_metrics,
                                   metric_colname,
                                   threshold,
                                   above = TRUE) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.numeric(threshold) && length(threshold) == 1 && threshold >= 0 && threshold <= 1)) {
    stop("`threshold` is not a numeric between 0 and 1.")
  }
  if (!is.logical(above)) {
    stop("`above` is not a logical (TRUE or FALSE).")
  }
  
  
  ## Exclude rows with NA values
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  if (above) {
    prevalence <- sum(grid_metrics[[metric_colname]] >= threshold) / nrow(grid_metrics) * 100
  }
  else {
    prevalence <- sum(grid_metrics[[metric_colname]] < threshold) / nrow(grid_metrics) * 100
  }
  
  return(prevalence)
}

calculate_spatial_autocorrelation3D <- function(grid_metrics,
                                                metric_colname,
                                                weight_method) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.character(weight_method) && weight_method %in% c("IDW", "rook", "queen"))) {
    stop("`weight_method` is not 'IDW', 'rook' or 'queen'.")
  }
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_metrics)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/3)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  z <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)))
  grid_prism_coords <- data.frame(x = x, y = y, z = z)
  
  ## Subset for non NA rows
  grid_prism_coords <- grid_prism_coords[!is.na(grid_metrics[[metric_colname]]), ]
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use rook method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within 1 unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "rook") {
    weight_matrix <- ifelse(weight_matrix > 1, 0, 1)
  }
  ## Use queen method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(3) unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "queen") {
    weight_matrix <- ifelse(weight_matrix > sqrt(3), 0, 1)
  }
  
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  n <- nrow(grid_metrics)
  
  # Center the data
  data_centered <- grid_metrics[[metric_colname]] - mean(grid_metrics[[metric_colname]])
  
  # Calculate numerator using matrix multiplication
  numerator <- sum(data_centered * (weight_matrix %*% data_centered))
  
  # Calculate denominator
  denominator <- sum(data_centered^2) * sum(weight_matrix)
  
  # Moran's I
  I <- (n * numerator) / denominator
  
  return(I)
}

get_spe_grid_metrics3D <- function(spe,
                                   n_splits,
                                   feature_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  spe_coords <- SpatialExperiment::spatialCoords(spe)
  
  ## Get dimensions of the window
  min_x <- min(spe_coords[ , "Cell.X.Position"])
  min_y <- min(spe_coords[ , "Cell.Y.Position"])
  min_z <- min(spe_coords[ , "Cell.Z.Position"])
  
  max_x <- max(spe_coords[ , "Cell.X.Position"])
  max_y <- max(spe_coords[ , "Cell.Y.Position"])
  max_z <- max(spe_coords[ , "Cell.Z.Position"])
  
  length <- (max_x - min_x)
  width  <- (max_y - min_y)
  height <- (max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  # Shift spe_coords so they begin at the origin
  spe_coords[, "Cell.X.Position"] <- spe_coords[, "Cell.X.Position"] - min_x
  spe_coords[, "Cell.Y.Position"] <- spe_coords[, "Cell.Y.Position"] - min_y
  spe_coords[, "Cell.Z.Position"] <- spe_coords[, "Cell.Z.Position"] - min_z
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$grid_prism_num <- floor(spe_coords[ , "Cell.X.Position"] / d_row) +
    floor(spe_coords[ , "Cell.Y.Position"] / d_col) * n_splits +
    floor(spe_coords[ , "Cell.Z.Position"] / d_lay) * n_splits^2 + 1
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spe[[feature_colname]], factor(spe$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + (min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + (min_y),
                                       z_coord = (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay + (min_z))
  
  spe@metadata[["grid_metrics"]] <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                                         "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(spe)
}

summarise_distances_between_cell_types3D <- function(distances_df) {
  
  # Just in case?
  pair <- distance <- NULL
  
  # Summarise the results
  distances_df_summarised <- distances_df %>%
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance),
                     min(distance),
                     max(distance),
                     stats::median(distance),
                     stats::sd(distance))
  
  distances_df_summarised <- data.frame(distances_df_summarised)
  
  colnames(distances_df_summarised) <- c("pair",
                                         "mean",
                                         "min",
                                         "max",
                                         "median",
                                         "std_dev")
  
  # Add columns for reference cell type and target cell type.
  for (i in seq(nrow(distances_df_summarised))) {
    # Get cell_types for each pair
    cell_types <- strsplit(distances_df_summarised[i,"pair"], "/")[[1]]
    
    distances_df_summarised[i, "reference"] <- cell_types[1]
    distances_df_summarised[i, "target"] <- cell_types[2]
  }
  
  return(distances_df_summarised)
}


# SPIAT2D ----
calculate_all_gradient_cc_metrics2D <- function(spe,
                                                reference_cell_type,
                                                target_cell_types,
                                                radii,
                                                feature_colname,
                                                plot_image = T) {
  
  # Define constants
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  mixing_score_df_colnames <- c("ref_cell_type",
                                "tar_cell_type",
                                "n_ref_cells",
                                "n_tar_cells",
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions",
                                "mixing_score",
                                "normalised_mixing_score")
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "neighbourhood_counts" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "neighbourhood_entropy" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["neighbourhood_counts"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["neighbourhood_entropy"]]) <- target_cell_types
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define individual data frames for mixing_score and cross_G
  for (target_cell_type in target_cell_types) {
    if (reference_cell_type != target_cell_type) {
      result[["mixing_score"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(mixing_score_df_colnames)))
      colnames(result[["mixing_score"]][[target_cell_type]]) <- mixing_score_df_colnames
    }
    result[["cross_G"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(cross_G_df_colnames)))
    colnames(result[["cross_G"]][[target_cell_type]]) <- cross_G_df_colnames
  }
  
  # Get gradient results for each metric
  for (i in seq(length(radii))) {
    df <- calculate_all_single_radius_cc_metrics2D(spe,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["neighbourhood_counts"]]$ref_cell_id <- NULL
    
    result[["neighbourhood_counts"]][i, ] <- apply(df[["neighbourhood_counts"]], 2, mean)
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["neighbourhood_entropy"]][i, ] <- apply(df[["neighbourhood_entropy"]][ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
    result[["cross_K"]][i, ] <- df[["cross_K"]]
    result[["cross_L"]][i, ] <- df[["cross_L"]]
    result[["co_occurrence"]][i, ] <- df[["co_occurrence"]]
    
    for (target_cell_type in names(df[["mixing_score"]])) {
      result[["mixing_score"]][[target_cell_type]][i, ] <- df[["mixing_score"]][[target_cell_type]]
    }
    for (target_cell_type in names(df[["cross_G"]])) {
      result[["cross_G"]][[target_cell_type]][i, ] <- df[["cross_G"]][[target_cell_type]]
    }
  }
  
  # Add radius column to each data frame
  result[["neighbourhood_counts"]]$radius <- radii
  result[["cells_in_neighbourhood"]]$radius <- radii
  result[["neighbourhood_entropy"]]$radius <- radii
  result[["cross_K"]]$radius <- radii
  result[["cross_L"]]$radius <- radii
  result[["co_occurrence"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  for (target_cell_type in names(df[["cross_G"]])) {
    result[["cross_G"]][[target_cell_type]]$radius <- radii
  }
  
  
  ## Plot
  if (plot_image) {
    fig_ANC <- plot_neighbourhood_counts_gradient2D(result[["neighbourhood_counts"]], reference_cell_type)
    methods::show(fig_ANC)
    
    fig_ACIN <- plot_cells_in_neighbourhood_gradient2D(result[["cells_in_neighbourhood"]], reference_cell_type)
    methods::show(fig_ACIN)
    
    fig_ANE <- plot_neighbourhood_entropy_gradient2D(result[["neighbourhood_entropy"]], reference_cell_type)
    methods::show(fig_ANE)
    
    for (target_cell_type in names(result[["mixing_score"]])) {
      fig_NMS <- plot_mixing_scores_gradient2D(result[["mixing_score"]][[target_cell_type]], "NMS")
      fig_MS <- plot_mixing_scores_gradient2D(result[["mixing_score"]][[target_cell_type]], "MS")
      fig_NMS_MS <- cowplot::plot_grid(fig_NMS, fig_MS, nrow = 2)
      methods::show(fig_NMS_MS)
    }
    fig_CK <- plot_cross_K_gradient2D(result[["cross_K"]])
    methods::show(fig_CK)
    
    fig_CL <- plot_cross_L_gradient2D(result[["cross_L"]])
    methods::show(fig_CL)
    
    for (target_cell_type in names(result[["cross_G"]])) {
      fig_CG <- plot_cross_G_gradient2D(result[["cross_G"]][[target_cell_type]], reference_cell_type, target_cell_type)
      methods::show(fig_CG)
    }
    
    fig_COO <- plot_co_occurrence_gradient2D(result[["co_occurrence"]])
    methods::show(fig_COO)
  }
  
  return(result)
}

calculate_all_single_radius_cc_metrics2D <- function(spe,
                                                     reference_cell_type,
                                                     target_cell_types,
                                                     radius,
                                                     feature_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop(paste(radius, " is not a positive numeric."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Define result
  result <- list("neighbourhood_counts" = list(),
                 "cells_in_neighbourhood" = list(),
                 "neighbourhood_entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list(),
                 "cross_L" = list(),
                 "cross_G" = list(),
                 "co_occurrence" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type",
                                "tar_cell_type",
                                "n_ref_cells",
                                "n_tar_cells",
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions",
                                "mixing_score",
                                "normalised_mixing_score")
  
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  # Get rough dimensions of window for cross_K
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  length <- (max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- (max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  
  # Get area of the window the cells are in
  area <- length * width
  
  # All single radius cc metrics stem from calculate_neighbourhood_entropy2D function
  neighbourhood_entropy_df <- calculate_neighbourhood_entropy2D(spe,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radius,
                                                                feature_colname)
  
  ## Neighbourhood counts ----------
  result[["neighbourhood_counts"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_prop", sep = ""))]
  
  ## Neighbourhood_entropy --------------
  result[["neighbourhood_entropy"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_entropy", sep = ""))]
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(neighbourhood_entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(neighbourhood_entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
  }
  
  ## Cross_K ---------------------
  cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
  colnames(cross_K_df) <- cross_K_df_colnames
  cross_K_df$reference <- reference_cell_type
  cross_K_df$expected <- pi * radius^2
  
  for (target_cell_type in target_cell_types) {
    cross_K_df[[target_cell_type]] <- (((area * sum(neighbourhood_entropy_df[[target_cell_type]])) / sum(spe[[feature_colname]] == reference_cell_type)) / sum(spe[[feature_colname]] == target_cell_type))
  }
  result[["cross_K"]] <- cross_K_df
  
  ## Cross_L ---------------------
  cross_L_df <- cross_K_df
  cross_L_df[ , c("expected", target_cell_types)] <- (cross_L_df[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  result[["cross_L"]] <- cross_L_df
  
  ## Cross_G ---------------------
  for (target_cell_type in target_cell_types) {
    cross_G_df <- data.frame(matrix(nrow = 1, ncol = length(cross_G_df_colnames)))
    colnames(cross_G_df) <- cross_G_df_colnames
    
    reference_target_interactions <- neighbourhood_entropy_df[[target_cell_type]]
    n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / area
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  }
  
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spe[[feature_colname]])
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  n_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe),
                                                  eps = radius,
                                                  query = SpatialExperiment::spatialCoords(spe),
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    
    # Account for case where target cell type is not present in data
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      co_occurrence_df[[target_cell_type]] <- NA
      next
    }
    
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe)[spe[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = SpatialExperiment::spatialCoords(spe),
                                                           sort = FALSE)
    n_target_cells_in_all_cell_type_radius <- sum(rapply(n_target_cells_in_all_cell_type_radius$id, base::length))
    
    target_cell_type_proportion_in_all_cell_type_radius <- n_target_cells_in_all_cell_type_radius / n_cells_in_all_cell_type_radius
    
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_all_cell_type_radius
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}

calculate_cell_proportion_grid_metrics2D <- function(spe,
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname,
                                                     plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics2D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^2
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  # Fill in the result data frame
  if (length(reference_cell_types) == 1) {
    result$reference <- grid_prism_cell_matrix[[reference_cell_types]]
  }
  else {
    result$reference <- rowSums(grid_prism_cell_matrix[ , reference_cell_types])
  }
  if (length(target_cell_types) == 1) {
    result$target <- grid_prism_cell_matrix[[target_cell_types]]
  }
  else {
    result$target <- rowSums(grid_prism_cell_matrix[ , target_cell_types])
  }
  result$total <- result$reference + result$target
  result$proportion <- result$target / result$total
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}

calculate_cell_proportions2D <- function(spe,
                                         cell_types_of_interest = NULL,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (ncol(spe) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spe[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Subset for cell types chosen by user
    cell_proportions <- cell_proportions[(cell_proportions$cell_type %in% cell_types_of_interest), ]
    
  }
  
  # Get frequency total for all cells
  cell_type_frequency_total <- sum(cell_proportions$frequency)
  
  # Get proportions and percentages
  cell_proportions$proportion <- cell_proportions$frequency / cell_type_frequency_total
  cell_proportions$percentage <- cell_proportions$proportion * 100
  
  # Order the cell types by proportion (highest cell proportion is first)
  cell_proportions <- cell_proportions[rev(order(cell_proportions$proportion)), ]
  rownames(cell_proportions) <- seq(nrow(cell_proportions))
  
  # Plot
  if (plot_image) {
    
    labels <- paste(round(cell_proportions$percentage, 1), "%", sep = "")
    
    fig <- ggplot(cell_proportions, aes(x = factor(cell_type, cell_type), y = percentage, fill = cell_type)) +
      geom_bar(stat='identity') +
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(fig)
  }
  
  return(cell_proportions)
}

calculate_cells_in_neighbourhood_gradient2D <- function(spe,
                                                        reference_cell_type,
                                                        target_cell_types,
                                                        radii,
                                                        feature_colname,
                                                        plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood2D(spe,
                                                                                        reference_cell_type,
                                                                                        target_cell_types,
                                                                                        radii[i],
                                                                                        feature_colname)
    
    if (is.null(cell_proportions_neighbourhood_proportions_df)) return(NULL)
    
    result[i, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_gradient2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cells_in_neighbourhood2D <- function(spe,
                                               reference_cell_type,
                                               target_cell_types,
                                               radius,
                                               feature_colname) {
  
  ## Get cells in neighbourhood df
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                              reference_cell_type,
                                                              c(reference_cell_type, target_cell_types),
                                                              radius,
                                                              feature_colname)
  
  if (is.null(neighbourhood_counts_df)) return(NULL)
  
  neighbourhood_counts_df[ , paste(target_cell_types, "_prop", sep = "")] <-
    neighbourhood_counts_df[ , target_cell_types] / (neighbourhood_counts_df[ , target_cell_types] + neighbourhood_counts_df[ , reference_cell_type])
  
  # If reference cell type is in target cell types, proportion should be 1
  if (reference_cell_type %in% target_cell_types) {
    neighbourhood_counts_df[neighbourhood_counts_df[[reference_cell_type]] != 0, paste(reference_cell_type, "_prop", sep = "")] <- 1
  }
  
  return(neighbourhood_counts_df)
}

calculate_co_occurrence_gradient2D <- function(spe,
                                               reference_cell_type,
                                               target_cell_types,
                                               radii,
                                               feature_colname,
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types) + 1))
  colnames(result) <- c("reference", target_cell_types)
  
  for (i in seq(length(radii))) {
    co_occurrence_df <- calculate_co_occurrence2D(spe,
                                                  reference_cell_type,
                                                  target_cell_types,
                                                  radii[i],
                                                  feature_colname)
    
    result[i, ] <- co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_co_occurrence_gradient2D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_co_occurrence2D <- function(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radius,
                                      feature_colname) {
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get all cell types in spe
  all_cell_types <- unique(spe[[feature_colname]])
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  # Get total number of cells in radius around all cell types
  n_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe),
                                                  eps = radius,
                                                  query = SpatialExperiment::spatialCoords(spe),
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    
    # Account for case where target cell type is not present in data
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      co_occurrence_df[[target_cell_type]] <- NA
      next
    }
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get total number of target cells in radius around all cell types
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe)[spe[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = SpatialExperiment::spatialCoords(spe),
                                                           sort = FALSE)
    n_target_cells_in_all_cell_type_radius <- sum(rapply(n_target_cells_in_all_cell_type_radius$id, base::length))
    
    # Get proportion of target cells in radius around all cell types
    target_cell_type_proportion_in_all_cell_type_radius <- n_target_cells_in_all_cell_type_radius / n_cells_in_all_cell_type_radius
    
    # Get co-occurence value for target cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_all_cell_type_radius
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}

calculate_cross_G_gradient2D <- function(spe,
                                         reference_cell_type,
                                         target_cell_type,
                                         radii,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check if radii input is valid
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 2))
  colnames(result) <- c("observed_cross_G",
                        "expected_cross_G")
  
  # Iterate through each radii
  for (i in seq(length(radii))) {
    cross_G_df <- calculate_cross_G2D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_G_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cross_G_gradient2D(result, reference_cell_type, target_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cross_G2D <- function(spe,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname) {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                              reference_cell_type,
                                                              target_cell_type,
                                                              radius,
                                                              feature_colname)
  
  reference_target_interactions <- neighbourhood_counts_df[[target_cell_type]]
  
  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  
  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  length <- (max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- (max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  
  # Get area of the window the cells are in
  area <- length * width
  
  # Get the number of target cells
  n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / area
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}

calculate_cross_K_gradient2D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check if radii input is valid
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  # Iterate through each radii
  for (i in seq(length(radii))) {
    cross_K_df <- calculate_cross_K2D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cross_K_gradient2D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cross_K2D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Check if inputs are valid
  if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  
  ## Get expected cross K-function
  expected_cross_K <- pi * radius^2
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  length <- (max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- (max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  
  ## Get area of the window the cells are in
  area <- length * width
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                              reference_cell_type,
                                                              target_cell_types,
                                                              radius,
                                                              feature_colname)
  
  # Calculate cross K-fucnction for each target cell type
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
    
    ## Get observed cross K-function
    if (n_tar_cells == 0) {
      observed_cross_K <- NA
    }
    else {
      observed_cross_K <- (area * n_ref_tar_interactions) / (n_ref_cells * n_tar_cells)
    }
    result[[target_cell_type]] <- observed_cross_K
  }
  
  return(result)
}

calculate_cross_L_gradient2D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_L_df <- calculate_cross_L2D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_L_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cross_L_gradient2D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cross_L2D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Calculate cross K first
  result <- calculate_cross_K2D(spe = spe,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  # Linearise cross K to get cross L
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  
  return(result)
}

calculate_entropy_background2D <- function(spe,
                                           cell_types_of_interest,
                                           feature_colname) {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions2D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy)
}

calculate_entropy_grid_metrics2D <- function(spe,
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname,
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spe[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics2D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^2
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  for (cell_type in cell_types_of_interest) {
    result[[cell_type]] <- grid_prism_cell_matrix[[cell_type]]
  }
  result$total <- rowSums(result)
  
  ## Get data frame containing proportions for cell_types_of_interest
  df_props <- result[ , cell_types_of_interest] / result$total
  
  ## Use proportion data frame to get entropy
  calculate_entropy <- function(x) {
    entropy <- -1 * sum(x * ifelse(is.infinite(log(x, 2)), 0, log(x, 2))) / log(length(x), 2)
    return(entropy)
  }
  result$entropy <- apply(df_props, 1, calculate_entropy)
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}

calculate_minimum_distances_between_cell_types2D <- function(spe,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname,
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Get spe coords
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spe[[feature_colname]] == cell_type1) == 0 || sum(spe[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spe_coords[spe[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spe_coords[spe[[feature_colname]] == cell_type2, ]
    
    # Find all of closest points
    # For each cell of cell_type1, find the closest cell of cell_type2
    if (cell_type1 != cell_type2) {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords,
                                      query = cell_type1_coords,
                                      k = 1)
    }
    # If we are comparing the same cell_type, and there is only one of this cell type, move on
    else if (nrow(cell_type1_coords) == 1) {
      warning("There is only 1 '", cell_type1, "' cell in your data. It has no nearest neighbour of the same cell type.", sep = "")
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    # If we are comparing the same cell_type, use the second closest neighbour
    else {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords,
                                      query = cell_type1_coords,
                                      k = 2)
      
      # Get cell IDs and distances
      nearest_neighbours[['nn.idx']] <- nearest_neighbours[['nn.idx']][ , 2]
      nearest_neighbours[['nn.dists']] <- nearest_neighbours[['nn.dists']][ , 2]
    }
    
    # Create the data frame containing the chosen cells and their ids, as well as the nearest cell to them and their ids, and the distance between
    df <- data.frame(
      ref_cell_id = cell_type_ids[[cell_type1]],
      ref_cell_type = cell_type1,
      nearest_cell_id = cell_type_ids[[cell_type2]][c(nearest_neighbours$nn.idx)],
      nearest_cell_type = cell_type2,
      distance = nearest_neighbours$nn.dists
    )
    result <- rbind(result, df)
  }
  
  # Add cell pair column
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(result))
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin2D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_mixing_scores_gradient2D <- function(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii,
                                               feature_colname,
                                               plot_image = TRUE) {
  
  # Check radii input
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 8))
  colnames(result) <- c("ref_cell_type",
                        "tar_cell_type",
                        "n_ref_cells",
                        "n_tar_cells",
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions",
                        "mixing_score",
                        "normalised_mixing_score")
  
  # Iterate through each radius
  for (i in seq(length(radii))) {
    mixing_scores <- calculate_mixing_scores2D(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii[i],
                                               feature_colname)
    
    result[i, ] <- mixing_scores
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_mixing_scores_gradient2D(result, "NMS")
    fig2 <- plot_mixing_scores_gradient2D(result, "MS")
    combined_fig <- cowplot::plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}

calculate_mixing_scores2D <- function(spe,
                                      reference_cell_types,
                                      target_cell_types,
                                      radius,
                                      feature_colname) {
  
  # Define result
  result <- data.frame()
  
  # Iterate through each reference cell type
  for (reference_cell_type in reference_cell_types) {
    
    # Iterate through each target cell type
    for (target_cell_type in target_cell_types) {
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Get number of reference cells and target cells
      n_ref <- sum(spe[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spe[[feature_colname]] == target_cell_type)
      
      
      # Can't get mixing scores if there are 0 or 1 reference cells
      if (n_ref == 0 || n_ref == 1) {
        result <-  rbind(result,
                         c(reference_cell_type,
                           target_cell_type,
                           n_ref,
                           n_tar,
                           0,
                           0,
                           NA,
                           NA))
      }
      
      
      ## Get cells in neighbourhood df
      neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                                  reference_cell_type,
                                                                  c(reference_cell_type, target_cell_type),
                                                                  radius,
                                                                  feature_colname)
      
      # Get number of ref-ref interactions
      # Halve it to avoid counting each ref-ref interaction twice
      n_ref_ref_interactions <- 0.5 * sum(neighbourhood_counts_df[[reference_cell_type]])
      
      # Get number of ref-tar interactions
      n_ref_tar_interactions <- sum(neighbourhood_counts_df[[target_cell_type]])
      
      
      # Can't get mixing scores if there are no target cells
      if (n_tar == 0) {
        
        result <-  rbind(result,
                         c(reference_cell_type,
                           target_cell_type,
                           n_ref,
                           0,
                           0,
                           n_ref_ref_interactions,
                           NA,
                           NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        
        if (n_ref_ref_interactions != 0) {
          mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
          normalised_mixing_score <- 0.5 * mixing_score * n_ref / n_tar
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        result <-  rbind(result,
                         c(reference_cell_type,
                           target_cell_type,
                           n_ref,
                           n_tar,
                           n_ref_tar_interactions,
                           n_ref_ref_interactions,
                           mixing_score,
                           normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  colnames(result) <- c("ref_cell_type",
                        "tar_cell_type",
                        "n_ref_cells",
                        "n_tar_cells",
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions",
                        "mixing_score",
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}

calculate_neighbourhood_counts_gradient2D <- function(spe,
                                                      reference_cell_type,
                                                      target_cell_types,
                                                      radii,
                                                      feature_colname,
                                                      plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radii[i],
                                                                feature_colname)
    
    if (is.null(neighbourhood_counts_df)) return(NULL)
    
    neighbourhood_counts_df$ref_cell_id <- NULL
    result[i, ] <- apply(neighbourhood_counts_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_neighbourhood_counts_gradient2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_neighbourhood_counts2D <- function(spe,
                                             reference_cell_type,
                                             target_cell_types,
                                             radius,
                                             feature_colname) {
  
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # Get spe coords
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spe_coords[spe[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(ref_cell_id = spe$Cell.ID[spe[[feature_colname]] == reference_cell_type])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spe_coords[spe[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells specified distance for each reference cell
    ref_tar_result <- dbscan::frNN(target_cell_type_coords,
                                   eps = radius,
                                   query = reference_cell_type_coords,
                                   sort = FALSE)
    
    n_targets <- rapply(ref_tar_result$id, length)
    
    
    # Don't want to include the reference cell as one of the target cells
    if (reference_cell_type == target_cell_type) n_targets <- n_targets - 1
    
    ## Add to data frame
    result[[target_cell_type]] <- n_targets
  }
  
  return(result)
}

calculate_neighbourhood_entropy_gradient2D <- function(spe,
                                                       reference_cell_type,
                                                       target_cell_types,
                                                       radii,
                                                       feature_colname,
                                                       plot_image = TRUE) {
  
  # Check radii input
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  # Iterate through each radii
  for (i in seq(length(radii))) {
    neighbourhood_entropy_df <- calculate_neighbourhood_entropy2D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radii[i],
                                                                  feature_colname)
    
    if (is.null(neighbourhood_entropy_df)) return(NULL)
    
    result[i, ] <- apply(neighbourhood_entropy_df[ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_neighbourhood_entropy_gradient2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_neighbourhood_entropy2D <- function(spe,
                                              reference_cell_type,
                                              target_cell_types,
                                              radius,
                                              feature_colname) {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  ## Get neighbourhood_entropy for target_cell_type
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_entropy", sep = "")] <-
    -1 *
    (cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] * log(cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")], 2) +
       (1 - cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")]) * log(1 - (cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")]), 2))
  
  return(cells_in_neighbourhood_df)
}

calculate_pairwise_distances_between_cell_types2D <- function(spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname,
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(SpatialExperiment::spatialCoords(spe))
  rownames(distance_matrix) <- spe$Cell.ID
  colnames(distance_matrix) <- spe$Cell.ID
  
  result <- data.frame()
  
  for (i in seq(length(cell_types_of_interest))) {
    
    for (j in i:length(cell_types_of_interest)) {
      
      # Get current cell types and cell ids
      cell_type1 <- names(cell_type_ids)[i]
      cell_type2 <- names(cell_type_ids)[j]
      
      cell_type1_ids <- cell_type_ids[[cell_type1]]
      cell_type2_ids <- cell_type_ids[[cell_type2]]
      
      ## Don't have a cell type, or the same cell type with only one cell
      if (length(cell_type1_ids) == 0 || length(cell_type2_ids) == 0) {
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      ## Same cell type only one cell
      if (cell_type1 == cell_type2 && length(cell_type1_ids) == 1) {
        warning("There is only 1 '", cell_type1, "' cell in your data. It has no pair of the same cell type.", sep = "")
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      # Subset distance_matrix for current cell types
      distance_matrix_subset <- distance_matrix[rownames(distance_matrix) %in% cell_type1_ids,
                                                colnames(distance_matrix) %in% cell_type2_ids]
      
      ## Different cell types, each only has one cell
      if (length(cell_type1_ids) == 1 && length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        rownames(distance_matrix_subset) <- cell_type1_ids
        colnames(distance_matrix_subset) <- cell_type2_ids
      }
      ## Different cell types, only one cell of cell_type1
      else if (length(cell_type1_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type1_ids
      }
      ## Different cell types, only one cell of cell_type2
      else if (length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type2_ids
      }
      ## Same cell type, only need part of the matrix (make irrelevant part of matrix equal to NA)
      if (cell_type1 == cell_type2) distance_matrix_subset[upper.tri(distance_matrix_subset, diag = TRUE)] <- NA
      
      # Convert distance_matrix_subset to a data frame
      df <- reshape2::melt(distance_matrix_subset, na.rm = TRUE)
      df$cell_type1 <- cell_type1
      df$cell_type2 <- cell_type2
      df$pair <- paste(cell_type1, cell_type2, sep="/")
      
      result <- rbind(result, df)
    }
  }
  
  # Rearrange columns
  colnames(result)[c(1, 2, 3)] <- c("cell_type1_id", "cell_type2_id", "distance")
  result <- result[ , c("cell_type1_id", "cell_type1", "cell_type2_id", "cell_type2", "distance", "pair")]
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(result))
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin2D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_prevalence_gradient2D <- function(grid_metrics,
                                            metric_colname,
                                            show_AUC = TRUE,
                                            plot_image = TRUE) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!is.logical(show_AUC)) {
    stop("`show_AUC` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  # Get prevalences for each threshold
  result$prevalence <- sapply(thresholds, function(threshold) {
    calculate_prevalence2D(grid_metrics, metric_colname, threshold)
  })
  
  # Show AUC of prevalence gradient graph
  if (show_AUC) {
    AUC_value <- sum(diff(result$threshold) * (head(result$prevalence, -1) + tail(result$prevalence, -1)) / 2)
    print(paste("AUC:", round(AUC_value, 2)))
  }
  
  # Plot
  if (plot_image) {
    fig <- ggplot(result, aes(threshold, prevalence)) +
      geom_line() +
      theme_bw() +
      labs(x = "Threshold",
           y = "Prevalence",
           title = paste("Prevalence vs Threshold (", metric_colname, ")", sep = "")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 100)
    methods::show(fig)
  }
  
  return(result)
}

calculate_prevalence2D <- function(grid_metrics,
                                   metric_colname,
                                   threshold,
                                   above = TRUE) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.numeric(threshold) && length(threshold) == 1 && threshold >= 0 && threshold <= 1)) {
    stop("`threshold` is not a numeric between 0 and 1.")
  }
  if (!is.logical(above)) {
    stop("`above` is not a logical (TRUE or FALSE).")
  }
  
  
  ## Exclude rows with NA values
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  if (above) {
    prevalence <- sum(grid_metrics[[metric_colname]] >= threshold) / nrow(grid_metrics) * 100
  }
  else {
    prevalence <- sum(grid_metrics[[metric_colname]] < threshold) / nrow(grid_metrics) * 100
  }
  
  return(prevalence)
}

calculate_spatial_autocorrelation2D <- function(grid_metrics,
                                                metric_colname,
                                                weight_method) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.character(weight_method) && weight_method %in% c("IDW", "rook", "queen"))) {
    stop("`weight_method` is not 'IDW', 'rook' or 'queen'.")
  }
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_metrics)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/2)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  grid_prism_coords <- data.frame(x = x, y = y)
  
  ## Subset for non NA rows
  grid_prism_coords <- grid_prism_coords[!is.na(grid_metrics[[metric_colname]]), ]
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use rook method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within 1 unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "rook") {
    weight_matrix <- ifelse(weight_matrix > 1, 0, 1)
  }
  ## Use queen method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(3) unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "queen") {
    weight_matrix <- ifelse(weight_matrix > sqrt(2), 0, 1)
  }
  
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  n <- nrow(grid_metrics)
  
  # Center the data
  data_centered <- grid_metrics[[metric_colname]] - mean(grid_metrics[[metric_colname]])
  
  # Calculate numerator using matrix multiplication
  numerator <- sum(data_centered * (weight_matrix %*% data_centered))
  
  # Calculate denominator
  denominator <- sum(data_centered^2) * sum(weight_matrix)
  
  # Moran's I
  I <- (n * numerator) / denominator
  
  return(I)
}

get_spe_grid_metrics2D <- function(spe,
                                   n_splits,
                                   feature_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  spe_coords <- SpatialExperiment::spatialCoords(spe)
  
  ## Get dimensions of the window
  min_x <- min(spe_coords[ , "Cell.X.Position"])
  min_y <- min(spe_coords[ , "Cell.Y.Position"])
  
  max_x <- max(spe_coords[ , "Cell.X.Position"])
  max_y <- max(spe_coords[ , "Cell.Y.Position"])
  
  length <- (max_x - min_x)
  width  <- (max_y - min_y)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  # Shift spe_coords so they begin at the origin
  spe_coords[, "Cell.X.Position"] <- spe_coords[, "Cell.X.Position"] - min_x
  spe_coords[, "Cell.Y.Position"] <- spe_coords[, "Cell.Y.Position"] - min_y
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$grid_prism_num <- floor(spe_coords[ , "Cell.X.Position"] / d_row) +
    floor(spe_coords[ , "Cell.Y.Position"] / d_col) * n_splits
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^2
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spe[[feature_colname]], factor(spe$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + (min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + (min_y))
  
  spe@metadata[["grid_metrics"]] <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                                         "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(spe)
}

summarise_distances_between_cell_types2D <- function(distances_df) {
  
  # Just in case?
  pair <- distance <- NULL
  
  # Summarise the results
  distances_df_summarised <- distances_df %>%
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance),
                     min(distance),
                     max(distance),
                     stats::median(distance),
                     stats::sd(distance))
  
  distances_df_summarised <- data.frame(distances_df_summarised)
  
  colnames(distances_df_summarised) <- c("pair",
                                         "mean",
                                         "min",
                                         "max",
                                         "median",
                                         "std_dev")
  
  # Add columns for reference cell type and target cell type.
  for (i in seq(nrow(distances_df_summarised))) {
    # Get cell_types for each pair
    cell_types <- strsplit(distances_df_summarised[i,"pair"], "/")[[1]]
    
    distances_df_summarised[i, "reference"] <- cell_types[1]
    distances_df_summarised[i, "target"] <- cell_types[2]
  }
  
  return(distances_df_summarised)
}


# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)
library(scales)
library(DescTools)
library(spaSim3D)

# Functions -----
# Generate data frame containing values of important parameters
generate_random_parameters <- function(
    n_simulations_per_parameter,
    bg_prop_A_range,
    bg_prop_B_range,
    E_radius_x_range,
    E_radius_y_range,
    E_radius_z_range,
    N_width_range,
    cluster_prop_A_range,
    ring_width_factor_range,
    cluster1_x_coord_range
) {
  
  parameters <- c("bg_prop_A", 
                  "bg_prop_B", 
                  "E_radius_x", # E_radius_x, E_radius_y, E_radius_z count as 1
                  "N_width", 
                  "cluster_prop_A", 
                  "ring_width_factor", 
                  "cluster1_x_coord")
  n_parameters <- length(parameters) 
  n_simulations <- n_simulations_per_parameter * n_parameters
  
  # Set parameter values to be the average of the min and max set in the parameter range, bg_props are an exception and will be equal to 0
  parameters_df <- data.frame(
    bg_prop_A = 0,
    bg_prop_B = 0,
    E_radius_x = mean(c(E_radius_x_range["min"], E_radius_x_range["max"])),
    E_radius_y = mean(c(E_radius_y_range["min"], E_radius_y_range["max"])),
    E_radius_z = mean(c(E_radius_z_range["min"], E_radius_z_range["max"])),
    N_width = mean(c(N_width_range["min"], N_width_range["max"])),
    cluster_prop_A = mean(c(cluster_prop_A_range["min"], cluster_prop_A_range["max"])),
    ring_width_factor = mean(c(ring_width_factor_range["min"], ring_width_factor_range["max"])),
    cluster1_x_coord = mean(c(cluster1_x_coord_range["min"], cluster1_x_coord_range["max"])),
    variable_parameter = rep(parameters, each = n_simulations_per_parameter)
  )
  
  # Put all parameter ranges into a list
  parameter_ranges <- list(
    bg_prop_A = bg_prop_A_range,
    bg_prop_B = bg_prop_B_range,
    E_radius_x = E_radius_x_range,
    E_radius_y = E_radius_y_range,
    E_radius_z = E_radius_z_range,
    N_width = N_width_range,
    cluster_prop_A = cluster_prop_A_range,
    ring_width_factor = ring_width_factor_range,
    cluster1_x_coord = cluster1_x_coord_range
  )
  
  # For the current variable parameter, vary the parameter...
  for (i in seq_len(nrow(parameters_df))) {
    
    curr_variable_parameter <- parameters_df[["variable_parameter"]][i]
    
    if (curr_variable_parameter == "E_radius_x") {
      parameters_df[i, "E_radius_x"] <- runif(1, E_radius_x_range[["min"]], E_radius_x_range[["max"]])
      parameters_df[i, "E_radius_y"] <- runif(1, E_radius_y_range[["min"]], E_radius_y_range[["max"]])
      parameters_df[i, "E_radius_z"] <- runif(1, E_radius_z_range[["min"]], E_radius_z_range[["max"]])
    }
    else {
      parameter_min <- parameter_ranges[[curr_variable_parameter]][["min"]]
      parameter_max  <- parameter_ranges[[curr_variable_parameter]][["max"]]
      parameters_df[i, curr_variable_parameter] <- runif(1, parameter_min, parameter_max)
    }
  }
  
  return(parameters_df)
}

# Generate updated parameters data frame to be specific for arrangement and shape
generate_updated_parameters_for_arrangement_and_shape <- function(
    n_simulations_per_parameter,
    bg_prop_A_range,
    bg_prop_B_range,
    E_radius_x_range,
    E_radius_y_range,
    E_radius_z_range,
    N_width_range,
    cluster_prop_A_range,
    ring_width_factor_range,
    cluster1_x_coord_range
) {
  
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  updated_parameter_df <- data.frame()
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      parameters_df <- generate_random_parameters(
        n_simulations_per_parameter = n_simulations_per_parameter,
        bg_prop_A_range = bg_prop_A_range,
        bg_prop_B_range = bg_prop_B_range,
        E_radius_x_range = E_radius_x_range,
        E_radius_y_range = E_radius_y_range,
        E_radius_z_range = E_radius_z_range,
        N_width_range = N_width_range,
        cluster_prop_A_range = cluster_prop_A_range,
        ring_width_factor_range = ring_width_factor_range,
        cluster1_x_coord_range = cluster1_x_coord_range
      )   
      
      parameters_df$arrangement <- arrangement
      parameters_df$shape <- shape
      
      if (arrangement == "mixed") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("ring_width_factor", "cluster1_x_coord")), ] 
      }
      else if (arrangement == "ringed") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("cluster_prop_A", "cluster1_x_coord")), ] 
      }
      else if (arrangement == "separated") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("cluster_prop_A", "ring_width_factor")), ] 
      }
      
      if (shape == "ellipsoid") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("N_width")), ] 
      }
      else if (shape == "network") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("E_radius_x")), ] 
      }
      updated_parameter_df <- rbind(updated_parameter_df, parameters_df)
    }
  }
  return(updated_parameter_df)
}


# Analyse SD2 simulations
analyse_SD2_simulations <- function(parameters_df) {
  
  # Set defined parameters/values
  cell_types <- c('A', 'B')
  n_simulations <- nrow(parameters_df)
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  n_splits <- 10
  thresholds <- seq(0.01, 1, 0.01)
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  # z-coords for slices
  bottom_z_coord_of_slices <- c(145, 175, 205)
  top_z_coord_of_slices <- bottom_z_coord_of_slices + 10
  n_slices <- length(bottom_z_coord_of_slices)
  
  # Generate simulation metadata using parameters data frame from 'generate_random_parameters' function
  generate_simulation_metadata <- function(parameters_df, simulation_index) {
    
    # Constant metadata for simulations
    simulation_metadata <- spe_metadata_background_template("random")
    simulation_metadata$background$n_cells <- 30000
    simulation_metadata$background$length <- 600
    simulation_metadata$background$width <- 600
    simulation_metadata$background$height <- 300
    simulation_metadata$background$minimum_distance_between_cells <- 10
    simulation_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later
    
    N_radius <- 125
    N_n_edges <- 20
    
    mixed_cluster_cell_types <- c("A", "B")
    mixed_cluster_centre_loc <- c(300, 300, 150)
    
    ringed_cluster_cell_type <- "A"
    ringed_cluster_cell_prop <- 1
    ringed_ring_cell_type <- "B"
    ringed_ring_cell_prop <- 1
    ringed_cluster_centre_loc <- c(300, 300, 150)
    
    separated_cluster1_cell_type <- "A"
    separated_cluster1_cell_prop <- 1
    separated_cluster1_y_z_centre_loc <- c(300, 150) # Append x-coord to it later
    separated_cluster2_shape <- "sphere"
    separated_cluster2_cell_type <- "B"
    separated_cluster2_cell_prop <- 1
    separated_cluster2_centre_loc <- c(450, 300, 150)
    separated_cluster2_radius <- 100
    
    # Alter background cell proportions
    simulation_metadata$background$cell_proportions <- c(parameters_df$bg_prop_A[simulation_index], 
                                                         parameters_df$bg_prop_B[simulation_index],
                                                         1 - parameters_df$bg_prop_A[simulation_index] - parameters_df$bg_prop_B[simulation_index]) # prop(O) = 1 - prop(A) - prop(B)
    
    shape <- parameters_df[simulation_index, "shape"]
    arrangement <- parameters_df[simulation_index, "arrangement"]
    
    # Determine what type of cluster to add to current metadata
    if (arrangement %in% c("mixed", "separated")) {
      simulation_metadata <- spe_metadata_cluster_template("regular", shape, simulation_metadata)
    }
    else if (arrangement == "ringed") {
      simulation_metadata <- spe_metadata_cluster_template("ring", shape, simulation_metadata)
    }
    
    # Modify shape parameters
    if (shape == "ellipsoid") {
      simulation_metadata$cluster_1$radii <- c(parameters_df$E_radius_x[simulation_index], parameters_df$E_radius_y[simulation_index], parameters_df$E_radius_z[simulation_index])
      simulation_metadata$cluster_1$axes_rotation <- c(runif(1, 0, 180), runif(1, 0, 180), runif(1, 0, 180)) # Choose random rotation angles
    }
    
    else if (shape == "network") {
      simulation_metadata$cluster_1$n_edges <- N_n_edges
      simulation_metadata$cluster_1$width <- parameters_df$N_width[simulation_index]
      simulation_metadata$cluster_1$radius <- N_radius
    }
    
    # Modify arrangemetn parameters
    if (arrangement == "mixed") {
      simulation_metadata$cluster_1$cluster_cell_types <- mixed_cluster_cell_types
      simulation_metadata$cluster_1$cluster_cell_proportions <- c(parameters_df$cluster_prop_A[simulation_index], 1 - parameters_df$cluster_prop_A[simulation_index])
      simulation_metadata$cluster_1$centre_loc <- mixed_cluster_centre_loc
    }
    else if (arrangement == "ringed") {
      simulation_metadata$cluster_1$cluster_cell_types <- ringed_cluster_cell_type
      simulation_metadata$cluster_1$cluster_cell_proportions <- ringed_cluster_cell_prop
      simulation_metadata$cluster_1$centre_loc <- ringed_cluster_centre_loc
      simulation_metadata$cluster_1$ring_cell_types <- ringed_ring_cell_type
      simulation_metadata$cluster_1$ring_cell_proportions <- ringed_ring_cell_prop
      simulation_metadata$cluster_1$ring_width <- parameters_df$ring_width_factor[simulation_index] * ifelse(shape == "network", 
                                                                                                             parameters_df$N_width[simulation_index],
                                                                                                             apply(parameters_df[, c("E_radius_x", "E_radius_y", "E_radius_z")], 1, mean))
    }
    else if (arrangement == "separated") {
      simulation_metadata$cluster_1$cluster_cell_types <- separated_cluster1_cell_type
      simulation_metadata$cluster_1$cluster_cell_proportions <- separated_cluster1_cell_prop
      simulation_metadata$cluster_1$centre_loc <- c(parameters_df$cluster1_x_coord[simulation_index], separated_cluster1_y_z_centre_loc)
      
      simulation_metadata <- spe_metadata_cluster_template("regular", "sphere", simulation_metadata)
      simulation_metadata$cluster_2$cluster_cell_types <- separated_cluster2_cell_type
      simulation_metadata$cluster_2$cluster_cell_proportions <- separated_cluster2_cell_prop
      simulation_metadata$cluster_2$centre_loc <- separated_cluster2_centre_loc
      simulation_metadata$cluster_2$radius <- separated_cluster2_radius
    }
    
    return(simulation_metadata)
  }
  
  # Function to get all slices from spe
  get_all_spe_slices <- function(spe, 
                                 bottom_z_coord_of_slices, 
                                 top_z_coord_of_slices) {
    
    spe_slices <- list()
    
    if (length(bottom_z_coord_of_slices) != length(top_z_coord_of_slices)) stop("Lengths of bottom_z_coords_of_slices and top_z_coords_of_slices should be equal.")
    number_of_slices <- length(bottom_z_coord_of_slices)
    
    for (i in seq(number_of_slices)) {
      bottom_z_coord_of_slice <- bottom_z_coord_of_slices[i]
      top_z_coord_of_slice <- top_z_coord_of_slices[i]
      z_coords_of_cells_in_spe <- spatialCoords(spe)[ , "Cell.Z.Position"]
      spe_for_slice <- spe[, bottom_z_coord_of_slice < z_coords_of_cells_in_spe & z_coords_of_cells_in_spe < top_z_coord_of_slice]
      spatialCoords(spe_for_slice) <- spatialCoords(spe_for_slice)[ , c("Cell.X.Position", "Cell.Y.Position")]
      
      spe_slices[[i]] <- spe_for_slice
    }
    
    return(spe_slices)
  }
  
  # Function to create empty metric df list
  create_empty_metric_df_list3D <- function(
    cell_types,
    n_simulations,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    # n_slices + 1, to include the 3D value
    AMD_df_colnames <- c("simulation", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("simulation", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("simulation", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("simulation", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("simulation", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           COO = COO_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  create_empty_metric_df_list2D <- function(
    cell_types,
    n_simulations,
    n_slices,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("simulation", "slice", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "slice", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("simulation", "slice", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("simulation", "slice", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("simulation", "slice", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("simulation", "slice", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           COO = COO_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  
  # Function analyse spe in 3D
  analyse_simulation3D <- function(spe, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   simulation_index, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    index <- n_cell_type_combinations * (simulation_index - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- simulation_index
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for 'pair' data which increments after each reference-target pair
    pair_index <- n_cell_type_combinations * (simulation_index - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics3D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F,
                                                           feature_colname = "Cell.Type")
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
        
        
        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CK"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CL"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CG"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- NA
          metric_df_list[["MS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- NA
        }
        else {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_counts"]][[target_cell_type]]
          metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
          metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
          metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
          
          if (reference_cell_type != target_cell_type) {
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
            metric_df_list[["MS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_entropy"]][[target_cell_type]]
          }
        }
        if (reference_cell_type == target_cell_type) {
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["MS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- Inf
        }
        
        # Spatial heterogeneity metrics
        metric_df_list[["PBSAC"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("simulation", "cell_types")] <- c(simulation_index, paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("simulation", "cell_types")] <- c(simulation_index, paste(reference_cell_type, target_cell_type, sep = ","))
        
        if (reference_cell_type != target_cell_type) {
          proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(spe, 
                                                                              n_splits,
                                                                              reference_cell_type, 
                                                                              target_cell_type,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
          
          if (is.null(proportion_grid_metrics)) {
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- NA
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            PBSAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                         "proportion",
                                                         weight_method = "queen")
            
            PBP_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                      "proportion",
                                                      show_AUC = F,
                                                      plot_image = F)
            
            
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
          } 
          
          entropy_grid_metrics <- calculate_entropy_grid_metrics3D(spe, 
                                                                   n_splits,
                                                                   c(reference_cell_type, target_cell_type), 
                                                                   plot_image = F,
                                                                   feature_colname = "Cell.Type")
          
          if (is.null(entropy_grid_metrics)) {
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            EBSAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                         "entropy",
                                                         weight_method = "queen")
            
            EBP_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                      "entropy",
                                                      show_AUC = F,
                                                      plot_image = F)
            
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- EBSAC
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- EBP_df$prevalence
          }    
        }
        else {
          metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- Inf
          metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- Inf
          metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- Inf
          metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- Inf
        }
        
        pair_index <- pair_index + 1
      }
    }
    return(metric_df_list)
  }
  
  # Function analyse spe in 2D
  analyse_simulation2D <- function(spe_slices, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   simulation_index, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    n_slices <- length(spe_slices)
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    for (slice_index in seq_len(n_slices)) {
      spe <- spe_slices[[slice_index]]
      
      index <- n_cell_type_combinations * (simulation_index - 1) * n_slices +  n_cell_type_combinations * (slice_index - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(spe,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- simulation_index
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- slice_index
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for 'pair' data which increments after each reference-target pair
      pair_index <- n_cell_type_combinations * (simulation_index - 1) * n_slices +  n_cell_type_combinations * (slice_index - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics2D(spe,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             plot_image = F,
                                                             feature_colname = "Cell.Type")
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ANC"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["ACIN"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CK"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CL"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CG"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["ANE"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, 
                                                                                                    paste(reference_cell_type, target_cell_type, sep = ","))
          
          
          if (is.null(gradient_data)) {
            metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
            metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CK"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CL"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CG"]][pair_index, radii_colnames] <- NA
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- NA
            metric_df_list[["MS"]][pair_index, radii_colnames] <- NA
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- NA
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- NA
          }
          else {
            metric_df_list[["ANC"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_counts"]][[target_cell_type]]
            metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
            metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
            metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
            metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
            
            if (reference_cell_type != target_cell_type) {
              metric_df_list[["ACIN"]][pair_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
              metric_df_list[["MS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
              metric_df_list[["NMS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
              metric_df_list[["ANE"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_entropy"]][[target_cell_type]]
            }
          }
          if (reference_cell_type == target_cell_type) {
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["MS"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- Inf
          }
          
          # Spatial heterogeneity metrics
          metric_df_list[["PBSAC"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][pair_index, c("simulation", "slice", "cell_types")] <- c(simulation_index, slice_index, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][pair_index, c("simulation", "slice", "cell_types")] <- c(simulation_index, slice_index, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(spe, 
                                                                                n_splits,
                                                                                reference_cell_type, 
                                                                                target_cell_type,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
            
            if (is.null(proportion_grid_metrics)) {
              metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- NA
              metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- NA
            }
            else {
              PBSAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                           "proportion",
                                                           weight_method = "queen")
              
              PBP_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                        "proportion",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              
              metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
              metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
            } 
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics2D(spe, 
                                                                     n_splits,
                                                                     c(reference_cell_type, target_cell_type), 
                                                                     plot_image = F,
                                                                     feature_colname = "Cell.Type")
            
            if (is.null(entropy_grid_metrics)) {
              metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
              metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
            }
            else {
              EBSAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                           "entropy",
                                                           weight_method = "queen")
              
              EBP_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                        "entropy",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- EBSAC
              metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- EBP_df$prevalence
            }    
          }
          else {
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- Inf
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- Inf
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- Inf
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- Inf
          }
          
          pair_index <- pair_index + 1
        }
      } 
    }
    return(metric_df_list)
  }
  
  # Define metric df lists
  metric_df_list3D <- create_empty_metric_df_list3D(cell_types, n_simulations, radii_colnames, thresholds_colnames)
  metric_df_list2D <- create_empty_metric_df_list2D(cell_types, n_simulations, n_slices, radii_colnames, thresholds_colnames)
  
  for (simulation_index in seq_len(n_simulations)) {
    print(simulation_index)
    # Simulate spes in 3D, and get 2D slices
    simulation_metadata <- generate_simulation_metadata(parameters_df, simulation_index)
    
    spe3D <- simulate_spe_metadata3D(simulation_metadata, plot_image = F)
    
    spe_slices2D <- get_all_spe_slices(spe3D, bottom_z_coord_of_slices, top_z_coord_of_slices)
    
    # Analyse spes
    metric_df_list3D <- analyse_simulation3D(spe3D,
                                             cell_types,
                                             radii,
                                             thresholds,
                                             n_splits,
                                             simulation_index,
                                             metric_df_list3D)
    
    metric_df_list2D <- analyse_simulation2D(spe_slices2D,
                                             cell_types,
                                             radii,
                                             thresholds,
                                             n_splits,
                                             simulation_index,
                                             metric_df_list2D)
  }
  
  ## Turn gradient metrics into AUC and add to metric_df list
  get_AUC_for_radii_gradient_metrics <- function(y) {
    x <- radii
    h <- diff(x)[1]
    n <- length(x)
    
    AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
    
    return(AUC)
  }
  
  add_AUC_for_radii_gradient_metrics_to_metric_df_list3D <- function(metric_df_list) {
    gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "COO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "reference", "target", metric_AUC_name)
      }
      else {
        stop("Unknown metric?")
      }
      
      df <- metric_df_list[[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      df <- df[ , subset_colnames]
      metric_df_list[[metric_AUC_name]] <- df
    }
    
    # PBP_AUC 3D
    PBP_df <- metric_df_list[["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    PBP_AUC_df <- PBP_df[ , c("simulation", "reference", "target", "PBP_AUC")]
    metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
    
    # EBP_AUC 3D
    EBP_df <- metric_df_list[["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    EBP_AUC_df <- EBP_df[ , c("simulation", "cell_types", "EBP_AUC")]
    metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
    
    return(metric_df_list)
  }
  add_AUC_for_radii_gradient_metrics_to_metric_df_list2D <- function(metric_df_list) {
    gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "COO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "slice", "reference", "target", metric_AUC_name)
      }
      else {
        stop("Unknown metric?")
      }
      
      df <- metric_df_list[[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      df <- df[ , subset_colnames]
      metric_df_list[[metric_AUC_name]] <- df
    }
    
    # PBP_AUC 3D
    PBP_df <- metric_df_list[["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    PBP_AUC_df <- PBP_df[ , c("simulation", "slice", "reference", "target", "PBP_AUC")]
    metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
    
    # EBP_AUC 3D
    EBP_df <- metric_df_list[["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    EBP_AUC_df <- EBP_df[ , c("simulation", "slice", "cell_types", "EBP_AUC")]
    metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
    
    return(metric_df_list)
  }
  
  metric_df_list3D <- add_AUC_for_radii_gradient_metrics_to_metric_df_list3D(metric_df_list3D)
  metric_df_list2D <- add_AUC_for_radii_gradient_metrics_to_metric_df_list2D(metric_df_list2D)
  
  # Combine 3D and 2D metric df lists
  metric_df_list_combined <- list()
  
  metrics <- names(metric_df_list3D)
  for (metric in metrics) {
    metric_df_list3D[[metric]][["slice"]] <- 0 # Set slice index for 3D to be 0
    
    metric_df_combined <- rbind(metric_df_list3D[[metric]], metric_df_list2D[[metric]])
    
    metric_df_list_combined[[metric]] <- metric_df_combined
  }
  
  return(metric_df_list_combined)
}


# Set seed for consistency ----
set.seed(678999821)
# Running the functions ------------
updated_parameters_df <- generate_updated_parameters_for_arrangement_and_shape(
  n_simulations_per_parameter = 1000,
  bg_prop_A_range = c("min" = 0, "max" = 0.10),
  bg_prop_B_range = c("min" = 0, "max" = 0.10),
  E_radius_x_range = c("min" = 75, "max" = 125),
  E_radius_y_range = c("min" = 75, "max" = 125),
  E_radius_z_range = c("min" = 75, "max" = 125),
  N_width_range = c("min" = 25, "max" = 35),
  cluster_prop_A_range = c("min" = 0.3, "max" = 0.7),
  ring_width_factor_range = c("min" = 0.1, "max" = 0.2) ,
  cluster1_x_coord_range = c("min" = 125, "max" = 175)
)

SD2_metric_df_list <- analyse_SD2_simulations(updated_parameters_df)


saveRDS(SD2_metric_df_list, "~/R/reproduce_all_code_data/SD2_metric_df_list.RDS")
saveRDS(updated_parameters_df, "~/R/reproduce_all_code_data/SD2_parameters_df.RDS")


# Read data and set up ----
setwd("~/R/reproduce_all_code_data")
metric_df_list <- readRDS("SD2_metric_df_list.RDS")
parameters_df <- readRDS("SD2_parameters_df.RDS")


# Functions for plotting-----
# The only change:
# arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  BECOMES
# arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.04, 1))    0.02 --> 0.04 

# Subset metric_df_list and parameters_df
add_pair_to_metric_df <- function(metric_df, metric) {
  
  if (metric %in% c("EBSAC", "EBP_AUC")) {
    # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
    metric_df$pair <- gsub(',', '/', metric_df$cell_types)
  }
  else if (metric %in% c("ANE_AUC")) {
    # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
    metric_df$pair <- gsub(',', '/', metric_df$target)
  }
  else {
    # Add reference-target column
    metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
  } 
  return(metric_df)
}

get_parameters <- function(arrangement, shape) {
  parameters <- c()
  
  # Add other parameters specific to arrangement and shape
  if (arrangement == "mixed") {
    parameters <- c(parameters, "cluster_prop_A")
  }
  else if (arrangement == "ringed") {
    parameters <- c(parameters, "ring_width_factor")        
  }
  else if (arrangement == "separated") {
    parameters <- c(parameters, "distance")
  }
  
  if (shape == "ellipsoid") {
    parameters <- c(parameters, "E_volume")
  }
  else if (shape == "network") {
    parameters <- c(parameters, "N_width")
  }
  
  # All plots include bg_prop_A and bg_prop_B as parameters
  parameters <- c(parameters, "bg_prop_A", "bg_prop_B")
  
  return(parameters)
}

combine_metric_and_parameters_df <- function(metric_df, parameters_df) {
  pairs <- unique(metric_df$pair)
  metric_and_parameters_df <- data.frame()
  for (slice in unique(metric_df$slice)) {
    metric_and_parameters_df <- rbind(metric_and_parameters_df, 
                                      cbind(metric_df[metric_df$slice == slice, ], 
                                            parameters_df[rep(1:nrow(parameters_df), each = length(pairs)), ]))
  }
  return(metric_and_parameters_df)
}

plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    x[!(x %in% range(x, na.rm = T))] <- NA
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- "A/B"
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Subset for A/B
  metric_df <- metric_df[metric_df$pair == "A/B", ]
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          fig <- ggplot(plot_df, aes_string(parameter, metric)) +
            geom_point(
              data = subset(plot_df, slice != "0"),
              aes(color = slice),
              size = 0.5,
              alpha = 0.5
            ) +
            geom_point(
              data = subset(plot_df, slice == "0"),
              color = "black",
              size = 0.5,
              alpha = 0.5
            ) +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x  = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              legend.position = "none"
            ) +
            scale_x_continuous(
              breaks = pretty_breaks(n = 3),
              labels = sci_clean_threshold
            ) +
            scale_y_continuous(
              breaks = pretty_breaks(n = 3),
              labels = sci_clean_threshold
            ) +
            scale_color_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50"
              )
            )
          
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        
        
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold", size = 12)
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.1, 1))  + 
        theme(plot.margin = margin(2.5, 2.5, 1.25, 1.25),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

# Running the functions ----
metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "COO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")

setwd("~/R/reproduce_all_code_plots")

pdf("fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_A_B_scatter_plot.pdf", width = 17.6, height = 4.4)

for (metric in metrics) {
  
  fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                                          parameters_df,
                                                                                                                                          metric)
  print(fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot)
  
}
dev.off()

