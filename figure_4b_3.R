# SPIAT3D ----
calculate_all_gradient_cc_metrics3D <- function(spe,
                                                reference_cell_type,
                                                target_cell_types,
                                                radii,
                                                feature_colname) {
  
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
  fast_co_occurrence_df_colnames <- c("reference",
                                      target_cell_types)
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "neighbourhood_counts" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "neighbourhood_entropy" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "fast_co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(fast_co_occurrence_df_colnames))))
  colnames(result[["neighbourhood_counts"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["neighbourhood_entropy"]]) <- target_cell_types
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["fast_co_occurrence"]]) <- fast_co_occurrence_df_colnames
  
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
    result[["fast_co_occurrence"]][i, ] <- df[["fast_co_occurrence"]]
    
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
  result[["fast_co_occurrence"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  for (target_cell_type in names(df[["cross_G"]])) {
    result[["cross_G"]][[target_cell_type]]$radius <- radii
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
                 "fast_co_occurrence" = list())
  
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
  fast_co_occurrence_df_colnames <- c("reference",
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
  
  
  ## Fast_co_occurrence ---------------
  all_cell_types <- unique(spe[[feature_colname]])
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  fast_co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(fast_co_occurrence_df_colnames)))
  colnames(fast_co_occurrence_df) <- fast_co_occurrence_df_colnames
  fast_co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  n_cells_in_spe <- length(spe[[feature_colname]])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spe[[feature_colname]] == target_cell_type) <= 1) {
      fast_co_occurrence_df[[target_cell_type]] <- NA
      next
    }
    
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    
    target_cell_type_fast_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    result[[target_cell_type]] <- target_cell_type_fast_co_occurrence
  }
  result[["fast_co_occurrence"]] <- fast_co_occurrence_df
  
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
    
    # Account for case where not enough target cells in the data
    if (sum(spe[[feature_colname]] == target_cell_type) <= 1) {
      result[[target_cell_type]] <- NA
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

calculate_fast_co_occurrence_gradient3D <- function(spe,
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
    fast_co_occurrence_df <- calculate_fast_co_occurrence3D(spe,
                                                            reference_cell_type,
                                                            target_cell_types,
                                                            radii[i],
                                                            feature_colname)
    
    result[i, ] <- fast_co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_fast_co_occurrence_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_fast_co_occurrence3D <- function(spe,
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
  
  # Get total number of cells in spe
  n_cells_in_spe <- length(spe[[feature_colname]])
  
  for (target_cell_type in target_cell_types) {
    
    # Account for case where not enough target cells in the data
    if (sum(spe[[feature_colname]] == target_cell_type) <= 1) {
      result[[target_cell_type]] <- NA
      next
    }
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get total number of target cells in spe
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    
    # Get proportion of target cells in spe
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    
    # Get fast co-occurence value for target cell type
    target_cell_type_fast_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_fast_co_occurrence
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
                                                feature_colname) {
  
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
  fast_co_occurrence_df_colnames <- c("reference",
                                      target_cell_types)
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "neighbourhood_counts" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "neighbourhood_entropy" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "fast_co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(fast_co_occurrence_df_colnames))))
  colnames(result[["neighbourhood_counts"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["neighbourhood_entropy"]]) <- target_cell_types
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["fast_co_occurrence"]]) <- fast_co_occurrence_df_colnames
  
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
    result[["fast_co_occurrence"]][i, ] <- df[["fast_co_occurrence"]]
    
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
  result[["fast_co_occurrence"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  for (target_cell_type in names(df[["cross_G"]])) {
    result[["cross_G"]][[target_cell_type]]$radius <- radii
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
                 "fast_co_occurrence" = list())
  
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
  fast_co_occurrence_df_colnames <- c("reference",
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
  
  
  ## Fast_co_occurrence ---------------
  all_cell_types <- unique(spe[[feature_colname]])
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  fast_co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(fast_co_occurrence_df_colnames)))
  colnames(fast_co_occurrence_df) <- fast_co_occurrence_df_colnames
  fast_co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  n_cells_in_spe <- length(spe[[feature_colname]])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spe[[feature_colname]] == target_cell_type) <= 1) {
      fast_co_occurrence_df[[target_cell_type]] <- NA
      next
    }
    
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    
    target_cell_type_fast_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    result[[target_cell_type]] <- target_cell_type_fast_co_occurrence
  }
  result[["fast_co_occurrence"]] <- fast_co_occurrence_df
  
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
    
    # Account for case where not enough target cells in the data
    if (sum(spe[[feature_colname]] == target_cell_type) <= 1) {
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

calculate_fast_co_occurrence_gradient2D <- function(spe,
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
    fast_co_occurrence_df <- calculate_fast_co_occurrence3D(spe,
                                                            reference_cell_type,
                                                            target_cell_types,
                                                            radii[i],
                                                            feature_colname)
    
    result[i, ] <- fast_co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_fast_co_occurrence_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_fast_co_occurrence2D <- function(spe,
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
  
  # Get total number of cells in spe
  n_cells_in_spe <- length(spe[[feature_colname]])
  
  for (target_cell_type in target_cell_types) {
    
    # Account for case where not enough target cells in the data
    if (sum(spe[[feature_colname]] == target_cell_type) <= 1) {
      result[[target_cell_type]] <- NA
      next
    }
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get total number of target cells in spe
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    
    # Get proportion of target cells in spe
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    
    # Get fast co-occurence value for target cell type
    target_cell_type_fast_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_fast_co_occurrence
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


# spaSim3D ----
library(SpatialExperiment)
library(dbscan)
library(alphashape3d)
library(apcluster)
library(plotly)
library(dplyr)
library(reshape2)
library(gtools)
library(cowplot)
library(Hmisc)

plot_cells3D <- function(spe,
                         plot_cell_types = NULL,
                         plot_colours = NULL,
                         feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.null(plot_cell_types) && !is.character(plot_cell_types)) {
    stop("`plot_cell_types` is not a character vector or NULL.")
  }
  if (!is.null(plot_colours) && !is.character(plot_colours)) {
    stop("`plot_colours` is not a character vector or NULL.")
  }
  if (is.character(plot_colours)) {
    non_colours <- plot_colours[which(!(sapply(plot_colours, function(X) {
      tryCatch(is.matrix(col2rgb(X)),
               error = function(e) FALSE)
    })))]
    if (length(non_colours) > 0) {
      stop(paste("The following plot_colours are not colours:\n   ",
                 paste(non_colours, collapse = ", ")))
    }
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste(feature_colname, "is not a valid column in your spe object."))
  }
  
  ## Convert spe object to data frame
  df <- data.frame(SpatialExperiment::spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    warning("plot_cell_types not specified, all cell types found in the spe object will be used.")
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If no colours inputted, use viridis (D) palette
  if (is.null(plot_colours)) {
    warning("plot_colours not specified, viridis (D) palette will be used.")
    plot_colours <- viridis::viridis(n = length(plot_cell_types), option = "D")
  }
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## If cell types have been chosen, check they are found in the spe object
  spe_cell_types <- unique(spe[[feature_colname]])
  unknown_cell_types <- setdiff(plot_cell_types, spe_cell_types)
  
  if (length(unknown_cell_types) == length(plot_cell_types)) {
    stop("None of the plot_cell_types are found in the spe object")
  }
  
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following plot_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    plot_colours <- plot_colours[which(plot_cell_types %in% spe_cell_types)]
    plot_cell_types <- intersect(plot_cell_types, spe_cell_types)
  }
  
  ## Factor for feature column
  df[, "Cell.Type"] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Plot
  fig <- plotly::plot_ly(df,
                         type = "scatter3d",
                         mode = 'markers',
                         x = ~Cell.X.Position,
                         y = ~Cell.Y.Position,
                         z = ~Cell.Z.Position,
                         color = ~Cell.Type,
                         colors = plot_colours,
                         marker = list(size = 2))
  
  fig <- plotly::layout(fig, scene = list(xaxis = list(title = 'x', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                       titlefont = list(size = 20), tickfont = list(size = 15)),
                                          yaxis = list(title = 'y', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                       titlefont = list(size = 20), tickfont = list(size = 15)),
                                          zaxis = list(title = 'z', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                       titlefont = list(size = 20), tickfont = list(size = 15))))
  
  return(fig)
}

simulate_background_interactive <- function() {
  
  ### Message strings
  message_background <- paste("Hello spaSim3D user, how do you want your background cells to look like?\n
          1. Random pattern\n
          2. Ordered pattern\n\n",
                              "In a random pattern, cells are placed randomly...\n",
                              "In a ordered pattern, cells follow a regularly spaced in a hexagonal grid\n",
                              "To choose, please enter 1 or 2.\n", sep = "")
  
  message_background_random <- paste("We will need a few parameters before we can obtain the simulation\n",
                                     "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                     "    Number of cells (e.g. 10000 cells)\n",
                                     "    Minimum distance between cells (e.g. minimum distance of 2)\n",
                                     "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_background_ordered <- paste("We will need a few parameters before we can obtain the simulation\n",
                                      "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                      "    Number of cells (e.g. 10000 cells)\n",
                                      "    Amount of jitter (choose to give a bit or a lot of randomness)\n",
                                      "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_mixing <- paste("Would you like to MIX the background cells with chosen cell types randomly?\n")
  
  # Supportive functions
  get_integer_input_from_options <- function(integer_options) {
    
    first_integers <- integer_options[1:(length(integer_options) - 1)]
    last_integer <- integer_options[length(integer_options)]
    integers_string <- paste(paste(first_integers, collapse = ", "), "or", last_integer)
    
    prompt <- paste("Enter either ", integers_string, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      if (interactive) {
        user_input <- readline(prompt = prompt)
      } else {
        user_input <- 1
      }
      
      # Try converting to numeric
      integer_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(integer_value)) {
        message("Invalid input. Please enter a numeric integer value.")
      }
      # Numeric but not integer
      else if (integer_value%%1 != 0) {
        message("Non-integer input. Please enter an integer value.")
      }
      # Check if conversion was successful and value is in integer_options
      else if (integer_value %in% integer_options) {
        valid_input <- TRUE
        message("Valid input received!")
      }
      else {
        message(paste("Invalid input. Please enter only", integers_string))
      }
    }
    
    return(integer_value)
  }
  
  get_non_negative_numeric_input <- function(parameter) {
    
    prompt <- paste("Enter a non-negative numeric value for the ", parameter, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      non_negative_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(non_negative_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Negative input
      else if (non_negative_value < 0) {
        message("Negative input. Please enter a non-negative number")
      }
      # Should be correct input
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }
    
    return(non_negative_value)
  }
  
  get_numeric_between_input <- function(parameter,
                                        lower,
                                        upper) {
    
    prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(numeric_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Out of bounds input
      else if (numeric_value < lower || numeric_value > upper) {
        message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
      }
      # Should be correct
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }
    
    return(numeric_value)
  }
  
  get_positive_numeric_input <- function(parameter) {
    
    prompt <- paste("Enter a positive numeric value for the ", parameter, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      positive_numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(positive_numeric_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Non-positive input
      else if (positive_numeric_value <= 0) {
        message("Non-positive input. Please enter a positive number")
      }
      # Should be correct input
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }
    
    return(positive_numeric_value)
  }
  
  get_y_or_n_input <- function() {
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = "Enter either y or n: ")
      
      if (user_input %in% c("y", "n")) {
        valid_input <- TRUE
        message("Valid input received!")
      }
      else {
        message("Invalid input. Please enter either y or n.")
      }
    }
    
    return(user_input)
  }
  
  display_parameters <- function(parameter_values) {
    
    message("Your current inputs are:\n")
    
    display_message <- ""
    
    for (i in seq(length(parameter_values))) {
      display_message <- paste(display_message, "    ", i, ". ", names(parameter_values)[i], ": ", parameter_values[[i]], '\n', sep = "")
    }
    message(display_message)
  }
  
  get_cell_types_and_proportions_for_mixing <- function(simulated_spe) {
    
    message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."
    
    # Supportive functions
    get_numeric_between_input <- function(parameter,
                                          lower,
                                          upper) {
      
      prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")
      
      valid_input <- FALSE
      while (!valid_input) {
        user_input <- readline(prompt = prompt)
        # Try converting to numeric
        numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
        
        # Non-numeric input
        if (is.na(numeric_value)) {
          message("Invalid input. Please enter a numeric value.")
        }
        # Out of bounds input
        else if (numeric_value < lower || numeric_value > upper) {
          message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
        }
        # Should be correct
        else {
          valid_input <- TRUE
          message("Valid input received!")
        }
      }
      
      return(numeric_value)
    }
    
    get_y_or_n_input <- function() {
      
      valid_input <- FALSE
      while (!valid_input) {
        user_input <- readline(prompt = "Enter either y or n: ")
        
        if (user_input %in% c("y", "n")) {
          valid_input <- TRUE
          message("Valid input received!")
        }
        else {
          message("Invalid input. Please enter either y or n.")
        }
      }
      
      return(user_input)
    }
    
    ## Get cell types from user
    cell_types <- c()
    user_input <- ""
    message(message_get_cell_types)
    while (user_input != "stop") {
      
      user_input <- readline(prompt = "Enter a cell type, or enter 'stop': ")
      
      ## Ignore if user enters a blank string
      if (user_input == "") {
        
      }
      ## Add inputted cell type to cell_types vector
      else if (user_input != "stop") {
        cell_types <- c(cell_types, user_input)
        message(paste("Cell type added:", user_input))
      }
      ## User wants to stop but hasn't entered any cell types
      else if (user_input == "stop" && length(cell_types) == 0) {
        message("You have not entered any cell types. Try again\n")
        user_input <- ""
      }
      ## User wants to stop
      else {
        message(paste("Your cell types chosen are:", paste(cell_types, collapse = ", ")))
        
        ## Allow user to re-choose cell types
        message("Would like to re-choose these cell types?\n")
        user_input_y_or_n <- get_y_or_n_input()
        if (user_input_y_or_n == "y") {
          cell_types <- c()
          message(message_get_cell_types)
          user_input <- ""
        }
      }
    }
    
    ## Get cell proportions from user
    cell_proportions <- c()
    max_proportion <- 1
    i <- 1
    message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
    while (i <= length(cell_types)) {
      
      ## For the last cell type, we can figure out what the cell proportion must be
      if (i == length(cell_types)) {
        cell_proportions <- c(cell_proportions, max_proportion)
        message("Cell proportion for ", cell_types[i], " must be ", round(max_proportion, 5))
      }
      ## Add inputted cell proportion to cell_proportions vector
      else {
        cell_proportion <- get_numeric_between_input(paste("cell proportion of", cell_types[i], "cells"), 0, max_proportion)
        cell_proportions <- c(cell_proportions, cell_proportion)
        max_proportion <- 1 - sum(cell_proportions)
        message("Cell proportion for ", cell_types[i], " is ", cell_proportion)
      }
      i <- i + 1
      
      if (i > length(cell_types)) {
        ## Generate simulation
        message("Generating simulation...")
        simulated_spe <- simulate_mixing3D(simulated_spe,
                                           cell_types,
                                           cell_proportions,
                                           plot_image = F)
        
        fig <- plot_cells3D(simulated_spe)
        print(fig)
        
        if (length(cell_types) == 1) break # If there is only one cell type, proportion is always 1
        
        ## Allow user to re-choose cell proportions
        message("Would like to re-choose these cell proportions?\n")
        user_input_y_or_n <- get_y_or_n_input()
        if (user_input_y_or_n == "y") {
          cell_proportions <- c()
          max_proportion <- 1
          i <- 1
          message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
        }
      }
    }
    
    return(simulated_spe)
  }
  
  # Ask if user wants a 'random' or 'ordered' patterned background
  message(message_background)
  user_input_background <- get_integer_input_from_options(c(1, 2))
  
  ### Simulate random pattern
  if (user_input_background == 1) {
    
    # Get required parameters for a random background from user
    message(message_background_random)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "minimum distance between cells" = get_non_negative_numeric_input("minimum distance between cells"))
    display_parameters(parameter_values)
    
    # Generate random background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                        parameter_values[["length"]],
                                                        parameter_values[["width"]],
                                                        parameter_values[["height"]],
                                                        parameter_values[["minimum distance between cells"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["minimum distance between cells"]] <- get_non_negative_numeric_input("minimum distance between cells")
      
      # Generate random background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                          parameter_values[["length"]],
                                                          parameter_values[["width"]],
                                                          parameter_values[["height"]],
                                                          parameter_values[["minimum distance between cells"]])
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Simulate ordered pattern
  else if (user_input_background == 2) {
    
    # Get required parameters for a ordered background from user
    message(message_background_ordered)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "amount of jitter" = get_numeric_between_input("amount of jitter", 0, 1))
    display_parameters(parameter_values)
    
    # Generate ordered background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
                                                         parameter_values[["length"]],
                                                         parameter_values[["width"]],
                                                         parameter_values[["height"]],
                                                         parameter_values[["amount of jitter"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your inputs?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values))) # 5 different parameters
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["amount of jitter"]] <- get_numeric_between_input("amount of jitter", 0, 1)
      
      # Generate ordered background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
                                                           parameter_values[["length"]],
                                                           parameter_values[["width"]],
                                                           parameter_values[["height"]],
                                                           parameter_values[["amount of jitter"]])
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  ### Simulate mixing
  message(message_mixing)
  choose_cell_types_y_or_n <- get_y_or_n_input()
  if (choose_cell_types_y_or_n == "y") {
    simulated_spe <- get_cell_types_and_proportions_for_mixing(simulated_spe)
  }
  message("All done!")
  
  return(simulated_spe)
}

simulate_cluster_interactive <- function(simulated_spe = NULL) {
  
  ### Message strings
  message_no_simulated_spe <- paste("Hello spaSim3D user. Please input your simulated spe object into this function.\n",
                                    "If you don't have any, you can use the spaSim3D_background_integrator function")
  
  message_shape_choice <- paste("Hello spaSim-3D user, hopefully you can see a plot of your current spe object. What type of shape do you want your cluster to be?\n
          1. Sphere\n
          2. Ellipsoid\n
          3. Cylinder\n
          4. Network\n\n",
                                "To choose, please enter 1, 2, 3 or 4.\n", sep = "")
  
  message_sphere_cluster <- paste("We will need a few parameters to generate a sphere cluster\n",
                                  "    Radius\n",
                                  "    Coordinates of sphere centre: x, y and z\n",
                                  "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_ellipsoid_cluster <- paste("We will need a few parameters to generate a ellipsoid cluster\n",
                                     "    Radii: x, y and z\n",
                                     "    Coordinates of ellipsoid centre: x, y and z\n",
                                     "    Angle of rotation in the x-axis, y-axis and z-axis\n",
                                     "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_cylinder_cluster <- paste("We will need a few parameters to generate a cylinder cluster\n",
                                    "    Radius\n",
                                    "    Coordinates of the cylinder start point: x, y and z\n",
                                    "    Coordinates of the cylinder end point: x, y and z\n",
                                    "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_network_cluster <- paste("We will need a few parameters to generate a network cluster\n",
                                   "    Number of branches\n",
                                   "    Width of each branch: x, y and z\n",
                                   "    Radius spanned by the whole network: x, y and z\n",
                                   "    Coordinates of network centre: x, y and z\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  
  message_cluster_choice <- paste("You can customise your cluster further if you'd like:\n
          1. Add a ring\n
          2. Add a double ring\n
          3. Continue\n\n",
                                  "To choose, please enter 1, 2 or 3.\n", sep = "")
  
  # Supporting functions
  get_integer_greater_than_or_equal_input <- function(parameter,
                                                      lower) {
    
    prompt <- paste("Enter an integer in value greater than or equal to ", lower, " for the ", parameter, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      integer_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(integer_value)) {
        message("Invalid input. Please enter a numeric integer value.")
      }
      # Numeric but not integer
      else if (integer_value%%1 != 0) {
        message("Non-integer input. Please enter an integer value.")
      }
      # Integer but below the lower bound
      else if (integer_value < lower) {
        message("Out of bounds input. Please a number greater than or equal to ", lower)
      }
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }
    
    return(integer_value)
  }
  
  get_integer_input_from_options <- function(integer_options) {
    
    first_integers <- integer_options[1:(length(integer_options) - 1)]
    last_integer <- integer_options[length(integer_options)]
    integers_string <- paste(paste(first_integers, collapse = ", "), "or", last_integer)
    
    prompt <- paste("Enter either ", integers_string, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      integer_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(integer_value)) {
        message("Invalid input. Please enter a numeric integer value.")
      }
      # Numeric but not integer
      else if (integer_value%%1 != 0) {
        message("Non-integer input. Please enter an integer value.")
      }
      # Check if conversion was successful and value is in integer_options
      else if (integer_value %in% integer_options) {
        valid_input <- TRUE
        message("Valid input received!")
      }
      else {
        message(paste("Invalid input. Please enter only", integers_string))
      }
    }
    
    return(integer_value)
  }
  
  get_non_negative_numeric_input <- function(parameter) {
    
    prompt <- paste("Enter a non-negative numeric value for the ", parameter, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      non_negative_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(non_negative_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Negative input
      else if (non_negative_value < 0) {
        message("Negative input. Please enter a non-negative number")
      }
      # Should be correct input
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }
    
    return(non_negative_value)
  }
  
  get_positive_numeric_input <- function(parameter) {
    
    prompt <- paste("Enter a positive numeric value for the ", parameter, ": ", sep = "")
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      positive_numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
      
      # Non-numeric input
      if (is.na(positive_numeric_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Non-positive input
      else if (positive_numeric_value <= 0) {
        message("Non-positive input. Please enter a positive number")
      }
      # Should be correct input
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }
    
    return(positive_numeric_value)
  }
  
  get_y_or_n_input <- function() {
    
    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = "Enter either y or n: ")
      
      if (user_input %in% c("y", "n")) {
        valid_input <- TRUE
        message("Valid input received!")
      }
      else {
        message("Invalid input. Please enter either y or n.")
      }
    }
    
    return(user_input)
  }
  
  display_parameters <- function(parameter_values) {
    
    message("Your current inputs are:\n")
    
    display_message <- ""
    
    for (i in seq(length(parameter_values))) {
      display_message <- paste(display_message, "    ", i, ". ", names(parameter_values)[i], ": ", parameter_values[[i]], '\n', sep = "")
    }
    message(display_message)
  }
  
  get_cell_types_and_proportions_for_clusters <- function(simulated_spe,
                                                          simulate_function,
                                                          properties,
                                                          cell_type_option,
                                                          cell_proportion_option,
                                                          temp_cell_type) {
    
    message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."
    
    # Supporting functions
    get_y_or_n_input <- function() {
      
      valid_input <- FALSE
      while (!valid_input) {
        user_input <- readline(prompt = "Enter either y or n: ")
        
        if (user_input %in% c("y", "n")) {
          valid_input <- TRUE
          message("Valid input received!")
        }
        else {
          message("Invalid input. Please enter either y or n.")
        }
      }
      
      return(user_input)
    }
    
    get_numeric_between_input <- function(parameter,
                                          lower,
                                          upper) {
      
      prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")
      
      valid_input <- FALSE
      while (!valid_input) {
        user_input <- readline(prompt = prompt)
        # Try converting to numeric
        numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
        
        # Non-numeric input
        if (is.na(numeric_value)) {
          message("Invalid input. Please enter a numeric value.")
        }
        # Out of bounds input
        else if (numeric_value < lower || numeric_value > upper) {
          message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
        }
        # Should be correct
        else {
          valid_input <- TRUE
          message("Valid input received!")
        }
      }
      
      return(numeric_value)
    }
    
    ## Display the cell types currently found in simulated_spe to the user
    current_cell_types <- setdiff(unique(simulated_spe[["Cell.Type"]]), temp_cell_type)
    message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
    
    ## Get cell types from user
    cell_types <- c()
    user_input <- ""
    message(message_get_cell_types)
    while (user_input != "stop") {
      
      user_input <- readline(prompt = "Enter a cell type, or enter 'stop': ")
      
      ## Ignore if user enters a blank string
      if (user_input == "") {
        
      }
      ## Add inputted cell type to cell_types vector
      else if (user_input != "stop") {
        cell_types <- c(cell_types, user_input)
        message(paste("Cell type added:", user_input))
      }
      ## User wants to stop but hasn't entered any cell types
      else if (user_input == "stop" && length(cell_types) == 0) {
        message("You have not entered any cell types. Try again\n")
        user_input <- ""
      }
      ## User wants to stop
      else {
        message(paste("Your cell types chosen are:", paste(cell_types, collapse = ", ")))
        
        ## Allow user to re-choose cell types
        message("Would like to re-choose these cell types?\n")
        user_input_y_or_n <- get_y_or_n_input()
        if (user_input_y_or_n == "y") {
          message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
          cell_types <- c()
          message(message_get_cell_types)
          user_input <- ""
        }
      }
    }
    properties[[1]][[cell_type_option]] <- cell_types
    
    ## Get cell proportions from user
    cell_proportions <- c()
    max_proportion <- 1
    i <- 1
    message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
    while (i <= length(cell_types)) {
      
      ## For the last cell type, we can figure out what the cell proportion must be
      if (i == length(cell_types)) {
        cell_proportions <- c(cell_proportions, max_proportion)
        message("Cell proportion for ", cell_types[i], " must be ", round(max_proportion, 5))
      }
      ## Add inputted cell proportion to cell_proportions vector
      else {
        cell_proportion <- get_numeric_between_input(paste("cell proportion of", cell_types[i], "cells"), 0, max_proportion)
        cell_proportions <- c(cell_proportions, cell_proportion)
        max_proportion <- 1 - sum(cell_proportions)
        message("Cell proportion for ", cell_types[i], " is ", cell_proportion)
      }
      i <- i + 1
      
      if (i > length(cell_types)) {
        properties[[1]][[cell_proportion_option]] <- cell_proportions
        
        ## Convert spe object to data frame
        df <- data.frame(spatialCoords(simulated_spe), "Cell.Type" = simulated_spe[["Cell.Type"]])
        
        ## Just change the cell type of the temp_cell_type, no need to actually re-simulate
        df[["Cell.Type"]] <- ifelse(df[["Cell.Type"]] == temp_cell_type,
                                    sample(cell_types, size = length(df[["Cell.Type"]]), replace = TRUE, prob = cell_proportions),
                                    df[["Cell.Type"]])
        
        # Add Cell.ID column to data frame
        df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
        
        # Update current meta data
        metadata <- simulated_spe@metadata
        metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_type_option]] <- cell_types
        metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_proportion_option]] <- cell_proportions
        
        # Convert data frame to spe object
        simulated_spe_new <- SpatialExperiment::SpatialExperiment(
          assay = matrix(data = NA, nrow = 0, ncol = nrow(df)),
          colData = df,
          spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
          metadata = metadata)
        
        ## Generate simulation
        message("Generating simulation...")
        fig <- plot_cells3D(simulated_spe_new)
        print(fig)
        
        if (length(cell_types) == 1) break # If there is only one cell type, proportion is always 1
        
        ## Allow user to re-choose cell proportions
        message("Would like to re-choose these cell proportions?\n")
        user_input_y_or_n <- get_y_or_n_input()
        if (user_input_y_or_n == "y") {
          cell_proportions <- c()
          max_proportion <- 1
          i <- 1
          message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
        }
      }
    }
    
    return(list(data = simulated_spe_new, properties = properties))
  }
  
  ## Start with checking if the user has inputted spe object
  if (class(simulated_spe) != "SpatialExperiment") {
    stop(message_no_simulated_spe)
  }
  
  ## Plot the user's data so they can see what they already have
  fig <- plot_cells3D(simulated_spe)
  methods::show(fig)
  
  ## Get user's choice for shape type (sphere, ellipsoid, cylinder or network)
  message(message_shape_choice)
  user_input_shape <- get_integer_input_from_options(1:4)
  
  ### Sphere
  if (user_input_shape == 1) {
    # Get required parameters for a sphere cluster from user
    message(message_sphere_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)
    
    # Generate sphere cluster simulation using these parameters
    cluster_properties <- list(list(shape = "sphere",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "sphere",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Ellipsoid
  else if (user_input_shape == 2) {
    # Get required parameters for an ellipsoid cluster from user
    message(message_ellipsoid_cluster)
    parameter_values <- list("x radius" = get_positive_numeric_input("x radius"),
                             "y radius" = get_positive_numeric_input("y radius"),
                             "z radius" = get_positive_numeric_input("z radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"),
                             "x-axis rotation angle" = get_non_negative_numeric_input("x-axis rotation angle"),
                             "y-axis rotation angle" = get_non_negative_numeric_input("y-axis rotation angle"),
                             "z-axis rotation angle" = get_non_negative_numeric_input("z-axis rotation angle"))
    display_parameters(parameter_values)
    
    # Generate ellipsoid cluster simulation using these parameters
    cluster_properties <- list(list(shape = "ellipsoid",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radii = c(parameter_values[["x radius"]],
                                              parameter_values[["y radius"]],
                                              parameter_values[["z radius"]]),
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]]),
                                    axes_rotation = c(parameter_values[["x-axis rotation angle"]],
                                                      parameter_values[["y-axis rotation angle"]],
                                                      parameter_values[["z-axis rotation angle"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["x radius"]] <- get_positive_numeric_input("x radius")
      if (user_input_parameter_choice == 2) parameter_values[["y radius"]] <- get_positive_numeric_input("y radius")
      if (user_input_parameter_choice == 3) parameter_values[["z radius"]] <- get_positive_numeric_input("z radius")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["x-axis rotation angle"]] <- get_non_negative_numeric_input("x-axis rotation angle")
      if (user_input_parameter_choice == 8) parameter_values[["y-axis rotation angle"]] <- get_non_negative_numeric_input("y-axis rotation angle")
      if (user_input_parameter_choice == 9) parameter_values[["z-axis rotation angle"]] <- get_non_negative_numeric_input("z-axis rotation angle")
      
      display_parameters(parameter_values)
      
      # Generate ellipsoid cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "ellipsoid",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radii = c(parameter_values[["x radius"]],
                                                parameter_values[["y radius"]],
                                                parameter_values[["z radius"]]),
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]]),
                                      axes_rotation = c(parameter_values[["x-axis rotation angle"]],
                                                        parameter_values[["y-axis rotation angle"]],
                                                        parameter_values[["z-axis rotation angle"]])))
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Cylinder
  else if (user_input_shape == 3) {
    # Get required parameters for a cylinder cluster from user
    message(message_cylinder_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "start x coordinate" = get_non_negative_numeric_input("start x coordinate"),
                             "start y coordinate" = get_non_negative_numeric_input("start y coordinate"),
                             "start z coordinate" = get_non_negative_numeric_input("start z coordinate"),
                             "end x coordinate" = get_non_negative_numeric_input("end x coordinate"),
                             "end y coordinate" = get_non_negative_numeric_input("end y coordinate"),
                             "end z coordinate" = get_non_negative_numeric_input("end z coordinate"))
    display_parameters(parameter_values)
    
    # Generate cylinder cluster simulation using these parameters
    cluster_properties <- list(list(shape = "cylinder",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    start_loc = c(parameter_values[["start x coordinate"]],
                                                  parameter_values[["start y coordinate"]],
                                                  parameter_values[["start z coordinate"]]),
                                    end_loc = c(parameter_values[["end x coordinate"]],
                                                parameter_values[["end y coordinate"]],
                                                parameter_values[["end z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["start x coordinate"]] <- get_non_negative_numeric_input("start x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["start y coordinate"]] <- get_non_negative_numeric_input("start y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["start z coordinate"]] <- get_non_negative_numeric_input("start z coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["end x coordinate"]] <- get_non_negative_numeric_input("end x coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["end y coordinate"]] <- get_non_negative_numeric_input("end y coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["end z coordinate"]] <- get_non_negative_numeric_input("end z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate cylinder cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "cylinder",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      start_loc = c(parameter_values[["start x coordinate"]],
                                                    parameter_values[["start y coordinate"]],
                                                    parameter_values[["start z coordinate"]]),
                                      end_loc = c(parameter_values[["end x coordinate"]],
                                                  parameter_values[["end y coordinate"]],
                                                  parameter_values[["end z coordinate"]])))
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Network
  else if (user_input_shape == 4) {
    # Get required parameters for a network cluster from user
    message(message_network_cluster)
    parameter_values <- list("number of branches" = get_integer_greater_than_or_equal_input("number of branches", 2),
                             "width of branch" = get_positive_numeric_input("width of branch"),
                             "radius spanned by network" = get_positive_numeric_input("radius spanned by network"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)
    
    # Generate network cluster simulation using these parameters
    cluster_properties <- list(list(shape = "network",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    n_edges = parameter_values[["number of branches"]],
                                    width = parameter_values[["width of branch"]],
                                    radius = parameter_values[["radius spanned by network"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["number of branches"]] <- get_integer_greater_than_or_equal_input("number of branches", 2)
      if (user_input_parameter_choice == 2) parameter_values[["width of branch"]] <- get_positive_numeric_input("width of branch")
      if (user_input_parameter_choice == 3) parameter_values[["radius spanned by network"]] <- get_positive_numeric_input("radius spanned by network")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "network",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      n_edges = parameter_values[["number of branches"]],
                                      width = parameter_values[["width of branch"]],
                                      radius = parameter_values[["radius spanned by network"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  # Allow user to change the cell composition of the cluster
  message("Let's change the cell composition of this cluster")
  simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                  simulate_clusters3D,
                                                                                  cluster_properties,
                                                                                  "cluster_cell_types",
                                                                                  "cluster_cell_proportions",
                                                                                  "Cluster")
  simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  cluster_properties <- simulated_spe_new_and_properties[["properties"]]
  
  
  
  ## Get user's choice for cluster type (ringed, double ringed or continue)
  message(message_cluster_choice)
  user_input_cluster <- get_integer_input_from_options(1:3)
  
  ### Ring
  if (user_input_cluster == 1) {
    # Get width of ring from user
    message("For a single ring, we needs its width.\n")
    
    # Generate cluster with ring simulation using this width
    cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
    cluster_properties[[1]][["ring_cell_types"]] <- c("Ring")
    cluster_properties[[1]][["ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_spe_new <- simulate_rings3D(simulated_spe,
                                          cluster_properties,
                                          plot_image = TRUE,
                                          plot_cell_types = NULL,
                                          plot_colours = NULL)
    
    # Allow user the option to change the ring width
    message("Would you like to change the ring width?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_rings3D(simulated_spe,
                                            cluster_properties,
                                            plot_image = TRUE,
                                            plot_cell_types = NULL,
                                            plot_colours = NULL)
      
      message("Would you like to the ring width?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the ring
    message("Let's change the cell composition of the ring")
    
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_rings3D,
                                                                                    cluster_properties,
                                                                                    "ring_cell_types",
                                                                                    "ring_cell_proportions",
                                                                                    "Ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  }
  ### Double ring
  else if (user_input_cluster == 2) {
    # Get width of inner and outer ring from user
    message("For a double ring, we needs the width of the inner and outer ring.\n")
    
    # Generate cluster with double ring simulation using both widths
    parameter_values <- list("inner ring width" = get_positive_numeric_input("inner ring width"),
                             "outer ring width" = get_positive_numeric_input("outer ring width"))
    display_parameters(parameter_values)
    
    cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
    cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
    cluster_properties[[1]][["inner_ring_cell_types"]] <- c("Inner ring")
    cluster_properties[[1]][["inner_ring_cell_proportions"]] <- 1
    cluster_properties[[1]][["outer_ring_cell_types"]] <- c("Outer ring")
    cluster_properties[[1]][["outer_ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                 cluster_properties,
                                                 plot_image = TRUE,
                                                 plot_cell_types = NULL,
                                                 plot_colours = NULL)
    
    # Allow user the option to change the widths of the inner or outer ring
    message("Would you like to change the widths of the inner or outer ring?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      if (user_input_parameter_choice == 1) parameter_values[["inner ring width"]] <- get_positive_numeric_input("inner ring width")
      if (user_input_parameter_choice == 2) parameter_values[["outer ring width"]] <- get_positive_numeric_input("outer ring width")
      
      cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
      cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
      
      display_parameters(parameter_values)
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                   cluster_properties,
                                                   plot_image = TRUE,
                                                   plot_cell_types = NULL,
                                                   plot_colours = NULL)
      
      message("Would you like to change the widths of the inner or outer ring?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the inner ring
    message("Let's change the cell composition of the inner ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "inner_ring_cell_types",
                                                                                    "inner_ring_cell_proportions",
                                                                                    "Inner ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
    cluster_properties <- simulated_spe_new_and_properties[["properties"]]
    
    # Allow user to change the cell composition of the outer ring
    message("Let's change the cell composition of the outer ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "outer_ring_cell_types",
                                                                                    "outer_ring_cell_proportions",
                                                                                    "Outer ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  }
  ### Continue
  else if (user_input_cluster == 3) {
    
  }
  
  message("All done!")
  return(simulated_spe_new)
}

simulate_clusters3D <- function(spe,
                                cluster_properties_list,
                                plot_image = TRUE,
                                plot_cell_types = NULL,
                                plot_colours = NULL) {
  
  # Check shape variable of cluster_properties
  shapes <- sapply(cluster_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`cluster_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(cluster_properties_list))) {
    
    shape <- shapes[[i]]
    
    ### Sphere shape
    if (shape == "sphere") {
      spe <- simulate_sphere_cluster(spe, cluster_properties_list[[i]])
    }
    
    ### Ellipsoid shape
    if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_cluster(spe, cluster_properties_list[[i]])
    }
    
    ### Cylinder shape
    if (shape == "cylinder") {
      spe <- simulate_cylinder_cluster(spe, cluster_properties_list[[i]])
    }
    
    ### Network shape
    if (shape == "network") {
      spe <- simulate_network_cluster(spe, cluster_properties_list[[i]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}

simulate_cylinder_cluster <- function(spe,
                                      cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get cylinder properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- SpatialExperiment::spatialCoords(spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}

simulate_cylinder_dr <- function(spe,
                                 dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get cylinder double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  start_loc <- dr_properties$start_loc
  end_loc <- dr_properties$end_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- SpatialExperiment::spatialCoords(spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  # Start with cells in outer ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width + outer_ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Start with cells in inner ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}

simulate_cylinder_ring <- function(spe,
                                   ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get cylinder ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  start_loc <- ring_properties$start_loc
  end_loc <- ring_properties$end_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- SpatialExperiment::spatialCoords(spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  # Start with cells in ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(spe)
}

simulate_double_rings3D <- function(spe,
                                    dr_properties_list,
                                    plot_image = TRUE,
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  # Check shape variable of dr_properties_list
  shapes <- sapply(dr_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`dr_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(dr_properties_list))) {
    
    shape <- shapes[[i]]
    
    ### Sphere shape with double ring
    if (shape == "sphere") {
      spe <- simulate_sphere_dr(spe, dr_properties_list[[i]])
    }
    
    ### Ellipsoid shape with double ring
    if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_dr(spe, dr_properties_list[[i]])
    }
    
    ### Cylinder shape with double ring
    if (shape == "cylinder") {
      spe <- simulate_cylinder_dr(spe, dr_properties_list[[i]])
    }
    
    ### Network shape with double ring
    if (shape == "network") {
      spe <- simulate_network_dr(spe, dr_properties_list[[i]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}

simulate_ellipsoid_cluster <- function(spe,
                                       cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get ellipsoid properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  x_radius <- cluster_properties$radii[1]
  y_radius <- cluster_properties$radii[2]
  z_radius <- cluster_properties$radii[3]
  centre_loc <- cluster_properties$centre_loc
  theta <- cluster_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- cluster_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- cluster_properties$axes_rotation[3] * (pi/180) # rotation in z-axis
  
  # Get rotation matrices for rotation in the y-z plane (T2), x-z plane (T3) and x-y plane (T4)
  T1 <- matrix(data = c(1, 0, 0,
                        0, cos(theta), -sin(theta),
                        0, sin(theta), cos(theta)), nrow = 3, ncol = 3, byrow = TRUE)
  T2 <- matrix(data = c(cos(alpha), 0, -sin(alpha),
                        0, 1, 0,
                        sin(alpha), 0, cos(alpha)), nrow = 3, ncol = 3, byrow = TRUE)
  T3 <- matrix(data = c(cos(beta), -sin(beta), 0,
                        sin(beta), cos(beta), 0,
                        0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Get translation matrix from ellipsoid centre (same as centre...)
  T4 <- centre_loc
  
  ## Change cell types in the ellipsoid cluster
  # Get spatial coords from spe (rows are x, y, z, columns are each cell)
  spe_coords <- t(SpatialExperiment::spatialCoords(spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}

simulate_ellipsoid_dr <- function(spe,
                                  dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get ellipsoid double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  x_radius <- dr_properties$radii[1]
  y_radius <- dr_properties$radii[2]
  z_radius <- dr_properties$radii[3]
  centre_loc <- dr_properties$centre_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  theta <- dr_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- dr_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- dr_properties$axes_rotation[3] * (pi/180) # rotation in z-axis
  
  # Get rotation matrices for rotation in the y-z plane (T2), x-z plane (T3) and x-y plane (T4)
  T1 <- matrix(data = c(1, 0, 0,
                        0, cos(theta), -sin(theta),
                        0, sin(theta), cos(theta)), nrow = 3, ncol = 3, byrow = TRUE)
  T2 <- matrix(data = c(cos(alpha), 0, -sin(alpha),
                        0, 1, 0,
                        sin(alpha), 0, cos(alpha)), nrow = 3, ncol = 3, byrow = TRUE)
  T3 <- matrix(data = c(cos(beta), -sin(beta), 0,
                        sin(beta), cos(beta), 0,
                        0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Get translation matrix from ellipsoid centre (same as centre...)
  T4 <- centre_loc
  
  ## Change cell types in the ellipsoid cluster
  # Get spatial coords from spe (rows are x, y, z, columns are each cell)
  spe_coords <- t(SpatialExperiment::spatialCoords(spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  
  # Start with cells in outer ring
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width + outer_ring_width))^2 +
                                 (y / (y_radius + inner_ring_width + outer_ring_width))^2 +
                                 (z / (z_radius + inner_ring_width + outer_ring_width))^2 <= 1,
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in inner ring
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width))^2 +
                                 (y / (y_radius + inner_ring_width))^2 +
                                 (z / (z_radius + inner_ring_width))^2 <= 1,
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}

simulate_ellipsoid_ring <- function(spe,
                                    ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get ellipsoid ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  x_radius <- ring_properties$radii[1]
  y_radius <- ring_properties$radii[2]
  z_radius <- ring_properties$radii[3]
  centre_loc <- ring_properties$centre_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  theta <- ring_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- ring_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- ring_properties$axes_rotation[3] * (pi/180) # rotation in z-axis
  
  # Get rotation matrices for rotation in the y-z plane (T2), x-z plane (T3) and x-y plane (T4)
  T1 <- matrix(data = c(1, 0, 0,
                        0, cos(theta), -sin(theta),
                        0, sin(theta), cos(theta)), nrow = 3, ncol = 3, byrow = TRUE)
  T2 <- matrix(data = c(cos(alpha), 0, -sin(alpha),
                        0, 1, 0,
                        sin(alpha), 0, cos(alpha)), nrow = 3, ncol = 3, byrow = TRUE)
  T3 <- matrix(data = c(cos(beta), -sin(beta), 0,
                        sin(beta), cos(beta), 0,
                        0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Get translation matrix from ellipsoid centre (same as centre...)
  T4 <- centre_loc
  
  ## Change cell types in the ellipsoid cluster
  # Get spatial coords from spe (rows are x, y, z, columns are each cell)
  spe_coords <- t(SpatialExperiment::spatialCoords(spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  # Start with cells in ring
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + ring_width))^2 +
                                 (y / (y_radius + ring_width))^2 +
                                 (z / (z_radius + ring_width))^2 <= 1,
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(spe)
}

simulate_mixing3D <- function(spe,
                              cell_types,
                              cell_proportions,
                              plot_image = TRUE,
                              plot_cell_types = NULL,
                              plot_colours = NULL) {
  
  # Check input parameters
  input_parameters <- list("spe" = spe,
                           "cell_types" = cell_types,
                           "cell_proportions" = cell_proportions,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Apply mixing
  spe[["Cell.Type"]] <- sample(cell_types, size = ncol(spe), replace = TRUE, prob = cell_proportions)
  
  spe@metadata[["simulation"]][["background"]][["cell_types"]] <- cell_types
  spe@metadata[["simulation"]][["background"]][["cell_proportions"]] <- cell_proportions
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}

simulate_network_cluster <- function(spe,
                                     cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  n_edges <- cluster_properties$n_edges
  width <- cluster_properties$width
  centre_loc <- cluster_properties$centre_loc
  radius <- cluster_properties$radius
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc.
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Lets stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated")
  random_coords <- data.frame(x = runif(1000, centre_loc[1] - radius, centre_loc[1] + radius),
                              y = runif(1000, centre_loc[2] - radius, centre_loc[2] + radius),
                              z = runif(1000, centre_loc[3] - radius, centre_loc[3] + radius))
  
  # Then subset points which are inside the sphere
  random_coords <- random_coords[(random_coords$x - centre_loc[1])^2 +
                                   (random_coords$y - centre_loc[2])^2 +
                                   (random_coords$z- centre_loc[3])^2 <= radius^2, ]
  
  ## Subset further and pick 'n_vertices' coords to represent the vertices
  random_coords <- dplyr::sample_n(random_coords, n_vertices)
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(random_coords)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  # Input is the adjacency matrix of the graph (i.e. output from -1 * apcluster::negDistMat(df of coords))
  prims_algorithm <- function(graph) {
    
    # Number of vertices is number of points
    num_vertices <- nrow(graph)
    
    # Start with no vertices selected except first
    selected <- rep(FALSE, num_vertices)
    selected[1] <- TRUE
    
    # Create tree_edge matrix. Currently zero, each row represents the two vertices the edge joins
    tree_edges <- matrix(0,
                         nrow = num_vertices - 1,
                         ncol = 2)
    
    # Iterate until we select enough edges (one less than the number of vertices for a MST)
    num_edges <- 0
    while (num_edges < num_vertices - 1) {
      # Set initial temp values for weight and vertex
      min_weight <- Inf
      min_vertex <- -1
      
      # Iterate through each currently selected vertex
      for (i in seq(num_vertices)) {
        
        # Found a currently selected vertex
        if (selected[i] == TRUE) {
          
          # Iterate through each unselected vertex and find the nearest one
          for (j in seq(num_vertices)) {
            if (!selected[j] && graph[i, j] < min_weight) {
              min_weight <- graph[i, j]
              min_vertex <- j
              curr_vertex <- i
            }
          }
        }
      }
      
      # Current edge connects the min_vertex and curr_vertex
      tree_edges[num_edges + 1, ] <- c(min_vertex, curr_vertex)
      selected[min_vertex] <- TRUE
      num_edges <- num_edges + 1
    }
    return(tree_edges)
  }
  
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  get_tree_depth <- function(tree_edges) {
    
    tree_edges <- data.frame(tree_edges)
    colnames(tree_edges) <- c("vertex1", "vertex2")
    
    # Set the initial depth of each tree_edge to be NA.
    tree_edges$depth <- NA
    
    # Get vertices on the 'outskirts' of MST (leaf_vertices which have a depth of 1)
    tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
    leaf_vertices <- as.numeric(names(table(tree_vertices))[table(tree_vertices) == 1])
    
    # Start with leaf_vertices
    curr_vertices <- leaf_vertices
    curr_depth <- 1
    
    while (NA %in% tree_edges$depth) {
      
      # New vertices will be those adjacent to the current vertices
      new_vertices <- c()
      
      # Check each current vertex
      for (vertex in curr_vertices) {
        # Start with vertex1
        curr_edges <- which(tree_edges$vertex1 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex2"])
        
        # Then vertex2
        curr_edges <- which(tree_edges$vertex2 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex1"])
        
        # Only keep unique vertices
        new_vertices <- unique(new_vertices)
      }
      
      curr_depth <- curr_depth + 1
      curr_vertices <- new_vertices
    }
    
    return(tree_edges)
  }
  
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_cluster_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (curr_width < 0) curr_width <- 0
    
    network_cluster_properties[[i]] <- list(shape = "cylinder",
                                            cluster_cell_types = cluster_cell_types,
                                            cluster_cell_proportions = cluster_cell_proportions,
                                            radius = curr_width,
                                            start_loc = start_loc,
                                            end_loc = end_loc)
  }
  
  network_spe <- simulate_clusters3D(spe,
                                     cluster_properties = network_cluster_properties,
                                     plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  cluster_properties[["cylinders"]] <- network_cluster_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- cluster_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}

simulate_network_dr <- function(spe,
                                dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  n_edges <- dr_properties$n_edges
  width <- dr_properties$width
  centre_loc <- dr_properties$centre_loc
  radius <- dr_properties$radius
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc.
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Let's stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated.")
  random_coords <- data.frame(x = runif(1000, centre_loc[1] - radius, centre_loc[1] + radius),
                              y = runif(1000, centre_loc[2] - radius, centre_loc[2] + radius),
                              z = runif(1000, centre_loc[3] - radius, centre_loc[3] + radius))
  
  # Then subset points which are inside the sphere
  random_coords <- random_coords[(random_coords$x - centre_loc[1])^2 +
                                   (random_coords$y - centre_loc[2])^2 +
                                   (random_coords$z- centre_loc[3])^2 <= radius^2, ]
  
  ## Subset further and pick 'n_vertices' coords to represent the vertices
  random_coords <- dplyr::sample_n(random_coords, n_vertices)
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(random_coords)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  # Input is the adjacency matrix of the graph (i.e. output from -1 * apcluster::negDistMat(df of coords))
  prims_algorithm <- function(graph) {
    
    # Number of vertices is number of points
    num_vertices <- nrow(graph)
    
    # Start with no vertices selected except first
    selected <- rep(FALSE, num_vertices)
    selected[1] <- TRUE
    
    # Create tree_edge matrix. Currently zero, each row represents the two vertices the edge joins
    tree_edges <- matrix(0,
                         nrow = num_vertices - 1,
                         ncol = 2)
    
    # Iterate until we select enough edges (one less than the number of vertices for a MST)
    num_edges <- 0
    while (num_edges < num_vertices - 1) {
      # Set initial temp values for weight and vertex
      min_weight <- Inf
      min_vertex <- -1
      
      # Iterate through each currently selected vertex
      for (i in seq(num_vertices)) {
        
        # Found a currently selected vertex
        if (selected[i] == TRUE) {
          
          # Iterate through each unselected vertex and find the nearest one
          for (j in seq(num_vertices)) {
            if (!selected[j] && graph[i, j] < min_weight) {
              min_weight <- graph[i, j]
              min_vertex <- j
              curr_vertex <- i
            }
          }
        }
      }
      
      # Current edge connects the min_vertex and curr_vertex
      tree_edges[num_edges + 1, ] <- c(min_vertex, curr_vertex)
      selected[min_vertex] <- TRUE
      num_edges <- num_edges + 1
    }
    return(tree_edges)
  }
  
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  get_tree_depth <- function(tree_edges) {
    
    tree_edges <- data.frame(tree_edges)
    colnames(tree_edges) <- c("vertex1", "vertex2")
    
    # Set the initial depth of each tree_edge to be NA.
    tree_edges$depth <- NA
    
    # Get vertices on the 'outskirts' of MST (leaf_vertices which have a depth of 1)
    tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
    leaf_vertices <- as.numeric(names(table(tree_vertices))[table(tree_vertices) == 1])
    
    # Start with leaf_vertices
    curr_vertices <- leaf_vertices
    curr_depth <- 1
    
    while (NA %in% tree_edges$depth) {
      
      # New vertices will be those adjacent to the current vertices
      new_vertices <- c()
      
      # Check each current vertex
      for (vertex in curr_vertices) {
        # Start with vertex1
        curr_edges <- which(tree_edges$vertex1 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex2"])
        
        # Then vertex2
        curr_edges <- which(tree_edges$vertex2 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex1"])
        
        # Only keep unique vertices
        new_vertices <- unique(new_vertices)
      }
      
      curr_depth <- curr_depth + 1
      curr_vertices <- new_vertices
    }
    
    return(tree_edges)
  }
  
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_dr_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_dr_properties[[i]] <- list(shape = "cylinder",
                                       cluster_cell_types = cluster_cell_types,
                                       cluster_cell_proportions = cluster_cell_proportions,
                                       radius = curr_width,
                                       start_loc = start_loc,
                                       end_loc = end_loc,
                                       inner_ring_cell_types = inner_ring_cell_types,
                                       inner_ring_cell_proportions = inner_ring_cell_proportions,
                                       inner_ring_width = inner_ring_width,
                                       outer_ring_cell_types = outer_ring_cell_types,
                                       outer_ring_cell_proportions = outer_ring_cell_proportions,
                                       outer_ring_width = outer_ring_width)
  }
  
  network_spe <- simulate_double_rings3D(spe,
                                         dr_properties = network_dr_properties,
                                         plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  dr_properties[["cylinders"]] <- network_dr_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- dr_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}

simulate_network_ring <- function(spe,
                                  ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  n_edges <- ring_properties$n_edges
  width <- ring_properties$width
  centre_loc <- ring_properties$centre_loc
  radius <- ring_properties$radius
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc.
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Lets stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated")
  random_coords <- data.frame(x = runif(1000, centre_loc[1] - radius, centre_loc[1] + radius),
                              y = runif(1000, centre_loc[2] - radius, centre_loc[2] + radius),
                              z = runif(1000, centre_loc[3] - radius, centre_loc[3] + radius))
  
  # Then subset points which are inside the sphere
  random_coords <- random_coords[(random_coords$x - centre_loc[1])^2 +
                                   (random_coords$y - centre_loc[2])^2 +
                                   (random_coords$z- centre_loc[3])^2 <= radius^2, ]
  
  ## Subset further and pick 'n_vertices' coords to represent the vertices
  random_coords <- dplyr::sample_n(random_coords, n_vertices)
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(random_coords)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  # Input is the adjacency matrix of the graph (i.e. output from -1 * apcluster::negDistMat(df of coords))
  prims_algorithm <- function(graph) {
    
    # Number of vertices is number of points
    num_vertices <- nrow(graph)
    
    # Start with no vertices selected except first
    selected <- rep(FALSE, num_vertices)
    selected[1] <- TRUE
    
    # Create tree_edge matrix. Currently zero, each row represents the two vertices the edge joins
    tree_edges <- matrix(0,
                         nrow = num_vertices - 1,
                         ncol = 2)
    
    # Iterate until we select enough edges (one less than the number of vertices for a MST)
    num_edges <- 0
    while (num_edges < num_vertices - 1) {
      # Set initial temp values for weight and vertex
      min_weight <- Inf
      min_vertex <- -1
      
      # Iterate through each currently selected vertex
      for (i in seq(num_vertices)) {
        
        # Found a currently selected vertex
        if (selected[i] == TRUE) {
          
          # Iterate through each unselected vertex and find the nearest one
          for (j in seq(num_vertices)) {
            if (!selected[j] && graph[i, j] < min_weight) {
              min_weight <- graph[i, j]
              min_vertex <- j
              curr_vertex <- i
            }
          }
        }
      }
      
      # Current edge connects the min_vertex and curr_vertex
      tree_edges[num_edges + 1, ] <- c(min_vertex, curr_vertex)
      selected[min_vertex] <- TRUE
      num_edges <- num_edges + 1
    }
    return(tree_edges)
  }
  
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  get_tree_depth <- function(tree_edges) {
    
    tree_edges <- data.frame(tree_edges)
    colnames(tree_edges) <- c("vertex1", "vertex2")
    
    # Set the initial depth of each tree_edge to be NA.
    tree_edges$depth <- NA
    
    # Get vertices on the 'outskirts' of MST (leaf_vertices which have a depth of 1)
    tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
    leaf_vertices <- as.numeric(names(table(tree_vertices))[table(tree_vertices) == 1])
    
    # Start with leaf_vertices
    curr_vertices <- leaf_vertices
    curr_depth <- 1
    
    while (NA %in% tree_edges$depth) {
      
      # New vertices will be those adjacent to the current vertices
      new_vertices <- c()
      
      # Check each current vertex
      for (vertex in curr_vertices) {
        # Start with vertex1
        curr_edges <- which(tree_edges$vertex1 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex2"])
        
        # Then vertex2
        curr_edges <- which(tree_edges$vertex2 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex1"])
        
        # Only keep unique vertices
        new_vertices <- unique(new_vertices)
      }
      
      curr_depth <- curr_depth + 1
      curr_vertices <- new_vertices
    }
    
    return(tree_edges)
  }
  
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_ring_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_ring_properties[[i]] <- list(shape = "cylinder",
                                         cluster_cell_types = cluster_cell_types,
                                         cluster_cell_proportions = cluster_cell_proportions,
                                         radius = curr_width,
                                         start_loc = start_loc,
                                         end_loc = end_loc,
                                         ring_cell_types = ring_cell_types,
                                         ring_cell_proportions = ring_cell_proportions,
                                         ring_width = ring_width)
  }
  
  network_spe <- simulate_rings3D(spe,
                                  ring_properties = network_ring_properties,
                                  plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  ring_properties[["cylinders"]] <- network_ring_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- ring_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}

simulate_ordered_background_cells3D <- function(n_cells,
                                                length,
                                                width,
                                                height,
                                                jitter_proportion = 0.25,
                                                background_cell_type = "Others",
                                                plot_image = TRUE) {
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "jitter_proportion" = jitter_proportion,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Obtain distance between each point using MAGIC formula
  d_cells <- ((sqrt(2) * length * width * height)/n_cells)^(1/3)
  
  # Get distance between rows, columns and layers using 'd_cells'
  d_rows <- d_cells
  d_cols <- (sqrt(3) / 2) * d_cells
  d_lays <- (sqrt(6) / 3) * d_cells
  
  # Get number of rows, columns and layers
  n_rows <- round(length / d_rows)
  n_cols <- round(width / d_cols)
  n_lays <- round(height / d_lays)
  
  # Step 0. Assume points are on a 3D rectangular grid
  rows <- rep(seq(n_rows), n_cols * n_lays) * d_rows
  cols <- rep(rep(seq(n_cols), each = n_rows), n_lays) * d_cols
  lays <- rep(seq(n_lays), each = n_rows * n_cols) * d_lays
  
  # Step 1. For every odd sheet, every even row shifts by d_cells/2 right
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2)
  }
  else {
    shift <- c(rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2), rep(0, n_rows))
  }
  rows <- rows + c(shift, rep(0, n_rows * n_cols)) # Shift each even row by d_cells/2 right
  
  # Step 2. For every even sheet, odd rows shift d_cells/2 right, all rows shift d_cells/(2*sqrt(3)) up
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2)
  }
  else {
    shift <- c(rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2), rep(d_cells/2, n_rows))
  }
  rows <- rows + c(rep(0, n_rows * n_cols), shift) # Shift each odd row by d_cells/2 right
  cols <- cols + rep(c(0, d_cells/(2 * sqrt(3))), each = n_rows * n_cols) # Shift all rows by d_cells/(2*sqrt(3)) up
  
  # Get total number of cells (should be roughly equal to n_cells)
  n_total <- n_rows * n_cols * n_lays
  
  # Add randomness to the location of the cells
  jitter <- jitter_proportion * d_cells # Jitter is proportional to distance between points in hexagonal grid
  jitter_row <- runif(n_total, -jitter, jitter)
  jitter_col <- runif(n_total, -jitter, jitter)
  jitter_lay <- runif(n_total, -jitter, jitter)
  
  rows <- rows + jitter_row
  cols <- cols + jitter_col
  lays <- lays + jitter_lay
  
  # Put data into data frame
  df <- data.frame("Cell.X.Position" = rows,
                   "Cell.Y.Position" = cols,
                   "Cell.Z.Position" = lays,
                   "Cell.Type" = background_cell_type)
  df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "normal",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "jitter_proportion" = jitter_proportion,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}

simulate_random_background_cells3D <- function(n_cells,
                                               length,
                                               width,
                                               height,
                                               minimum_distance_between_cells,
                                               background_cell_type = "Others",
                                               plot_image = TRUE) {
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "minimum_distance_between_cells" = minimum_distance_between_cells,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Need to over-sample as cells which are too close will be removed later
  n_cells_inflated <- n_cells * 2
  
  spatial_df <- data.frame("Cell.X.Position" = runif(n_cells_inflated, 0, length),
                           "Cell.Y.Position" = runif(n_cells_inflated, 0, width),
                           "Cell.Z.Position" = runif(n_cells_inflated, 0, height))
  
  # Add integer rownames to data frame - each cell is labelled by an integer
  rownames(spatial_df) <- seq(nrow(spatial_df))
  
  ### Check if all other cells are to close to the current cell
  # Use frNN function: for each point, get all points within min_d of it
  spatial_df_distances <- dbscan::frNN(spatial_df,
                                       eps = minimum_distance_between_cells,
                                       query = NULL,
                                       sort = FALSE)
  
  # For each cell, get all other cells which were within 'minimum_distance_between_cells'
  spatial_df_distances_ids <- spatial_df_distances$id
  
  # Filter out zero length cells
  spatial_df_distances_ids <- Filter(function(x) length(x) != 0, spatial_df_distances_ids)
  
  # Get integer labels for the remaining cells
  spatial_df_distances_ids_names <- as.integer(names(spatial_df_distances_ids))
  
  # Determine which cells should be chosen from spatial_df
  cells_chosen <- rep(T, nrow(spatial_df))
  for (i in seq_len(length(spatial_df_distances_ids))) {
    cells_too_close <- spatial_df_distances_ids[[i]]
    
    if (cells_chosen[spatial_df_distances_ids_names[i]]) cells_chosen[cells_too_close] <- F
  }
  
  spatial_df <- spatial_df[cells_chosen, ]
  
  # If number of cells remaining is still higher than n_cells, randomly subset n_cells
  if (nrow(spatial_df) > n_cells) {
    spatial_df <- dplyr::sample_n(spatial_df, n_cells)
  }
  else {
    warning("The number of cells in your object could not reach the desired number of cells, consider lowering the number of cells, the minimum distance between cells, or increasing length, width or height.")
  }
  
  # Add Cell.Type and Cell.ID
  spatial_df$Cell.Type <- background_cell_type
  spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "random",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "minimum_distance_between_cells" = minimum_distance_between_cells,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(spatial_df)),
    colData = spatial_df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}

simulate_rings3D <- function(spe,
                             ring_properties_list,
                             plot_image = TRUE,
                             plot_cell_types = NULL,
                             plot_colours = NULL) {
  
  # Check shape variable of ring_properties_list
  shapes <- sapply(ring_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`ring_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(ring_properties_list))) {
    
    shape <- shapes[[i]]
    
    ### Sphere shape with ring
    if (shape == "sphere") {
      spe <- simulate_sphere_ring(spe, ring_properties_list[[i]])
    }
    
    ### Ellipsoid shape with ring
    else if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_ring(spe, ring_properties_list[[i]])
    }
    
    ### Cylinder shape with ring
    else if (shape == "cylinder") {
      spe <- simulate_cylinder_ring(spe, ring_properties_list[[i]])
    }
    
    ### Network shape with ring
    else if (shape == "network") {
      spe <- simulate_network_ring(spe, ring_properties_list[[i]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}

simulate_spe_metadata3D <- function(spe_metadata,
                                    plot_image = TRUE,
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  # First element should contain background metadata
  bg_metadata <- spe_metadata[[1]]
  if (!(is.character(bg_metadata$background_type) && length(bg_metadata$background_type) == 1)) {
    stop("background_type parameter found in the metadata background list is not a character.")
  }
  
  if (bg_metadata$background_type == "random") {
    spe <- simulate_random_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$minimum_distance_between_cells,
                                              plot_image = F)
  }
  else if (bg_metadata$background_type == "ordered") {
    spe <- simulate_ordered_background_cells3D(bg_metadata$n_cells,
                                               bg_metadata$length,
                                               bg_metadata$width,
                                               bg_metadata$height,
                                               bg_metadata$jitter_proportion,
                                               plot_image = F)
  }
  else {
    stop("background_type parameter found in the first list must be 'random' or 'ordered'.")
  }
  # Apply background mixing
  spe <- simulate_mixing3D(spe,
                           bg_metadata$cell_types,
                           bg_metadata$cell_proportions,
                           plot_image = F)
  
  ### If there is only background metadata, we are done
  if (length(spe_metadata) == 1) {
    
    # Plot
    if (plot_image) {
      fig <- plot_cells3D(spe,
                          plot_cell_types,
                          plot_colours)
      methods::show(fig)
    }
    
    return(spe)
  }
  
  ### All other elements should help to simulate clusters
  for (i in 2:length(spe_metadata)) {
    cluster_metadata <- spe_metadata[[i]]
    
    if (!(is.character(cluster_metadata$cluster_type) && length(cluster_metadata$cluster_type) == 1)) {
      stop(paste("cluster_type parameter found in the metadata cluster list", i,"is not a character."))
    }
    
    if (cluster_metadata$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(cluster_metadata), plot_image = F)
    }
    else if (cluster_metadata$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(cluster_metadata), plot_image = F)
    }
    else if (cluster_metadata$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(cluster_metadata), plot_image = F)
    }
    else {
      stop("cluster_type parameter must be either 'regular', 'ring' or 'double ring'.")
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}

simulate_sphere_cluster <- function(spe,
                                    cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  centre_loc <- cluster_properties$centre_loc
  
  # Change cell types in the sphere cluster
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}

simulate_sphere_dr <- function(spe,
                               dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  # Start with cells in outer ring
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + inner_ring_width + outer_ring_width)^2,
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in inner ring
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + inner_ring_width)^2,
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}

simulate_sphere_ring <- function(spe,
                                 ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  centre_loc <- ring_properties$centre_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  
  # Start with cells in ring
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + ring_width)^2,
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(spe)
}

spe_metadata_background_template <- function(background_type,
                                             original_spe_metadata = NULL) {
  
  if (background_type == "random") {
    background_metadata <- list(background = list(background_type = "random",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 600,
                                                  height = 300,
                                                  minimum_distance_between_cells = 10,
                                                  cell_types = c("Tumour", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "ordered") {
    background_metadata <- list(background = list(background_type = "ordered",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 300,
                                                  height = 300,
                                                  jitter_proportion = 0.25,
                                                  cell_types = c("Immune", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else {
    stop("background_type parameter must be 'random' or 'ordered'.")
  }
  
  
  # If original_spe_metadata input is not null, replace its background metadata with new background metadata
  if (!is.null(original_spe_metadata) && !is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[["background"]] <- background_metadata
    return(original_spe_metadata)
  }
  else if (!is.null(original_spe_metadata) && is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata <- c(background_metadata, original_spe_metadata)
    return(original_spe_metadata)
  }
  
  # Else, just return the background_metadata
  return(background_metadata)
}

spe_metadata_cluster_template <- function(cluster_type,
                                          shape,
                                          original_spe_metadata = NULL) {
  
  ### Get template for different shapes
  if (shape == "sphere") {
    cluster_metadata <- list(shape = "sphere",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radius = 100,
                             centre_loc = c(200, 150, 200))
  }
  else if (shape == "ellipsoid") {
    cluster_metadata <- list(shape = "ellipsoid",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radii = c(75, 100, 125),
                             centre_loc = c(450, 300, 100),
                             axes_rotation = c(0, 45, 0))
  }
  else if (shape == "cylinder") {
    cluster_metadata <- list(shape = "cylinder",
                             cluster_cell_types = c("Endothelial", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             radius = 40,
                             start_loc = c(400, 0, 0),
                             end_loc   = c(600, 400, 200))
  }
  else if (shape == "network") {
    cluster_metadata <- list(shape = "network",
                             cluster_cell_types = c("Immune", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             n_edges = 20,
                             width = 30,
                             centre_loc = c(200, 400, 150),
                             radius = 200)
  }
  else {
    stop("shape parameter must be 'sphere', 'ellipsoid', 'cylinder' or 'network'")
  }
  
  ### Add extra metadata for different cluster types
  if (cluster_type == "regular") {
    cluster_metadata <- append(list(cluster_type = "regular"), cluster_metadata)
  }
  else if (cluster_type == "ring") {
    cluster_metadata <- append(list(cluster_type = "ring"), cluster_metadata)
    cluster_metadata$ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$ring_width <- 12
  }
  else if (cluster_type == "double ring") {
    cluster_metadata <- append(list(cluster_type = "double ring"), cluster_metadata)
    cluster_metadata$inner_ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$inner_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$inner_ring_width <- 10
    cluster_metadata$outer_ring_cell_types <- c("Immune2", "Others")
    cluster_metadata$outer_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$outer_ring_width <- 10
  }
  else {
    stop("cluster_type parameter must be 'regular', 'ring' or 'double ring'")
  }
  
  # If original_spe_metadata input is not null, add new cluster_metadata to it
  if (!is.null(original_spe_metadata) && !is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[[paste("cluster", length(original_spe_metadata), sep="_")]] <- cluster_metadata
    return(original_spe_metadata)
  }
  else if (!is.null(original_spe_metadata) && is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[[paste("cluster", length(original_spe_metadata) + 1, sep="_")]] <- cluster_metadata
    return(original_spe_metadata)
  }
  
  # Else, just return the new cluster_metadata
  return(list("cluster_1" = cluster_metadata))
}

check_input_parameters <- function(input_parameters) {
  
  input_parameter_names <- names(input_parameters)
  
  check_value <- 0
  
  is_equal_with_tolerance <- function(x, y, tolerance = 1e-6) {
    abs(x - y) <= tolerance
  }
  
  for (input_parameter_name in input_parameter_names) {
    input_parameter <- input_parameters[[input_parameter_name]]
    
    # spe
    if (input_parameter_name == "spe" && class(input_parameter) != "SpatialExperiment") {
      check_value <- 1
      break
    }
    # Positive integer
    if (input_parameter_name %in% c("n_cells", "n_edges") && !(is.integer(input_parameter) && length(input_parameter) == 1 || (is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0 && input_parameter%%1 == 0))) {
      check_value <- 2
      break
    }
    # Positive numeric
    if (input_parameter_name %in% c("length", "width", "height", "radius", "x_radius", "y_radius", "z_radius", "ring_width", "inner_ring_width", "outer_ring_width") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0)) {
      check_value <- 3
      break
    }
    # Non-negative numeric
    if (input_parameter_name %in% c("minimum_distance_between_cells") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0)) {
      check_value <- 4
      break
    }
    # Numeric between 0 and 1
    if (input_parameter_name %in% c("jitter_proportion") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0 && input_parameter <= 1)) {
      check_value <- 5
      break
    }
    # Character
    if (input_parameter_name %in% c("background_cell_type") && !(is.character(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 6
      break
    }
    # Logical
    if (input_parameter_name %in% c("plot_image") && !(is.logical(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 7
      break
    }
    # Character vector
    if (input_parameter_name %in% c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types") &&
        !(is.character(input_parameter))) {
      check_value <- 8
      break
    }
    # Numeric vector
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") &&
        !(is.numeric(input_parameter))) {
      check_value <- 9
      break
    }
    # Numeric vector contains values between 0 and 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") &&
        sum(input_parameter < 0 | input_parameter > 1) != 0) {
      check_value <- 10
      break
    }
    # Numeric vector contains values that sum to 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") &&
        !is_equal_with_tolerance(sum(input_parameter), 1)) {
      check_value <- 11
      break
    }
    # Numeric vector of length 3
    if (input_parameter_name %in% c("centre_loc", "start_loc", "end_loc") && !(is.numeric(input_parameter) && length(input_parameter) == 3)) {
      check_value <- 12
      break
    }
    # Numeric
    if (input_parameter_name %in% c("y_z_rotation", "x_z_rotation", "x_y_rotation") && !(is.numeric(input_parameter) && length(input_parameter) == 1)) {
      check_value <- 13
      break
    }
  }
  
  # Two vectors match in length
  if (check_value == 0) {
    pairs <- data.frame(name1 = c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types"),
                        name2 = c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions"))
    
    for (i in seq(nrow(pairs))) {
      name1 <- pairs[["name1"]][i]
      name2 <- pairs[["name2"]][i]
      if (name1 %in% input_parameter_names && name2 %in% input_parameter_names && length(input_parameters[[name1]]) != length(input_parameters[[name2]])) {
        check_value <- 14
        input_parameter_name <- c(name1, name2)
        break
      }
    }
  }
  
  # If check_value equals 0, all inputs are valid.
  if (check_value == 0) {
    return(TRUE)
  }
  # At least one input is not valid, return the first invalid input.
  else {
    return(list(input_parameter_name = input_parameter_name, check_value = check_value))
  }
}

input_parameter_error_message <- function(input_parameter_check_value) {
  
  input_parameter_name <- input_parameter_check_value[[1]]
  check_value <- input_parameter_check_value[[2]]
  
  error_message <- switch(check_value,
                          "1" = paste(input_parameter_name, "is not a SpatialExperiment object."),
                          "2" = paste(input_parameter_name, "is not a positive integer."),
                          "3" = paste(input_parameter_name, "is not a positive numeric."),
                          "4" = paste(input_parameter_name, "is not a non-negative numeric."),
                          "5" = paste(input_parameter_name, "is not a numeric between 0 and 1."),
                          "6" = paste(input_parameter_name, "is not a character."),
                          "7" = paste(input_parameter_name, "is not a logical (TRUE or FALSE)."),
                          "8" = paste(input_parameter_name, "is not a character vector."),
                          "9" = paste(input_parameter_name, "is not a numeric vector."),
                          "10" = paste(input_parameter_name, "cannot be negative or greater than 1."),
                          "11" = paste(input_parameter_name, "does not sum to 1."),
                          "12" = paste(input_parameter_name, "is not a numeric vector of length 3."),
                          "13" = paste(input_parameter_name, "is not a numeric."),
                          "14" = paste(input_parameter_name[1], "and", input_parameter_name[2], "do not match in length."))
  
  return(error_message)
}

# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)
library(scales)
library(DescTools)
library(SpatialExperiment)
# Set seed for consistency ----
set.seed(678999821)


# NORMAL RANGES --------------------------------------------------- 
# ellipsoid_x_radii_range <- c(min = 75, max = 125)
# ellipsoid_y_radii_range <- c(min = 75, max = 125)
# ellipsoid_z_radii_range <- c(min = 75, max = 125)
# network_width_range <- c(min = 25, max = 35)
# 
# mixed_cell_type_A_proportion_range <- c(min = 0.3, max = 0.7)
# ringed_ring_width_factor_range <- c(min = 0.1, max = 0.2)
# separated_cluster1_x_coordinate_range <- c(min = 125, max = 175)


# OUTSIDE OF NORMAL RANGES ---------------------------------------- 
# ellipsoid_x_radii_range_small <- c(min = 50, max = 100)
# ellipsoid_x_radii_range_large <- c(min = 100, max = 150)
# network_width_range_small <- c(min = 20, max = 30)
# network_width_range_large <- c(min = 30, max = 40)
#
# mixed_cell_type_A_proportion_range_small <- c(min = 0.1, max = 0.5)
# mixed_cell_type_A_proportion_range_large <- c(min = 0.5, max = 0.9)
# ringed_ring_width_factor_range_small <- c(min = 0.05, max = 0.15)
# ringed_ring_width_factor_range_large <- c(min = 0.15, max = 0.25)
# separated_cluster1_x_coordinate_range_small <- c(min = 100, max = 150)
# separated_cluster1_x_coordinate_range_large <- c(min = 150, max = 200)


# OTHER PARAMETERS ------------------------------------------------ 

# BACKGROUND ***
# number_of_cells = 30000
# length = 600
# width = 600
# height = 300
# minimum_distance_between_cells = 10\

# SEPARATED SPHERE CLUSTER ***
# radius = 100
# centre_loc <- c(450, 300, 150)

# NETWORK CLUSTER ***
# edges = 20
# containment sphere radius = 125


# FUNCTIONS ----
# Code to perform simulated collection analysis (SD3 analysis)

# Set seed for consistency
set.seed(678999821)


# NORMAL RANGES --------------------------------------------------- #
# ellipsoid_x_radii_range <- c(min = 75, max = 125)
# ellipsoid_y_radii_range <- c(min = 75, max = 125)
# ellipsoid_z_radii_range <- c(min = 75, max = 125)
# network_width_range <- c(min = 25, max = 35)
# 
# mixed_cell_type_A_proportion_range <- c(min = 0.3, max = 0.7)
# ringed_ring_width_factor_range <- c(min = 0.1, max = 0.2)
# separated_cluster1_x_coordinate_range <- c(min = 125, max = 175)


# OUTSIDE OF NORMAL RANGES ---------------------------------------- #
# ellipsoid_x_radii_range_small <- c(min = 50, max = 100)
# ellipsoid_x_radii_range_large <- c(min = 100, max = 150)
# network_width_range_small <- c(min = 20, max = 30)
# network_width_range_large <- c(min = 30, max = 40)
#
# mixed_cell_type_A_proportion_range_small <- c(min = 0.1, max = 0.5)
# mixed_cell_type_A_proportion_range_large <- c(min = 0.5, max = 0.9)
# ringed_ring_width_factor_range_small <- c(min = 0.05, max = 0.15)
# ringed_ring_width_factor_range_large <- c(min = 0.15, max = 0.25)
# separated_cluster1_x_coordinate_range_small <- c(min = 100, max = 150)
# separated_cluster1_x_coordinate_range_large <- c(min = 150, max = 200)


# OTHER PARAMETERS ------------------------------------------------ #

# BACKGROUND ***
# number_of_cells = 30000
# length = 600
# width = 600
# height = 300
# minimum_distance_between_cells = 10\

# SEPARATED SPHERE CLUSTER ***
# radius = 100
# centre_loc <- c(450, 300, 150)

# NETWORK CLUSTER ***
# edges = 20
# containment sphere radius = 125


# FUNCTIONS ------------------------------------------------------- #
# Output is a 'parameters' data frame with 'n_simulations' rows, and columns 
# containing parameters specific to the chosen shape and arrangement.
# Values are randomly chosen from within the 'normal ranges' specified above.
create_random_parameters_for_simulated_group <- function(n_simulations, 
                                                         shape, 
                                                         arrangement) {
  
  # Create data frame for parameters
  parameters_df <- data.frame(number = seq_len(n_simulations))
  
  # Set random parameter ranges here
  ellipsoid_x_radii_range <- c(min = 75, max = 125)
  ellipsoid_y_radii_range <- c(min = 75, max = 125)
  ellipsoid_z_radii_range <- c(min = 75, max = 125)
  network_width_range <- c(min = 25, max = 35)
  
  mixed_cell_type_A_proportion_range <- c(min = 0.3, max = 0.7)
  ringed_ring_width_factor_range <- c(min = 0.1, max = 0.2)
  separated_cluster1_x_coordinate_range <- c(min = 125, max = 175)
  
  
  ### Shapes
  if (shape == "ellipsoid") {
    ellipsoid_x_radii <- runif(n_simulations, 
                               min = ellipsoid_x_radii_range["min"], 
                               max = ellipsoid_x_radii_range["max"])
    
    ellipsoid_y_radii <- runif(n_simulations, 
                               min = ellipsoid_y_radii_range["min"], 
                               max = ellipsoid_y_radii_range["max"])
    
    ellipsoid_z_radii <- runif(n_simulations, 
                               min = ellipsoid_z_radii_range["min"], 
                               max = ellipsoid_z_radii_range["max"])
    
    parameters_df$ellipsoid_x_radius <- ellipsoid_x_radii
    parameters_df$ellipsoid_y_radius <- ellipsoid_y_radii
    parameters_df$ellipsoid_z_radius <- ellipsoid_z_radii
    
  }
  else if (shape == "network") {
    network_widths <- runif(n_simulations,
                            min = network_width_range["min"],
                            max = network_width_range["max"])
    
    parameters_df$network_width <- network_widths
  }
  else {
    stop(paste(shape, "is not a valid shape", 
               "Should be 'ellipsoid', 'network'."))
  }
  
  
  
  ### Arrangements
  if (arrangement == "mixed") {
    mixed_cell_type_A_proportions <- runif(n_simulations, 
                                           min = mixed_cell_type_A_proportion_range["min"], 
                                           max = mixed_cell_type_A_proportion_range["max"])
    
    parameters_df$mixed_cell_type_A_proportion <- mixed_cell_type_A_proportions
    
  }
  else if (arrangement == "ringed") {
    ringed_ring_width_factors <- runif(n_simulations,
                                       min = ringed_ring_width_factor_range["min"],
                                       max = ringed_ring_width_factor_range["max"])
    
    # Get ring width as well, depending on if it is an ellipsoid or network shape
    if (shape == "ellipsoid") {
      ringed_ring_widths <- ringed_ring_width_factors * 
        (ellipsoid_x_radii + ellipsoid_y_radii + ellipsoid_z_radii) / 3
    }
    else if (shape == "network") {
      ringed_ring_widths <- ringed_ring_width_factors * network_widths
    }
    parameters_df$ringed_ring_width_factor <- ringed_ring_width_factors
    parameters_df$ringed_ring_width <- ringed_ring_widths
    
  }
  else if (arrangement == "separated") {
    separated_cluster1_x_coordinates <- runif(n_simulations,
                                              min = separated_cluster1_x_coordinate_range["min"],
                                              max = separated_cluster1_x_coordinate_range["max"])
    parameters_df$separated_cluster1_x_coordinate <- separated_cluster1_x_coordinates
  }
  else {
    stop(paste(arrangement, "is not a valid arrangement.", 
               "Should be 'mixed', 'ringed', 'separated'."))
  }  
  
  return(parameters_df)
}

# Output is a 'parameters' data frame
# The 'chosen_parameter' column is altered to include values randomly chosen from the 'range_for_chosen_parameter'
alter_chosen_parameter_for_simulated_group <- function(parameters, 
                                                       chosen_parameter, 
                                                       range_for_chosen_parameter) {
  
  if (is.null(parameters[[chosen_parameter]])) {
    stop(paste(chosen_parameter, "is not a parameter in 'parameters'."))
  }
  
  n_simulations <- nrow(parameters)
  
  parameters[[chosen_parameter]] <- runif(n_simulations, 
                                          min = range_for_chosen_parameter["min"], 
                                          max = range_for_chosen_parameter["max"])
  
  if (chosen_parameter == "ringed_ring_width_factor") {
    if ("network_width" %in% colnames(parameters)) {
      parameters[["ringed_ring_width"]] <- parameters[["ringed_ring_width_factor"]] * parameters[["network_width"]]
    }
    else {
      parameters[["ringed_ring_width"]] <- parameters[["ringed_ring_width_factor"]] * 
        (parameters[["ellipsoid_x_radius"]] + parameters[["ellipsoid_y_radius"]] + parameters[["ellipsoid_z_radius"]]) / 3
    }
  }
  
  return(parameters)
}

# Output is a list of simulation metadata
# Simulation metadata is made using 'parameters' data frame as input
create_simulation_metadata_for_simulated_group <- function(parameters, 
                                                           shape, 
                                                           arrangement) {
  
  # Define fixed background metadata parameters
  background_parameters <- list(
    number_of_cells = 30000,
    length = 600, # Units: micrometers (um)
    width = 600,
    height = 300,
    minimum_distance_between_cells = 10,
    cell_types = c("A", "B", "O"),
    cell_proportions = c(0, 0, 1)
  )
  
  spe_metadata_background <- spe_metadata_background_template("random", NULL)
  spe_metadata_background$background$length <- background_parameters$length
  spe_metadata_background$background$width <- background_parameters$width
  spe_metadata_background$background$height <- background_parameters$height
  spe_metadata_background$background$n_cells <- background_parameters$number_of_cells
  spe_metadata_background$background$minimum_distance_between_cells <- background_parameters$minimum_distance_between_cells
  spe_metadata_background$background$cell_types <- background_parameters$cell_types
  spe_metadata_background$background$cell_proportions <- background_parameters$cell_proportions
  
  # Determine cluster_type
  if (arrangement == "mixed") {
    cluster_type <- "regular"
  }
  else if (arrangement == "ringed") {
    cluster_type <- "ring"
  }
  else if (arrangement == "separated") {
    cluster_type <- "regular"
  }
  else {
    stop(paste(arrangement, "is not a valid arrangement.", 
               "Should be 'mixed', 'ringed', 'separated'."))
  }
  
  # Add cluster metadata
  spe_metadata_cluster <- spe_metadata_cluster_template(cluster_type, shape, spe_metadata_background)
  spe_metadata_cluster$cluster_1$cluster_cell_types <- c('A', 'B')
  spe_metadata_cluster$cluster_1$centre_loc <- c(300, 300, 150)
  
  # Add sphere if arrangement is separated
  if (arrangement == "separated") {
    spe_metadata_cluster <- spe_metadata_cluster_template(cluster_type, "sphere", spe_metadata_cluster)
    spe_metadata_cluster$cluster_2$cluster_cell_types <- c('A', 'B')
    spe_metadata_cluster$cluster_2$cluster_cell_proportions <- c(0, 1)
    spe_metadata_cluster$cluster_2$radius <- 100
    spe_metadata_cluster$cluster_2$centre_loc <- c(450, 300, 150)
    
  }
  
  simulation_metadata <- list()
  
  for (index in seq_len(nrow(parameters))) {
    simulation_metadata[[index]] <- spe_metadata_cluster
    
    if (shape == "ellipsoid") {
      simulation_metadata[[index]][["cluster_1"]][["x_radius"]] <- parameters$ellipsoid_x_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["y_radius"]] <- parameters$ellipsoid_y_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["z_radius"]] <- parameters$ellipsoid_z_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["y_z_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
      simulation_metadata[[index]][["cluster_1"]][["x_z_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
      simulation_metadata[[index]][["cluster_1"]][["x_y_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
    }
    else if (shape == "network") {
      simulation_metadata[[index]][["cluster_1"]][["n_edges"]] <- 20
      simulation_metadata[[index]][["cluster_1"]][["width"]] <- parameters$network_width[index]
      simulation_metadata[[index]][["cluster_1"]][["radius"]] <- 125
    }
    else {
      stop(paste(shape, "is not a valid shape", 
                 "Should be 'ellipsoid', 'network'."))
    }
    
    
    ### Arrangements
    if (arrangement == "mixed") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- 
        c(parameters$mixed_cell_type_A_proportion[index], 1 - parameters$mixed_cell_type_A_proportion[index])
      
    }
    else if (arrangement == "ringed") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- c(1, 0)
      simulation_metadata[[index]][["cluster_1"]][["ring_cell_types"]] <- c('A', 'B')
      simulation_metadata[[index]][["cluster_1"]][["ring_cell_proportions"]] <- c(0, 1)
      simulation_metadata[[index]][["cluster_1"]][["ring_width"]] <- parameters$ringed_ring_width[index]
    }
    else if (arrangement == "separated") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- c(1, 0)
      simulation_metadata[[index]][["cluster_1"]][["centre_loc"]] <- 
        c(parameters$separated_cluster1_x_coordinate[index], 300, 150)
    }
    else {
      stop(paste(arrangement, "is not a valid arrangement.", 
                 "Should be 'mixed', 'ringed', 'separated'."))
    }  
  }
  return(simulation_metadata)
}


# Input are two lists of simulation metadata for two groups
# Output is complete analysis of these two groups (i.e. one set)
analyse_simulated_set <- function(simulation_metadata_for_simulated_group1,
                                  simulation_metadata_for_simulated_group2) {
  
  # Check
  if (length(simulation_metadata_for_simulated_group1) != length(simulation_metadata_for_simulated_group2)) {
    stop("Number of simulations in set 1 and set 2 are different.")
  }
  
  # Set defined parameters/values
  cell_types <- c('A', 'B')
  n_simulations <- length(simulation_metadata_for_simulated_group1)
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  n_splits <- 10
  thresholds <- seq(0.01, 1, 0.01)
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  # z-coords for slices
  bottom_z_coord_of_slices <- c(85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205)
  top_z_coord_of_slices <- bottom_z_coord_of_slices + 10
  n_slices <- length(bottom_z_coord_of_slices)
  
  # Function to get a random slice from spe
  get_random_slice_from_spe <- function(spe, 
                                        bottom_z_coord_of_slices, 
                                        top_z_coord_of_slices) {
    
    slices_from_spe <- list()
    
    if (length(bottom_z_coord_of_slices) != length(top_z_coord_of_slices)) stop("Lengths of bottom_z_coords_of_slices and top_z_coords_of_slices should be equal.")
    number_of_slices <- length(bottom_z_coord_of_slices)
    
    i <- sample(length(bottom_z_coord_of_slices), 1) # random index to choose from slice coordinates
    bottom_z_coord_of_random_slice <- bottom_z_coord_of_slices[i]
    top_z_coord_of_random_slice <- top_z_coord_of_slices[i]
    z_coords_of_cells_in_spe <- spatialCoords(spe)[ , "Cell.Z.Position"]
    random_slice_from_spe <- spe[ , bottom_z_coord_of_random_slice < z_coords_of_cells_in_spe & z_coords_of_cells_in_spe < top_z_coord_of_random_slice]
    spatialCoords(random_slice_from_spe) <- spatialCoords(random_slice_from_spe)[ , c("Cell.X.Position", "Cell.Y.Position")]
    
    return(random_slice_from_spe)
  }
  
  # Function to create empty metric df list
  create_empty_metric_df_list <- function(
    cell_types,
    n_simulations,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("simulation", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, FCO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- FCO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
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
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           FCO = FCO_df,
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
                                   i, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    index <- n_cell_type_combinations * (i - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- i
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for gradient-based data which increments after each target cell type
    pair_index <- n_cell_type_combinations * (i - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics3D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           feature_colname = "Cell.Type")
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["FCO"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
        
        
        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["FCO"]][pair_index, radii_colnames] <- NA
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
          metric_df_list[["FCO"]][pair_index, radii_colnames] <- gradient_data[["fast_co_occurrence"]][[target_cell_type]]
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
        metric_df_list[["PBSAC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        
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
  analyse_simulation2D <- function(spe, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   i, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    index <- n_cell_type_combinations * (i - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- i
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for gradient-based data which increments after each target cell type
    pair_index <- n_cell_type_combinations * (i - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics2D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           feature_colname = "Cell.Type")
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["FCO"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
        
        
        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["FCO"]][pair_index, radii_colnames] <- NA
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
          metric_df_list[["FCO"]][pair_index, radii_colnames] <- gradient_data[["fast_co_occurrence"]][[target_cell_type]]
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
        metric_df_list[["PBSAC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        
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
    return(metric_df_list)
  }
  
  
  # Define metric df lists
  group1_metric_df_list3D <- group1_metric_df_list2D <- group2_metric_df_list3D <- group2_metric_df_list2D <-
    create_empty_metric_df_list(cell_types, n_simulations, radii_colnames, thresholds_colnames)
  
  for (i in seq(n_simulations)) {
    print(i)
    # Simulate spes in 3D, and get 2D slices
    group1_spe3D <- simulate_spe_metadata3D(simulation_metadata_for_simulated_group1[[i]], plot_image = F)
    group2_spe3D <- simulate_spe_metadata3D(simulation_metadata_for_simulated_group2[[i]], plot_image = F)
    
    group1_spe2D <- get_random_slice_from_spe(group1_spe3D, bottom_z_coord_of_slices, top_z_coord_of_slices)
    group2_spe2D <- get_random_slice_from_spe(group2_spe3D, bottom_z_coord_of_slices, top_z_coord_of_slices)
    
    # Analyse spes
    group1_metric_df_list3D <- analyse_simulation3D(group1_spe3D,
                                                    cell_types,
                                                    radii,
                                                    thresholds,
                                                    n_splits,
                                                    i,
                                                    group1_metric_df_list3D)
    
    group2_metric_df_list3D <- analyse_simulation3D(group2_spe3D,
                                                    cell_types,
                                                    radii,
                                                    thresholds,
                                                    n_splits,
                                                    i,
                                                    group2_metric_df_list3D)
    
    group1_metric_df_list2D <- analyse_simulation2D(group1_spe2D,
                                                    cell_types,
                                                    radii,
                                                    thresholds,
                                                    n_splits,
                                                    i,
                                                    group1_metric_df_list2D)
    
    group2_metric_df_list2D <- analyse_simulation2D(group2_spe2D,
                                                    cell_types,
                                                    radii,
                                                    thresholds,
                                                    n_splits,
                                                    i,
                                                    group2_metric_df_list2D)
  }
  
  simulated_set_analysis <- list(
    group1_metric_df_list3D = group1_metric_df_list3D,
    group2_metric_df_list3D = group2_metric_df_list3D,
    group1_metric_df_list2D = group1_metric_df_list2D,
    group2_metric_df_list2D = group2_metric_df_list2D
  )
  
  ## Turn gradient metrics into AUC and add to metric_df list
  get_AUC_for_radii_gradient_metrics <- function(y) {
    x <- radii
    h <- diff(x)[1]
    n <- length(x)
    
    AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
    
    return(AUC)
  }
  
  gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "FCO", "CK", "CL", "CG")
  
  for (simulated_group_analysis in names(simulated_set_analysis)) {
    metric_df_list <- simulated_set_analysis[[simulated_group_analysis]]
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "FCO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "reference", "target", metric_AUC_name)
      }
      else {
        subset_colnames <- c("simulation", "reference", metric_AUC_name)
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
    
    simulated_set_analysis[[simulated_group_analysis]] <- metric_df_list
  }
  
  return(simulated_set_analysis)
}


# Output is p-values obtained from analysis of sets
# group1_chosen_parameters is a list of the form: list("mixed_cell_type_A_proportion" = c(min = 0.1, max = 0.5), ...)
analyse_simulated_collection <- function(n_simulations_in_group,
                                         n_sets,
                                         group1_shape,
                                         group1_arrangement,
                                         group1_chosen_parameters = NULL,
                                         group2_shape,
                                         group2_arrangement,
                                         group2_chosen_parameters)  {
  
  # Define constants
  metrics <- c("AMD",
               "MS_AUC", "NMS_AUC",
               "ACIN_AUC", "ANE_AUC", "ANC_AUC", "FCO_AUC", "CK_AUC", "CL_AUC", "CG_AUC",
               "PBSAC", "EBSAC", "PBP_AUC", "EBP_AUC")
  
  cell_types <- c("A", "B")
  
  # Define result
  simulated_collection_p_values <- data.frame()
  
  # Subset metric_df
  subset_metric_df <- function(metric_df, reference_cell_type, target_cell_type, metric) {
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      metric_df <- metric_df[
        metric_df[["cell_types"]] == paste(reference_cell_type, target_cell_type, sep = ","), 
      ]
    }
    else if (metric %in% c("ANE_AUC")) {
      metric_df <- metric_df[
        metric_df[["target"]] == paste(reference_cell_type, target_cell_type, sep = ","), 
      ]
    }
    else {
      metric_df <- metric_df[
        metric_df[["reference"]] == reference_cell_type & metric_df[["target"]] == target_cell_type, 
      ]
    }
    return(metric_df)
  }
  
  # Calculate p-values from a single set
  calculate_p_values_from_simulated_set <- function(simulated_set_analysis,
                                                    metrics,
                                                    cell_types) {
    
    multiple_testing_correction_factor <- length(metrics) * n_sets
    
    group1_metric_df_list3D <- simulated_set_analysis[["group1_metric_df_list3D"]]
    group2_metric_df_list3D <- simulated_set_analysis[["group2_metric_df_list3D"]]
    group1_metric_df_list2D <- simulated_set_analysis[["group1_metric_df_list2D"]]
    group2_metric_df_list2D <- simulated_set_analysis[["group2_metric_df_list2D"]]
    
    combined_p_value_df <- data.frame()
    
    for (metric in metrics) {
      metric_p_value_df <- data.frame()
      
      group1_metric_df3D <- group1_metric_df_list3D[[metric]]
      group2_metric_df3D <- group2_metric_df_list3D[[metric]]
      group1_metric_df2D <- group1_metric_df_list2D[[metric]]
      group2_metric_df2D <- group2_metric_df_list2D[[metric]]
      
      for (reference_cell_type in cell_types) {
        for (target_cell_type in cell_types) {
          
          group1_metric_df3D_subset <- subset_metric_df(group1_metric_df3D, reference_cell_type, target_cell_type, metric)
          group2_metric_df3D_subset <- subset_metric_df(group2_metric_df3D, reference_cell_type, target_cell_type, metric)
          group1_metric_df2D_subset <- subset_metric_df(group1_metric_df2D, reference_cell_type, target_cell_type, metric)
          group2_metric_df2D_subset <- subset_metric_df(group2_metric_df2D, reference_cell_type, target_cell_type, metric)
          
          p_value3D <- wilcox.test(group1_metric_df3D_subset[[metric]], group2_metric_df3D_subset[[metric]])$p.value
          p_value2D <- wilcox.test(group1_metric_df2D_subset[[metric]], group2_metric_df2D_subset[[metric]])$p.value
          
          # Multiple testing correction
          p_value3D <- p_value3D * multiple_testing_correction_factor
          p_value2D <- p_value2D * multiple_testing_correction_factor
          
          # Add rows for 3D and 2D results
          metric_p_value_df <- rbind(
            metric_p_value_df,
            data.frame(
              reference = reference_cell_type,
              target    = target_cell_type,
              metric    = metric,
              dimension = "3D",
              p_value   = p_value3D
            ),
            data.frame(
              reference = reference_cell_type,
              target    = target_cell_type,
              metric    = metric,
              dimension = "2D",
              p_value   = p_value2D
            )
          )
        }
      }
      combined_p_value_df <- rbind(combined_p_value_df, metric_p_value_df)
    }
    return(combined_p_value_df)
  }
  
  
  # Iterate for each set
  for (set_index in seq_len(n_sets)) {
    print(paste("Set:", set_index))
    group1_parameters <- create_random_parameters_for_simulated_group(n_simulations = n_simulations_in_group,
                                                                      shape = group1_shape,
                                                                      arrangement = group1_arrangement)
    
    group2_parameters <- create_random_parameters_for_simulated_group(n_simulations = n_simulations_in_group,
                                                                      shape = group2_shape,
                                                                      arrangement = group2_arrangement)
    
    if (!is.null(group1_chosen_parameters)) {
      for (i in seq_len(length(group1_chosen_parameters))) {
        group1_parameters <- alter_chosen_parameter_for_simulated_group(parameters = group1_parameters,
                                                                        chosen_parameter = names(group1_chosen_parameters)[i],
                                                                        range_for_chosen_parameter = group1_chosen_parameters[i])
      }
    }
    
    if (!is.null(group2_chosen_parameters)) {
      for (i in seq_len(length(group2_chosen_parameters))) {
        group2_parameters <- alter_chosen_parameter_for_simulated_group(parameters = group2_parameters,
                                                                        chosen_parameter = names(group2_chosen_parameters)[i],
                                                                        range_for_chosen_parameter = group2_chosen_parameters[[i]])
      }
    }
    
    group1_simulation_metadata <- create_simulation_metadata_for_simulated_group(parameters = group1_parameters,
                                                                                 shape = group1_shape,
                                                                                 arrangement = group1_arrangement)
    
    group2_simulation_metadata <- create_simulation_metadata_for_simulated_group(parameters = group2_parameters,
                                                                                 shape = group2_shape,
                                                                                 arrangement = group2_arrangement)  
    
    simulated_set_analysis <- analyse_simulated_set(group1_simulation_metadata,
                                                    group2_simulation_metadata)
    
    simulated_set_p_values <- calculate_p_values_from_simulated_set(simulated_set_analysis,
                                                                    metrics,
                                                                    cell_types)
    
    simulated_set_p_values$set <- set_index
    
    simulated_collection_p_values <- rbind(simulated_collection_p_values, simulated_set_p_values)
    
  }
  return(simulated_collection_p_values)
}



# MN vs MN ----
MN_vs_MN_collection <- analyse_simulated_collection(n_simulations_in_group = 20,
                                                    n_sets = 200,
                                                    group1_shape = "network",
                                                    group1_arrangement = "mixed",
                                                    group1_chosen_parameters = NULL,
                                                    group2_shape = "network",
                                                    group2_arrangement = "mixed",
                                                    group2_chosen_parameters = NULL)
saveRDS(MN_vs_MN_collection, "~/R/reproduce_all_code_data/MN_vs_MN_collection.RDS")



print("DONE")


