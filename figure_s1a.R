### SPIAT3D ----
calculate_all_gradient_cc_metrics3D <- function(spatial_df,
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
    df <- calculate_all_single_radius_cc_metrics3D(spatial_df,
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

calculate_all_single_radius_cc_metrics3D <- function(spatial_df,
                                                     reference_cell_type,
                                                     target_cell_types,
                                                     radius,
                                                     feature_colname) {
  
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  spatial_df_coords <- spatial_df
  length <- (max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- (max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  height <- (max(spatial_df_coords$Cell.Z.Position) - min(spatial_df_coords$Cell.Z.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # All single radius cc metrics stem from calculate_neighbourhood_entropy3D function
  neighbourhood_entropy_df <- calculate_neighbourhood_entropy3D(spatial_df,
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
      mixing_score_df$n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
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
    cross_K_df[[target_cell_type]] <- (((volume * sum(neighbourhood_entropy_df[[target_cell_type]])) / sum(spatial_df[[feature_colname]] == reference_cell_type)) / sum(spatial_df[[feature_colname]] == target_cell_type))
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
    n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / volume
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  }
  
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spatial_df[[feature_colname]])
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spatial_df,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  n_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")],
                                                  eps = radius,
                                                  query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")],
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")][spatial_df[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")],
                                                           sort = FALSE)
    n_target_cells_in_all_cell_type_radius <- sum(rapply(n_target_cells_in_all_cell_type_radius$id, base::length))
    
    target_cell_type_proportion_in_all_cell_type_radius <- n_target_cells_in_all_cell_type_radius / n_cells_in_all_cell_type_radius
    
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_all_cell_type_radius
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}

calculate_cell_proportion_grid_metrics3D <- function(spatial_df,
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname,
                                                     plot_image = TRUE) {
  
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(reference_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  spatial_df_coords <- spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get dimensions of the window
  min_x <- min(spatial_df_coords[ , "Cell.X.Position"])
  min_y <- min(spatial_df_coords[ , "Cell.Y.Position"])
  min_z <- min(spatial_df_coords[ , "Cell.Z.Position"])
  
  max_x <- max(spatial_df_coords[ , "Cell.X.Position"])
  max_y <- max(spatial_df_coords[ , "Cell.Y.Position"])
  max_z <- max(spatial_df_coords[ , "Cell.Z.Position"])
  
  length <- (max_x - min_x)
  width  <- (max_y - min_y)
  height <- (max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  # Shift spatial_df_coords so they begin at the origin
  spatial_df_coords[, "Cell.X.Position"] <- spatial_df_coords[, "Cell.X.Position"] - min_x
  spatial_df_coords[, "Cell.Y.Position"] <- spatial_df_coords[, "Cell.Y.Position"] - min_y
  spatial_df_coords[, "Cell.Z.Position"] <- spatial_df_coords[, "Cell.Z.Position"] - min_z
  
  ## Figure out which 'grid prism number' each cell is inside
  spatial_df$grid_prism_num <- floor(spatial_df_coords[ , "Cell.X.Position"] / d_row) +
    floor(spatial_df_coords[ , "Cell.Y.Position"] / d_col) * n_splits +
    floor(spatial_df_coords[ , "Cell.Z.Position"] / d_lay) * n_splits^2 + 1
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spatial_df[[feature_colname]], factor(spatial_df$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + (min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + (min_y),
                                       z_coord = (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay + (min_z))
  
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
  result <- cbind(result, grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}

calculate_cell_proportions3D <- function(spatial_df,
                                         cell_types_of_interest = NULL,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (nrow(spatial_df) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spatial_df[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
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

calculate_cells_in_neighbourhood_gradient3D <- function(spatial_df,
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
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood3D(spatial_df,
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

calculate_cells_in_neighbourhood3D <- function(spatial_df,
                                               reference_cell_type,
                                               target_cell_types,
                                               radius,
                                               feature_colname) {
  
  ## Get cells in neighbourhood df
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spatial_df,
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

calculate_co_occurrence_gradient3D <- function(spatial_df,
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
    co_occurrence_df <- calculate_co_occurrence3D(spatial_df,
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

calculate_co_occurrence3D <- function(spatial_df,
                                      reference_cell_type,
                                      target_cell_types,
                                      radius,
                                      feature_colname) {
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get all cell types in spatial_df
  all_cell_types <- unique(spatial_df[[feature_colname]])
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spatial_df,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  # Get total number of cells in radius around all cell types
  n_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")],
                                                  eps = radius,
                                                  query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")],
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get total number of target cells in radius around all cell types
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")][spatial_df[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")],
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

calculate_cross_G_gradient3D <- function(spatial_df,
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
    cross_G_df <- calculate_cross_G3D(spatial_df,
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

calculate_cross_G3D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname) {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spatial_df,
                                                              reference_cell_type,
                                                              target_cell_type,
                                                              radius,
                                                              feature_colname)
  
  reference_target_interactions <- neighbourhood_counts_df[[target_cell_type]]
  
  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  
  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")])
  
  length <- (max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- (max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  height <- (max(spatial_df_coords$Cell.Z.Position) - min(spatial_df_coords$Cell.Z.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # Get the number of target cells
  n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}

calculate_cross_K_gradient3D <- function(spatial_df,
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
    cross_K_df <- calculate_cross_K3D(spatial_df,
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

calculate_cross_K3D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Check if inputs are valid
  if (is.null(spatial_df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spatial_df object"))
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  
  ## Get expected cross K-function
  expected_cross_K <- (4/3) * pi * radius^3
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")])
  
  length <- (max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- (max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  height <- (max(spatial_df_coords$Cell.Z.Position) - min(spatial_df_coords$Cell.Z.Position))
  
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spatial_df,
                                                              reference_cell_type,
                                                              target_cell_types,
                                                              radius,
                                                              feature_colname)
  
  # Calculate cross K-fucnction for each target cell type
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    
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

calculate_cross_L_gradient3D <- function(spatial_df,
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
    cross_L_df <- calculate_cross_L3D(spatial_df,
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

calculate_cross_L3D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Calculate cross K first
  result <- calculate_cross_K3D(spatial_df = spatial_df,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  # Linearise cross K to get cross L
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  
  return(result)
}

calculate_entropy_background3D <- function(spatial_df,
                                           cell_types_of_interest,
                                           feature_colname) {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions3D(spatial_df, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy)
}

calculate_entropy_grid_metrics3D <- function(spatial_df,
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname,
                                             plot_image = TRUE) {
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spatial_df[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  spatial_df_coords <- spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get dimensions of the window
  min_x <- min(spatial_df_coords[ , "Cell.X.Position"])
  min_y <- min(spatial_df_coords[ , "Cell.Y.Position"])
  min_z <- min(spatial_df_coords[ , "Cell.Z.Position"])
  
  max_x <- max(spatial_df_coords[ , "Cell.X.Position"])
  max_y <- max(spatial_df_coords[ , "Cell.Y.Position"])
  max_z <- max(spatial_df_coords[ , "Cell.Z.Position"])
  
  length <- (max_x - min_x)
  width  <- (max_y - min_y)
  height <- (max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  # Shift spatial_df_coords so they begin at the origin
  spatial_df_coords[, "Cell.X.Position"] <- spatial_df_coords[, "Cell.X.Position"] - min_x
  spatial_df_coords[, "Cell.Y.Position"] <- spatial_df_coords[, "Cell.Y.Position"] - min_y
  spatial_df_coords[, "Cell.Z.Position"] <- spatial_df_coords[, "Cell.Z.Position"] - min_z
  
  ## Figure out which 'grid prism number' each cell is inside
  spatial_df$grid_prism_num <- floor(spatial_df_coords[ , "Cell.X.Position"] / d_row) +
    floor(spatial_df_coords[ , "Cell.Y.Position"] / d_col) * n_splits +
    floor(spatial_df_coords[ , "Cell.Z.Position"] / d_lay) * n_splits^2 + 1
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spatial_df[[feature_colname]], factor(spatial_df$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + (min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + (min_y),
                                       z_coord = (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay + (min_z))
  
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
  result <- cbind(result, grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}

calculate_minimum_distances_between_cell_types3D <- function(spatial_df,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname,
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (nrow(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Get spatial_df coords
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")])
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spatial_df[[feature_colname]] == cell_type1) == 0 || sum(spatial_df[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type2, ]
    
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

calculate_mixing_scores_gradient3D <- function(spatial_df,
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
    mixing_scores <- calculate_mixing_scores3D(spatial_df,
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

calculate_mixing_scores3D <- function(spatial_df,
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
      n_ref <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spatial_df[[feature_colname]] == target_cell_type)
      
      
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
      neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spatial_df,
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

calculate_neighbourhood_counts_gradient3D <- function(spatial_df,
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
    neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spatial_df,
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

calculate_neighbourhood_counts3D <- function(spatial_df,
                                             reference_cell_type,
                                             target_cell_types,
                                             radius,
                                             feature_colname) {
  
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # Get spatial_df coords
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")])
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(ref_cell_id = spatial_df$Cell.ID[spatial_df[[feature_colname]] == reference_cell_type])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spatial_df[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == target_cell_type, ]
    
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

calculate_neighbourhood_entropy_gradient3D <- function(spatial_df,
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
    neighbourhood_entropy_df <- calculate_neighbourhood_entropy3D(spatial_df,
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

calculate_neighbourhood_entropy3D <- function(spatial_df,
                                              reference_cell_type,
                                              target_cell_types,
                                              radius,
                                              feature_colname) {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
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

calculate_pairwise_distances_between_cell_types3D <- function(spatial_df,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname,
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (nrow(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatial_df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")])
  rownames(distance_matrix) <- spatial_df$Cell.ID
  colnames(distance_matrix) <- spatial_df$Cell.ID
  
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



### SPIAT2D ----
calculate_all_gradient_cc_metrics2D <- function(spatial_df,
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
    df <- calculate_all_single_radius_cc_metrics2D(spatial_df,
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

calculate_all_single_radius_cc_metrics2D <- function(spatial_df,
                                                     reference_cell_type,
                                                     target_cell_types,
                                                     radius,
                                                     feature_colname) {
  
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  spatial_df_coords <- spatial_df
  length <- (max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- (max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  
  # Get area of the window the cells are in
  area <- length * width
  
  # All single radius cc metrics stem from calculate_neighbourhood_entropy2D function
  neighbourhood_entropy_df <- calculate_neighbourhood_entropy2D(spatial_df,
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
      mixing_score_df$n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
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
    cross_K_df[[target_cell_type]] <- (((area * sum(neighbourhood_entropy_df[[target_cell_type]])) / sum(spatial_df[[feature_colname]] == reference_cell_type)) / sum(spatial_df[[feature_colname]] == target_cell_type))
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
    n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / area
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  }
  
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spatial_df[[feature_colname]])
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spatial_df,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  n_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")],
                                                  eps = radius,
                                                  query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position")],
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")][spatial_df[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position")],
                                                           sort = FALSE)
    n_target_cells_in_all_cell_type_radius <- sum(rapply(n_target_cells_in_all_cell_type_radius$id, base::length))
    
    target_cell_type_proportion_in_all_cell_type_radius <- n_target_cells_in_all_cell_type_radius / n_cells_in_all_cell_type_radius
    
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_all_cell_type_radius
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}

calculate_cell_proportion_grid_metrics2D <- function(spatial_df,
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname,
                                                     plot_image = TRUE) {
  
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(reference_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  spatial_df_coords <- spatial_df[, c("Cell.X.Position", "Cell.Y.Position")]
  
  ## Get dimensions of the window
  min_x <- min(spatial_df_coords[ , "Cell.X.Position"])
  min_y <- min(spatial_df_coords[ , "Cell.Y.Position"])
  
  max_x <- max(spatial_df_coords[ , "Cell.X.Position"])
  max_y <- max(spatial_df_coords[ , "Cell.Y.Position"])
  
  length <- (max_x - min_x)
  width  <- (max_y - min_y)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  # Shift spatial_df_coords so they begin at the origin
  spatial_df_coords[, "Cell.X.Position"] <- spatial_df_coords[, "Cell.X.Position"] - min_x
  spatial_df_coords[, "Cell.Y.Position"] <- spatial_df_coords[, "Cell.Y.Position"] - min_y
  
  ## Figure out which 'grid prism number' each cell is inside
  spatial_df$grid_prism_num <- floor(spatial_df_coords[ , "Cell.X.Position"] / d_row) +
    floor(spatial_df_coords[ , "Cell.Y.Position"] / d_col) * n_splits
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^2
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spatial_df[[feature_colname]], factor(spatial_df$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + (min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + (min_y))
  
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
  result <- cbind(result, grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}

calculate_cell_proportions2D <- function(spatial_df,
                                         cell_types_of_interest = NULL,
                                         feature_colname,
                                         plot_image = TRUE) {
  
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (nrow(spatial_df) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spatial_df[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
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

calculate_cells_in_neighbourhood_gradient2D <- function(spatial_df,
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
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood2D(spatial_df,
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

calculate_cells_in_neighbourhood2D <- function(spatial_df,
                                               reference_cell_type,
                                               target_cell_types,
                                               radius,
                                               feature_colname) {
  
  ## Get cells in neighbourhood df
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spatial_df,
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

calculate_co_occurrence_gradient2D <- function(spatial_df,
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
    co_occurrence_df <- calculate_co_occurrence2D(spatial_df,
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

calculate_co_occurrence2D <- function(spatial_df,
                                      reference_cell_type,
                                      target_cell_types,
                                      radius,
                                      feature_colname) {
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get all cell types in spatial_df
  all_cell_types <- unique(spatial_df[[feature_colname]])
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spatial_df,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  # Get total number of cells in radius around all cell types
  n_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")],
                                                  eps = radius,
                                                  query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position")],
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get total number of target cells in radius around all cell types
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")][spatial_df[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = spatial_df[, c("Cell.X.Position", "Cell.Y.Position")],
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

calculate_cross_G_gradient2D <- function(spatial_df,
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
    cross_G_df <- calculate_cross_G2D(spatial_df,
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

calculate_cross_G2D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname) {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spatial_df,
                                                              reference_cell_type,
                                                              target_cell_type,
                                                              radius,
                                                              feature_colname)
  
  reference_target_interactions <- neighbourhood_counts_df[[target_cell_type]]
  
  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  
  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")])
  
  length <- (max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- (max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  
  # Get area of the window the cells are in
  area <- length * width
  
  # Get the number of target cells
  n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / area
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}

calculate_cross_K_gradient2D <- function(spatial_df,
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
    cross_K_df <- calculate_cross_K2D(spatial_df,
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

calculate_cross_K2D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Check if inputs are valid
  if (is.null(spatial_df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spatial_df object"))
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  
  ## Get expected cross K-function
  expected_cross_K <- pi * radius^2
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")])
  
  length <- (max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- (max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  
  ## Get area of the window the cells are in
  area <- length * width
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spatial_df,
                                                              reference_cell_type,
                                                              target_cell_types,
                                                              radius,
                                                              feature_colname)
  
  # Calculate cross K-fucnction for each target cell type
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    
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

calculate_cross_L_gradient2D <- function(spatial_df,
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
    cross_L_df <- calculate_cross_L2D(spatial_df,
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

calculate_cross_L2D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {
  
  # Calculate cross K first
  result <- calculate_cross_K2D(spatial_df = spatial_df,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  # Linearise cross K to get cross L
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  
  return(result)
}

calculate_entropy_background2D <- function(spatial_df,
                                           cell_types_of_interest,
                                           feature_colname) {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions2D(spatial_df, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy)
}

calculate_entropy_grid_metrics2D <- function(spatial_df,
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname,
                                             plot_image = TRUE) {
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spatial_df[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  spatial_df_coords <- spatial_df[, c("Cell.X.Position", "Cell.Y.Position")]
  
  ## Get dimensions of the window
  min_x <- min(spatial_df_coords[ , "Cell.X.Position"])
  min_y <- min(spatial_df_coords[ , "Cell.Y.Position"])
  
  max_x <- max(spatial_df_coords[ , "Cell.X.Position"])
  max_y <- max(spatial_df_coords[ , "Cell.Y.Position"])
  
  length <- (max_x - min_x)
  width  <- (max_y - min_y)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  # Shift spatial_df_coords so they begin at the origin
  spatial_df_coords[, "Cell.X.Position"] <- spatial_df_coords[, "Cell.X.Position"] - min_x
  spatial_df_coords[, "Cell.Y.Position"] <- spatial_df_coords[, "Cell.Y.Position"] - min_y
  
  ## Figure out which 'grid prism number' each cell is inside
  spatial_df$grid_prism_num <- floor(spatial_df_coords[ , "Cell.X.Position"] / d_row) +
    floor(spatial_df_coords[ , "Cell.Y.Position"] / d_col) * n_splits
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^2
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spatial_df[[feature_colname]], factor(spatial_df$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + (min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + (min_y))
  
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
  result <- cbind(result, grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}

calculate_minimum_distances_between_cell_types2D <- function(spatial_df,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname,
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (nrow(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Get spatial_df coords
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")])
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spatial_df[[feature_colname]] == cell_type1) == 0 || sum(spatial_df[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type2, ]
    
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

calculate_mixing_scores_gradient2D <- function(spatial_df,
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
    mixing_scores <- calculate_mixing_scores2D(spatial_df,
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

calculate_mixing_scores2D <- function(spatial_df,
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
      n_ref <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spatial_df[[feature_colname]] == target_cell_type)
      
      
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
      neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spatial_df,
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

calculate_neighbourhood_counts_gradient2D <- function(spatial_df,
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
    neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spatial_df,
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

calculate_neighbourhood_counts2D <- function(spatial_df,
                                             reference_cell_type,
                                             target_cell_types,
                                             radius,
                                             feature_colname) {
  
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # Get spatial_df coords
  spatial_df_coords <- data.frame(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")])
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(ref_cell_id = spatial_df$Cell.ID[spatial_df[[feature_colname]] == reference_cell_type])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spatial_df[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == target_cell_type, ]
    
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

calculate_neighbourhood_entropy_gradient2D <- function(spatial_df,
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
    neighbourhood_entropy_df <- calculate_neighbourhood_entropy2D(spatial_df,
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

calculate_neighbourhood_entropy2D <- function(spatial_df,
                                              reference_cell_type,
                                              target_cell_types,
                                              radius,
                                              feature_colname) {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
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

calculate_pairwise_distances_between_cell_types2D <- function(spatial_df,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname,
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spatial_df[[feature_colname]][trimws(spatial_df[[feature_colname]]) == ""]) > 0) {
    stop("spatial_df cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (nrow(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatial_df[, c("Cell.X.Position", "Cell.Y.Position")])
  rownames(distance_matrix) <- spatial_df$Cell.ID
  colnames(distance_matrix) <- spatial_df$Cell.ID
  
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


### Libraries -----
library(SpatialExperiment)
library(dbscan)
library(apcluster)
library(plotly)
library(dplyr)
library(reshape2)
library(gtools)
library(cowplot)
library(Hmisc)

### cycif human colorectal cancer dataset -----
setwd("~/R/data3D/CyCIF_colorectal_cancer")
data3D <- read.csv("colorectal_cancer_df.csv")
data3D <- data3D[, -1]
colnames(data3D)[1:3] <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")
colnames(data3D)[4:5] <- c("Cell.Type.Specific", "Cell.Type")
data3D$Cell.Type[data3D$Cell.Type == "Tumor/Epi"] <- "Tumour"

cell_types <- unique(data3D$Cell.Type)

### Code to analyse a public 3D dataset made up of multiple slices. -----

# ******** READ THIS *********
# This function assumes the following:
# data3D is a data frame with column names: "Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position", ""Cell.Type"
# "Cell.Z.Position" contains discrete values.
# cell_types consists of only cells in the "Cell.Type" column, and at least two cells.

analyse_3D_data_with_2D_slices <- function(
    data3D,
    cell_types,
    radii = seq(20, 100, 10),
    n_splits = 10,
    thresholds = seq(0.01, 1, 0.01)
) {
  
  slice_z_coords <- unique(data3D$Cell.Z.Position)
  n_slices <- length(unique(data3D$Cell.Z.Position))
  n_cell_type_combinations <- length(cell_types)^2
  radii_colnames <- paste("r", radii, sep = "")
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  create_empty_metric_df_list <- function(
    cell_types,
    n_slices,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("slice", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("slice", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("slice", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("slice", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("slice", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("slice", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
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
  
  metric_df_list <- create_empty_metric_df_list(
    cell_types,
    n_slices,
    radii_colnames,
    thresholds_colnames
  )
  
  for (i in seq(n_slices + 1, 1)) {
    print(i)
    
    # i represents the current slice index
    if (i == n_slices + 1) {
      df <- data3D
    }
    # if i is the last index, analyse in 3D instead
    else {
      df <- data3D[data3D$Cell.Z.Position == slice_z_coords[i], ]
    }
    
    ### 3D analysis -----------------------------
    if (i == n_slices + 1) {
      
      index <- n_cell_type_combinations * (i - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(df,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- i
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for gradient-based data which increments after each target cell type
      pair_index <- n_cell_type_combinations * (i - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics3D(df,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             plot_image = F,
                                                             feature_colname = "Cell.Type")
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ANC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ACIN"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CK"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CL"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CG"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ANE"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, 
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
          metric_df_list[["PBSAC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(df, 
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
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics3D(df, 
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
    }
    ### 2D analysis -----------------------------
    else {
      
      index <- n_cell_type_combinations * (i - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(df,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- i
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for gradient-based data which increments after each target cell type
      pair_index <- n_cell_type_combinations * (i - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics2D(df,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             plot_image = F,
                                                             feature_colname = "Cell.Type")
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ANC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ACIN"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CK"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CL"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CG"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ANE"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, 
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
          metric_df_list[["PBSAC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(df, 
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
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics2D(df, 
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
  }
  return(metric_df_list)
}


human_colorectal_cancer_metric_df_list <- 
  analyse_3D_data_with_2D_slices(
    data3D = data3D,
    cell_types = cell_types,
    radii = seq(20, 100, 10),
    n_splits = 10,
    thresholds = seq(0.01, 1, 0.01)
  )


saveRDS(human_colorectal_cancer_metric_df_list, "~/R/reproduce_all_code_data/human_colorectal_cancer_metric_df_list.RDS")

# Code to plot public 3D dataset analysis.

### Read data and set file name and save directory ------
metric_df_list <- readRDS("~/R/reproduce_all_code_data/human_colorectal_cancer_metric_df_list.RDS") # e.g. "~/R/public_data_analysis/CyCIF_metric_df_list.RDS"

file_name_prefix <- "human_colorectal_cancer" # e.g. "colorectal_cancer"
save_directory <- "~/R/reproduce_all_code_plots" # e.g. "~/R/plots/public_data/colorectal_cancer"

### Format data -----

get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "CK", "CL", "CG")) {
    return("radius")
  }
  else if (metric %in% c("PBP", "EBP")) {
    return("threshold")  
  }
  else {
    stop("Invalid metric. Must be gradient-based")
  }
}


## Turn gradient radii metrics into AUC and add to metric_df list
get_AUC_for_radii_gradient_metrics <- function(y) {
  x <- radii
  h <- diff(x)[1]
  n <- length(x)
  
  AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
  
  return(AUC)
}


radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "CKR", "CLR", "COO", "CGR", "CK", "CL", "CG")


for (metric in gradient_radii_metrics) {
  metric_AUC_name <- paste(metric, "AUC", sep = "_")
  
  if (metric %in% c("MS", "NMS", "ANC", "CKR", "CLR", "COO", "CGR", "CK", "CL", "CG")) {
    subset_colnames <- c("slice", "reference", "target", metric_AUC_name)
  }
  else {
    subset_colnames <- c("slice", "reference", metric_AUC_name)
  }
  
  df <- metric_df_list[[metric]]
  df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
  metric_df_list[[metric_AUC_name]] <- df
  
}

## Turn threshold radii metrics into AUC and add to metric_df list
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

# PBP_AUC 3D
PBP_df <- metric_df_list[["PBP"]]
PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
  sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
})
PBP_AUC_df <- PBP_df[ , c("slice", "reference", "target", "PBP_AUC")]
metric_df_list[["PBP_AUC"]] <- PBP_AUC_df

# EBP_AUC 3D
EBP_df <- metric_df_list[["EBP"]]
EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
  sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
})
EBP_AUC_df <- EBP_df[ , c("slice", "cell_types", "EBP_AUC")]
metric_df_list[["EBP_AUC"]] <- EBP_AUC_df


### Plotting functions -----

## Functions to plot

# Utility function to add pair column
add_pair_column_to_metric_df <- function(metric_df,
                                         metric) {
  
  # Add pair column
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

plot_median_percentage_difference_of_each_slice_vs_metrics_box_plot <- function(metric_df_list,
                                                                                metrics) {
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
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
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    # Factor slice column so it is from 1, 2, 3, ...
    metric_df$slice <- factor(metric_df$slice, as.character(sort(as.integer(unique(metric_df$slice)))))
    
    
    median_df <- metric_df %>%
      group_by(slice) %>%
      dplyr::summarize(median_error = median(error, na.rm = TRUE), .groups = "drop")
    
    median_df$metric <- metric  
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, median_df)
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Factor metrics
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  summary_df <- plot_df %>%
    group_by(metric) %>%
    summarise(
      min     = min(median_error, na.rm = TRUE),
      max     = max(median_error, na.rm = TRUE),
      range   = max - min,
      median  = median(median_error, na.rm = TRUE),
      average = mean(median_error, na.rm = TRUE)
    )
  write.csv(summary_df, paste("~/R/values_from_figures/", file_name_prefix, "_summary_df.csv", sep = ""))
  
  fig <- ggplot(plot_df, aes(x = metric, y = median_error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(
      # title = "Box plots showing median percentage difference between 2D and 3D metrics of each slice",
      x = "Metric",
      y = "Percentage difference (%)"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Axis tick labels 
      axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 15), 
      
      # Axis titles 
      axis.title.x = element_text(size = 15), 
      axis.title.y = element_text(size = 15), 
      
      # Plot title 
      plot.title = element_text(size = 15)
    ) +
    scale_y_continuous(labels = sci_clean_threshold)
  
  return(fig)
}


## Get plot
metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "COO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")


### Plotting and upload ------
setwd(save_directory)

fig_median_percentage_difference_of_each_slice_vs_metrics_box_plot <- 
  plot_median_percentage_difference_of_each_slice_vs_metrics_box_plot(metric_df_list,
                                                                      metrics)
pdf(paste(file_name_prefix, "_fig_median_percentage_difference_of_each_slice_vs_metrics_box_plot.pdf", sep = ""), width = 9, height = 5)
print(fig_median_percentage_difference_of_each_slice_vs_metrics_box_plot)
dev.off()


print("DONE")




