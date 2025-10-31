#' Manual spatial cell selection within a Seurat object
#'
#' @importFrom Seurat ImageDimPlot
#' @importFrom utils head
#' @importFrom methods is
#'
#' @param seurat_obj A Seurat object containing spatial data.
#' @param meta_id Name of the metadata column identifying the FOVs (e.g. "fov", "Region", "sample_id")
#' @param fovs Names of the images if `meta_id` is not used (e.g. in Xenium objects)
#'            or character vector of FOV identifiers to analyze (must match the values in `meta_id`)
#' @param images Optional character vector of image names (same length as `fovs`).
#'               If NULL, the function will automatically attempt to use `fovs` as image names.
#' @param metadata_col Name of the metadata column to store selection results (default "ManualSelection").
#' @param positive_label Label for selected cells (default "Selected").
#' @param negative_label Label for unselected cells (default "Not_Selected").
#' @param incremental Logical. If TRUE, keep existing labels in `metadata_col` and add new `positive_label`
#'                    instead of resetting the column (default FALSE).
#' @param point_size Numeric. Size of points in the interactive selection plot (default: 3).
#' @param invert_coordinates Logical. If TRUE, swaps x/y coordinates in interactive selection (default: FALSE).
#' @param cols Optional color vector for idents passed to ImageDimPlot.
#'
#' @return The modified Seurat object with updated metadata column.
#' @export
run_spatial_selection <- function(seurat_obj,
                                  meta_id = NULL,
                                  fovs,
                                  images = NULL,
                                  metadata_col = "ManualSelection",
                                  positive_label = "Selected",
                                  negative_label = "Not_Selected",
                                  incremental = FALSE,
                                  point_size = 3,
                                  invert_coordinates = FALSE,
                                  cols = NULL) {
  # ---- Input checks ----
  if (!methods::is(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.")
  if (missing(fovs) || length(fovs) < 1) stop("Provide at least one fov or image name in 'fovs'.")
  if (!is.null(images) && length(images) != length(fovs)) {
    stop("If 'images' is provided, it must have the same length as 'fovs'.")
  }


  # ---- Initialize or preserve metadata ----
  if (!(metadata_col %in% colnames(seurat_obj@meta.data))) {
    message("Creating new metadata column: ", metadata_col)
    seurat_obj@meta.data[[metadata_col]] <- negative_label
  } else {
    # coerce to character to avoid factor issues
    seurat_obj@meta.data[[metadata_col]] <- as.character(seurat_obj@meta.data[[metadata_col]])

    # reset only if not incremental
    if (!incremental) {
      message("Resetting metadata column: ", metadata_col)
      seurat_obj@meta.data[[metadata_col]] <- rep(negative_label, ncol(seurat_obj))
    } else {
      message("Incremental mode enabled: existing labels in ", metadata_col, " will be preserved.")
    }
  }


  # ---- Main loop ----
  for (i in seq_along(fovs)) {
    fov <- fovs[i]
    message("Processing FOV/Image: ", fov)

    # If meta_id provided : subset by metadata
    if (!is.null(meta_id) && meta_id %in% colnames(seurat_obj@meta.data)) {
      subset_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[meta_id]] == fov]

      if (length(subset_cells) == 0) {
        warning("No cells found for FOV: ", fov, ". Skipping.")
        next
      }

      # determine image name to pass to ImageDimPlot (if images provided)
      image_name <- NULL
      if (!is.null(images)) image_name <- images[i]

      # build the plot
      subset_plot <- NULL

      # if image_name is NULL, do not pass fov (ImageDimPlot will use default)
      subset_plot <- tryCatch(
        {
          if (is.null(image_name)) {
            Seurat::ImageDimPlot(seurat_obj, cells = subset_cells,  cols = cols, boundaries = "centroids")
          } else {
            Seurat::ImageDimPlot(seurat_obj, cells = subset_cells, fov = image_name, cols = cols, boundaries = "centroids")
          }
        },
        error = function(e) {
          warning("ImageDimPlot failed for FOV ", fov, ": ", conditionMessage(e))
          return(NULL)
        }
      )

    # If meta_id NOT provided : subset by IMAGE layer
    } else {
      # Otherwise, assume fov is directly the image name (for Xenium-like objects)

      # build the plot
      subset_plot <- NULL

      # if image_name is NULL, do not pass fov (ImageDimPlot will use default)
      subset_plot <- tryCatch(
          Seurat::ImageDimPlot(seurat_obj, fov = fov, cols = cols, boundaries = "centroids"),
        error = function(e) {
          warning("ImageDimPlot failed for FOV ", fov, ": ", conditionMessage(e))
          return(NULL)
        }
      )

    }


    if (is.null(subset_plot)) next

    # ---- Interactive selection ----
    selectedPoints <- tryCatch(
      interactive_point_selection(subset_plot,
                                  invert_coordinates = invert_coordinates,
                                  point_size = point_size),
      error = function(e) {
        warning("Interactive selection failed for ", fov, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(selectedPoints) || nrow(selectedPoints) == 0) {
      message("No points selected for FOV/Image: ", fov)
      next
    }

    # ---- Map selected points back ----
    plot_data <- subset_plot$data
    merged_selected <- tryCatch(
      merge(plot_data, selectedPoints, by = c("x", "y")),
      error = function(e) {
        warning("Merge of plot data and selected points failed for ", fov, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(merged_selected) || nrow(merged_selected) == 0) {
      warning("No mapped cells after merging selected points for FOV ", fov, ". Skipping.")
      next
    }


    # clean cell names and deduplicate
    selected_subset <- unique(sub("^centroids_", "", merged_selected$cell))
    # keep only those that truly exist in the Seurat object
    selected_subset <- selected_subset[selected_subset %in% colnames(seurat_obj)]
    if (length(selected_subset) == 0) {
      warning("No valid cell names found after mapping for FOV ", fov, ". Skipping.")
      next
    }

    # ---- Update metadata safely by rows ----
    mask <- colnames(seurat_obj) %in% selected_subset
    seurat_obj@meta.data[mask, metadata_col] <- positive_label

    message("Number of selected cells: ", length(selected_subset))
  }

  return(seurat_obj)
}
