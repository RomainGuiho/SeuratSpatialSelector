
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SeuratSpatialSelector

<!-- badges: start -->

<!-- badges: end -->

The goal of SeuratSpatialSelector is to …

## Installation

You can install the development version of SeuratSpatialSelector like
so:

``` r
devtools::install_github("RomainGuiho/SeuratSpatialSelector")
```

## Basic example

This is a basic example:

``` r
library(SeuratSpatialSelector)

## Create fov metadata column in the Seurat object

### Nanostring CosMx object 
CosMx_seurat_object@meta.data[["fovs"]] <- gsub("_[^_]+_", "_", CosMx_seurat_object@images[["FOV"]]@boundaries[["centroids"]]@cells)

### 10X Xenium object 
Xenium_seurat_object@meta.data[["fovs"]] <- gsub("_[^_]+_", "_", Xenium_seurat_object@images[["FOV"]]@boundaries[["centroids"]]@cells)


## basic example code

list_fov <- c(1,3,5)

spatial_seurat_object <- run_spatial_selection(
  seurat_obj   = spatial_seurat_object,
  meta_id      = "slide_fov",
  fovs         = list_fov,
  metadata_col = "ManualSelection",
  positive_label = "Selected",
  negative_label = "Not_Selected",
  cols = NULL
)
```

## Multi-images example

### Nanostring CosMx object created with LoadNanostring()

``` r

## Create slide_fov column in the Seurat object

spatial_seurat_object@meta.data[["slide_fov"]] <- c(
  gsub("_[^_]+_", "_", spatial_seurat_object@images[["slide_1"]]@boundaries[["centroids"]]@cells),
  gsub("_[^_]+_", "_", spatial_seurat_object@images[["slide_2"]]@boundaries[["centroids"]]@cells)
  )
```

### 10X Xenium object created with LoadXenium()

``` r

spatial_seurat_object@meta.data[["slide_fov"]] <- c(
  gsub("_[^_]+_", "_", spatial_seurat_object@images[["slide_1"]]@boundaries[["centroids"]]@cells),
  gsub("_[^_]+_", "_", spatial_seurat_object@images[["slide_2"]]@boundaries[["centroids"]]@cells)
  )
```
