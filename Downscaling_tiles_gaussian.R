# ============================================================
# R Script for Tiles Downscaling (Gaussian interpolation on tiles)
#
# Input:
#   - .dat files containing point data (x, y, valeur)
#   - Tile polygons with buffer (for interpolation window)
#   - Tile polygons without buffer (for masking final raster)
#
# Output:
#   - One GeoTIFF per tile and per .dat file (e.g. monthly), saved
#     under a directory based on year_month extracted from the .dat filename.
#
# Notes:
#   - Interpolation is done with a Gaussian kernel (C++/Rcpp function).
#   - Parallelization is used over tiles (parLapply + cluster).
# ============================================================

# --- Libraries ---
library(sf)
library(dplyr)
library(parallel)
library(Rcpp)
library(raster)
library(pryr)   # optional, for memory usage monitoring
library(sp)

# --- USER PARAMETERS (edit as needed) ---

# Directory containing .dat files
dat_files_directory <- "PATH/TO/DAT/FILES/"

# Directories for tiles (with and without buffer)
tiles_directory_with_buffer    <- "PATH/TO/TILES/with_buffer/"
tiles_directory_without_buffer <- "PATH/TO/TILES/without_buffer/"

# Base output directory for GeoTIFF results
output_base_directory <- "PATH/TO/OUTPUT/RESULTS/"

# Pixel resolution in meters
resolution <- 50

# Coordinate reference system (e.g. EPSG:2154 = Lambert 93)
crs_proj <- 2154

# Gaussian kernel bandwidth (in the same units as coordinates, e.g. meters)
sigma_value <- 500

# Maximum number of nearest points used for interpolation
nmax_points <- 12

# Number of cores used in parallel processing
num_cores <- 6

# ======================================================================
# PREPARATION
# ======================================================================

# List all .dat files
dat_files <- list.files(
  dat_files_directory,
  pattern = "\\.dat$",
  full.names = TRUE
)

# Create cluster
cl <- makeCluster(num_cores)

# Export needed objects to the cluster
clusterExport(
  cl,
  list(
    "tiles_directory_with_buffer",
    "tiles_directory_without_buffer",
    "resolution",
    "output_base_directory",
    "crs_proj",
    "sigma_value",
    "nmax_points"
  )
)

# Define interpolation function (Rcpp) inside cluster
clusterEvalQ(cl, {
  library(Rcpp)
  
  cppFunction('
    NumericVector smooth_gaussian_cpp(
      NumericVector x, NumericVector y, NumericVector z,
      NumericVector xi, NumericVector yi,
      double sigma, int nmax
    ) {
      int n = xi.size();
      NumericVector zi(n);

      for (int k = 0; k < n; ++k) {
        double x0 = xi[k];
        double y0 = yi[k];

        std::vector<std::pair<double, int>> distances;

        // Collect candidate points within 3 * sigma
        for (int i = 0; i < x.size(); ++i) {
          double dist = sqrt(pow(x0 - x[i], 2) + pow(y0 - y[i], 2));
          if (dist < sigma * 3) {
            distances.push_back(std::make_pair(dist, i));
          }
        }

        // Sort by distance
        std::sort(distances.begin(), distances.end());

        // Keep at most nmax nearest points
        if (distances.size() > nmax) {
          distances.resize(nmax);
        }

        double num = 0.0;
        double den = 0.0;

        for (int j = 0; j < distances.size(); ++j) {
          int i = distances[j].second;
          double dist = distances[j].first;
          double w = exp(-pow(dist, 2) / (2 * pow(sigma, 2)));
          num += w * z[i];
          den += w;
        }

        zi[k] = (den > 0) ? num / den : NA_REAL;
      }
      return zi;
    }
  ')
})

# ======================================================================
# FUNCTIONS
# ======================================================================

# Function to process a single .dat file
process_dat_file <- function(dat_file) {
  tryCatch({
    
    dat_base_name <- tools::file_path_sans_ext(basename(dat_file))
    
    # Extract year_month (e.g. 2010_01) from filename
    year_month <- sub(".*(\\d{4}_\\d{2}).*", "\\1", dat_base_name)
    
    # Output directory for this .dat file
    output_dir <- file.path(output_base_directory, year_month)
    
    if (dir.exists(output_dir)) {
      cat("Output exists, skipping:", output_dir, "\n")
      return(NULL)
    } else {
      dir.create(output_dir, recursive = TRUE)
      cat("Created directory:", output_dir, "\n")
    }
    
    # Read .dat file
    data <- read.table(
      dat_file,
      header = TRUE,
      sep = "",
      stringsAsFactors = FALSE
    )
    data$valeur <- as.numeric(data$valeur)
    
    # Convert to sf object
    coordinates_sf <- st_as_sf(data, coords = c("x", "y"), crs = crs_proj)
    coords <- st_coordinates(coordinates_sf)
    coordinates_sf$x <- coords[, 1]
    coordinates_sf$y <- coords[, 2]
    
    # List tile files (with and without buffer)
    tiles_files_with_buffer <- list.files(
      tiles_directory_with_buffer,
      pattern = "tile_.*_with_buffer\\.gpkg$",
      full.names = TRUE
    )
    tiles_files_without_buffer <- list.files(
      tiles_directory_without_buffer,
      pattern = "tile_.*_without_buffer\\.gpkg$",
      full.names = TRUE
    )
    
    cat("Tiles to process:", length(tiles_files_with_buffer), "\n")
    
    # Function to process a single tile (index i)
    process_tile <- function(i) {
      library(sf)
      library(raster)
      library(sp)
      library(dplyr)
      
      # Read tile with buffer
      tile_with_buffer <- st_read(tiles_files_with_buffer[i], quiet = TRUE)
      bbox_tile <- st_bbox(tile_with_buffer)
      
      xmin <- bbox_tile["xmin"]
      xmax <- bbox_tile["xmax"]
      ymin <- bbox_tile["ymin"]
      ymax <- bbox_tile["ymax"]
      
      # Filter points within tile + buffer
      tile_data <- coordinates_sf %>%
        dplyr::filter(
          x >= xmin & x <= xmax &
            y >= ymin & y <= ymax
        )
      
      if (nrow(tile_data) == 0) {
        cat("No data in tile with buffer", i, "\n")
        return(NULL)
      }
      
      # Create interpolation grid
      x_seq <- seq(xmin, xmax, by = resolution)
      y_seq <- seq(ymin, ymax, by = resolution)
      grid <- expand.grid(x = x_seq, y = y_seq)
      
      # Gaussian interpolation
      interp_vals <- smooth_gaussian_cpp(
        x  = tile_data$x,
        y  = tile_data$y,
        z  = tile_data$valeur,
        xi = grid$x,
        yi = grid$y,
        sigma = sigma_value,
        nmax  = nmax_points
      )
      
      grid$valeur <- interp_vals
      coordinates(grid) <- ~ x + y
      gridded(grid) <- TRUE
      proj4string(grid) <- CRS(paste0("+init=epsg:", crs_proj))
      
      r <- raster(grid["valeur"])
      
      # Masking with tile without buffer
      tile_without_buffer <- st_read(
        tiles_files_without_buffer[i],
        quiet = TRUE
      )
      tile_mask <- rasterize(tile_without_buffer, r)
      r_masked <- mask(r, tile_mask)
      
      # Output filename (here hard-coded as "pr_..." – change if needed)
      output_tif_path <- file.path(
        output_dir,
        paste0(
          "pr_",
          year_month,
          "_tile_",
          i,
          "_Gaussian_kernel_",
          resolution,
          "m.tif"
        )
      )
      
      if (!file.exists(output_tif_path)) {
        writeRaster(r_masked, output_tif_path, format = "GTiff", overwrite = TRUE)
        cat("Saved:", output_tif_path, "\n")
      } else {
        cat("Already exists:", output_tif_path, "\n")
      }
    }
    
    # Parallel loop over tiles
    parLapply(cl, 1:length(tiles_files_with_buffer), process_tile)
    
  }, error = function(e) {
    cat("Error processing:", dat_file, "\nMessage:", e$message, "\n")
  })
}

# ======================================================================
# MAIN LOOP
# ======================================================================

for (dat_file in dat_files) {
  cat("Processing file:", dat_file, "\n")
  process_dat_file(dat_file)
}

stopCluster(cl)
cat("All .dat files processed.\n")
