##02_prep_optim_funtions

# Packages: sf, terra
# External (wp4_functions_utils.R): min_max_normalize(), zero_one_scale()

fill_path <- function(template, site) gsub("{site}", site, template, fixed = TRUE)

## ---- 1. Feature layers -----------------------------------------------------
load_es_means <- function(site, es_dir, es_groups, na_rm = FALSE, normalize = TRUE) {
  files <- list.files(es_dir, pattern = "\\.tif$", full.names = TRUE)
  if (!length(files)) stop("No .tif files in ", es_dir)
  r <- terra::rast(files)
  names(r) <- sub(paste0("^", site, "_(.*)_pred$"), "\\1", names(r))
  
  if (anyDuplicated(names(r))) stop("Duplicated ES layer names: ", paste(names(r), collapse = ", "))
  missing <- setdiff(unlist(es_groups), names(r))
  if (length(missing)) {
    stop("ES layers not found: ", paste(missing, collapse = ", "),
         "\nAvailable: ", paste(names(r), collapse = ", "))
  }
  if (normalize) stopifnot(is.function(min_max_normalize))
  
  lapply(es_groups, function(nm) {
    m <- mean(r[[nm]], na.rm = na_rm)
    if (normalize) min_max_normalize(m) else m
  })
}

load_feature_layers <- function(site, es_dir, es_groups, layer_files,
                                es_na_rm = FALSE, normalize_es = TRUE, skip_missing = FALSE) {
  layers <- load_es_means(site, fill_path(es_dir, site), es_groups, es_na_rm, normalize_es)
  for (nm in names(layer_files)) {
    f <- fill_path(layer_files[[nm]], site)
    if (!file.exists(f)) {
      if (skip_missing) { message("Skipping '", nm, "': file not found (", f, ")"); next }
      stop("File not found for '", nm, "': ", f)
    }
    layers[[nm]] <- terra::rast(f)
  }
  layers
}

## ---- 2. Sample rasters per planning unit -----------------------------------
sample_rasters <- function(grid, rasters, prefix = "sampled_", fun = mean,
                           fallback_point = TRUE) {
  geom <- sf::st_geometry(grid)
  pts  <- sf::st_point_on_surface(geom)
  
  for (nm in names(rasters)) {
    r <- rasters[[nm]]
    if (terra::nlyr(r) > 1) {
      warning("'", nm, "' has ", terra::nlyr(r), " layers; using the first")
      r <- r[[1]]
    }
    rc <- terra::crs(r)
    to_r <- function(g) if (nzchar(rc)) sf::st_transform(g, rc) else g   # sample in the raster's CRS
    
    v <- terra::extract(r, terra::vect(sf::st_sf(geometry = to_r(geom))),
                        fun = fun, na.rm = TRUE)[, 2]
    v[is.nan(v)] <- NA
    
    if (fallback_point && anyNA(v)) {            # polygons without a cell centre inside
      na <- which(is.na(v))
      v[na] <- terra::extract(r, terra::vect(sf::st_sf(geometry = to_r(pts[na]))))[, 2]
    }
    grid[[paste0(prefix, nm)]] <- v
    message(sprintf("%-14s %6d NA cells (%.1f%%)", nm, sum(is.na(v)), 100 * mean(is.na(v))))
  }
  grid
}

## ---- 3. Distance to existing core PA ---------------------------------------
add_distance_to_core <- function(grid, core_col = "existing_corePA") {
  core <- grid[[core_col]] %in% TRUE
  if (!any(core)) {
    warning("No core PA cells: min_distance set to NA")
    grid$min_distance <- NA_real_
    return(grid)
  }
  grid$min_distance <- 0
  if (any(!core)) {
    out <- grid[!core, ]
    pa  <- grid[core, ]
    nearest <- sf::st_nearest_feature(out, pa)
    grid$min_distance[!core] <- as.numeric(sf::st_distance(out, pa[nearest, ], by_element = TRUE))
  }
  grid                                          # original row order is preserved
}

## ---- 4. NA handling -------------------------------------------------------------
fill_na_neighbours <- function(grid, cols, max_iter = 100) {
  nb <- sf::st_touches(grid)
  for (col in cols) {
    s <- grid[[col]]
    for (it in seq_len(max_iter)) {
      miss <- which(is.na(s))
      if (!length(miss)) break
      new <- vapply(miss, function(i) {
        v <- s[nb[[i]]]
        if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
      }, numeric(1))
      if (all(is.na(new))) break
      s[miss] <- new
    }
    grid[[col]] <- s
  }
  grid
}

clean_features <- function(grid, cols, na_action = c("drop", "neighbour_mean"),
                           keep_core = TRUE, core_col = "existing_corePA") {
  na_action <- match.arg(na_action)
  if (na_action == "neighbour_mean") grid <- fill_na_neighbours(grid, cols)
  
  bad  <- Reduce(`|`, lapply(cols, function(cc) is.na(grid[[cc]])))
  core <- if (core_col %in% names(grid)) grid[[core_col]] %in% TRUE else rep(FALSE, nrow(grid))
  drop <- bad & !(keep_core & core)
  
  if (any(drop)) {
    d <- sf::st_drop_geometry(grid)[drop, , drop = FALSE]
    message(sprintf("Dropping %d of %d cells (%.1f km2) with NA in feature columns",
                    sum(drop), nrow(grid), sum(d$area, na.rm = TRUE) / 1e6))
    if ("lulc_name" %in% names(d)) {
      tab <- table(d$lulc_name)
      message("  by LULC: ", paste(names(tab), tab, sep = " = ", collapse = ", "))
    }
  }
  
  # Existing core PAs are kept; their NAs are median-filled (features of locked-in
  # cells are zeroed in pa_optim() anyway) and flagged
  grid$na_core_filled <- FALSE
  fill <- bad & core & keep_core
  if (any(fill)) {
    for (cc in cols) {
      idx <- fill & is.na(grid[[cc]])
      if (any(idx)) grid[[cc]][idx] <- stats::median(grid[[cc]], na.rm = TRUE)
    }
    grid$na_core_filled <- fill
    message(sum(fill), " existing core PA cells had NA features: kept, median-filled")
  }
  
  out <- grid[!drop, ]
  attr(out, "dropped_ids") <- if ("id" %in% names(grid)) grid$id[drop] else which(drop)
  out
}

#scale
zero_one_scale <- function(sf_df, cols = NULL, na.rm = TRUE) {
  
  # If no columns specified: use all numeric columns
  if (is.null(cols)) {
    cols <- sf_df %>%
      st_drop_geometry() %>%
      select(where(is.numeric)) %>%
      names()
  }
  
  sf_df %>%
    mutate(
      across(
        all_of(cols),
        .fns = function(x) {
          
          min_x <- min(x, na.rm = na.rm)
          max_x <- max(x, na.rm = na.rm)
          rng   <- max_x - min_x
          
          if (rng == 0) {
            return(rep(0, length(x)))
          } else {
            return((x - min_x) / rng)
          }
        },
        .names = "{.col}_scaled"
      )
    )
}

#normalize
min_max_normalize <- function(r) {
  # Normalize each layer individually
  norm_r <- lapp(r, function(x) {
    r_min <- min(x, na.rm = TRUE)
    r_max <- max(x, na.rm = TRUE)
    (x - r_min) / (r_max - r_min)
  })
  return(norm_r)
  
}

## ---- 5. Output ---------------------------------------------------------------------
write_grid <- function(grid, dir, file_stub, format = "json") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, paste0(file_stub, ".", format))
  if (format %in% c("json", "geojson")) {
    sf::st_write(grid, f, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  } else {
    sf::st_write(grid, f, delete_dsn = TRUE, quiet = TRUE)
  }
  invisible(f)
}

## ---- 6. Pipeline ------------------------------------------------------------------
prepare_optim_grid <- function(
    target_site,
    grid         = NULL,                      # sf from build_pa_grid(); if NULL it is read from grid_dir
    grid_dir     = "outputs/WP4/01_PA_analysis",
    out_dir      = "outputs/WP4/02_optim",
    es_dir       = "data/WP2/{site}/es_all",
    es_groups    = list(                      # names become layer names (es_prov, ...)
      es_prov = c("wild_plant", "wild_hunt", "farm", "mat"),
      es_reg  = c("erosion", "flood", "habitat"),
      es_cult = c("aest", "recr", "sense")
    ),
    layer_files  = c(
      eco_cond    = "data/WP2/{site}/int.tif",
      cost_policy = "data/WP3/pol_resistance_{site}.tif",
      ghm         = "data/WP4/{site}_GHM.tif"
    ),
    es_na_rm       = FALSE,
    normalize_es   = TRUE,
    skip_missing   = FALSE,
    fun            = mean,
    fallback_point = TRUE,
    na_action      = c("drop", "neighbour_mean"),
    keep_core      = TRUE,
    scale_cols     = NULL,                    # default: all sampled_* columns
    write          = TRUE,
    out_format     = "json"
) {
  na_action <- match.arg(na_action)
  stopifnot(is.function(zero_one_scale))
  
  if (is.null(grid)) {
    grid <- sf::st_read(file.path(grid_dir, paste0(target_site, "_input_grid.json")), quiet = TRUE)
  }
  stopifnot(inherits(grid, "sf"), "existing_corePA" %in% names(grid))
  
  layers    <- load_feature_layers(target_site, es_dir, es_groups, layer_files,
                                   es_na_rm, normalize_es, skip_missing)
  grid      <- sample_rasters(grid, layers, fun = fun, fallback_point = fallback_point)
  feat_cols <- paste0("sampled_", names(layers))
  na_summary <- vapply(feat_cols, function(cc) sum(is.na(grid[[cc]])), integer(1))
  
  grid <- add_distance_to_core(grid)
  
  grid    <- clean_features(grid, feat_cols, na_action, keep_core)
  dropped <- attr(grid, "dropped_ids")
  
  if (is.null(scale_cols)) scale_cols <- feat_cols
  grid <- zero_one_scale(grid, cols = scale_cols)
  
  if (write) write_grid(grid, out_dir, paste0(target_site, "_input_final_grid"), out_format)
  
  list(grid = grid, layers = names(layers), na_summary = na_summary, dropped_ids = dropped)
}