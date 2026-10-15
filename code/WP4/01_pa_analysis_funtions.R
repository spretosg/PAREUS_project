############################################################################
# Pareus WP4 - Protected area analysis functions
# Date: 10.10.2026
# Author: Reto Spielhofer (Norwegian Institute for Nature Research), code revisons by Claude Code 
############################################################################

# Packages: sf, terra, dplyr, tidyr, ggplot2, scales




## ---- 1. Inputs ----------------------------------------------------------------
site_crs <- function(site, lookup = c(TRD = 25833, FRA = 2154), default = 5514) {
  if (site %in% names(lookup)) unname(lookup[[site]]) else default
}

read_site_inputs <- function(site, shared_dir = "data/shared", pa_dir = "data/WP4",
                             crs = site_crs(site)) {
  stud_area <- sf::read_sf(file.path(shared_dir, paste0(site, ".gpkg")))
  PA        <- sf::st_read(file.path(pa_dir, paste0("WDPA_", site, ".gpkg")), quiet = TRUE)
  lulc      <- terra::rast(file.path(shared_dir, paste0(site, "_lulc.tif")))
  list(stud_area = sf::st_transform(stud_area, crs),
       PA        = sf::st_transform(PA, crs),
       lulc      = lulc,                       # kept in its own CRS on purpose
       crs       = crs)
}

clip_pa <- function(PA, stud_area) {
  stopifnot("IUCN_CAT" %in% names(PA))
  area_geom <- sf::st_union(sf::st_geometry(stud_area))
  PA <- sf::st_make_valid(PA)                                   # validate BEFORE intersecting
  PA <- suppressWarnings(sf::st_intersection(PA, area_geom))
  PA <- sf::st_collection_extract(PA, "POLYGON")
  PA[!sf::st_is_empty(PA), ]
}

## ---- 2. Planning units -------------------------------------------------------
build_planning_units <- function(stud_area, cell_size = 600, hexagonal = TRUE) {
  area_geom <- sf::st_union(sf::st_geometry(stud_area))
  g    <- sf::st_make_grid(area_geom, cellsize = cell_size, square = !hexagonal)
  grid <- sf::st_sf(id = seq_along(g), geometry = g)
  grid <- suppressWarnings(sf::st_intersection(grid, area_geom))
  grid <- sf::st_collection_extract(grid, "POLYGON")
  grid <- grid[!sf::st_is_empty(grid), ]
  stopifnot(!anyDuplicated(grid$id))
  grid$area     <- as.numeric(sf::st_area(grid))   # m2, as used by protection_gap()
  grid$area_km2 <- grid$area / 1e6
  grid
}

## ---- 3. LULC -------------------------------------------------------------------
assign_lulc <- function(grid, lulc, divisor = 100, lulc_names = default_lulc_names,
                        method = c("centroid", "majority")) {
  method   <- match.arg(method)
  lulc_crs <- terra::crs(lulc)
  stopifnot(nzchar(lulc_crs))
  
  if (method == "centroid") {
    pts  <- sf::st_transform(sf::st_point_on_surface(sf::st_geometry(grid)), lulc_crs)
    code <- floor(terra::extract(lulc, terra::vect(sf::st_sf(geometry = pts)))[, 2] / divisor)
  } else {
    polys <- sf::st_transform(sf::st_geometry(grid), lulc_crs)
    code  <- terra::extract(
      lulc, terra::vect(sf::st_sf(geometry = polys)),
      fun = function(x) {
        x <- floor(x[!is.na(x)] / divisor)
        if (!length(x)) NA_real_ else as.numeric(names(which.max(table(x))))
      }
    )[, 2]
  }
  
  grid$sampled_habitat <- code
  grid <- grid[!is.na(grid$sampled_habitat) & grid$sampled_habitat != 0, ]
  grid$lulc_name <- unname(lulc_names[as.character(grid$sampled_habitat)])
  if (anyNA(grid$lulc_name)) {
    warning("LULC codes without a name in `lulc_names`: ",
            paste(unique(grid$sampled_habitat[is.na(grid$lulc_name)]), collapse = ", "))
  }
  grid
}

lulc_area_stats <- function(grid) {
  grid |>
    sf::st_drop_geometry() |>
    dplyr::group_by(sampled_habitat, lulc_name) |>
    dplyr::summarise(area_km2 = sum(area) / 1e6, .groups = "drop")
}

## ---- 4. PA status ----------------------------------------------------------------
assign_pa_status <- function(grid, PA, core_pa = c("Ia", "II"),
                             iucn_rank = default_iucn_rank) {
  if (nrow(PA) == 0) {
    grid$n_pa <- 0L; grid$n_pa_poly <- 0L
    grid$IUCN_CAT <- NA_character_; grid$prot_class <- NA_real_
    grid$existing_corePA <- FALSE
    return(grid)
  }
  unknown <- setdiff(unique(PA$IUCN_CAT), names(iucn_rank))
  if (length(unknown)) {
    warning("IUCN categories not in `iucn_rank` (rank 0): ", paste(unknown, collapse = ", "))
  }
  PA$prot_class <- unname(iucn_rank[PA$IUCN_CAT])
  PA$prot_class[is.na(PA$prot_class)] <- 0
  
  ints <- sf::st_intersects(grid, PA)
  
  grid$n_pa      <- vapply(ints, function(i) length(unique(PA$IUCN_CAT[i])), integer(1))  # as in original
  grid$n_pa_poly <- lengths(ints)                                                          # number of PA polygons
  
  top <- vapply(ints, function(i) if (length(i)) i[which.max(PA$prot_class[i])] else NA_integer_,
                integer(1))
  grid$IUCN_CAT   <- PA$IUCN_CAT[top]       # highest-ranked category per cell
  grid$prot_class <- PA$prot_class[top]
  
  grid$existing_corePA <- vapply(ints, function(i) any(PA$IUCN_CAT[i] %in% core_pa), logical(1))
  grid
}

## ---- 5. Gap analysis -----------------------------------------------------------------
compute_gap <- function(grid, global_target, lulc_targets, lockout = "built-up",
                        lulc_col = "lulc_name", corePA_col = "existing_corePA") {
  gap_glob <- protection_gap(grid, lulc_col = lulc_col, corePA_col = corePA_col,
                             lockout = lockout, mode = "global", target = global_target)
  gap_lulc <- protection_gap(grid, lulc_col = lulc_col, corePA_col = corePA_col,
                             lockout = lockout, mode = "class", target = lulc_targets)
  gap <- rbind(gap_lulc, gap_glob)
  gap[!is.na(gap$lulc_name), ]
}


protection_gap <- function(x,
                           mode = c("class", "global"),
                           target,
                           lulc_col = "lulc_name",
                           corePA_col = "core_PA",
                           lockout = NULL,
                           area_col = NULL) {
  
  mode <- match.arg(mode)
  
  dat <- x
  
  # Remove lock-out LULC classes
  if (!is.null(lockout)) {
    dat <- dat %>%
      filter(!(.data[[lulc_col]] %in% lockout))
  }
  
  # area in m²
  if (is.null(area_col)) {
    dat$area <- as.numeric(st_area(dat))
  } else {
    dat$area <- dat[[area_col]]
  }
  
  
  
  if (mode == "class") {
    
    # Summary by LULC
    gap <- dat %>%
      st_drop_geometry() %>%
      group_by(.data[[lulc_col]]) %>%
      summarise(
        total_area = sum(area),
        protected_area = sum(area[.data[[corePA_col]]]),
        .groups = "drop"
      )
    stopifnot(!is.null(names(target)))
    
    gap <- gap %>%
      mutate(
        target_fraction = target[as.character(.data[[lulc_col]])],
        target_area = target_fraction * total_area
      )
    
  } else {
    gap <- dat %>%
      st_drop_geometry() %>%
      summarise(
        lulc_name = paste0("global_",target,"_prot"),
        total_area = sum(area),
        protected_area = sum(area[.data[[corePA_col]]])
      )
    stopifnot(length(target) == 1)
    #excluding lockout
    total_area <- sum(gap$total_area)
    
    gap <- gap %>%
      mutate(
        target_fraction = target,
        target_area = total_area * target
      )
  }
  
  gap %>%
    mutate(
      gap_area = pmax(target_area - protected_area, 0),
      surplus_area = pmax(protected_area - target_area, 0),
      protected_fraction = protected_area / total_area,
      gap_fraction = gap_area / total_area
    )
}
## ---- 6. Plots ------------------------------------------------------------------------
plot_n_pa <- function(grid, stud_area, text_size = 20) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = grid[grid$n_pa > 0, ], ggplot2::aes(fill = .data$n_pa), colour = NA) +
    ggplot2::scale_fill_gradient(low = "#FFC2C2", high = "#C50202",
                                 breaks = scales::breaks_width(1),
                                 labels = function(x) as.integer(x), na.value = "grey50") +
    ggplot2::geom_sf(data = stud_area, fill = NA, colour = "black") +
    ggplot2::labs(fill = "n PA") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = text_size))
}

plot_iucn <- function(grid, stud_area, iucn_colours = default_iucn_colours,
                      built_up = "built-up") {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = grid, ggplot2::aes(fill = .data$IUCN_CAT), colour = NA) +
    ggplot2::scale_fill_manual(values = iucn_colours, na.value = NA, name = "IUCN") +
    ggplot2::geom_sf(data = grid[grid$lulc_name %in% built_up, ], fill = "grey", colour = NA) +
    ggplot2::geom_sf(data = stud_area, fill = NA, colour = "black") +
    ggplot2::theme_minimal()
}

plot_iucn_area <- function(grid, iucn_rank = default_iucn_rank, text_size = 20) {
  d <- sf::st_drop_geometry(grid)
  d <- d[!is.na(d$IUCN_CAT), ]
  d <- stats::aggregate(area ~ IUCN_CAT, data = d, FUN = sum)
  d$IUCN_CAT <- factor(d$IUCN_CAT,
                       levels = union(names(sort(iucn_rank, decreasing = TRUE)), d$IUCN_CAT))
  ggplot2::ggplot(d, ggplot2::aes(IUCN_CAT, area / 1e6)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::labs(x = "IUCN", y = expression(Area~(km^2))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = text_size))
}

plot_gap <- function(gap) {
  d <- gap |>
    dplyr::select(lulc_name, protected_area, gap_area) |>
    tidyr::pivot_longer(c(protected_area, gap_area), names_to = "status", values_to = "area")
  ggplot2::ggplot(d, ggplot2::aes(y = stats::reorder(lulc_name, area), x = area / 1e6, fill = status)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(values = c(protected_area = "#00A300", gap_area = "#D95F02"),
                               breaks = c("gap_area", "protected_area"),
                               labels = c("Gap", "Existing core PA")) +
    ggplot2::labs(x = expression("Area (km"^2*")"), y = NULL, fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top", panel.grid.major.y = ggplot2::element_blank())
}

## ---- 7. Pipeline: slow part (grid) ---------------------------------------------------
build_pa_grid <- function(
    target_site,
    core_pa      = c("Ia", "II"),      # WDPA categories counted as strictly protected
    pu_size      = 600,                # hexagon flat-to-flat distance in map units (m)
    hexagonal    = TRUE,
    lulc_divisor = 100,
    lulc_names   = default_lulc_names,
    lulc_method  = c("centroid", "majority"),
    iucn_rank    = default_iucn_rank,
    crs          = site_crs(target_site),
    shared_dir   = "data/shared",
    pa_dir       = "data/WP4",
    out_dir      = NULL,               # set to write the grid
    grid_format  = "json"
) {
  lulc_method <- match.arg(lulc_method)
  inp  <- read_site_inputs(target_site, shared_dir, pa_dir, crs)
  PA   <- clip_pa(inp$PA, inp$stud_area)
  
  grid <- build_planning_units(inp$stud_area, pu_size, hexagonal) |>
    assign_lulc(inp$lulc, lulc_divisor, lulc_names, lulc_method) |>
    assign_pa_status(PA, core_pa, iucn_rank)
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(out_dir, paste0(target_site, "_input_grid.", grid_format))
    if (grid_format %in% c("json", "geojson")) {
      sf::st_write(grid, f, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
    } else {
      sf::st_write(grid, f, delete_dsn = TRUE, quiet = TRUE)
    }
  }
  
  list(site = target_site, grid = grid, stud_area = inp$stud_area, PA = PA,
       lulc_stats = lulc_area_stats(grid), core_pa = core_pa, iucn_rank = iucn_rank)
}

## ---- 8. Pipeline: fast part (gap + plots) --------------------------------------------
analyse_pa_gap <- function(
    x,                                   # result of build_pa_grid()
    global_target = 0.1,
    lulc_targets  = c(forest = 0.1, water = 0.1, wetland = 0.1, agricultural = 0),
    lockout       = "built-up",
    out_dir       = NULL
) {
  gap <- compute_gap(x$grid, global_target, lulc_targets, lockout)
  
  plots <- list(
    n_pa      = plot_n_pa(x$grid, x$stud_area),
    iucn      = plot_iucn(x$grid, x$stud_area),
    iucn_area = plot_iucn_area(x$grid, x$iucn_rank),
    gap       = plot_gap(gap)
  )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(gap, file.path(out_dir, paste0(x$site, "_gap_analysis.csv")), row.names = FALSE)
    suffix <- c(n_pa = "_n_pa.png", iucn = "_IUCN.png",
                iucn_area = "_area_stats_IUCN.png", gap = "_gap.png")
    for (nm in names(plots)) {
      ggplot2::ggsave(file.path(out_dir, paste0(x$site, suffix[[nm]])), plots[[nm]],
                      width = 8, height = 6, dpi = 300)
    }
  }
  list(gap = gap, plots = plots)
}