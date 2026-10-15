# Packages: sf, prioritizr, dplyr, ggplot2
# Uses write_grid() from the previous step
# Gap table columns: lulc_name, target_area (m2), protected_area, gap_area

## ---- 0. Defaults --------------------------------------------------------------
pa_levels <- c("existing core pa", "new core PA", "proposed upgrade existing PA",
               "other protected areas", "not protected")
pa_cols <- c("existing core pa" = "#00A300", "new core PA" = "#bf8bff",
             "proposed upgrade existing PA" = "#e5d0ff",
             "other protected areas" = "#ADD8E6", "not protected" = "#1A1A1A")



## ---- 1. Optimiser ----------------------------------------------------------------
pa_optim <- function(pu, features, area_budget,
                     lockin_col = NULL, target_lulc = NULL, lockout = NULL,
                     boundary_penalty = 5e-4, edge_factor = 0.5,
                     boundary_data = NULL,   # boundary_matrix() of the INPUT pu (same row order)
                     mip_gap = 0.01, time_limit = 600) {
  stopifnot(
    inherits(pu, "sf"),
    is.numeric(area_budget), length(area_budget) == 1, is.finite(area_budget), area_budget > 0,
    all(c("area_km2", "lulc_name", features) %in% names(pu))
  )
  if (anyNA(sf::st_drop_geometry(pu)[features])) stop("NA in feature columns")
  
  n0 <- nrow(pu)
  if (!is.null(lockin_col)) {
    if (!lockin_col %in% names(pu)) stop("lockin_col not found in pu")
    if (anyNA(pu[[lockin_col]]))    stop("lockin_col contains NA")
    locked <- as.logical(pu[[lockin_col]])
  } else {
    locked <- rep(FALSE, n0)
  }
  
  # Filters restrict candidates only; existing PAs stay in as neutral context
  cand <- rep(TRUE, n0)
  if (!is.null(lockout))     cand <- cand & !(pu$lulc_name %in% lockout)
  if (!is.null(target_lulc)) cand <- cand &  (pu$lulc_name %in% target_lulc)
  
  keep   <- cand | locked
  pu     <- pu[keep, ]
  locked <- locked[keep]
  if (!any(!locked)) return(NULL)
  
  # Locked-in cells: no cost, no feature contribution (budget = ADDITIONAL area)
  pu$cost_km2 <- ifelse(locked, 0, pu$area_km2)
  for (f in features) pu[[f]][locked] <- 0
  pu$.is_lockin <- locked
  
  avail <- sum(pu$cost_km2)
  diag  <- list(budget_km2 = area_budget, available_km2 = avail)
  if (area_budget >= avail) {
    warning(sprintf("area_budget (%.1f km2) >= available candidate area (%.1f km2): selecting all candidates",
                    area_budget, avail), call. = FALSE)
    out <- pu[!locked, ]
    out$solution_1 <- 1
    attr(out, "diag") <- c(diag, list(status = "all_candidates", runtime = 0,
                                      selected_km2 = sum(out$area_km2)))
    return(out)
  }
  
  p <- prioritizr::problem(pu, features = features, cost_column = "cost_km2") |>
    prioritizr::add_max_utility_objective(budget = area_budget) |>
    prioritizr::add_binary_decisions() |>
    prioritizr::add_default_solver(gap = mip_gap, time_limit = time_limit, verbose = FALSE)
  
  if (boundary_penalty > 0) {
    bd <- if (is.null(boundary_data)) NULL else boundary_data[keep, keep]
    p <- prioritizr::add_boundary_penalties(p, penalty = boundary_penalty,
                                            edge_factor = edge_factor, data = bd)
  }
  if (any(locked)) p <- prioritizr::add_locked_in_constraints(p, ".is_lockin")
  
  s   <- solve(p)
  st  <- list(status = attr(s, "status"), runtime = as.numeric(attr(s, "runtime")),
              objective = attr(s, "objective"))
  out <- s[s$solution_1 == 1 & !s$.is_lockin, ]            # new picks only
  attr(out, "diag") <- c(diag, st, list(selected_km2 = sum(out$area_km2)))
  out
}

## ---- 2. Features ------------------------------------------------------------------
scale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

add_distance_scaled <- function(pu, transform = c("linear", "log", "rank"), col = "min_distance") {
  transform <- match.arg(transform)
  stopifnot(col %in% names(pu))
  d <- pu[[col]]
  pu[[paste0(col, "_scaled")]] <- switch(
    transform,
    linear = scale01(d),
    log    = scale01(log1p(d)),
    rank   = (rank(d, ties.method = "average") - 1) / (length(d) - 1)
  )
  pu
}

drop_partial_pus <- function(pu, min_area_km2 = 0.075, keep_core = TRUE,
                             core_col = "existing_corePA") {
  core  <- if (core_col %in% names(pu)) pu[[core_col]] %in% TRUE else rep(FALSE, nrow(pu))
  small <- pu$area_km2 <= min_area_km2
  drop  <- small & !(keep_core & core)
  message(sprintf("Removing %d partial PUs (<= %.3f km2; %.1f km2 in total); %d small core PA cells kept",
                  sum(drop), min_area_km2, sum(pu$area_km2[drop]), sum(small & core & keep_core)))
  pu[!drop, ]
}

# spec: data.frame(column, direction = "max"/"min", weight)
build_features <- function(pu, spec, area_weighted = FALSE) {
  stopifnot(is.data.frame(spec), all(c("column", "direction") %in% names(spec)),
            all(spec$direction %in% c("max", "min")))
  if (is.null(spec$weight)) spec$weight <- 1
  missing <- setdiff(spec$column, names(pu))
  if (length(missing)) stop("Feature columns not found: ", paste(missing, collapse = ", "))
  
  w_area <- if (area_weighted) pu$area_km2 / max(pu$area_km2) else 1
  feats  <- character(nrow(spec))
  for (i in seq_len(nrow(spec))) {
    x <- pu[[spec$column[i]]]
    if (anyNA(x)) stop("NA in ", spec$column[i])
    if (any(x < -1e-9 | x > 1 + 1e-9)) stop(spec$column[i], " is not on a 0-1 scale")
    if (spec$direction[i] == "min") x <- 1 - x
    feats[i] <- paste0("opt_", spec$column[i])
    pu[[feats[i]]] <- x * spec$weight[i] * w_area
  }
  list(pu = pu, features = feats)
}

## ---- 3. Runs ------------------------------------------------------------------------
get_budget_km2 <- function(gap, name, budget_col = "target_area") {
  v <- gap[[budget_col]][gap$lulc_name == name]
  if (length(v) != 1 || !is.finite(v) || v <= 0) {
    stop(sprintf("No single positive '%s' for '%s' in the gap table", budget_col, name))
  }
  v / 1e6
}

run_optimizations <- function(
    pu, gap, features,
    lulcs            = NULL,            # default: all classes with a budget > 0 (no global row)
    run_global       = TRUE,
    global_name      = NULL,            # default: the single 'global*' row
    global_lockout   = c("built-up", "agricultural"),
    budget_col       = "target_area",
    lockin_col       = "existing_corePA",
    boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
    mip_gap = 0.01, time_limit = 600
) {
  is_global <- grepl("^global", gap$lulc_name)
  if (is.null(lulcs)) {
    ok    <- !is.na(gap$lulc_name) & !is_global & gap[[budget_col]] > 0
    lulcs <- setdiff(gap$lulc_name[ok], "built-up")
  }
  if (run_global && is.null(global_name)) {
    gl <- gap$lulc_name[is_global]
    if (length(gl) != 1) stop("Expected one 'global*' row in the gap table, found: ",
                              paste(gl, collapse = ", "))
    global_name <- gl
  }
  
  bd <- if (share_boundary && boundary_penalty > 0) prioritizr::boundary_matrix(pu) else NULL
  
  one <- function(run, budget, ...) {
    t0  <- Sys.time()
    err <- NA_character_
    res <- tryCatch(
      pa_optim(pu, features = features, area_budget = budget, lockin_col = lockin_col,
               boundary_penalty = boundary_penalty, edge_factor = edge_factor,
               boundary_data = bd, mip_gap = mip_gap, time_limit = time_limit, ...),
      error = function(e) { err <<- conditionMessage(e); NULL }
    )
    d <- attr(res, "diag")
    list(res = res, log = data.frame(
      run = run, budget_km2 = budget,
      available_km2 = if (is.null(d$available_km2)) NA_real_ else d$available_km2,
      selected_km2  = if (is.null(res)) 0 else sum(res$area_km2),
      n_cells       = if (is.null(res)) 0L else nrow(res),
      status        = if (is.null(d$status)) NA_character_ else as.character(d$status),
      runtime_s     = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      error         = err))
  }
  
  lulc_runs  <- lapply(lulcs, function(x) one(x, get_budget_km2(gap, x, budget_col), target_lulc = x))
  global_run <- if (run_global) one(global_name, get_budget_km2(gap, global_name, budget_col),
                                    lockout = global_lockout)
  
  res_l <- Filter(Negate(is.null), lapply(lulc_runs, `[[`, "res"))
  list(
    lulc   = if (length(res_l)) dplyr::bind_rows(res_l) else NULL,
    global = if (run_global) global_run$res else NULL,
    log    = dplyr::bind_rows(c(lapply(lulc_runs, `[[`, "log"),
                                if (run_global) list(global_run$log)))
  )
}

## ---- 4. Merge and classify ----------------------------------------------------------------
classify_pa <- function(pu, new_ids, id_col = "id") {
  new <- pu[[id_col]] %in% new_ids
  factor(
    dplyr::case_when(
      pu$existing_corePA %in% TRUE ~ "existing core pa",
      new & pu$n_pa == 0           ~ "new core PA",
      new & pu$n_pa > 0            ~ "proposed upgrade existing PA",
      pu$n_pa > 0                  ~ "other protected areas",
      TRUE                         ~ "not protected"
    ),
    levels = pa_levels
  )
}

add_solutions <- function(pu, runs, id_col = "id") {
  stopifnot(!anyDuplicated(pu[[id_col]]))
  lid <- if (is.null(runs$lulc))   integer(0) else runs$lulc[[id_col]]
  gid <- if (is.null(runs$global)) integer(0) else runs$global[[id_col]]
  pu$optim_PA_single_lulc <- pu[[id_col]] %in% lid
  pu$optim_PA_global      <- pu[[id_col]] %in% gid
  pu$core_pa_lulc         <- classify_pa(pu, lid, id_col)
  pu$core_pa_global       <- classify_pa(pu, gid, id_col)
  pu$is_core_prot_lulc    <- pu$core_pa_lulc   %in% pa_levels[1:3]
  pu$is_core_prot_global  <- pu$core_pa_global %in% pa_levels[1:3]
  pu
}

## ---- 5. Statistics and plots -----------------------------------------------------------------
summarise_optim <- function(pu,
                            scenario_cols = c(LULC = "core_pa_lulc", GLOBAL = "core_pa_global"),
                            stat_cols = NULL, group_col = "lulc_name") {
  d <- sf::st_drop_geometry(pu)                                  # drop geometry BEFORE grouping
  if (is.null(stat_cols)) {
    stat_cols <- c(grep("^sampled_.*_scaled$", names(d), value = TRUE), "min_distance")
  }
  stat_cols <- intersect(stat_cols, names(d))
  
  dplyr::bind_rows(lapply(names(scenario_cols), function(sc) {
    d$core_pa <- d[[scenario_cols[[sc]]]]
    d |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col)), core_pa) |>
      dplyr::summarise(
        n_cells = dplyr::n(),
        km2     = sum(area_km2),
        dplyr::across(dplyr::all_of(stat_cols), list(mean = mean, sd = stats::sd),
                      .names = "{.col}_{.fn}"),
        .groups = "drop"
      ) |>
      dplyr::mutate(scenario = sc)
  }))
}

plot_pa_map <- function(pu, class_col, stud_area, title = NULL,
                        built_up = "built-up", text_size = 15) {
  show <- c(pa_cols[1:3], "other protected areas" = "transparent", "not protected" = "transparent")
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = pu, ggplot2::aes(fill = .data[[class_col]]), colour = NA) +
    ggplot2::geom_sf(data = pu[pu$lulc_name %in% built_up, ], fill = "grey", colour = NA) +
    ggplot2::geom_sf(data = stud_area, fill = NA, colour = "black") +
    ggplot2::scale_fill_manual(values = show, breaks = names(pa_cols)[1:3], name = NULL) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", text = ggplot2::element_text(size = text_size))
}

plot_pa_area <- function(stats, exclude_lulc = c("built-up", "agricultural"),
                         classes = pa_levels[1:3]) {
  d <- stats[stats$core_pa %in% classes & !is.na(stats$lulc_name) &
               !stats$lulc_name %in% exclude_lulc, ]
  ggplot2::ggplot(d, ggplot2::aes(x = scenario, y = km2, fill = core_pa)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::facet_wrap(~ lulc_name, scales = "free_y") +
    ggplot2::scale_fill_manual(values = pa_cols, breaks = classes) +
    ggplot2::labs(x = NULL, y = expression("Area (km"^2*")"), fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top", strip.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

## ---- 6. Pipeline --------------------------------------------------------------------------------
run_pa_optimization <- function(
    target_site,
    scenario         = "nature",
    feature_sets     = default_feature_sets,
    pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
    gap_dir          = "outputs/WP4/01_PA_analysis",
    pu_dir           = "outputs/WP4/02_optim",
    shared_dir       = "data/shared",
    out_dir          = "outputs/WP4/02_optim",
    budget_col       = "target_area",           # or "gap_area"
    dist_transform   = c("linear", "log", "rank"),
    min_area_km2     = 0.075, keep_core = TRUE,
    area_weighted    = FALSE,
    lulcs            = NULL, run_global = TRUE,
    global_lockout   = c("built-up", "agricultural"),
    boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
    mip_gap = 0.01, time_limit = 600,
    write            = TRUE
) {
  dist_transform <- match.arg(dist_transform)
  if (!scenario %in% names(feature_sets)) {
    stop("Unknown scenario '", scenario, "'. Available: ", paste(names(feature_sets), collapse = ", "))
  }
  
  if (is.null(gap)) gap <- utils::read.csv(file.path(gap_dir, paste0(target_site, "_gap_analysis.csv")))
  if (is.null(pu))  pu  <- sf::st_read(file.path(pu_dir, paste0(target_site, "_input_final_grid.json")), quiet = TRUE)
  if (is.null(stud_area)) stud_area <- sf::read_sf(file.path(shared_dir, paste0(target_site, ".gpkg")))
  stud_area <- sf::st_transform(stud_area, sf::st_crs(pu))
  if (!"area_km2" %in% names(pu)) pu$area_km2 <- pu$area / 1e6
  stopifnot(!anyDuplicated(pu$id), all(c("existing_corePA", "n_pa", "lulc_name") %in% names(pu)))
  
  pu <- add_distance_scaled(pu, dist_transform)       # scaled on the full grid, before any subsetting
  pu <- drop_partial_pus(pu, min_area_km2, keep_core)
  fb <- build_features(pu, feature_sets[[scenario]], area_weighted)
  pu <- fb$pu
  
  runs <- run_optimizations(
    pu, gap, fb$features, lulcs = lulcs, run_global = run_global,
    global_lockout = global_lockout, budget_col = budget_col,
    boundary_penalty = boundary_penalty, edge_factor = edge_factor,
    share_boundary = share_boundary, mip_gap = mip_gap, time_limit = time_limit
  )
  pu    <- add_solutions(pu, runs)
  stats <- summarise_optim(pu)
  
  plots <- list(
    map_lulc   = plot_pa_map(pu, "core_pa_lulc",   stud_area, "LULC-specific targets"),
    map_global = plot_pa_map(pu, "core_pa_global", stud_area, "Global target"),
    area       = plot_pa_area(stats)
  )
  
  if (write) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stub <- file.path(out_dir, paste0(target_site, "_", scenario))
    write_grid(pu, out_dir, paste0("PA_optim_grid_", target_site, "_", scenario), "geojson")
    utils::write.csv(stats,    paste0(stub, "_optim_stats.csv"), row.names = FALSE)
    utils::write.csv(runs$log, paste0(stub, "_optim_log.csv"),   row.names = FALSE)
    ggplot2::ggsave(paste0(stub, "_pa_lulc_optim.png"),     plots$map_lulc,   width = 18, height = 10, dpi = 300)
    ggplot2::ggsave(paste0(stub, "_pa_global_optim.png"),   plots$map_global, width = 18, height = 10, dpi = 300)
    ggplot2::ggsave(paste0(stub, "_pa_area_statsProt.png"), plots$area,       width = 13, height = 10, dpi = 300)
  }
  
  list(pu = pu, runs = runs, log = runs$log, stats = stats, plots = plots, features = fb$features)
}