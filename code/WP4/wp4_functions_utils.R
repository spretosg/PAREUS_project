## calculate core protection gap





# pa_optim <- function(
#     pu,
#     target_lulc = NULL,
#     lockout = NULL,
#     lockin_col,
#     features,
#     area_budget
# ) {
#   
#   # Filter out lockout classes
#   if (!is.null(lockout)) {
#     pu <- pu %>%
#       filter(!(.data$lulc_name %in% lockout))
#   }
#   
#   # Filter target LULC
#   if (!is.null(target_lulc)) {
#     pu <- pu %>%
#       filter(.data$lulc_name == target_lulc)
#   }
#   
#   # Nothing to optimize
#   if (nrow(pu) == 0) {
#     return(NULL)
#   }
#   
#   pa_problem <- problem(
#     pu,
#     features = features,
#     cost_column = "area_km2"
#   ) %>%
#     add_max_utility_objective(budget = area_budget) %>%
#     add_boundary_penalties(penalty = 0.0005) %>%
#     add_binary_decisions()
#   
#   if (lockin_col %in% names(pu) &&
#       any(pu[[lockin_col]], na.rm = TRUE)) {
#     pa_problem <- pa_problem %>%
#       add_locked_in_constraints(lockin_col)
#   }
#   
#   solve(pa_problem) %>%
#     filter(solution_1 == 1)
# }

####V2 optim
pa_optim <- function(
    pu,
    features,
    area_budget,
    lockin_col       = NULL,
    target_lulc      = NULL,
    lockout          = NULL,
    boundary_penalty,
    gap              = 0.01,
    time_limit       = 600
) {
  stopifnot(
    inherits(pu, "sf"),
    is.numeric(area_budget), area_budget > 0,
    all(c("area_km2", "lulc_name", features) %in% names(pu))
  )
  if (anyNA(sf::st_drop_geometry(pu)[features])) stop("NA in feature columns")
  
  # Lock-in flag (existing PAs)
  if (!is.null(lockin_col)) {
    if (!lockin_col %in% names(pu)) stop("lockin_col not found in pu")
    if (anyNA(pu[[lockin_col]]))    stop("lockin_col contains NA")
    locked <- as.logical(pu[[lockin_col]])
  } else {
    locked <- rep(FALSE, nrow(pu))
  }
  
  # Filters apply to candidate cells only; existing PAs are kept
  # so that boundary penalties still see them as neighbours
  cand <- rep(TRUE, nrow(pu))
  if (!is.null(lockout))     cand <- cand & !(pu$lulc_name %in% lockout)
  if (!is.null(target_lulc)) cand <- cand &  (pu$lulc_name %in% target_lulc)
  
  keep   <- cand | locked
  pu     <- pu[keep, ]
  locked <- locked[keep]
  
  if (!any(!locked)) return(NULL)  # nothing left to optimise
  
  # Neutralise existing PAs: no cost, no feature contribution
  pu$cost_km2 <- ifelse(locked, 0, pu$area_km2)
  for (f in features) pu[[f]][locked] <- 0
  pu$.is_lockin <- locked
  
  p <- prioritizr::problem(pu, features = features, cost_column = "cost_km2") |>
    prioritizr::add_max_utility_objective(budget = area_budget) |>
    prioritizr::add_binary_decisions() |>
    prioritizr::add_default_solver(gap = gap, time_limit = time_limit,
                                   verbose = FALSE)
  
  if (boundary_penalty > 0) {
    p <- prioritizr::add_boundary_penalties(p, penalty = boundary_penalty)
  }
  if (any(locked)) {
    p <- prioritizr::add_locked_in_constraints(p, ".is_lockin")
  }
  
  s <- solve(p)
  
  # Return only NEW selections
  s[s$solution_1 == 1 & !s$.is_lockin, ]
}

###
classify_pa <- function(pu, new_ids) {
  new <- pu$id %in% new_ids
  factor(
    dplyr::case_when(
      pu$existing_corePA     ~ "existing core pa",
      new & pu$n_pa == 0     ~ "new core PA",
      new & pu$n_pa  > 0     ~ "proposed upgrade existing PA",
      pu$n_pa > 0            ~ "other protected areas",
      TRUE                   ~ "not protected"
    ),
    levels = pa_levels
  )
}








