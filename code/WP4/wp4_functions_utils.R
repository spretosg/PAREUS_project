## calculate core protection gap


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

min_max_normalize <- function(r) {
  # Normalize each layer individually
  norm_r <- lapp(r, function(x) {
    r_min <- min(x, na.rm = TRUE)
    r_max <- max(x, na.rm = TRUE)
    (x - r_min) / (r_max - r_min)
  })
  return(norm_r)
  
}


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

oecm_lin_w<-function(cult_es,prov_es,ec,conn,w_cult,w_prov,w_ec,w_connect){
  weighted_sum <- cult_es * w_cult +
    prov_es * w_prov +
    ec * w_ec +
    conn * w_connect
  
  weight_total <- w_cult + w_prov + w_ec + w_connect
  
  return(weighted_sum / weight_total)
  
}


smooth_suitability <- function(dat, suitability_col, w = 0.5, order = 1) {
  nb <- st_touches(dat)                      # neighbours among eligible PUs
  s  <- dat[[suitability_col]]
  for (k in seq_len(order)) {                # order > 1 = wider smoothing
    nb_mean <- vapply(seq_along(nb), function(i)
      if (length(nb[[i]])) mean(s[nb[[i]]]) else s[i], numeric(1))
    s <- (1 - w) * s + w * nb_mean
  }
  s
}

## oecm selection function either global % highest oecm suitability, or per LULC % of highest suitability
select_oecm <- function(pu, mode = c("class", "global"), coverage,
                        lulc_col = "lulc",
                        suitability_col = "oecm_suitability",
                        corePA_col = "core_pa_lulc",
                        search_oecm_in = c("not protected", "other protected areas"),
                        w = 0.7, order = 1) {
  mode <- match.arg(mode)
  
  dat <- pu %>% filter(.data[[corePA_col]] %in% search_oecm_in)
  dat$.score <- smooth_suitability(dat, suitability_col, w, order)
  
  pick <- function(d, cov) {
    if (nrow(d) == 0) return(d)
    target <- sum(d$area) * cov
    d <- d %>% arrange(desc(.score)) %>% mutate(cum_area = cumsum(area))
    d %>% slice(seq_len(which(cum_area >= target)[1])) %>% select(-cum_area)
  }
  
  if (mode == "global") {
    stopifnot(length(coverage) == 1)
    return(pick(dat, coverage) %>% select(-.score))
  }
  
  stopifnot(!is.null(names(coverage)))
  map_dfr(names(coverage), function(cls)
    pick(filter(dat, .data[[lulc_col]] == cls), coverage[[cls]])) %>%
    select(-.score)
}


plot_pca_map<-function(pu,
                       corePA_col = "core_pa_lulc",
                       oecm_df,
                       stud_area,
                       scen,
                       save_output = F
){
  dat <-pu
  oecm<-oecm_df
  core_pa <- dat %>% filter(!.data[[corePA_col]] %in% c("not protected","other protected areas"))%>%mutate(pca = .data[[corePA_col]])
  
  other_pa<-dat %>% filter(.data[[corePA_col]]=="other protected areas")%>%mutate(pca = "Other PA (IUCN III-VI)")

  # oecm$pca <- ifelse(oecm$n_pa > 0, "Other PA (IUCN III-VI) high suitability for OECM", "High OECM suitability not protected - potential OECM")
  oecm$pca <- ifelse(
    oecm[[corePA_col]] == "other protected areas",
    "Other PA (IUCN III-VI) high suitability for OECM",
    "High OECM suitability not protected - potential OECM"
  )

  ## other pa, remove ids from OECM
  no_oecm_pa <- other_pa[!other_pa$id %in% oecm$id, ]
  
  other<-dat %>% filter(.data[[corePA_col]] == "not protected")%>%mutate(pca = "not protected")
  #also here remove OECMs
  other <- other[!other$id %in% oecm$id, ]

  
  pu_pca<-rbind(other,oecm,no_oecm_pa,core_pa)
  pu_pca<-pu_pca
  
  if(save_output == T){
    st_write(pu_pca,paste0("outputs/WP4/03_pca_landscape/PCA_PU_",target_site,"_",scen,".geojson"))
  }
  
  stats<-pu_pca%>%group_by(pca,lulc_name)%>%summarise(area = sum(area))%>%st_drop_geometry()%>%mutate(scenario = scen)

  total<-dat%>%group_by(lulc_name)%>%summarise(area = sum(area))%>%st_drop_geometry()%>%mutate(scenario = scen)
  
  ############
  p<-ggplot() +
    # Layer 0: other pa
    geom_sf(data = pu_pca%>%filter(pca != "not protected"),
            aes(fill = pca),
            color = NA) +
    scale_fill_manual(
      values = c(    "existing core pa" = "#00A300",
                     "new core PA" = "#bf8bff",
                     "proposed upgrade existing PA" = "#e5d0ff",
                     "Other PA (IUCN III-VI) high suitability for OECM" = "#44a6c6",
                 "High OECM suitability not protected - potential OECM" = "#194553",
                 "Other PA (IUCN III-VI)" = "#ADD8E6"),
      name = NULL
    )  +
    geom_sf(
      data = pu%>%filter(lulc_name == "built-up"),
      aes(fill = lulc_name),
      color = "grey", fill="grey"  )+
    
    new_scale_fill() +
    # stud area
    geom_sf(data = stud_area, fill = NA, color = "black") +
    theme_minimal()+
    
    theme(legend.position="bottom")+
    ggtitle(scen)
  
  list(
    plot = p,
    stats = stats,
    total_area = total,
    pu = pu_pca
  )
  
}
