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

### LCP on OECM suitability
lcp_connect <- function(
    pu,                         # sf: ALL planning units, with a suitability column
    sel,                        # logical vector (length nrow(pu)): selected cells to connect
    suit_col,                   # suitability column (high = easy to move through)
    topology         = c("mst", "nearest", "all"),
    k                = 1,       # neighbours per patch if topology = "nearest"
    transform        = c("linear", "inverse", "exp"),
    max_resistance   = 50,      # resistance of the least suitable cell (best cell = 1)
    contiguity       = c("queen", "rook"),
    exclude          = NULL,    # logical vector: cells that cannot be crossed (e.g. built-up)
    patch_resistance = NULL,    # e.g. 1 to make paths prefer running through selected cells
    min_patch_km2    = 0,
    id_col           = "id"
) {
  topology   <- match.arg(topology)
  transform  <- match.arg(transform)
  contiguity <- match.arg(contiguity)
  n <- nrow(pu)
  
  stopifnot(
    inherits(pu, "sf"), !sf::st_is_longlat(pu),
    is.logical(sel), length(sel) == n, !anyNA(sel), sum(sel) > 1,
    suit_col %in% names(pu), is.numeric(pu[[suit_col]]), 
    is.numeric(max_resistance), max_resistance > 1
  )
  if (is.null(exclude)) exclude <- rep(FALSE, n)
  stopifnot(is.logical(exclude), length(exclude) == n, !anyNA(exclude))
  exclude <- exclude & !sel          # selected cells are always traversable
  
  # 1. Suitability -> resistance, rescaled over ALL pu so it does not depend on the selection
  s <- pu[[suit_col]]
  if (anyNA(s[!exclude])) {
    stop(sum(is.na(s[!exclude])), " non-excluded cells have NA in '", suit_col,
         "'. Fill them (fill_suit_na) or pass them in `exclude`.")
  }
  s[is.na(s)] <- median(s, na.rm = TRUE)
  if (diff(range(s)) == 0) stop("suit_col is constant")
  s01 <- (s - min(s)) / (max(s) - min(s))
  res <- switch(
    transform,
    linear  = 1 + (1 - s01) * (max_resistance - 1),
    inverse = 1 / (s01 * (1 - 1 / max_resistance) + 1 / max_resistance),
    exp     = exp((1 - s01) * log(max_resistance))
  )
  if (!is.null(patch_resistance)) res[sel] <- patch_resistance
  
  # 2. Cell graph: nodes = pu rows, edges = neighbours, weight = distance x mean resistance
  pts      <- sf::st_point_on_surface(sf::st_geometry(pu))
  xy       <- sf::st_coordinates(pts)
  area_km2 <- as.numeric(sf::st_area(pu)) / 1e6
  
  nb <- if (contiguity == "queen") sf::st_touches(pu) else
    sf::st_relate(pu, pu, pattern = "F***1****")
  
  from <- rep(seq_along(nb), lengths(nb))
  to   <- unlist(nb)
  ok   <- from < to & !exclude[from] & !exclude[to]
  from <- from[ok]; to <- to[ok]
  
  d <- sqrt(rowSums((xy[from, , drop = FALSE] - xy[to, , drop = FALSE])^2)) / 1000  # km
  w <- d * (res[from] + res[to]) / 2
  
  g <- igraph::make_undirected_graph(edges = c(rbind(from, to)), n = n)
  igraph::E(g)$w <- w
  w_all <- igraph::E(g)$w
  
  # 3. Patches = connected components of the selected cells
  both  <- sel[from] & sel[to]
  gs    <- igraph::make_undirected_graph(edges = c(rbind(from[both], to[both])), n = n)
  comp  <- igraph::components(gs)$membership
  patch_of <- ifelse(sel, comp, NA_integer_)
  pid   <- sort(unique(patch_of[sel]))
  
  patch_area <- vapply(pid, function(p) sum(area_km2[which(patch_of == p)]), numeric(1))
  # Representative cell: the selected cell closest to the area-weighted patch centre,
  # so it always lies inside the patch (even for concave patches)
  start <- vapply(pid, function(p) {
    idx <- which(patch_of == p)
    a   <- area_km2[idx]
    ctr <- c(sum(xy[idx, 1] * a), sum(xy[idx, 2] * a)) / sum(a)
    idx[which.min((xy[idx, 1] - ctr[1])^2 + (xy[idx, 2] - ctr[2])^2)]
  }, integer(1))
  
  keep  <- patch_area >= min_patch_km2
  pid <- pid[keep]; patch_area <- patch_area[keep]; start <- start[keep]
  if (length(start) < 2) stop("Need at least two patches to connect")
  
  # 4. Least-cost distances between patch representatives
  D <- igraph::distances(g, v = start, to = start, weights = w_all)
  D[!is.finite(D)] <- NA
  
  edges <- switch(
    topology,
    all = {
      ij <- which(upper.tri(D), arr.ind = TRUE)
      ij[!is.na(D[ij]), , drop = FALSE]
    },
    nearest = {
      do.call(rbind, lapply(seq_len(nrow(D)), function(i) {
        di <- D[i, ]; di[i] <- NA
        j  <- order(di)[seq_len(min(k, sum(!is.na(di))))]
        if (length(j)) cbind(i, j)
      }))
    },
    mst = {
      W <- D; W[is.na(W)] <- 0; diag(W) <- 0          # 0 = no edge
      g2 <- igraph::graph_from_adjacency_matrix(W, mode = "max", weighted = TRUE, diag = FALSE)
      igraph::as_edgelist(igraph::mst(g2), names = FALSE)
    }
  )
  if (is.null(edges) || !nrow(edges)) stop("No connectable patch pairs (check `exclude`)")
  edges <- unique(t(apply(edges, 1, sort)))
  if (!is.matrix(edges)) edges <- matrix(edges, ncol = 2)
  
  # 5. Paths as cell sequences -> sf lines
  paths <- lapply(seq_len(nrow(edges)), function(e) {
    as.integer(igraph::shortest_paths(g, from = start[edges[e, 1]], to = start[edges[e, 2]],
                                      weights = w_all, output = "vpath")$vpath[[1]])
  })
  
  geom <- sf::st_sfc(lapply(paths, function(ids) sf::st_linestring(xy[ids, , drop = FALSE])),
                     crs = sf::st_crs(pu))
  lines <- sf::st_sf(
    from    = pid[edges[, 1]],
    to      = pid[edges[, 2]],
    cost    = D[edges],
    n_cells = lengths(paths),
    geometry = geom
  )
  lines$length_km <- as.numeric(sf::st_length(lines)) / 1000
  
  # Non-selected cells the paths run through (the corridor itself)
  used     <- unique(unlist(paths))
  corridor <- used[!sel[used]]
  corridor_ids <- if (id_col %in% names(pu)) pu[[id_col]][corridor] else corridor
  
  nodes <- sf::st_sf(patch_id = pid, area_km2 = patch_area, geometry = pts[start])
  
  list(paths = lines, nodes = nodes, corridor_ids = corridor_ids,
       cost_matrix = D, resistance = res)
}


### helper function to impute na in oecm suitability
fill_suit_na <- function(pu, suit_col,
                         method   = c("neighbours", "idw"),
                         k        = 8,      # idw: number of nearest valid cells
                         power    = 2,      # idw: distance exponent
                         max_iter = 100,    # neighbours: max passes
                         new_col  = paste0(suit_col, "_filled")) {
  method <- match.arg(method)
  s  <- pu[[suit_col]]
  na <- is.na(s)
  
  if (any(na)) {
    xy <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(pu)))
    
    if (method == "idw") {
      # inverse-distance weighting from the k nearest valid cells (needs RANN)
      nn   <- RANN::nn2(xy[!na, , drop = FALSE], xy[na, , drop = FALSE],
                        k = min(k, sum(!na)))
      wts  <- 1 / pmax(nn$nn.dists, 1)^power                # distances in metres
      vals <- matrix(s[!na][nn$nn.idx], nrow = nrow(nn$nn.idx))
      s[na] <- rowSums(wts * vals) / rowSums(wts)
      
    } else {
      # iterative neighbour mean: each pass fills cells that touch a valid cell
      nb <- sf::st_touches(pu)
      for (it in seq_len(max_iter)) {
        miss <- which(is.na(s))
        if (!length(miss)) break
        new <- vapply(miss, function(i) {
          v <- s[nb[[i]]]
          if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
        }, numeric(1))
        s[miss] <- new                                       # simultaneous update
      }
      if (anyNA(s)) {
        warning(sum(is.na(s)), " cells still NA after ", max_iter,
                " passes (isolated from valid cells); filled with the median")
        s[is.na(s)] <- median(s, na.rm = TRUE)
      }
    }
  }
  
  pu[[new_col]]       <- s
  pu[[paste0(new_col, "_imputed")]] <- na
  pu
}
