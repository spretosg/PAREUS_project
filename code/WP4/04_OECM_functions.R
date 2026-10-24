### OECM functions

oecm_lin_w<-function(cult_es,prov_es,ec,conn,w_cult,w_prov,w_ec,w_connect){
  weighted_sum <- cult_es * w_cult +
    prov_es * w_prov +
    ec * w_ec +
    conn * w_connect
  
  weight_total <- w_cult + w_prov + w_ec + w_connect
  
  return(weighted_sum / weight_total)
  
}

plot_area_pies <- function(
    stats,
    slice_col,                       # column defining the slices (e.g. "pca" or "lulc_name")
    value_col       = "AREA_SQKM",
    scenario_col    = "scenario",
    facet_col       = NULL,          # optional: one row of pies per level (e.g. "lulc_name")
    keep_slices     = NULL,          # optional: show only these slice values
    slice_levels    = NULL,          # slice order (default: largest first)
    scenario_levels = NULL,
    fill_values     = NULL,          # named colours; default viridis
    label           = c("km2", "pct", "both"),
    min_label_pct   = 1,             # no label for slices below this share
    digits          = 1,
    unit            = "km²",
    legend_title    = NULL,
    subtitle        = "Pie area proportional to total area",
    base_size       = 14,
    pad             = 2.2            # horizontal room for labels
) {
  label <- match.arg(label)
  d <- as.data.frame(stats)
  d$.scn   <- d[[scenario_col]]
  d$.slice <- d[[slice_col]]
  d$.val   <- d[[value_col]]
  d$.facet <- if (is.null(facet_col)) "" else d[[facet_col]]
  
  if (!is.null(keep_slices)) d <- d[d$.slice %in% keep_slices, ]
  d <- d[!is.na(d$.slice) & !is.na(d$.val) & d$.val > 0, ]
  if (!nrow(d)) stop("Nothing to plot after filtering")
  
  d <- d |>
    dplyr::group_by(.scn, .facet, .slice) |>
    dplyr::summarise(.val = sum(.val), .groups = "drop")
  
  if (is.null(slice_levels)) {
    slice_levels <- d |>
      dplyr::group_by(.slice) |>
      dplyr::summarise(t = sum(.val), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(t)) |>
      dplyr::pull(.slice) |> as.character()
  }
  d$.slice <- factor(d$.slice, levels = slice_levels)
  d <- d[!is.na(d$.slice), ]
  d$.scn <- factor(d$.scn, levels = if (is.null(scenario_levels)) unique(d$.scn) else scenario_levels)
  
  # angles per pie
  d <- d |>
    dplyr::arrange(.scn, .facet, .slice) |>
    dplyr::group_by(.scn, .facet) |>
    dplyr::mutate(total = sum(.val),
                  frac  = .val / total,
                  end   = 2 * pi * cumsum(frac),
                  start = end - 2 * pi * frac,
                  mid   = (start + end) / 2) |>
    dplyr::ungroup()
  d$end <- ifelse(d$frac > 0.9999, d$start + 2 * pi - 1e-4, d$end)   # single-slice pie
  d$r   <- sqrt(d$total / max(d$total))                              # same scale in all panels
  
  # strip labels: the total goes under the facet name
  fmt_total <- function(x) sprintf("%s %s total", format(round(x, 0), big.mark = ",", trim = TRUE), unit)
  
  if (is.null(facet_col)) {
    # one strip per scenario: "<scenario>\n<total>"
    lab <- paste0(d$.scn, "\n", fmt_total(d$total))
    d$.col_lab <- factor(lab, levels = unique(lab[order(d$.scn)]))
  } else {
    # columns = scenarios; rows = facet values with their total (one number if the
    # scenarios share it, otherwise "a | b | ...")
    d$.col_lab <- d$.scn
    row_tot <- d |>
      dplyr::distinct(.scn, .facet, total) |>
      dplyr::arrange(.facet, .scn) |>
      dplyr::group_by(.facet) |>
      dplyr::summarise(
        row_total = {
          u <- unique(round(total, 6))
          if (length(u) == 1) fmt_total(u) else paste(fmt_total(total), collapse = " | ")
        },
        .groups = "drop")
    d   <- dplyr::left_join(d, row_tot, by = ".facet")
    lab <- paste0(d$.facet, "\n", d$row_total)
    d$.row_lab <- factor(lab, levels = unique(lab[order(d$.facet)]))
  }
  
  # slice labels
  lab_df <- d[d$frac >= min_label_pct / 100, ]
  num <- format(round(lab_df$.val, digits), nsmall = digits, trim = TRUE)
  pct <- sprintf("%.0f%%", 100 * lab_df$frac)
  lab_df$lab <- switch(label,
                       km2  = paste0(num, " ", unit),
                       pct  = pct,
                       both = paste0(num, " ", unit, "\n(", pct, ")"))
  lab_df$lx <- 1.15 * lab_df$r * sin(lab_df$mid)
  lab_df$ly <- 1.15 * lab_df$r * cos(lab_df$mid)
  lab_df$hj <- ifelse(sin(lab_df$mid) >= 0, 0, 1)
  
  p <- ggplot2::ggplot() +
    ggforce::geom_arc_bar(
      data = d,
      ggplot2::aes(x0 = 0, y0 = 0, r0 = 0, r = r, start = start, end = end, fill = .slice),
      colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(data = lab_df, ggplot2::aes(lx, ly, label = lab, hjust = hj),
                       size = 3.5, lineheight = 0.9) +
    (if (is.null(fill_values)) ggplot2::scale_fill_viridis_d(drop = FALSE)
     else ggplot2::scale_fill_manual(values = fill_values, drop = FALSE)) +
    ggplot2::coord_fixed(xlim = c(-pad, pad), ylim = c(-1.4, 1.4), clip = "off") +
    ggplot2::labs(fill = legend_title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text   = ggplot2::element_text(face = "bold", lineheight = 0.95),
                   strip.text.y = ggplot2::element_text(face = "bold", lineheight = 0.95, angle = 0),
                   plot.margin  = ggplot2::margin(10, 20, 10, 20))
  
  if (is.null(facet_col)) {
    p + ggplot2::facet_wrap(ggplot2::vars(.col_lab))
  } else {
    p + ggplot2::facet_grid(rows = ggplot2::vars(.row_lab), cols = ggplot2::vars(.col_lab))
  }
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
  
  if(save_output == T){
    st_write(pu_pca,paste0("outputs/WP4/03_pca_landscape/PCA_PU_",stud_site,"_",scen,".geojson"))
    ggsave(paste0("outputs/WP4/03_pca_landscape/PCA_PU_",stud_site,"_",scen,".png"),p)
  }
  
  list(
    plot = p,
    stats = stats,
    total_area = total,
    pu = pu_pca
  )
  
}


# pu-level data -> long area table (generalises corridor_lulc_stats)
summarise_area <- function(pu, group_col = "lulc_name", flag_col = NULL,
                           scenario, area_col = "area_km2") {
  d <- sf::st_drop_geometry(pu)
  if (!is.null(flag_col)) d <- d[d[[flag_col]] %in% TRUE, , drop = FALSE]
  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::summarise(km2 = sum(.data[[area_col]]), n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(pct = 100 * km2 / sum(km2), scenario = scenario)
}

# comparison table: classes x scenarios, with a TOTAL row
area_table <- function(stats, slice_col, value_col = "km2", scenario_col = "scenario") {
  wide <- stats |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(slice_col, scenario_col)))) |>
    dplyr::summarise(v = sum(.data[[value_col]]), .groups = "drop") |>
    tidyr::pivot_wider(names_from = dplyr::all_of(scenario_col), values_from = v,
                       values_fill = 0)
  tot <- wide |> dplyr::summarise(dplyr::across(where(is.numeric), sum))
  tot[[slice_col]] <- "TOTAL"
  dplyr::bind_rows(wide, tot)
}

# pie chart: one pie per scenario (optionally one row per facet_col),
# pie AREA proportional to the total area of the displayed slices
# Needs: sample_rasters() (step 2), fill_suit_na(), lcp_connect() (patched version),
#        summarise_area(), write_grid(), zero_one_scale(), oecm_lin_w()

corridor_map_cols <- c(
  "existing core pa"             = "#00A300",
  "new core PA"                  = "#bf8bff",
  "proposed upgrade existing PA" = "#e5d0ff",
  "least-cost corridor"          = "orange",
  "other protected areas"        = "#F1EE8D",
  "not protected"                = "transparent"
)

plot_corridor_map <- function(pu, stud_area, lcp = NULL, map_col = ".map_class",
                              built_up = "built-up", show_paths = FALSE, show_nodes = TRUE,
                              cols = corridor_map_cols, title = NULL, text_size = 14) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = pu, ggplot2::aes(fill = .data[[map_col]]), colour = NA) +
    ggplot2::geom_sf(data = pu[pu$lulc_name %in% built_up, ], fill = "grey", colour = NA)
  if (show_paths && !is.null(lcp)) p <- p + ggplot2::geom_sf(data = lcp$paths, colour = "red", linewidth = 0.6)
  if (show_nodes && !is.null(lcp)) p <- p + ggplot2::geom_sf(data = lcp$nodes, size = 1.5)
  p +
    ggplot2::geom_sf(data = stud_area, fill = NA, colour = "black") +
    ggplot2::scale_fill_manual(values = cols, breaks = setdiff(names(cols), "not protected"), name = NULL) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(legend.position = "top")
}

run_pca_corridors <- function(
    site,
    scenario        = c("agri", "global"),
    weights         = list(w_cult = 0.4, w_prov = 0.4, w_reg = 0.1, w_connect = 0.1),  # human centric
    pu              = NULL,                       # optional: sf instead of reading the file
    stud_area       = NULL,
    optim_dir       = "outputs/WP4/02_optim",
    shared_dir      = "data/shared",
    out_dir         = "outputs/WP4/03_pca_landscape",
    conn_layer      = "normalized_current",
    class_col       = "core_pa_global",           # same column name in the agri and global files
    core_classes    = c("existing core pa", "new core PA", "proposed upgrade existing PA"),
    barrier_lulc    = "built-up",
    transform       = "linear",
    max_resistance  = 50,
    topology        = "mst",
    scenario_labels = c(agri   = "Agricultural areas as PA possible",
                        global = "No PA in agricultural areas"),
    show_paths      = FALSE,
    write           = TRUE,
    write_pu        = FALSE,
    ...                                           # further arguments to lcp_connect()
) {
  scenario <- match.arg(scenario)
  stopifnot(is.function(zero_one_scale), is.function(oecm_lin_w), is.function(fill_suit_na))
  if (!setequal(names(weights), c("w_cult", "w_prov", "w_reg", "w_connect"))) {
    stop("weights needs w_cult, w_prov, w_reg and w_connect")
  }
  if (abs(sum(unlist(weights)) - 1) > 1e-6) warning("weights do not sum to 1")
  
  # 1. Inputs
  if (is.null(pu)) {
    pu <- sf::st_read(file.path(optim_dir, sprintf("PA_optim_grid_%s_%s.geojson", site, scenario)),
                      quiet = TRUE)
  }
  if (is.null(stud_area)) stud_area <- sf::read_sf(file.path(shared_dir, paste0(site, ".gpkg")))
  stud_area <- sf::st_transform(stud_area, sf::st_crs(pu))
  if (!"area_km2" %in% names(pu)) pu$area_km2 <- pu$area / 1e6
  stopifnot(all(c("id", "lulc_name", class_col) %in% names(pu)), !anyDuplicated(pu$id))
  
  # 2. Connectivity -> PUs (sampled in the raster's own CRS), scaled 0-1
  conn <- terra::rast(file.path(optim_dir, paste0("conn_", site, ".tif")))
  if (!conn_layer %in% names(conn)) {
    stop("Layer '", conn_layer, "' not in connectivity raster. Available: ", paste(names(conn), collapse = ", "))
  }
  pu <- sample_rasters(pu, list(connectivity = conn[[conn_layer]]), prefix = "")
  pu <- zero_one_scale(pu, cols = "connectivity")            # -> connectivity_scaled
  
  # 3. OECM suitability, NA gaps filled from neighbours
  req <- c("sampled_es_cult_scaled", "sampled_es_prov_scaled", "sampled_es_reg_scaled", "connectivity_scaled")
  if (!all(req %in% names(pu))) stop("Missing columns: ", paste(setdiff(req, names(pu)), collapse = ", "))
  pu$oecm_suit <- oecm_lin_w(pu$sampled_es_cult_scaled, pu$sampled_es_prov_scaled,
                             pu$sampled_es_reg_scaled, pu$connectivity_scaled,
                             weights$w_cult, weights$w_prov, weights$w_reg, weights$w_connect)
  pu <- fill_suit_na(pu, "oecm_suit", method = "neighbours")  # -> oecm_suit_filled
  
  # 4. Least-cost corridors between the optimised core PAs
  sel <- pu[[class_col]] %in% core_classes
  lcp <- lcp_connect(pu, sel = sel, suit_col = "oecm_suit_filled",
                     exclude = pu$lulc_name %in% barrier_lulc,
                     transform = transform, max_resistance = max_resistance,
                     topology = topology, ...)
  
  corridor_col <- paste0("is_corridor_", scenario)
  pu[[corridor_col]] <- pu$id %in% lcp$corridor_ids
  
  # 5. Map classes (core stays core, corridor takes precedence over the rest)
  cls <- as.character(pu[[class_col]])
  cls[pu[[corridor_col]] & !(cls %in% core_classes)] <- "least-cost corridor"
  pu$.map_class <- factor(cls, levels = names(corridor_map_cols))
  
  label <- unname(scenario_labels[[scenario]])
  map <- plot_corridor_map(pu, stud_area, lcp, show_paths = show_paths, title = label)
  
  # 6. Statistics
  stats <- summarise_area(pu, "lulc_name", corridor_col, label)
  
  if (write) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(out_dir, sprintf("%s_PCA_%s.png", site, scenario)), map,
                    width = 14, height = 12, dpi = 300)
    if (write_pu) write_grid(pu, out_dir, sprintf("PCA_corridor_grid_%s_%s", site, scenario), "geojson")
  }
  
  list(pu = pu, lcp = lcp, stats = stats, map = map, scenario = scenario, label = label,
       corridor_km2 = sum(pu$area_km2[pu[[corridor_col]]]),
       path_km = sum(lcp$paths$length_km))
}
