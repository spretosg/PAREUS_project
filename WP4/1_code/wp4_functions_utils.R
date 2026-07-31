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


pa_optim<-function(pu,
                   target_lulc = NULL,
                   lockout = NULL,
                   lockin_col,
                   features,
                   area_budget
){
  #filter out lockout e.g. built-up
  if (!is.null(lockout)) {
    pu <- pu %>%
      filter(!(.data[["lulc_name"]] %in%lockout))
  }
  #subset
  if(!is.null(target_lulc)){
    pu<-pu %>% filter(.data[["lulc_name"]] == target_lulc)
  }
  
  
  pa_problem<- problem(
    pu,
    features = features,
    cost_column = "area_km2" # Cost is set to area to enforce the area budget
  ) %>%
    # Maximizes feature values strictly within your 30% area budget
    add_max_utility_objective(budget = area_budget) %>%
    # add_locked_out_constraints("lock_out")%>%
    add_boundary_penalties(penalty = 0.0005) %>% 
    
    # Ensures planning units are either completely selected (1) or not (0)
    add_binary_decisions()
  
  #in case already protected area present, add lockin
  if(nrow(pu%>%filter(.data[[lockin_col]]==T)) !=0){
    
    pa_problem<-pa_problem %>% add_locked_in_constraints(lockin_col)
  }
  
  optim_pa <- solve(pa_problem) 
  # plot(
  #   st_as_sf(optim_pa[, "solution_1"]), main = "Prioritization",
  #   pal = c("grey90", "darkgreen")
  # )

  return(optim_pa%>%filter(solution_1==1))
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


## oecm selection function either global % highest oecm suitability, or per LULC % of highest suitability
select_oecm <- function(pu,
                        mode = c("class", "global"),
                        coverage,
                        lulc_col = "lulc",
                        suitability_col = "oecm_suitability",
                        corePA_col = "core_pa_lulc",
                        search_oecm_in = c("not protected","other protected areas")) {
  
  mode <- match.arg(mode)
  
  dat <- pu
  

  dat <- dat %>% filter(.data[[corePA_col]] %in% search_oecm_in )
  
  ## ---------- Global selection ----------
  if (mode == "global") {
    
    stopifnot(length(coverage) == 1)
    
    target_area <- sum(dat$area) * coverage
    
    return(
      dat %>%
        arrange(desc(.data[[suitability_col]])) %>%
        mutate(cum_area = cumsum(area)) %>%
        slice(seq_len(which(cum_area >= target_area)[1])) %>%
        dplyr::select(-cum_area)
    )
  }
  
  ## ---------- Class-specific selection ----------
  stopifnot(!is.null(names(coverage)))
  
  map_dfr(names(coverage), function(cls) {
    
    dat_cls <- dat %>%
      filter(.data[[lulc_col]] == cls) %>%
      arrange(desc(.data[[suitability_col]]))
    
    if (nrow(dat_cls) == 0)
      return(dat_cls)
    
    target_area <- sum(dat_cls$area) * coverage[[cls]]
    
    dat_cls %>%
      mutate(cum_area = cumsum(area)) %>%
      slice(seq_len(which(cum_area >= target_area)[1]))
  }) %>%
    dplyr::select( -cum_area)
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
  pu_pca<-pu_pca%>%select(id,lulc_name, sampled_es_reg_scaled,sampled_es_prov_scaled,sampled_es_cult_scaled,sampled_grand_mean_es_scaled, sampled_eco_cond_scaled, sampled_cost_policy,connectivity_scaled,area,pca)
  
  if(save_output == T){
    st_write(pu_pca,paste0("WP4/2_output/03_pca_landscape/PCA_PU_",siteID,"_",scen,".geojson"))
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
    
    new_scale_fill() +

    # stud area
    geom_sf(data = stud_area, fill = NA, color = "black") +
    theme_minimal()+
    ggtitle(scen)
  
  list(
    plot = p,
    stats = stats,
    total_area = total,
    pu = pu_pca
  )
  
}
