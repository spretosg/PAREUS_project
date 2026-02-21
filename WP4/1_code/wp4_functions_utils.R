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