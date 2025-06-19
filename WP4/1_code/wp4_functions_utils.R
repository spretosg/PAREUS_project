min_max_normalize <- function(r) {
  # Normalize each layer individually
  norm_r <- lapp(r, function(x) {
    r_min <- min(x, na.rm = TRUE)
    r_max <- max(x, na.rm = TRUE)
    (x - r_min) / (r_max - r_min)
  })
  return(norm_r)
}