#AHP functions
convert_saaty <- function(x) {
  if(is.na(x)) return(NA)
  if(x < 0) return(1 / abs(x))
  return(x)
}

build_ahp_matrix <- function(df_row, atts, negconvert = TRUE) {
  n <- length(atts)
  mat <- matrix(1, n, n)
  colnames(mat) <- atts
  rownames(mat) <- atts
  
  for(col in colnames(df_row)) {
    val <- as.numeric(df_row[[col]])
    pair <- strsplit(col, "_")[[1]]
    
    i <- pair[1]
    j <- pair[2]
    
    if(i %in% atts & j %in% atts & !is.na(val)) {
      
      if(negconvert) val <- convert_saaty(val)
      if(val == 0) val <- 1
      
      mat[i, j] <- val
      mat[j, i] <- 1 / val
    }
  }
  
  return(mat)
}

ahp_eigen <- function(mat) {
  ev <- eigen(mat)
  lambda_max <- Re(ev$values[1])
  w <- Re(ev$vectors[,1])
  w <- w / sum(w)
  
  return(list(weights = w, lambda = lambda_max))
}


ahp_cr <- function(lambda_max, n) {
  CI <- (lambda_max - n) / (n - 1)
  
  RI_table <- c(
    0.00, # n=1
    0.00, # n=2
    0.58, # n=3
    0.90, # n=4
    1.12, # n=5
    1.24, # n=6
    1.32, # n=7
    1.41, # n=8
    1.45, # n=9
    1.49  # n=10
  )
  
  RI <- RI_table[n]
  CR <- ifelse(RI == 0, 0, CI / RI)
  
  return(CR)
}


aggregate_geom_mean <- function(weight_list) {
  W <- do.call(rbind, weight_list)
  
  # geometric mean across respondents
  gm <- apply(W, 2, function(x) exp(mean(log(x))))
  
  # normalize
  gm <- gm / sum(gm)
  
  return(gm)
}