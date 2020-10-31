

# Outer Products ----------------------------------------------------------

get_outer_product_vectors <- function(m) {
  if (qr(m)$rank != 1)
    stop("Matrix not seperable")
  
  svd_out <- svd(m)
  
  outerProd_h <- matrix(svd_out$v[, 1] * sqrt(svd_out$d[1]),
                        nrow = 1,
                        ncol = ncol(m))
  
  outerProd_v <- matrix(svd_out$u[, 1] * sqrt(svd_out$d[1]),
                        nrow = nrow(m),
                        ncol = 1)
  
  list(x = outerProd_h, y = outerProd_v)
}

# Gaussian Weights --------------------------------------------------------

gaussian_weights <- function(radius, scaling, type = "vector") {
  sigma <- ceiling(radius * scaling)
  n <- sigma * 3
  n_seq <- seq(-n, n, 1)
  n_weights <- dnorm(n_seq, sd = sigma)
  if (type == "matrix") {
    x <- matrix(n_weights, ncol = length(n_weights))
    y <- matrix(n_weights, ncol = 1)
    m <- y %*% x
    return(m)
  } else if (type == "vector") {
    return(n_weights)
  } else if (type == "outer vectors") {
    x <- matrix(n_weights, ncol = length(n_weights))
    y <- matrix(n_weights, ncol = 1)
    return(list(x = x, y = y))
  } else {
    print("type must be one of vector, matrix, outer vectors")
  }
}


# Seperated Focal ---------------------------------------------------------

focal_seperable <- function(r, outer_vectors) {
  return(focal(focal(r, outer_vectors$y),  outer_vectors$x))
  
}
