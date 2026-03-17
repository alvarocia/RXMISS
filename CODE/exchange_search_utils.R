############################################################
# exchange_search_utils.R
#
# Utilities to optimize a subsample using an exchange
# algorithm and a KDE-based Monte Carlo criterion.
#
# Main function:
#   exchange_search(X, n, J_hat, ...)
############################################################


############################################################
# Scott bandwidth matrix
############################################################
H_scott <- function(data, delta = 1) {
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  S <- cov(data)
  
  h <- n^(-1/(p+4))
  H <- (delta^2)*(h^2)*S
  H <- H + diag(1e-10, p)
  return(H)
}



############################################################
# KDE evaluation
############################################################

kde_eval_point <- function(x, data, H) {
  
  dens <- mean(
    apply(
      data,
      1,
      function(mu) mvtnorm::dmvnorm(x, mean = mu, sigma = H)
    )
  )
  
  return(dens)
}



kde_eval <- function(Xq, data, H) {
  
  return(apply(
    Xq,
    1,
    kde_eval_point,
    data = data,
    H = H
  )
  )
}



############################################################
# Sampling from Gaussian KDE
############################################################

rkde_gaussian <- function(m, data, H) {
  
  data <- as.matrix(data)
  
  n <- nrow(data)
  p <- ncol(data)
  
  centers <- data[
    sample.int(n, size = m, replace = TRUE),
    ,
    drop = FALSE
  ]
  
  noise <- mvtnorm::rmvnorm(
    m,
    mean = rep(0, p),
    sigma = H
  )
  
  return(centers + noise)
  
}



############################################################
# Trace helper
############################################################

trace_mat <- function(A) {
  
  sum(diag(A))
  
}



############################################################
# Criterion function g(X*)
############################################################

g_of_Xstar <- function(Xstar, J, ridge = 1e-8) {
  
  XtX <- crossprod(Xstar)
  
  XtX_r <- XtX + diag(ridge, ncol(XtX))
  
  inv <- solve(XtX_r)
  
  trace_mat(J %*% inv) * trace_mat(XtX)
  
}



############################################################
# Monte Carlo estimation of the integral
############################################################

estimate_integral_MC <- function(n,
                                 data_for_kde,
                                 H,
                                 J_hat,
                                 B = 2000,
                                 ridge = 1e-8,
                                 seed = 999) {
  
  set.seed(seed)
  
  vals <- numeric(B)
  
  for (b in 1:B) {
    
    Xstar <- rkde_gaussian(
      m = n,
      data = data_for_kde,
      H = H
    )
    
    Xstar <- cbind(1, Xstar)
    
    vals[b] <- g_of_Xstar(
      Xstar,
      J = J_hat,
      ridge = ridge
    )
    
  }
  
  list(
    estimate = mean(vals),
    sd = sd(vals),
    B = B,
    vals = vals
  )
  
}



############################################################
# Criterion for a given subsample
############################################################

criterion_subsample <- function(idx,
                                X,
                                n,
                                J_hat,
                                B = 1000,
                                ridge = 1e-8,
                                seed = 1,
                                delta_sub = 1) {
  
  X <- as.matrix(X)
  
  X_star <- X[idx, , drop = FALSE]
  
  H_sub <- H_scott(
    X_star,
    delta = delta_sub
  )
  
  res <- estimate_integral_MC(
    n = n,
    data_for_kde = X_star,
    H = H_sub,
    J_hat = J_hat,
    B = B,
    ridge = ridge,
    seed = seed
  )
  
  res$estimate
  
}



############################################################
# One exchange step
############################################################

exchange_step <- function(idx,
                          X,
                          n,
                          J_hat,
                          B,
                          ridge,
                          seed,
                          delta_sub) {
  
  X <- as.matrix(X)
  
  N <- nrow(X)
  
  improved_any <- FALSE
  
  current_val <- criterion_subsample(
    idx,
    X,
    n,
    J_hat,
    B,
    ridge,
    seed,
    delta_sub
  )
  
  for (k in seq_along(idx)) {
    
    in_idx <- idx[k]
    
    outside <- setdiff(seq_len(N), idx)
    
    best_local_val <- current_val
    
    best_out_idx <- in_idx
    
    for (out_idx in outside) {
      
      idx_new <- idx
      
      idx_new[k] <- out_idx
      
      val_new <- criterion_subsample(
        idx_new,
        X,
        n,
        J_hat,
        B,
        ridge,
        seed + 1,
        delta_sub
      )
      
      if (val_new < best_local_val) {
        
        best_local_val <- val_new
        
        best_out_idx <- out_idx
        
      }
      
    }
    
    if (best_out_idx != in_idx) {
      
      idx[k] <- best_out_idx
      
      current_val <- best_local_val
      
      improved_any <- TRUE
      
    }
    
  }
  
  list(
    idx = idx,
    value = current_val,
    improved = improved_any
  )
  
}



############################################################
# Main exchange search
############################################################

exchange_search <- function(X,
                            n,
                            J_hat,
                            B = 1000,
                            ridge = 1e-8,
                            max_iter = 30,
                            seed = 123,
                            verbose = TRUE,
                            delta_sub = 1,
                            init_idx = NULL) {
  
  X <- as.matrix(X)
  
  N <- nrow(X)
  
  set.seed(seed)
  
  if (is.null(init_idx)) {
    
    idx <- sample.int(N, n, replace = FALSE)
    
  } else {
    
    idx <- init_idx
    
  }
  
  best_val <- criterion_subsample(
    idx,
    X,
    n,
    J_hat,
    B,
    ridge,
    seed,
    delta_sub
  )
  
  if (verbose)
    cat("Initial criterion:", best_val, "\n")
  
  for (it in seq_len(max_iter)) {
    
    step <- exchange_step(
      idx,
      X,
      n,
      J_hat,
      B,
      ridge,
      seed + it,
      delta_sub
    )
    
    if (step$improved) {
      
      idx <- step$idx
      
      best_val <- step$value
      
      if (verbose)
        cat("Iter", it, "- improved:", best_val, "\n")
      
    } else {
      
      if (verbose)
        cat("Iter", it, "- no improvement\n")
      
    }
    
  }
  
  list(
    idx = idx,
    value = best_val
  )
  
}