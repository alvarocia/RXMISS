############################################################
# Kernel density (multivariate Gaussian KDE) estimation for a subsample
# + Monte Carlo approx of:
#  ∫ tr(J(X*^T X*)^{-1}) tr(X*^T X*)  ∏_{i=1}^n h_n(x_i*) dx_1*...dx_n*
# + exchange algorithm
############################################################

# ---------- Packages ----------
pkgs <- c("mvtnorm")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
library(mvtnorm)

source("CODE/exchange_search_utils.R")

my_seed <- 12345678
set.seed(my_seed)
X <- rnorm(500)
# X <- rt(100, df=2)
# X <- runif(100, min=-3, max=3)
# Sigma <- matrix(c(1, 0.7,
#                   0.7, 1), nrow = 2)
# 
# X <- MASS::mvrnorm(100, mu = c(0, 0), Sigma = Sigma)
X <- scale(X)  # optional standardization
N <- nrow(X)
p <- ncol(X)
n <- 100  # subsample size
delta_full <- 1


# Estimate matrix J with full data
Xmatrix <- cbind(1, X)   # n x (p+1)
J_hat <- crossprod(Xmatrix) / nrow(Xmatrix)  # (p+1) x (p+1)


res_opt <- exchange_search(
  X = X,
  n = n,
  J_hat = J_hat,
  B = 50,
  ridge = 1e-8,
  max_iter = 2,
  seed = my_seed+1,
  verbose = TRUE,
  delta_sub = 0.7
)

idx_opt <- res_opt$idx
X_star_data <- X[idx_opt, , drop = FALSE]



H_full <- H_scott(X, delta = delta_full)
H_sub  <- H_scott(X_star_data, delta = 0.7)

# B <- 300
# res_I <- estimate_integral_MC(
#   n = n,
#   data_for_kde = X_star_data,
#   H = H_sub,
#   B = B,
#   ridge = 1e-8,
#   seed = my_seed + 1
# )



if (p==2){
# 2D:
x_seq <- seq(-3, 3, by = 0.1)
y_seq <- seq(-3, 3, by = 0.1)

grid <- expand.grid(x_seq, y_seq)
colnames(grid) <- c("x1", "x2")

z_full <- kde_eval(as.matrix(grid), X, H_full)
z_sub  <- kde_eval(as.matrix(grid), X_star_data, H_sub)

Z_full_mat <- matrix(z_full, nrow = length(x_seq), ncol = length(y_seq))
Z_sub_mat  <- matrix(z_sub,  nrow = length(x_seq), ncol = length(y_seq))

Z_full_scaled <- (N / n) * Z_full_mat

  
# Heatmap de Z_sub_mat
image(
  x_seq, y_seq, Z_sub_mat,
  col = heat.colors(50),
  xlab = "X1",
  ylab = "X2",
  main = "Full (contours) vs Sub (heatmap)"
)

# Curvas de nivel de Z_full_scaled encima
contour(
  x_seq, y_seq, Z_full_scaled,
  add = TRUE,
  drawlabels = FALSE,
  col = "black",
  lwd = 1
)
}


# ONE DIMENSION:

if (p==1){
  

# plot(seq(-3,3,by=0.01),kde_eval(as.matrix(seq(-3,3,by=0.01)), X, H_full) )
# plot(seq(-3,3,by=0.01),kde_eval(as.matrix(seq(-3,3,by=0.01)), X_star_data, H_sub) )

xx <- seq(-3, 3, by = 0.05)
y1 <- N/n * kde_eval(as.matrix(xx), X, H_full)
y2 <- kde_eval(as.matrix(xx), X_star_data, H_sub)

plot(
  xx, y1,
  type = "p",          # points; use "l" for lines
  col = "blue",
  pch = 16,
  xlab = "x",
  ylab = "Density"
)

points(
  xx, y2,
  col = "red",
  pch = 17
)

legend(
  "topright",
  legend = c("Full KDE", "Sub KDE"),
  col = c("blue", "red"),
  pch = c(16, 17),
  bty = "n"
)

}