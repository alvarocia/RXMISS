list.of.packages <- c("crmReg", "MASS", "gss", "foreach", "doParallel", "randtoolbox", "nabor", "mvtnorm", "future", "progressr", "doFuture")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(crmReg)
library(MASS)
library(gss)
library(randtoolbox)
library(foreach)
library(doParallel)
library(nabor)

library(doFuture)
library(progressr)
library(future)



source("CODE/exchange_search_utils.R")


x_standardize = function(x) {
  x = (x - min(x)) / (max(x) - min(x))
  return(x)
}

mse_beta = function(x) {
  mean(apply(x - bt.true, 1, crossprod))
}

mspe_beta <- function(beta, x, y) {
  X <- cbind(1, as.matrix(x))
  
  mspe_each <- apply(beta, 1, function(b) {
    y_pred <- X %*% b
    mean((y - y_pred)^2)
  })
  
  mean(mspe_each)
}

# Setting working directory to source file location
# setwd(getSrcDirectory()[1]) #For R
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # For R studio



################################
#First dataset
################################
mat_raw = read.csv("DATA/soil.csv")
mat = mat_raw[, 3580:3600]
# select the variable: "CTI", "ELEV", "RELI", "TMAP", "TMF"
x = (apply(mat[, c(4, 5, 13, 14, 15)], 2, x_standardize) - 0.5) * 2
y = scale(mat[, 21], scale = T)


################################
#Second dataset
################################
# mat = read.csv("DATA/diamonds.csv")
# x = (apply(mat[, c(2, 6, 7)], 2, x_standardize)-0.5)*2
# y = scale(mat[, 8], scale=T)



N = length(y)
p = dim(x)[2]
nloop = 50
theta = 1

J_hat <- crossprod(cbind(1, as.matrix(x))) / nrow(as.matrix(x))



#calculate the sampling probabilities for leveraging methods
svdf <- svd(x)
U <- svdf$u
PP <- apply(U, 1, crossprod)
prob.lev <- PP / p
prob.shrink9 = 0.9 * prob.lev + 0.1 * 1 / N
#calculate the design space for Lowcon, using the parameter "theta"
lhd_data = apply(x, 2, x_standardize)
lowcon_space <-
  apply(lhd_data, 2, quantile, probs = c(theta / 100, 1 - theta / 100))



######################################################################
# SSANOVA
######################################################################
source("CODE/smoothing_spline_utility.R")
set.seed(1000)
sam.size = min(ceiling(200 * (length(y) ^ 0.25)), length(y))
esps.re <- esps(y ~ ., data = data.frame(y, x), sam.size, r = 3)
pro.gcv.lamb <- esps.re$lambda * ((N / sam.size) ^ (-3 / (1 * 3 + 1)))
lambda.pro <- log10(N * pro.gcv.lamb)
theta <- esps.re$theta
fit_ss <-
  ssa(y ~ .,
      data = as.data.frame(x),
      lambda = lambda.pro,
      theta = theta)

######################################################################
# OLS
######################################################################
fit_ols = lm(y ~ x)

######################################################################
# M-estimator
######################################################################
rfit <- rlm(y ~ x, method = "M")

######################################################################
# CRM-estimator
######################################################################
if (length(y) > 10000) {
  crmfit_coef = rep(NA, p + 1)
} else{
  crmfit <- crm(y ~ ., data = data.frame(y, x))
  crmfit_coef = crmfit$coefficient
}



# numCores <- parallel::detectCores()-1
# cl <- makeCluster(numCores)
# registerDoParallel(cl)

plan(multisession, workers = 6)
registerDoFuture()
handlers(global = TRUE)

for (cc in c(10)) {
  subsize = cc * p
  
  set.seed(1000)
  IBOSS.ind = rep(0, subsize)
  for (i in 1:p) {
    IBOSS.ind[(i - 1) * cc + 1:(floor(cc / 2))] <-
      which(rank(x[, i], ties.method = "random") %in% c(1:floor(cc / 2)))
    IBOSS.ind[i * cc - 0:(ceiling(cc / 2) - 1)] <-
      which(rank(x[, i], ties.method = "random") %in% c(N + 1 - 1:ceiling(cc /
                                                                            2)))
  }
  
  
  ptm = proc.time()
  
  with_progress({
    prog <- progressor(steps = nloop)
    
    res_parallel <- foreach(
      i = 1:nloop,
      .packages = c("nabor", "randtoolbox"),
      .combine = "cbind"
    ) %dopar% {
      on.exit(prog(sprintf("iter %d", i)), add = TRUE)
  
    tryCatch({
      setseed = 2517 + i * 821
      set.seed(setseed)
      
      ###################
      # UNIF
      
      unif.ind = sample.int(N,
                            size = subsize,
                            replace = TRUE,
                            prob = rep(1, N) / N)
      yy = as.matrix(y[unif.ind])
      xx = as.matrix(x[unif.ind, ])
      datt = data.frame(yy = yy, xx = xx)
      lm.unif = lm(yy ~ xx, data = datt)
      
      
      ###################
      # BLEV & LEVUNW
      
      lev.ind  = sample.int(N, size = subsize, replace = TRUE, prob.lev)
      yy = as.matrix(y[lev.ind])
      xx = as.matrix(x[lev.ind, ])
      wgt = 1 / prob.lev[lev.ind]
      datt = data.frame(yy = yy, xx = xx, wgt = wgt)
      lm.lev = lm(yy ~ xx, weights = wgt, data = datt)
      lm.levunwt = lm(yy ~ xx, data = datt)
      
      
      ###################
      # SLEV
      
      slev.ind = sample.int(N,
                            size = subsize,
                            replace = TRUE,
                            prob = prob.shrink9)
      yy = as.matrix(y[slev.ind])
      xx = as.matrix(x[slev.ind, ])
      wgt = 1 / prob.shrink9[slev.ind]
      datt = data.frame(yy = yy, xx = xx, wgt = wgt)
      lm.shrink9 = lm(yy ~ xx, weights = wgt, data = datt)
      
      
      ###################
      # LowCon
      design = sobol(
        n = subsize ,
        dim = p,
        init = T,
        scrambling = 1,
        seed = setseed
      )
      
      for (dd in 1:p) {
        design[, dd] = (design[, dd]) * (lowcon_space[2, dd] - lowcon_space[1, dd]) +
          lowcon_space[1, dd]
      }
      
      lhd.ind <- knn(lhd_data, design, k = 1)$nn.idx[, 1]
      if (length(unique(lhd.ind)) < subsize) {
        lhd.ind <-
          c(unique(lhd.ind), sample(N, subsize - length(unique(lhd.ind))))
      }
      
      yy = as.matrix(y[lhd.ind])
      xx = as.matrix(x[lhd.ind, ])
      datt = data.frame(yy = yy, xx = xx)
      lm.LowCon = lm(yy ~ xx, data = datt)
      
      
      ###################
      # MISSPECIFICATION B
      
      
      ok <- apply(lhd_data, 1, function(r) all(r >= lowcon_space[1, ] & r <= lowcon_space[2, ]))
      
      res <- exchange_search(
        X = as.matrix(x[ok, ]),
        n = subsize,
        J_hat = J_hat,
        B = 10,
        max_iter = 1,
        seed = setseed,
        delta_sub = 0.7,
        verbose = TRUE
      )
      
      MISS.ind <- which(ok)[res$idx]
      
      # res <- exchange_search(
      #   X = as.matrix(x),
      #   n = subsize,
      #   J_hat = J_hat,
      #   B = 10,
      #   max_iter = 5,
      #   seed = setseed,
      #   delta_sub = 0.7,
      #   verbose = TRUE
      # )
      # 
      # MISS.ind <- res$idx
      
      yy = as.matrix(y[MISS.ind])
      xx = as.matrix(x[MISS.ind, ])
      datt = data.frame(yy = yy, xx = xx)
      lm.MISS = lm(yy ~ xx, data = datt)
      
      ###################
      # IBOSS
      
      yy = as.matrix(y[IBOSS.ind])
      xx = as.matrix(x[IBOSS.ind, ])
      datt = data.frame(yy = yy, xx = xx)
      lm.IBOSS = lm(yy ~ xx, data = datt)
      
      
      rbind(
        fit_ols$coefficients,
        lm.unif$coefficients,
        lm.lev$coefficients,
        lm.levunwt$coefficients,
        lm.shrink9$coefficients,
        lm.IBOSS$coefficients,
        lm.LowCon$coefficients,
        lm.MISS$coefficients
      )
      
    }, error = function(e) {
      cat(i, "ERROR :", conditionMessage(e), "\n")
      NULL
    })
    }
  })
  
  mse.mat = matrix(NA, 4, 7)
  colnames(mse.mat) = c("unif", "lev", "shrink9", "unwtlev", "IBOSS", "LowCon", "MISS")
  rownames(mse.mat) = c("Full_OLS", "M-estimator", "CRM", "Spline")
  
  mspe.mat = matrix(NA, 1, 7)
  colnames(mspe.mat) = c("unif", "lev", "shrink9", "unwtlev", "IBOSS", "LowCon", "MISS")
  rownames(mspe.mat) = c("MSPE")
  
  
  for (i in 1:4) {
    if (i == 1)
      bt.temp <- fit_ols$coefficients
    if (i == 2)
      bt.temp <- rfit$coefficients
    if (i == 3)
      bt.temp <- crmfit_coef
    if (i == 4)
      bt.temp <- fit_ss$d
    
    
    bt.true = t(matrix(rep(bt.temp, nloop), ncol = nloop))
    bt.unif = t(matrix(res_parallel[2, ], ncol = nloop))
    bt.lev = t(matrix(res_parallel[3, ], ncol = nloop))
    bt.unwt.lev = t(matrix(res_parallel[4, ], ncol = nloop))
    bt.shrink9 = t(matrix(res_parallel[5, ], ncol = nloop))
    bt.IBOSS = t(matrix(res_parallel[6, ], ncol = nloop))
    bt.LowCon = t(matrix(res_parallel[7, ], ncol = nloop))
    bt.MISS = t(matrix(res_parallel[8, ], ncol = nloop))
    
    mse.unif = mse_beta(bt.unif)
    mse.lev = mse_beta(bt.lev)
    mse.unwt.lev = mse_beta(bt.unwt.lev)
    mse.shrink9 = mse_beta(bt.shrink9)
    mse.LowCon = mse_beta(bt.LowCon)
    mse.IBOSS = mse_beta(bt.IBOSS)
    mse.MISS = mse_beta(bt.MISS)
    
    mse.mat[i, ] = c(mse.unif,
                     mse.lev,
                     mse.shrink9,
                     mse.unwt.lev,
                     mse.IBOSS,
                     mse.LowCon,
                     mse.MISS)
    if (i==1){
    
    mspe.unif = mspe_beta(bt.unif, x, y)
    mspe.lev = mspe_beta(bt.lev, x ,y)
    mspe.unwt.lev = mspe_beta(bt.unwt.lev, x ,y)
    mspe.shrink9 = mspe_beta(bt.shrink9, x ,y)
    mspe.LowCon = mspe_beta(bt.LowCon, x ,y)
    mspe.IBOSS = mspe_beta(bt.IBOSS, x ,y)
    mspe.MISS = mspe_beta(bt.MISS, x ,y)
    
    mspe.mat[i, ] = c(mspe.unif,
                     mspe.lev,
                     mspe.shrink9,
                     mspe.unwt.lev,
                     mspe.IBOSS,
                     mspe.LowCon,
                     mspe.MISS)
    }
  }
  cat("#################  r/p=",cc,"#######################","\n")
  print(round(mse.mat, digits = 2))
  cat("#################  MSPE  #######################","\n")
  print(round(mspe.mat, digits = 2))
}
# stopCluster(cl)
