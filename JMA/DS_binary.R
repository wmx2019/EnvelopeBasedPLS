################################################################################
# Generalized Functional Envelope Partial Least Squares (GFEPLS) Simulation Study 
# for Binary Response
#
# Purpose: Compare GFEPLS performance against GLM, GLMNET, and GPLS methods
#          for binary response prediction with functional predictors
# 
# Paper Section: Section 6.2
################################################################################


rm(list = ls())
# Required packages
source("FEPLS.R")
library(doFuture)
library(doRNG)
library(glmnet)
library(fda.usc)
library(xtable)
library(mixOmics)
Sys.setenv(LANG="en")

################################################################################
# CONFIGURATION
################################################################################

# Number of replications
n_sim <- 100

# Set the sample sizes for sample set and test (validation) set
n.test <- 10000

# Set sample sizes
n.v <- 40*2^(0:4)

# Set seeds for reproducibility
set.seed(6)

# Setup parallel processing and the number of cores
plan(multisession, workers = min(20, parallel::detectCores() - 1,n_sim))

################################################################################
# HELPER FUNCTIONS
################################################################################

psi <- function(x){
   #' Fourier Basis Function
   #'
   #' Generates (true) Fourier basis functions (sine and cosine terms)
   #'
   #' @param x Evaluation point(s) in [0,1]
   #' @return Vector of basis function values
   r <- 1
   for (i in 1:6) {
      r <- c(r,sqrt(2)*sin(2*i*pi*x),sqrt(2)*cos(2*i*pi*x))
   }
   return(r)
}   


DG_binary <- function(B_list,n,t){
   #' Generate binary Response Data
   #'
   #' Simulates functional predictors and binary responses for the GFEPLS study
   #'
   #' @param B_list List of coefficient matrices (one per predictor)
   #' @param tau_list List of eigenvalue vectors (one per predictor)
   #' @param n Sample size
   #' @param t Time points vector
   #' @return List with X_list (functional predictors) and Y (response vector)
   #' @details
   #'   See Section 6.2.
   
   
   # Validate inputs
   stopifnot(
      "n must be positive integer" = n > 0 && n == round(n),
      "t must be numeric vector" = is.numeric(t) && length(t) > 0,
      "B_list cannot be empty" = length(B_list) > 0,
      "tau_list length must match B_list" = length(tau_list) == length(B_list)
   )
   
   nt <- length(t)
   X_list <- list()
   xi_list <- list()
   p <- length(B_list)
   
   
   # Generate functional predictors
   for (l in 1:p) {
      x <- matrix(0,n,nt)
      xi <- matrix(rnorm(13*n),n)
      xi_list <- append(xi_list,list(xi))
      for (i in 1:n) {
         for (j in 1:nt) {
            x[i,j] <- sum(tau_list[[l]]^{1/2}*xi[i,]*psi(t[j]))
         }
      }
      X_list <- append(X_list,list(x))
   }
   
   X_list <- lapply(X_list, scale)
   
   # Generate responses
   eta <- matrix(0,n,1)
   for (l in 1:p) {
      xi <- xi_list[[l]]
      B <- B_list[[l]]
      for (i in 1:n) {
         eta[i] <- eta[i] + sum((tau_list[[l]]^{1/2}*xi[i,])%*%t(B))
      }
   }
   lp <- exp(eta)/(1+exp(eta))
   Y <- rbinom(n,1,lp)
   
   return(list(X_list=X_list,Y=Y,lp=lp))
}



################################################################################
# MODEL SPECIFICATION
################################################################################

# Coordinate of operator B given the Fourier (TRUE) basis
B1 <- matrix(0,nrow = 1,ncol=13)
B1[1,2] <- -0.8
B1[1,3] <- 4.4
B_list <- list(B1=B1)


# Eigenvalues of covariance matrix given the Fourier (TRUE) basis
tau1 <- c(8,5,0.9,7.3,2,runif(8,0.1,0.2))
tau_list <- list(tau1)


# Basis configuration for GFEPLS
knots <- seq(0,1,1/4)
order <- 4
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

ux <- 2
t <- seq(0,1.0,1/14)
tx.list <- list(t1=t)



################################################################################
# SIMULATION STUDY
################################################################################


mmr.l <- replicate(length(n.v),NA, simplify = FALSE)
results.l <- list()
names(mmr.l) <- n.v

# Suppress the warning caused by GLM due to imbalanced data
options(warn=-1)

for (i in 1:length(n.v))  {
   n_samples <- n.v[i]
   set.seed(6)
   mmr <- foreach(i=1:n_sim,.combine = "rbind", .options.future = list(seed = TRUE)) %dofuture% {
      
      # Generate training data
      train_data <- DG_binary(B_list,n=n_samples,t)
      X_train_list <- train_data$X_list
      Y_train <- train_data$Y
      
      # X_train_cord <- get_coord_dir_sy_sp(X_train_list,tx.list,x.knots.list,x.order.list)$X.cord
      # dx <- ncol(X_train_cord)
      
      # Generate test data
      test_data <- DG_binary(B_list,n=n.test,t)
      X_test_list <- test_data$X_list
      X_test <- X_test_list[[1]]
      Y_test <- test_data$Y
      
      # Fit GFEPLS and competing methods
      res_dirc <- cfelmdir(X_train_list,Y_train,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T,do.gpls =T)
      
      # Estimated coordinates of test data given selected basis functions
      X_test_dir <- bvfitl(X_test_list,res_dirc$basis.value.X)
      
      acc_env <- pred_acc(res_dirc$beta,res_dirc$alpha,X_test_dir,Y_test)
      acc_env_full <- pred_acc(res_dirc$betafull,res_dirc$alphafull,X_test_dir,Y_test)
      acc_glmnet <- pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_test_dir,Y_test)
      acc_gpls <- pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_test_dir,Y_test)
      
      c(1-c(acc_env,
            acc_env_full, 
            acc_glmnet,
            acc_gpls))
      
   }
   colnames(mmr) <- c("env","fglm","fglmnet","fgpls")
   
   results.l[[i]] <- mmr
   mmr.l[[i]] <- signif(colMeans(mmr),3)
}

plan(sequential)

if(!file.exists("results_DS_binary.rds")){
   saveRDS(list(results=results.l,mmr=mmr.l),file="results_DS_binary.rds")
} 


mmr.m <- do.call(cbind,mmr.l)
xtable(mmr.m,digits = 3)
print(mmr.m)
