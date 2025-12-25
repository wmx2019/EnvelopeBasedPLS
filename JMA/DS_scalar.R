################################################################################
# Functional Envelope Partial Least Squares (FEPLS) Simulation Study for Scalar 
# Response
#
# Purpose: Compare FEPLS performance against OLS, PCR, and PLS methods
#          for scalar response prediction with functional predictors
# 
# Paper Section: Appendix D
################################################################################

rm(list = ls())
# Required packages
library(doFuture)
library(doRNG)
library(fda.usc)
source("FEPLS.R")

Sys.setenv(LANG = "en")

################################################################################
# CONFIGURATION
################################################################################

# Number of replications
n_sim <- 100

# Set the sample sizes for sample set and test (validation) set
n_validation <- 1000

# Set seeds for reproducibility
set.seed(8)
registerDoRNG(144)

# Setup parallel processing and the number of cores
registerDoFuture()
plan(multisession,workers=min(20, parallel::detectCores() - 1,n_sim))


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

DG_s <- function(B_list,n,t){
   #' Generate Scalar Response Data
   #'
   #' Simulates functional predictors and scalar responses for the FEPLS study
   #'
   #' @param B_list List of coefficient matrices (one per predictor)
   #' @param tau_list List of eigenvalue vectors (one per predictor)
   #' @param n Sample size
   #' @param t Time points vector
   #' @return List with X_list (functional predictors) and Y (response matrix)
   #' @details
   #'   See Appendix D.
   
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
   
   r <- dim(B_list[[1]])[1]
   
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
   
   
   # Generate responses
   Y <- matrix(0,n,r)
   e <- matrix(rnorm(n*r,sd=1),n,r)
   for (l in 1:p) {
      xi <- xi_list[[l]]
      B <- B_list[[l]]
      for (i in 1:n) {
         Y[i,] <- Y[i,] + sum((tau_list[[l]]^{1/2}*xi[i,])%*%t(B))
      }
   }
   
   # Add Gaussian error
   Y <- Y+0.8*e
   
   return(list(X_list=X_list,Y=Y))
}

################################################################################
# MODEL SPECIFICATION
################################################################################

# Coordinate of operator B given the Fourier (TRUE) basis
B1 <- matrix(0, nrow = 4, ncol = 13)
B1[1, 2:3] <- c(-1.2, -0.4)
B1[2, 2:3] <- c(0.4, 0.6)
B1[3, 2:3] <- c(0.3, -0.4)
B1[4, 2:3] <- c(-1.2, 0.3)
B_list <- list(B1 = B1)

# Eigenvalues of covariance matrix given the Fourier (TRUE) basis
tau1 <- c(8,2,0.08,12,5.5,runif(8,2,3)) 
tau_list <- list(tau1)



# Basis configuration for FEPLS
knots <- seq(0,1,1/5)
order <- 4
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

ux <- 2
t <- seq(0,1.0,1/12)
tx.list <- list(t1=t)

################################################################################
# SIMULATION STUDY
################################################################################
mse <- NULL
for (i in 1:5) {
   n_samples <- 2^(i-1)*25
   acc <- foreach(i=1:n_sim) %dopar% {
      # Generate training data
      train_data <- DG_s(B_list,n=n_samples,t)
      X_list <- train_data$X_list
      X.list <- X_list
      Y <- train_data$Y
      
      # Generate test data
      test_data <- DG_s(B_list,n=n_validation,t)
      X_list_test <- test_data$X_list
      Y_test <- test_data$Y
      
      # Fit FEPLS and competing methods
      res_dirs <- sfpedir(X_list,Y,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
      
      # Estimated coordinates of test data given selected basis functions
      X_test_dir <- bvfitl(X_list_test,res_dirs$basis.value.X)
      
      c(MSEFs(res_dirs$beta,res_dirs$alpha,X_test_dir,Y_test),
        MSEFs(res_dirs$betafull,res_dirs$alphafull,X_test_dir,Y_test),
        MSEPs(res_dirs$mpcr,res_dirs$ncomp.pcr,X_test_dir,Y_test),
        MSEPs(res_dirs$mpls,res_dirs$ncomp.pls,X_test_dir,Y_test))
   }
   acc_M_matrix <- matrix(unlist(acc),byrow = TRUE,ncol = 4)
   
   colnames(acc_M_matrix) <- c("FEPLS","OLS","PCR","PLS")
   mse <- cbind(mse,signif(colMeans(acc_M_matrix),digits = 3))
}

plan(sequential)

if(!file.exists("results_DS_scalar.rds")){
   saveRDS(list(mse=mse),file="results_DS_scalar.rds")
} 

colnames(mse) <- 25*2^(0:4)
print(mse)








