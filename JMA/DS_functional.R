################################################################################
# Functional Envelope Partial Least Squares (FEPLS) Simulation Study for 
# Functional Response
#
# Purpose: Compare FEPLS performance against FFFR, PCR, and PLS methods
#          for functional response prediction with functional predictors
# 
# Paper Section: Section 6.1
################################################################################

rm(list=ls())
# Required packages
source("FEPLS.R")
library(doFuture)
library(doRNG)
library(fda.usc)
library(xtable)
Sys.setenv(LANG="en")

################################################################################
# CONFIGURATION
################################################################################

# Number of replications
n_sim <- 100

# Set the sample sizes for sample set and test (validation) set
n_validation <- 10000

# Set sample sizes
n.v <- 50*2^(0:4)

# Set seeds for reproducibility
set.seed(6)

# Setup parallel processing and the number of cores
plan(multisession, workers = min(20, parallel::detectCores() - 1,n_sim))

################################################################################
# HELPER FUNCTIONS
################################################################################

psi <- function(x){
   #' Fourier Basis Function for predictor
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

chi <- psi # Fourier Basis Function for response

DG_f <- function(B_list,n,t){
   #' Generate Functional Response Data
   #'
   #' Simulates functional predictors and functional responses for the FEPLS study
   #'
   #' @param B_list List of coefficient matrices (one per predictor)
   #' @param tau_list List of eigenvalue vectors (one per predictor)
   #' @param n Sample size
   #' @param t Time points vector
   #' @return List with X_list (functional predictors) and Y (response matrix)
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
      B <- B_list[[l]]
      xi <- matrix(rnorm(dim(B)[2]*n),n)
      xi_list <- append(xi_list,list(xi))
      for (i in 1:n) {
         for (j in 1:nt) {
            x[i,j] <- sum(tau_list[[l]]^{1/2}*xi[i,]*psi(t[j]))
         }
      }
      X_list <- append(X_list,list(x))
   }
   
   # Generate error terms
   e <- matrix(0,n,nt)
   nu <- matrix(rnorm(13*n),nrow = n)
   for (i in 1:n) {
      for (j in 1:nt) {
         e[i,j] <- 0.2*sum(rho^{1/2}*nu[i,]*chi(t[j]))
      }
   }
   
   # Generate responses
   Y <- matrix(0,n,nt)
   for (l in 1:p) {
      xi <- xi_list[[l]]
      B <- B_list[[l]]
      for (i in 1:n) {
         for (j in 1:nt){
            Y[i,j] <- Y[i,j] + sum((tau_list[[l]]^{1/2}*xi[i,])%*%t(B)%*%chi(t[j]))
         }
      }
   }
   Y <- Y+e
   
   
   return(list(X_list=X_list,Y=Y))
}


################################################################################
# MODEL SPECIFICATION
################################################################################

# Coordinate of operator B given the Fourier (TRUE) basis
B1 <- matrix(0,nrow = 13,ncol=13)
B1[2,2] <- -1.2
B1[4,3] <- 2.4
B2 <- matrix(0,nrow = 13,ncol=13)
B2[2,3] <- -0.04
B2[1,3] <- -0.05
B3 <- matrix(0,nrow = 13,ncol=13)
B3[3,3] <- 0.03
B3[4,3] <- -0.01

B_list <- list(B1,B2,B3)


# Basis configuration for FEPLS
x.knots.list <- list(knots1=seq(0,1,1/4),knots2=seq(0,1,1/4),knots3=seq(0,1,1/4))
x.order.list <- list(order1=4,order2=4,order3=4)
y.knots <- seq(0,1,1/5)
y.order <- 4

ux <- 4
t <- seq(0,1.0,1/15)
tx.list <- list(t1=t,t2=t,t3=t)
ty <- t

# Eigenvalues of covariance matrix given the Fourier (TRUE) basis
tau1 <- c(8,4,1.6,7.3,5.5,runif(8,0.2,0.3))
tau2 <- c(8,2,1,7.3,5.5,runif(8,0.2,0.3))
tau3 <- c(8,2,0.5,3,5.5,runif(8,0.2,0.3))
tau_list <- list(tau1,tau2,tau3) # Predictor
rho <- rep(1,13) # Response


################################################################################
# SIMULATION STUDY
################################################################################

results.l <- replicate(length(n.v), NA, simplify = FALSE)
mse.l <- replicate(length(n.v),NA, simplify = FALSE)

for (i in 1:length(n.v)) {
   n_samples <- n.v[i]
   set.seed(6)
   result_mpse <- foreach(i=1:n_sim,.combine = "rbind", .options.future = list(seed = TRUE)) %dofuture% {
      
      # Generate training data
      train_data <- DG_f(B_list,n=n_samples,t)
      X_train_list <- train_data$X_list
      Y_train <- train_data$Y
      
      # X_train_cord <- get_coord_dir_sy_sp(X_train_list,
      #                               tx.list,
      #                               x.knots.list,
      #                               x.order.list)$X.cord
      # dx <- ncol(X_train_cord)
      
      
      # Generate test data
      test_data <- DG_f(B_list,n=n_validation,t)
      X_list_test <- test_data$X_list
      Y_test <- test_data$Y
      
      
      
      # Fit FEPLS and competing methods
      res_dir <- mfpedir(X_train_list,Y_train,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                         ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)
      
      
      # Estimated coordinates of test data given selected basis functions
      X_test_dir <- bvfitl(X_list_test,res_dir$basis.value.X)
      
      mse_env <- MSEF(res_dir$beta,res_dir$alpha,res_dir$basis.value.Y,X_test_dir,Y_test)
      mse_full <- MSEF(res_dir$betafull,res_dir$alphafull,res_dir$basis.value.Y,X_test_dir,Y_test)
      mse_pcr <- MSEP(res_dir$mpcr,NULL,res_dir$basis.value.Y,X_test_dir,Y_test)
      mse_pls <- MSEP(res_dir$mpls,NULL,res_dir$basis.value.Y,X_test_dir,Y_test)
      
      c(mse_env,mse_full,mse_pcr,mse_pls)
      
   }
   colnames(result_mpse) <- c("FEPLS","FFFR","PCR","PLS")
   mpse <- colMeans(result_mpse)
   
   mse.l[[i]] <- mpse
   results.l[[i]] <- result_mpse
}

plan(sequential)

if(!file.exists("results_DS_functional.rds")){
   saveRDS(list(results=results.l,mse=mse.l),file="results_DS_functional.rds")
} 


mse.m <- do.call(rbind,mse.l)
colnames(mse.m) <- c("FEPLS","FFFR","PCR","PLS")
rownames(mse.m) <- n.v
xtable(t(mse.m))
print(xtable(t(mse.m)))













