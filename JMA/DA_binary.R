################################################################################
# Generalized Functional Envelope Partial Least Squares (GFEPLS) Empirical Study
# for Binary Response
#
# Purpose: Compare GFEPLS performance against GLM, GLMNET, and GPLS methods
#          for varying numbers of components using cross-validation.
# 
# Paper Section: Section 7.2
################################################################################

remove(list = ls())
Sys.setenv(LANG="en")
# Required packages
library(refund)
library(doFuture)
library(doRNG)
library(glmnet)
source("FEPLS.R")
options(warn = -1)


################################################################################
# CONFIGURATION
################################################################################

# Set seeds for reproducibility
set.seed(4534)

# Number of replications
n_sim <- 100

# Cross-validation setup
n_fold <- 5


################################################################################
# HELPER FUNCTIONS
################################################################################
read_data_binary <- function(){
   #' Read and Preprocess Binary Classification Data
   #'
   #' @return List containing Y (binary response), X_list (functional predictors),
   #'   X.cord (basis coordinates), dx (number of basis functions), tx.list,
   #'   x.knots.list, x.order.list (basis specifications), and t (time grid)
   
   # Load gasoline NIR spectroscopy data
   data("gasoline")
   X <- gasoline$NIR
   Y <- as.matrix(gasoline$octane)
   
   # Create binary response: octane > 88
   Y <- Y>88
   
   # Set up parameters for spline basis
   t <- seq(0,1,1/400)
   X_list <- list(X=X)
   
   # Set the parameters for the spline basis
   knots <- seq(0,1,1/8)
   order <- 4
   tx.list <- list(t1=t)
   x.knots.list <- list(knots1=knots)
   x.order.list <- list(order1=order)
   
   # Estimate coordinates w.r.t the given basis functions
   res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
   X.cord <- res_cord$X.cord
   dx <- ncol(X.cord)
   
   
   return(list(Y=Y,X_list=X_list,X.cord=X.cord,dx=dx,tx.list=tx.list,
               x.knots.list=x.knots.list,x.order.list=x.order.list,t=t))
}


################################################################################
# EMPIRICAL STUDY
################################################################################
# Load and prepare data
data <- read_data_binary()
Y <- data$Y
X_list <- data$X_list

n <- length(Y)
n_val <- as.integer(n/n_fold)

dx <- data$dx
tx.list <- data$tx.list
x.knots.list <- data$x.knots.list
x.order.list <- data$x.order.list
t <- data$t
knots <- x.knots.list$knots1
order <- x.order.list$order1

# Parallel cross-validation: compare methods across all component dimensions
plan(multisession,workers=min(parallel::detectCores()-1,20,n_sim))
mmr.l <- foreach(i=1:n_sim, .options.future = list(seed = TRUE)) %dofuture% {
   options(warn = -1)
   
   # Create random folds
   folds <- split(sample(1:n),sample(rep(1:n_fold, length.out = n)))
   tmp <- array(0, dim=c(dx, 6))
   
   # Loop over all possible numbers of components
   for (ux in 1:dx) {
      for (k in 1:n_fold) {
         
         # Split data into training and validation sets
         idx_val <- folds[[k]]
         idx_train <- setdiff(1:n,idx_val)
         Y_train <- Y[idx_train]
         Y_val <- as.matrix(Y[idx_val])
         X_train_list <- list()
         X_val_list <- list()
         for (X in X_list) {
            X_train <- X[idx_train,]
            X_val <- X[idx_val,]
            X_train_list <- append(X_train_list,list(X_train))
            X_val_list <- append(X_val_list,list(X_val))
         }
         
         # Fit GFEPLS and competing methods
         res_dirc <- cfelmdir(X_train_list,Y_train,ux,tx.list=tx.list,
                              x.knots.list=x.knots.list,x.order.list=x.order.list,
                              u.selected=F,do.gpls = T)
         
         
         Y_val <- factor(Y_val, levels = c(FALSE, TRUE))
         
         # Estimated coordinates of test data given selected basis functions
         X_val_dir <- bvfitl(X_val_list,res_dirc$basis.value.X)
         
         acc_env <- pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val)
         acc_full <- pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val)
         acc_mglmnet <- pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val)
         acc_mgpls <- pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_val_dir,Y_val)
         
         # Compute MMR for all three methods
         tmp[ux,] <-  tmp[ux,] + c(1-c(acc_env,acc_full,acc_mglmnet,acc_mgpls),res_dirc$ux,res_dirc$ncomp.gpls)
      }
      tmp[ux,] <- tmp[ux,]/n_fold
   }
   tmp
}
plan(sequential)


mmr.m <- Reduce("+", mmr.l)/n_sim


colnames(mmr.m) <- c("GFEPLS","GLM","GLMNET","GPLS","u_env",
                     "u_gpls")

minimal_mmr <- signif(c(mmr.m[which.min(mmr.m[,1]),1],
                        mmr.m[which.min(mmr.m[,2]),2],
                        mmr.m[which.min(mmr.m[,3]),3],
                        mmr.m[which.min(mmr.m[,4]),4],
                        mmr.m[which.min(mmr.m[,1]),5],
                        mmr.m[which.min(mmr.m[,4]),6]),3)
print(minimal_mmr)

if (!file.exists("results_DA_binary.Rds")) {
   saveRDS(list(mmr.m=mmr.m),file = "results_DA_binary.Rds")
}

################################################################################
# VISUALIZATION
################################################################################
library(ggplot2)
data <- data.frame(mmr.m[,c(1:4)])
data$x <- 1:dx

p1 <- ggplot(data) +
   geom_point(aes(x = x, y = GFEPLS, color = "GFEPLS"), size = 1) +
   geom_line(aes(x = x, y = GFEPLS, color = "GFEPLS"), linewidth = 1) +
   geom_point(aes(x = x, y = GLM, color = "GLM"), size = 1) +
   geom_line(aes(x = x, y = GLM, color = "GLM"), linewidth = 1) +
   geom_point(aes(x = x, y = GLMNET, color = "GLMNET"), size = 1) +
   geom_line(aes(x = x, y = GLMNET, color = "GLMNET"), linewidth = 1) +
   geom_point(aes(x = x, y = GPLS, color = "GPLS"), size = 1) +
   geom_line(aes(x = x, y = GPLS, color = "GPLS"), linewidth = 1) +
   labs( x = "Number of components", y = "Mean misclassification rate", color = "Method") +
   scale_color_manual(values = c("GFEPLS" = 1, "GLM" = 2, "GLMNET" = 3, "GPLS" = 4)) +
   theme_minimal()

p1
if (!file.exists("DA_binary.pdf")) {
   ggsave("DA_binary.pdf",p1,width = 6,height = 4)
}


