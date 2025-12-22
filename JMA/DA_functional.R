################################################################################
# Functional Envelope Partial Least Squares (FEPLS) Empirical Study for Functional 
# Response
#
# Purpose: Compare FEPLS performance against PCR and PLS methods
#          for varying numbers of components using cross-validation.
# 
# Paper Section: Section 7.1
################################################################################

remove(list = ls())
# Required packages
library(doFuture)
library(doRNG)
Sys.setenv(LANG="en")
source("FEPLS.R")

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
read_data_func <- function(){
   #' Load and Prepare Ocean Salinity Data
   #'
   #' Loads ocean salinity data and converts functional predictors and response
   #' to spline basis coefficients.
   #'
   #' @return List containing:
   #'   \item{Y}{Response matrix (salinity)}
   #'   \item{X_list}{List of predictor matrices}
   #'   \item{X.cord}{Predictor coordinates w.r.t. spline basis}
   #'   \item{dx}{Number of basis dimensions for predictors}
   #'   \item{tx.list, x.knots.list, x.order.list}{Predictor basis parameters}
   #'   \item{ty, y.knots, y.order}{Response basis parameters}
   #'

   # ocean.RData in the Functional_Response_Data folder
   load("Functional_Response_Data/ocean.RData")
   Y <- ocean$Salinity
   X_list <- ocean[names(ocean)!="Salinity"]
   

   # Set the parameters for the spline basis
   x.knots.list <- list(knots1=seq(0,1,1/5),knots2=seq(0,1,1/5),knots3=seq(0,1,1/5),knots4=seq(0,1,1/5))
   x.order.list <- list(order1=4,order2=4,order3=4,order4=4)
   t <- seq(0,1.0,1/(dim(X_list[[1]])[2]-1))
   tx.list <- list(t1=t,t2=t,t3=t,t4=t)
   
   y.knots <- seq(0,1,1/4)
   y.order <- 4
   ty <- seq(0,1.0,1/(dim(Y)[2]-1))
   
   # Estimate coordinates w.r.t the given basis functions
   res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
   X.cord <- res_cord$X.cord
   dx <- ncol(X.cord)
   return(list(Y=Y,X_list=X_list,X.cord=X.cord,dx=dx,tx.list=tx.list,
               x.knots.list=x.knots.list,x.order.list=x.order.list,y.knots=y.knots,y.order=y.order,ty=ty))
}

pred_error_cv_func_all <- function(X_list, Y, ux=NULL, n_fold=5,n_sim=50){
   ' Cross-Validation Prediction Error for Functional Methods
#'
   #' Estimates prediction error via repeated cross-validation for FEPLS, PCR, 
   #' and PLS across all possible numbers of components (1 to dx).
   #'
   #' @param X_list List of predictor matrices
   #' @param Y Response matrix
   #' @param ux Not used (included for compatibility)
   #' @param n_fold Number of cross-validation folds (default: 5)
   #' @param n_sim Number of replications (default: 50)
   #'
   #' @return List of length n_sim, each element is a dx x 6 matrix with columns:
   #'   FEPLS, PCR, PLS prediction errors, and selected dimensions for each method
   
   n <- dim(Y)[1]
   
   # registerDoFuture()
   
   # Set up parallel processing
   plan(multisession,workers=min(parallel::detectCores()-1,20,n_sim))
   
   
   # Parallel cross-validation: compare methods across all component dimensions
   pse <- foreach(i=1:n_sim, .options.future = list(seed = TRUE)) %dofuture%  {
      
      # Create random folds
      folds <- split(sample(1:n),sample(rep(1:n_fold, length.out = n)))
      tmp_res <- array(0, dim=c(dx, 6))
      
      # Loop over all possible numbers of components
      for (ux in 1:dx) {
         for (j in 1:n_fold) {
            
            # Split data into training and validation sets
            idx_val <- folds[[j]]
            idx_train <- setdiff(1:n,idx_val)
            Y_train <- Y[idx_train,]
            Y_val <- Y[idx_val,]
            X_list_train <- list()
            X_list_val <- list()
            for (X in X_list) {
               X_train <- X[idx_train,]
               X_val <- X[idx_val,]
               X_list_train <- append(X_list_train,list(X_train))
               X_list_val <- append(X_list_val,list(X_val))
            }
            
            # Fit FEPLS and competing methods
            res_dir <- mfpedir(X_list_train,Y_train,ux, tx.list=tx.list, ty=ty, x.knots.list=x.knots.list
                               ,y.knots=y.knots, x.order.list=x.order.list, y.order = y.order,u.selected=F)
            
            # Estimated coordinates of test data given selected basis functions
            X_val_dir <- bvfitl(X_list_val,res_dir$basis.value.X)
            
            # Compute prediction errors for all three methods
            tmp_res[ux,] <- tmp_res[ux,] + c(MSEF(res_dir$beta,res_dir$alpha,res_dir$basis.value.Y,X_val_dir,Y_val),
                                             MSEP(res_dir$mpcr,NULL,res_dir$basis.value.Y,X_val_dir,Y_val),
                                             MSEP(res_dir$mpls,NULL,res_dir$basis.value.Y,X_val_dir,Y_val),
                                             res_dir$ux,res_dir$mpcr$ncomp,res_dir$mpls$ncomp)
         }
         tmp_res[ux,] <- tmp_res[ux,]/n_fold
      }
      
      
      tmp_res
   }
   
   plan(sequential)
   return(pse)
}



################################################################################
# EMPIRICAL STUDY
################################################################################

# Load and prepare data
data <- read_data_func()
Y <- data$Y
X_list <- data$X_list

dx <- data$dx
tx.list <- data$tx.list
x.knots.list <- data$x.knots.list
x.order.list <- data$x.order.list
y.knots <- data$y.knots
y.order <- data$y.order
ty <- data$ty



pse <- pred_error_cv_func_all(X_list,Y,n_fold=n_fold,n_sim=n_sim)

# Average results across replications
pse.m <- Reduce("+",pse)/length(pse)
colnames(pse.m) <- c("FEPLS","PCR","PLS","u_env_dir","u_pcr","u_pls")
# Save results
if(!file.exists("results_DA_functional.rds")){
   saveRDS(list(pse=pse,pse.m=pse.m),file = "results_DA_functional.rds")
}


################################################################################
# VISUALIZATION
################################################################################
library(ggplot2)
data <- data.frame((pse.m[,1:3]))
data$x <- 1:dx
p1 <- ggplot(data) +
   geom_point(aes(x = x, y = FEPLS, color = "FEPLS", group = 1), size = 2) +  
   geom_line(aes(x = x, y = FEPLS, color = "FEPLS", group = 1), linewidth = 1) +    
   geom_point(aes(x = x, y = PCR, color = "PCR", group = 2), size = 2) +  
   geom_line(aes(x = x, y = PCR, color = "PCR", group = 2), linewidth = 1) +    
   geom_point(aes(x = x, y = PLS, color = "PLS", group = 3), size = 2) +  
   geom_line(aes(x = x, y = PLS, color = "PLS", group = 3), linewidth = 1) +    
   labs(x = "Number of components", y = "Mean squared prediction error", color = "Method") +
   scale_color_manual(values = c("FEPLS" = "black", "PCR" = "red", "PLS" = "green")) +
   theme_minimal()

p1

if(!file.exists("DA_functional.pdf")){
   ggsave("DA_functional.pdf",p1,width = 6,height = 4)
}


