################################################################################
# Basic utility functions for FEPLS.R
################################################################################
cbindl <- function(X.list){
   # Column-bind all matrices in list
   X <- NULL
   for (x in X.list) {
      X <- cbind(X,x)
   }
   return(X)
}

dot <- function(a,b){
   # Dot product of vectors a and b
   return(sum(a*b))
}

normv <- function(z){
   # L2 norm of vector z
   return(sqrt(sum(z^2)))
}

gauss_elimination <- function(A,tol=1e-12){
   #' Gaussian elimination with partial pivoting
   #'
   #' @param A Numeric matrix (r x u, full column rank, r >= u)
   #' @param tol Tolerance for rank deficiency check (default: 1e-12)
   #' @return Integer vector of row indices where first rank(A) are pivot rows
   #' @examples
   #' A <- matrix(c(1,2,3,4,5,6), nrow=3)
   #' gauss_elimination(A)
   r <- dim(A)[1]
   u <- dim(A)[2]
   idx <- rep(0,u)
   ridx <- 1:r
   
   for (i in 1:u) {
      # Find pivot row (largest absolute value in column i)
      sort_idx <- sort(abs(A[,i]),decreasing = TRUE,index.return=TRUE)$ix
      max_idx <- setdiff(sort_idx,idx)[1]
      
      if(abs(A[max_idx,i])<tol){
         stop("The matrix is not full rank!")
      }
      
      idx[i] <- max_idx
      ridx <- setdiff(ridx,max_idx)
      # Eliminate below pivot
      for (j in 1:(r-i)) {
         A[ridx[j],] <- A[ridx[j],]-A[max_idx,]*A[ridx[j],i]/A[max_idx,i]
      }
   }
   c(idx,ridx)
}

################################################################################
# Functions for estimating the coordinates of the function
################################################################################

bvfit <- function(X,basis.value){
   #' OLS projection onto basis functions
   #'
   #' @param X Numeric matrix (n x p1) to project
   #' @param basis.value Basis matrix (p1 x p2, p1 > p2)
   #' @return Matrix of coordinates (n x p2)
   X <- as.matrix(X)
   basis.value <- as.matrix(basis.value)
   X.cord <- X%*%basis.value%*%solve(t(basis.value)%*%basis.value)
   return(X.cord)
}

bvfitl <- function(X.list,basis.value,list=TRUE){
   #' OLS projection for list of matrices
   #'
   #' @param X.list List of matrices to project
   #' @param basis.value Basis matrix or list of basis matrices
   #' @param list Logical: if TRUE, use basis.value[[i]] for X.list[[i]]; 
   #'   if FALSE, use same basis.value for all
   #' @return Matrix of concatenated coordinates
   X.cord <- NULL
   if(list){
      n <- length(X.list)
      basis.value.list <- basis.value
      if(n!=length(basis.value.list))
         stop("X.list and basis.value.list must have the same length.")
      for (i in 1:n) {
         x <- as.matrix(X.list[[i]])
         basis.value <- as.matrix(basis.value.list[[i]])
         x.cord <- bvfit(x,basis.value)
         X.cord <- cbind(X.cord,x.cord)
      }
   }
   else{
      n <- length(X.list)
      for (i in 1:n) {
         x <- as.matrix(X.list[[i]])
         x.cord <- bvfit(x,basis.value)
         X.cord <- cbind(X.cord,x.cord)
      }
   }
   
   return(X.cord)
}


get_coord_dir_sy_sp <- function(X.list,tx.list,x.knots.list,x.order.list){
   #' Get coordinates with respect to orthogonal spline basis
   #'
   #' @param X.list List of data matrices to project
   #' @param tx.list List of observation points for each predictor
   #' @param x.knots.list List of knot locations for spline basis
   #' @param x.order.list List of spline orders
   #' @return List containing:
   #'   \item{X.cord}{Estimated coordinates}
   #'   \item{OX}{Block-diagonal orthogonalization matrix}
   #'   \item{basis.value.x.list}{List of basis function matrices}
   
   p <- length(X.list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X.list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Coordinate of X
   basis.value.x.list <- NULL
   Ox.list <- NULL
   for (i in 1:p) {
      # Obtain orthogonal basis for predictor i
      knots <- x.knots.list[[i]]
      tx <- tx.list[[i]]
      x.order <- x.order.list[[i]]
      
      if (min(knots) < 0) stop("Location of the knots cannot be negative.")
      if (is.null(knots)) {
         knots <- orthogonalsplinebasis::expand.knots(c(0, 0.25, 0.5, 0.75, 1), order = x.order)
      } else {
         knots <- knots / max(knots)
         knots <- orthogonalsplinebasis::expand.knots(knots, order = x.order)
      }
      
      spbasis <- orthogonalsplinebasis::SplineBasis(knots,order=x.order)
      Gx <- orthogonalsplinebasis::GramMatrix(spbasis)
      eigX <- eigen(Gx)
      Gx12 <- eigX$vectors %*% diag(sqrt(eigX$values)) %*% t(eigX$vectors)
      Gx12inv <- eigX$vectors %*% diag(1/sqrt(eigX$values)) %*% t(eigX$vectors)
      # spbasis.x.list <- append(spbasis.x.list,spbasis)
      basis.value.x <- orthogonalsplinebasis::evaluate(spbasis,tx) %*% Gx12inv
      Ox.list <- append(Ox.list,list(Gx12inv))
      basis.value.x.list <- append(basis.value.x.list,list(basis.value.x))
   }
   
   X.cord <- bvfitl(X.list, basis.value.x.list)
   
   OX <- Matrix::bdiag(Ox.list)
   return(list(X.cord=X.cord, OX=OX, basis.value.x.list=basis.value.x.list))
}

generate_basis_value <- function(knots,norder,t,Oy){
   #' Generate orthogonal spline basis values
   #'
   #' @param knots Knot locations
   #' @param norder Spline order
   #' @param t Evaluation points
   #' @param Oy Orthogonalization matrix
   #' @return Matrix of basis function values at points t
   
   if (min(knots) < 0) stop("Location of the knots cannot be negative.")
   if (is.null(knots)) {
      knots <- orthogonalsplinebasis::expand.knots(c(0, 0.25, 0.5, 0.75, 1), order = norder)
   } else {
      knots <- knots / max(knots)
      knots <- orthogonalsplinebasis::expand.knots(knots, order = norder)
   }
   
   spbasis <- orthogonalsplinebasis::SplineBasis(knots,order=norder)
   basis.value <- orthogonalsplinebasis::evaluate(spbasis,t) %*% Oy
   return(basis.value)
}

################################################################################
# Cross-validation functions
################################################################################
CV_pcr <- function(X,Y,n_fold=5){
   #' Cross-validation for principal component regression
   #'
   #' @param X Predictor matrix (n x p)
   #' @param Y Response matrix (n x r)
   #' @param n_fold Number of CV folds (default: 5)
   #' @return Optimal number of components
   n <- dim(X)[1]
   p <- dim(X)[2]
   n_test <- n/n_fold
   n_train <- n-n_test
   
   folds <- split(sample(1:n),sample(rep(1:n_fold, length.out = n)))
   m <- min(n_train-1,p)
   mse <- rep(0,m)
   
   for (i in 1:n_fold) {
      idx_test <- folds[[i]]
      idx_train <- setdiff(1:n,idx_test)
      X_train <- X[idx_train,]
      Y_train <- as.matrix(Y[idx_train,])
      X_test <- X[idx_test,]
      Y_test <- as.matrix(Y[idx_test,])
      for (j in 1:m) {
         reg <- pls::pcr(Y_train~X_train,ncomp=j)
         mse[j] <- mse[j] + MSEPs(reg,j,X_test,Y_test)
      }
   }
   ncomp <- which.min(mse)
   return(ncomp)
   
}


CV_pls <- function(X,Y,n_fold=5){
   #' Cross-validation for partial least squares regression
   #'
   #' @param X Predictor matrix (n x p)
   #' @param Y Response matrix (n x r)
   #' @param n_fold Number of CV folds (default: 5)
   #' @return Optimal number of components
   n <- dim(X)[1]
   p <- dim(X)[2]
   n_test <- n/n_fold
   n_train <- n-n_test
   
   folds <- split(sample(1:n),sample(rep(1:n_fold, length.out = n)))
   m <- min(n_train-1,p)
   mse <- rep(0,m)
   
   for (i in 1:n_fold) {
      idx_test <- folds[[i]]
      idx_train <- setdiff(1:n,idx_test)
      X_train <- X[idx_train,]
      Y_train <- as.matrix(Y[idx_train,])
      X_test <- X[idx_test,]
      Y_test <- as.matrix(Y[idx_test,])
      for (j in 1:m) {
         reg <- pls::plsr(Y_train~X_train,ncomp=j,method = "simpls")
         mse[j] <- mse[j] + MSEPs(reg,j,X_test,Y_test)
      }
   }
   ncomp <- which.min(mse)
   return(ncomp)
}

u_gpls <- function(X,Y,n_fold=5){
   #' Cross-validation for generalized PLS (classification)
   #'
   #' @param X Predictor matrix (n x p)
   #' @param Y Response vector (n x 1, class labels)
   #' @param n_fold Number of CV folds (default: 5)
   #' @return Optimal number of components
   #' 
   X <- as.matrix(X)
   Y <- as.matrix(Y)
   n <- dim(X)[1]
   p <- dim(X)[2]
   n_test <- n/n_fold
   n_train <- n-n_test
   m <- min(n_train,p)
   
   folds <- split(sample(1:n),sample(rep(1:n_fold, length.out = n)))
   
   mse <- rep(0,m)
   
   for (i in 1:n_fold) {
      idx_test <- folds[[i]]
      idx_train <- setdiff(1:n,idx_test)
      X_train <- X[idx_train,]
      Y_train <- Y[idx_train,]
      X_test <- X[idx_test,]
      Y_test <- Y[idx_test,]
      for (j in 1:m) {
         gpls1 <- gpls::gpls(X_train,Y_train,j)
         # gpls(Y_train~X_train,K.prov = i)
         y_test_pred <- pred_logit(coef(gpls1)[-1],coef(gpls1)[1],X_test)
         mse[j] <- mse[j] + mean(Y_test!=y_test_pred)
      }
   }
   ncomp <- which.min(mse)
   return(ncomp)
}


################################################################################
# BIC Based Dimension Selection
################################################################################

# Select the dimension by BIC for generalized envelope linear model (GMELM)
u_env_glm <- function(X,Y,init=NULL,type="logit",tol_inc=3){
   
   #' BIC-based dimension selection for generalized envelope linear model
   #'
   #' @param X Predictor matrix (n x p)
   #' @param Y Response matrix (n x r)
   #' @param init Initial values
   #' @param type Model type (default: "logit")
   #' @param tol_inc Tolerance for consecutive BIC increases before stopping (default: 3)
   #' @return List containing:
   #'   \item{u}{Optimal dimension}
   #'   \item{bic}{Minimum BIC value}
   #'   \item{bic_seq}{Sequence of BIC values}
   
   # Check arguments
   X <- as.matrix(X)
   Y <- as.matrix(Y)
   n <- dim(X)[1]
   p <- dim(X)[2]
   
   bic_min <- Inf
   num_inc <- 0
   bic_seq <- NULL
   
   u <- NULL
   for (i in 0:p) {
      res <- env_glm(X,Y,i,asy = asy, init = init)
      # BIC = penalty term + 2*objective
      bic <- log(n)*(p*(p+1)/2+i+p)+2*res$obj
      bic_seq <- c(bic_seq,bic)
      if(bic<bic_min){
         bic_min <- bic
         u <- i
      } else {
         num_inc <- num_inc + 1
      }
      if(num_inc>tol_inc-1){
         break
      }
   }
   return(list(u=u,bic=bic_min,bic_seq=bic_seq))
}


################################################################################
# Prediction and error calculation
################################################################################
MSEP <- function(model,ncomp,basis.value,X_new,Y_new){
   #' Mean squared error for PCR/PLS with functional response
   #'
   #' @param model Fitted PCR or PLS model
   #' @param ncomp Number of components
   #' @param basis.value Basis function matrix
   #' @param X_new New predictor matrix
   #' @param Y_new New response matrix
   #' @return Total mean squared error
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   if(is.null(ncomp))
      ncomp <- model$ncomp
   Y_hat <- predict(model,X_new,ncomp=ncomp)[,,1]%*%t(basis.value)
   return(sum(colMeans((Y_hat-Y_new)^2)))
}


MSEF <- function(beta,mu,basis.value,X_new,Y_new){
   #' Mean squared error for functional regression with functional response
   #'
   #' @param beta Coefficient matrix
   #' @param mu Intercept vector
   #' @param basis.value Basis function matrix
   #' @param X_new New predictor matrix
   #' @param Y_new New response matrix
   #' @return Total mean squared error
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   Y_hat <- (X_new %*% beta+matrix(1,n,1)%*%t(mu))%*%t(basis.value)
   return(sum(colMeans((Y_hat-Y_new)^2)))
}


MSEPs <- function(model,ncomp,X_new,Y_new){
   #' Mean squared error for PCR/PLS with scalar response
   #'
   #' @param model Fitted PCR or PLS model
   #' @param ncomp Number of components
   #' @param X_new New predictor matrix
   #' @param Y_new New response matrix
   #' @return Total mean squared error
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   
   if(is.null(ncomp))
      ncomp <- model$ncomp
   Y_hat <- predict(model,X_new,ncomp=ncomp)[,,1]
   return(sum(colMeans((Y_hat-Y_new)^2)))
}


MSEFs <- function(beta,mu,X_new,Y_new){
   #' Mean squared error for functional regression with scalar response
   #'
   #' @param beta Coefficient matrix
   #' @param mu Intercept scalar or vector
   #' @param X_new New predictor matrix
   #' @param Y_new New response matrix
   #' @return Total mean squared error
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   Y_hat <- (X_new %*% beta+matrix(1,n,1)%*%t(mu))
   return(sum(colMeans((Y_hat-Y_new)^2)))
} 


pred <- function(beta,mu,X_new){
   #' Predict response coefficients with respect to spline basis
   #'
   #' @param beta Coefficient matrix
   #' @param mu Intercept vector
   #' @param X_new New predictor matrix
   #' @return Matrix of predicted response coefficients
   n <- dim(X_new)[1]
   Y_hat.coef <- (X_new %*% beta+matrix(1,n,1)%*%t(mu))
   return(Y_hat.coef)
}

MMAC_F <- function(model,X_new,Y_new){
   #' Misclassification rate for classification models
   #'
   #' @param model Fitted classification model
   #' @param X_new New predictor matrix
   #' @param Y_new New response vector (class labels)
   #' @return Misclassification rate
   y_pred <- predict(model,newdata=X_new)$class
   return(mean(Y_new!=y_pred))
}

pred_logit_p <- function(beta,mu,X_new){
   #' Predicted probabilities for logistic regression
   #'
   #' @param beta Coefficient vector
   #' @param mu Intercept
   #' @param X_new New predictor matrix
   #' @return Vector of predicted probabilities
   n <- dim(X_new)[1]
   if(is.null(beta)){
      eta <- rep(mu,n)
   } else{
      eta <- X_new%*%beta+mu
   }
   plot(1:n,1/(1+exp(-eta)))
   return(1/(1+exp(-eta)))
}

pred_logit <- function(beta,mu,X_new){
   #' Predicted classes for logistic regression
   #'
   #' @param beta Coefficient vector
   #' @param mu Intercept
   #' @param X_new New predictor matrix
   #' @return Logical vector of predicted classes (threshold = 0.5)
   n <- dim(X_new)[1]
   if(is.null(beta)){
      eta <- rep(mu,n)
   } else{
      eta <- X_new%*%beta+mu
   }
   # plot(1:n,1/(1+exp(-eta)))
   return(1/(1+exp(-eta))>0.5)
}

pred_acc <- function(beta,mu,X_new,Y_new){
   #' Prediction accuracy for logistic regression
   #'
   #' @param beta Coefficient vector
   #' @param mu Intercept
   #' @param X_new New predictor matrix
   #' @param Y_new New response vector (true class labels)
   #' @return Prediction accuracy (proportion correct)
   Y_pred <- pred_logit(beta,mu,X_new)
   Y_new <- as.logical(Y_new)
   return(mean(Y_new==Y_pred))
}