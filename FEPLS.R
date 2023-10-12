bvfit <- function(X,basis.value){
   # OLS of X on basis.value 
   # basis.value is p1xp2 and X is nxp1
   # p1>p2
   X <- as.matrix(X)
   basis.value <- as.matrix(basis.value)
   X.cord <- X%*%basis.value%*%solve(t(basis.value)%*%basis.value)
   return(X.cord)
}

bvfitr <- function(x,bv,lambda=0.1){
   # L2 Regularization LS of X on bv (basis.value)
   # basis.value is p1xp2 and X is nxp1
   X <- as.matrix(X)
   basis.value <- as.matrix(basis.value)
   p2 <- dim(basis.value)[2]
   X.cord <- X%*%basis.value%*%solve(t(basis.value)%*%basis.value+lambda*diag(p2))
   return(X.cord)
}

bvfitrl <- function(X_list,basis.value, list=TRUE){
   # If list is TRUE, L2 Regularization LS   of X_list on basis.value(list) respectively
   # If list is FALSE, L2 Regularization LS of X_list on basis.value
   X.cord <- NULL
   if(list){
      n <- length(X_list)
      basis.value.list <- basis.value
      if(n!=length(basis.value.list))
         stop("X_list and basis.value.list must have the same length.")
      for (i in 1:n) {
         x <- as.matrix(X_list[[i]])
         basis.value <- as.matrix(basis.value.list[[i]])
         x.cord <- bvfitr(x,basis.value)
         X.cord <- cbind(X.cord,x.cord)
      }
   }
   else{
      n <- length(X_list)
      for (i in 1:n) {
         x <- as.matrix(X_list[[i]])
         x.cord <- bvfitr(x,basis.value)
         X.cord <- cbind(X.cord,x.cord)
      }
   }
   
   return(X.cord)
}

bvfitl <- function(X_list,basis.value,list=TRUE){
   # If list is TRUE, OLS of X_list on basis.value(list) respectively
   # If list is FALSE, OLS of X_list on basis.value
   X.cord <- NULL
   if(list){
      n <- length(X_list)
      basis.value.list <- basis.value
      if(n!=length(basis.value.list))
         stop("X_list and basis.value.list must have the same length.")
      for (i in 1:n) {
         x <- as.matrix(X_list[[i]])
         basis.value <- as.matrix(basis.value.list[[i]])
         x.cord <- bvfit(x,basis.value)
         X.cord <- cbind(X.cord,x.cord)
      }
   }
   else{
      n <- length(X_list)
      for (i in 1:n) {
         x <- as.matrix(X_list[[i]])
         x.cord <- bvfit(x,basis.value)
         X.cord <- cbind(X.cord,x.cord)
      }
   }
   
   return(X.cord)
}

cbindl <- function(X_list){
   # Extension of cbind on list
   X <- NULL
   for (x in X_list) {
      X <- cbind(X,x)
   }
   return(X)
}

dot <- function(a,b){
   return(sum(a*b))
}

normv <- function(z){
   # L2 norm for vector z
   return(sqrt(sum(z^2)))
}

gauss_elimination <- function(A,tol=1e-12){
   # Gaussian elimination
   #
   # Args:
   #   A: r x u full colmun rank matrix, where r is not less than u 
   #
   # Returns:
   #   The index: the first rank(A) indexes are indexes of non zero row after Guassian elimination
   r <- dim(A)[1]
   u <- dim(A)[2]
   idx <- rep(0,u)
   ridx <- 1:r
   
   for (i in 1:u) {
      sort_idx <- sort(abs(A[,i]),decreasing = TRUE,index.return=TRUE)$ix
      max_idx <- setdiff(sort_idx,idx)[1]
      
      if(abs(A[max_idx,i])<tol){
         stop("The matrix is not full rank!")
      }
      
      idx[i] <- max_idx
      ridx <- setdiff(ridx,max_idx)
      for (j in 1:(r-i)) {
         A[ridx[j],] <- A[ridx[j],]-A[max_idx,]*A[ridx[j],i]/A[max_idx,i]
      }
   }
   c(idx,ridx)
}

get_coord_dir_sy_sp <- function(X_list,tx.list,x.knots.list,x.order.list){
   # Get the estimated coordinate w.r.t. spline basis
   # tx.list, x.knots.list, and x.order.list are the arguments to help specify 
   # the spline baiss
   
   # Returns: 
   #   The estimated coordinate
   #   The orthogonalization matrix.
   #   The matrix of basis function at the give observation points tx.list
   p <- length(X_list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X_list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Coordinate of X
   basis.value.x.list <- NULL
   Ox.list <- NULL
   for (i in 1:p) {
      # Obtain an orthogonal basis list for x
      knots <- x.knots.list[[i]]
      tx <- tx.list[[i]]
      x.order <- x.order.list[[i]]
      
      # TODO
      # a check of order 
      
      
      # wrap the below as a function #
      # TODO 
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
   
   X.cord <- bvfitl(X_list, basis.value.x.list)
   
   OX <- Matrix::bdiag(Ox.list)
   return(list(X.cord=X.cord, OX=OX, basis.value.x.list=basis.value.x.list))
}

# # Dimension select via BIC
# u_stenv <- functionfunction(X,Y,asy=TRUE,init=NULL,tol_inc=3){
#    X <- as.matrix(X)
#    Y <- as.matrix(Y)
#    n <- dim(X)[1]
#    p <- dim(X)[2]
#    r <- dim(Y)[2]
#    bic_min <- Inf
#    num_inc <- 0
#    # tol_inc <- 3
#    uy <- NULL
#    ux <- NULL
#    
#    for (i in 1:r) {
#       for (j in 1:p) {
#          
#       }
#    }
# }

CV_pcr <- function(X,Y,n_fold=5){
   # Cross validation for choosing the number of components for pcr
   # Args:
   #   X: n x p predictor matrix
   #   Y: n x r response matrix
   #   n_fold: number of folds
   n <- dim(X)[1]
   p <- dim(X)[2]
   n_test <- n/n_fold
   n_train <- n-n_test
   
   m <- min(n_train-1,p)
   mse <- rep(0,m)
   
   for (j in 1:n_fold) {
      idx_test <- sample(1:n,n/n_test)
      idx_train <- setdiff(1:n,idx_test)
      X_train <- X[idx_train,]
      Y_train <- as.matrix(Y[idx_train,])
      X_test <- X[idx_test,]
      Y_test <- as.matrix(Y[idx_test,])
      for (i in 1:m) {
         reg <- pls::pcr(Y_train~X_train,ncomp=i)
         mse[i] <- mse[i] + MSEPs(reg,i,X_test,Y_test)
      }
   }
   ncomp <- which.min(mse)
   return(ncomp)
}

CV_pls <- function(X,Y,n_fold=5){
   # Cross validation for choosing the number of components for pls
   # Args:
   #   X: n x p predictor matrix
   #   Y: n x r response matrix
   #   n_fold: number of folds
   n <- dim(X)[1]
   p <- dim(X)[2]
   
   n_test <- n/n_fold
   n_train <- n-n_test
   
   m <- min(n_train-1,p)
   mse <- rep(0,m)
   
   for (j in 1:n_fold) {
      idx_test <- sample(1:n,n_test)
      idx_train <- setdiff(1:n,idx_test)
      X_train <- X[idx_train,]
      Y_train <- as.matrix(Y[idx_train,])
      X_test <- X[idx_test,]
      Y_test <- as.matrix(Y[idx_test,])
      for (i in 1:m) {
         reg <- pls::plsr(Y_train~X_train,ncomp=i,method = "simpls")
         mse[i] <- mse[i] + MSEPs(reg,i,X_test,Y_test)
      }
   }
   
   ncomp <- which.min(mse)
   return(ncomp)
}

# Multivariate Functional Envelope Linear Model with spline basis (Direct Method)
mfelmdir <- function(X_list,Y,ux,uy,tx.list,ty,x.knots.list,y.knots,x.order.list,y.order,u.selected=FALSE,alpha=.01){
   # Direct methods with spline basis
   # X_i and Y use the same spline bases
   
   # Arguments Check 
   p <- length(X_list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X_list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Coordinate of X
   tmp <- get_coord_dir_sy_sp(X_list=X_list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
   X.cord <- tmp$X.cord
   OX <- tmp$OX
   basis.value.x.list <- tmp$basis.value.x.list
   
   
   # Coordinate of Y
   Y <- as.matrix(Y)
   tmp <- get_coord_dir_sy_sp(list(Y),list(ty),list(y.knots),list(y.order))
   Y.cord <- tmp$X.cord
   Oy <- tmp$OX
   basis.value.Y <- tmp$basis.value.x.list[[1]]
   
   
   dx <- dim(X.cord)[2]
   dy <- dim(Y.cord)[2]
   
   if(u.selected){
      
      u.mat <- Renvlp::u.stenv(X.cord,Y.cord,alpha = alpha)
      ux <- u.mat$u.bic[1]
      uy <- u.mat$u.bic[2]
      
      
      ncomp.pcr <- CV_pcr(X.cord,Y.cord)
      mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
      
      
      ncomp.pls <- CV_pls(X.cord,Y.cord)
      mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp.pls,validation="CV")
   }
   else{
      ncomp <- ux
      mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp)
      mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp,validation="CV")
   }
   
   
   menv <- Renvlp::stenv(X.cord,Y.cord,ux,uy)
   mfull <- Renvlp::stenv(X.cord,Y.cord,dx,dy)
   
   
   # basis.value <- Matrix::bdiag(replicate(p,basis.value.x,simplify = FALSE))
   
   
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$mu,alphafull=mfull$mu,
               mpcr=mpcr,mpls=mpls,basis.value.Y=basis.value.Y,
               basis.value.X=basis.value.x.list,ux=ux,uy=uy,OX=OX,Oy=Oy))
}

# Multivariate Functional Predictor Envelope Linear Model with spline basis (Direct Method)
mfpedir <- function(X_list,Y,ux,tx.list,ty,x.knots.list,y.knots,x.order.list,y.order,u.selected=FALSE,alpha=.01){
   # Direct methods with spline basis
   # X_i and Y use the same spline bases
   
   # Arguments Check 
   p <- length(X_list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X_list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Coordinate of X
   tmp <- get_coord_dir_sy_sp(X_list=X_list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
   X.cord <- tmp$X.cord
   OX <- tmp$OX
   basis.value.x.list <- tmp$basis.value.x.list
   
   
   # Coordinate of Y
   Y <- as.matrix(Y)
   tmp <- get_coord_dir_sy_sp(list(Y),list(ty),list(y.knots),list(y.order))
   Y.cord <- tmp$X.cord
   Oy <- tmp$OX
   basis.value.Y <- tmp$basis.value.x.list[[1]]
   
   
   dx <- dim(X.cord)[2]
   dy <- dim(Y.cord)[2]
   
   if(u.selected){
      
      u.mat <- Renvlp::u.xenv(X.cord,Y.cord,alpha = alpha)
      ux <- u.mat$u.bic
      
      
      ncomp.pcr <- CV_pcr(X.cord,Y.cord)
      mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
      
      
      ncomp.pls <- CV_pls(X.cord,Y.cord)
      mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp.pls,validation="CV")
   }
   else{
      ncomp <- ux
      mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp)
      mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp,validation="CV")
   }
   
   
   menv <- Renvlp::xenv(X.cord,Y.cord,ux)
   mfull <- Renvlp::xenv(X.cord,Y.cord,dx)
   
   
   # basis.value <- Matrix::bdiag(replicate(p,basis.value.x,simplify = FALSE))
   
   
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$mu,alphafull=mfull$mu,
               mpcr=mpcr,mpls=mpls,basis.value.Y=basis.value.Y,
               basis.value.X=basis.value.x.list,ux=ux,OX=OX,Oy=Oy))
}



# Dimension selection via BIC #
u_xenv <- function(X,Y,asy=TRUE,init=NULL,tol_inc=3){
   # The asy and init is not available for now.
   
   # Check arguments
   # TODO
   
   
   X <- as.matrix(X)
   Y <- as.matrix(Y)
   n <- dim(X)[1]
   p <- dim(X)[2]
   r <- dim(Y)[2]
   
   # Check arguments
   if(dim(Y)[1]!=n)
      stop("The sample sizes of predictors and responses should be the same!(u_xenv)")
   bic_min <- Inf
   num_inc <- 0
   bic_seq <- NULL
   # tol_inc <- 3
   u <- NULL
   for (i in 0:p) {
      res <- Renvlp::xenv(X,Y,i,asy = asy, init = init)
      # BIC
      bic <- log(n)*(r+p+r*i+p*(p+1)/2+r*(r+1)/2) - 2*res$loglik
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


sfpedir <- function(X_list,Y,ux,tx.list,x.knots.list,x.order.list,u.selected=FALSE,alpha=0.01){
   # Functional predictor envelope with scalar (vector) response. Direct methods with spline basis
   # Args:
   #   
   
   # Arguments check
   p <- length(X_list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X_list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Coordinate of X
   tmp <- get_coord_dir_sy_sp(X_list=X_list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
   X.cord <- tmp$X.cord
   OX <- tmp$OX
   basis.value.x.list <- tmp$basis.value.x.list
   
   # Coordinate of Y
   Y <- as.matrix(Y)
   
   Y.cord <- Y
   # return(list(X.cord,Y.cord))
   dx <- dim(X.cord)[2]
   dy <- dim(Y.cord)[2]
   
   if(u.selected){
      u.mat <- Renvlp::u.xenv(X.cord,Y.cord)
      ux <- u.mat$u.bic[1]
      
      # u.mat <- u_xenv(X.cord,Y.cord)
      # ux <- u.mat$u
      
      #####
      # # Old method works for one dimensional scalar response
      # mpcr <- pls::pcr(Y.cord~X.cord,validation="CV")
      # ncomp.pcr <- pls::selectNcomp(mpcr,"onesigma")
      # mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
      # 
      # mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",validation="CV")
      # ncomp.pls <- pls::selectNcomp(mpls,"onesigma")
      # mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp.pls,validation="CV")
      #####
      if(dim(Y.cord)[2]==1){
         # Old method works for one dimensional scalar response
         mpcr <- pls::pcr(Y.cord~X.cord,validation="CV")
         ncomp.pcr <- pls::selectNcomp(mpcr,"onesigma")
         
         if(is.null(ncomp.pcr)|ncomp.pcr==0){
            ncomp.pcr <- CV_pcr(X.cord,Y.cord)
         }
         ncomp.pcr <- CV_pcr(X.cord,Y.cord)
         mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
         
         mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",validation="CV")
         ncomp.pls <- pls::selectNcomp(mpls,"onesigma")
         if(is.null(ncomp.pls)|ncomp.pls==0){
            ncomp.pls <- CV_pls(X.cord,Y.cord)
         }
         ncomp.pls <- CV_pls(X.cord,Y.cord)
         mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp.pls,validation="CV")
      }
      else{
         ncomp.pcr <- CV_pcr(X.cord,Y.cord)
         mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
         
         ncomp.pls <- CV_pls(X.cord,Y.cord)
         mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp.pls,validation="CV")
      }
   }
   else{
      ncomp <- ux
      mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp)
      mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp,validation="CV")
   }
   
   menv <- Renvlp::xenv(X.cord,Y.cord,ux)
   mfull <- Renvlp::xenv(X.cord,Y.cord,dx)
   
   # basis.value <- Matrix::bdiag(replicate(p,basis.value.x,simplify = FALSE))
   
   
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$mu,alphafull=mfull$mu,
               mpcr=mpcr,mpls=mpls,
               basis.value.X=basis.value.x.list,ux=ux,OX=OX,ncomp.pcr=mpcr$ncomp,ncomp.pls=mpls$ncomp,X.cord=X.cord,Y.cord=Y.cord))
}





# The envelope method for logistics regression
env_glm <- function(X,Y,u,type="logit",asy=TRUE,init=NULL){
   # The asy and init is not available for now
   
   X <- as.matrix(X)
   Y <- as.matrix(Y)
   n <- dim(X)[1]
   p <- dim(X)[2]
   
   SX <- stats::cov(X)*(n-1)/n
   SX_inv <- chol2inv(chol(SX))
   
   normv <- function(z){
      return(sqrt(sum(z^2)))
   }
   
   # Check arguments
   # TODO
   
   if(type=="logit"){
      
      # b function for Bernoulli
      b_function <- function(z){
         tmp <- function(z1){
            if(is.infinite(exp(z1))){
               return(log(1+exp(-z1))+z1)
            } else{
               return(log(1+exp(z1)))
            }
         }
         return(vapply(z, tmp, FUN.VALUE = numeric(1)))
      }
      if (u==p){
         Gamma=diag(p)
         Gamma0=NULL
         
         # Given Gamma, optimize it over alpha and eta
         glm1 <- glm(as.logical(Y)~X%*%Gamma,family=binomial(link='logit'))
         
         
         alpha <- glm1$coefficients[1]
         eta <- glm1$coefficients[2:(u+1)]
         
         
         # Objective function
         obj_func_G <- function(G){
            # Objective function with given eta and alpha.
            theta <- alpha+X%*%G%*%eta
            eig_values1 <- eigen(t(G)%*%SX%*%G)$values
            eig_values2 <- eigen(t(G)%*%SX_inv%*%G)$values
            eig_values3 <- eigen(SX)$values
            Lu <- -2/n*(sum(Y*theta-b_function(theta))) + sum(log(eig_values1)) +
               sum(log(eig_values2)) + sum(log(eig_values3))
            return(Lu)
         }
         beta <- Gamma%*%eta
         return(list(Gamma=Gamma,Gamma0=Gamma0,obj=obj_func_G(Gamma)*n/2,alpha=alpha,eta=eta,beta=beta))
      } else if (u==0){
         Gamma=NULL
         Gamma0=diag(p)
         # Given Gamma, optimize it over alpha and eta
         glm1 <- glm(as.logical(Y)~1,family=binomial(link='logit'))
         alpha <- glm1$coefficients[1]
         
         # Objective function
         obj_func_G <- function(G){
            theta <- alpha
            eig_values <- eigen(SX)$values
            Lu <- -2/n*(sum(Y*theta-b_function(theta))) + sum(log(eig_values)) 
            return(Lu)
         }
         return(list(Gamma=Gamma,Gamma0=Gamma0,obj=obj_func_G(Gamma)*n/2,alpha=alpha,eta=NULL,beta=NULL))
      } else if (u==1) {
         
         # Starting value via glmnet
         cv.reg <- cv.glmnet(x=X,y=Y,alpha=0,family="binomial",nfolds = 5)
         reg <- glmnet(x=X,y=Y,alpha=0,family="binomial",lambda = cv.reg$lambda.1se)
         beta_0 <- coef(reg)[2:(p+1)]
         if(is.null(init)){
            eigv <- eigen(SX)$vectors
            Gamma_init <- eigv[,sort(abs(beta_0%*%eigv),decreasing=TRUE,index.return=TRUE)$ix[1:u]]
            Gamma_init <- matrix(Gamma_init,nrow = p)
            Gamma <- Gamma_init
         }
         else{
            Gamma_init <- init
            Gamma <-  Gamma_init
         }
         obj_func <- function(G,alpha,eta){
            # Objective function with given eta and alpha.
            theta <- alpha+X%*%G%*%eta
            eig_values1 <- eigen(t(G)%*%SX%*%G)$values
            eig_values2 <- eigen(t(G)%*%SX_inv%*%G)$values
            eig_values3 <- eigen(SX)$values
            Lu <- -2/n*(sum(Y*theta-b_function(theta))) + sum(log(eig_values1)) +
               sum(log(eig_values2)) + sum(log(eig_values3))
            return(Lu)
         }
         
         obj_func_G <- function(G){
            # Objective function with given eta and alpha.
            theta <- alpha+X%*%G%*%eta
            eig_values1 <- eigen(t(G)%*%SX%*%G)$values
            eig_values2 <- eigen(t(G)%*%SX_inv%*%G)$values
            eig_values3 <- eigen(SX)$values
            Lu <- -2/n*(sum(Y*theta-b_function(theta))) + sum(log(eig_values1)) +
               sum(log(eig_values2)) + sum(log(eig_values3))
            return(Lu)
         }
         
         
         # initial obj value
         glm <- glm(as.logical(Y)~X%*%Gamma_init,family=binomial(link = "logit"))
         alpha <- glm$coefficients[1]
         eta <- glm$coefficients[2:(u+1)]
         obj_init <- obj_func(Gamma_init,alpha,eta)
         
         
         iter <- 0
         max_iter <- 30
         
         while (iter<max_iter) {
            iter <- iter + 1
            G <- Gamma
            
            # Give Gamma, optimize it over eta and alpha
            glm <- glm(as.logical(Y)~X%*%G,family=binomial(link = "logit"))
            alpha <- glm$coefficients[1]
            eta <- glm$coefficients[2:(u+1)]
            
            if(normv(eta)>1e7|abs(alpha)>1e7){
               break
            }
            
            # Give eta and alpha, optimize it over Gamma via reparameterization
            obj_func_col <- function(g){
               g <- as.matrix(g)
               theta <- alpha+X%*%g*eta/normv(g)
               Lu <- -2/n*(sum(Y*theta-b_function(theta))) + log(dot(g, SX%*%g)) +
                  log(dot(g, SX_inv%*%g)) -2*log(dot(g,g))
               return(Lu)
            }
            
            obj_func_col_grad <- function(g){
               g <- as.matrix(g)
               theta <- alpha+X%*%g*eta/normv(g)
               theta_grad <- eta*(X/normv(g)-X%*%g%*%t(g)/normv(g)^3)
               Lu_grad <- -2/n*(t(theta_grad)%*%Y-t(theta_grad)%*%(1/(1+exp(-theta)))) + 2*SX%*%g/dot(g, SX%*%g) +
                  2*SX_inv%*%g/dot(g, SX_inv%*%g)-4*g/dot(g,g)
               return(Lu_grad)
            }
            ##### DEBUG ####
            numeric_grad <- function(f,x,k,dx=1e-8){
               d <- length(x)
               e <- rep(0,d)
               e[k] <- 1
               return((f(x+dx*e)-f(x))/dx)
            }
            numeric_gradv <- function(f,x,dx=1e-8){
               d <- length(x)
               grad <- rep(0,d)
               for (i in 1:d) {
                  grad[i] <- numeric_grad(f,x,i,dx)
               }
               return(matrix(grad,ncol=1))
            }
            # for (j in 1:10) {
            #    a <- runif(p,-5,5)
            #    if(normv(numeric_gradv(obj_func_col,a,dx=1e-9)-obj_func_col_grad(a))>1e-4*normv(obj_func_col_grad(a))){
            #       # print(normv(numeric_gradv(obj_func_col,a,dx=1e-9)-obj_func_col_grad(a)))
            #       # print(1e-4*normv(obj_func_col_grad(a)))
            #       stop("Gradient does not match the numerical gradient!!")
            #    }
            # }
            #####
            res1 <- optim(G,obj_func_col,obj_func_col_grad,control=list(maxit=300),method = "L-BFGS")
            if(res1$convergence==1){
               print(u)
               print("Warning! Row optimization does not converge.")
            }
            if(sqrt(sum(res1$par^2))>1e4){
               print("Much")
            }
            if(obj_func_G(res1$par/normv(res1$par))-obj_init>1e-9){
               print(u)
               cat("1",obj_func_G(res1$par/normv(res1$par))-obj_init,"\n")
            }
            
            # glm <- glm(as.logical(Y)~X%*%G,family=binomial(link='logit'))
            # alpha <- glm$coefficients[1]
            # eta <- glm$coefficients[2:(u+1)]
            
            # print(obj_func(G,alpha,eta)-obj_func(Gamma_t,alpha_t,eta_t))
            
            G <- as.matrix(res1$par)/normv(res1$par)
            if(abs(obj_func_G(G)-obj_func_G(Gamma))<1e-3*abs(obj_func_G(Gamma))){
               Gamma <- G
               break
            } else{
               Gamma <- G
            }
         }
         
         glm1 <- glm(as.logical(Y)~X%*%Gamma,family=binomial(link='logit'))
         alpha <- glm1$coefficients[1]
         eta <- glm1$coefficients[2:(u+1)]
         
         Gamma0 <- lmreg::compbasis(Gamma)
         beta <- Gamma%*%eta
         return(list(Gamma=Gamma,Gamma0=Gamma0,obj=obj_func_G(Gamma)*n/2,alpha=alpha,eta=eta,beta=beta))
      } else {
         # Starting value via glmnet
         cv.reg <- cv.glmnet(x=X,y=Y,alpha=0,family="binomial",nfolds = 5)
         reg <- glmnet(x=X,y=Y,alpha=0,family="binomial",lambda = cv.reg$lambda.1se)
         beta_0 <- coef(reg)[2:(p+1)]
         if(is.null(init)){
            eigv <- eigen(SX)$vectors
            Gamma_init <- eigv[,sort(abs(beta_0%*%eigv),decreasing=TRUE,index.return=TRUE)$ix[1:u]]
            Gamma_init <- matrix(Gamma_init,nrow = p)
            Gamma <- Gamma_init
         }
         else{
            Gamma_init <- init
            Gamma <-  Gamma_init
         }
         obj_func <- function(G,alpha,eta){
            # Objective function with given eta and alpha.
            theta <- alpha+X%*%G%*%eta
            eig_values1 <- eigen(t(G)%*%SX%*%G)$values
            eig_values2 <- eigen(t(G)%*%SX_inv%*%G)$values
            eig_values3 <- eigen(SX)$values
            Lu <- -2/n*(sum(Y*theta-b_function(theta))) + sum(log(eig_values1)) +
               sum(log(eig_values2)) + sum(log(eig_values3))
            return(Lu)
         }
         
         obj_func_G <- function(G){
            # Objective function with given eta and alpha.
            theta <- alpha+X%*%G%*%eta
            eig_values1 <- eigen(t(G)%*%SX%*%G)$values
            eig_values2 <- eigen(t(G)%*%SX_inv%*%G)$values
            eig_values3 <- eigen(SX)$values
            Lu <- -2/n*(sum(Y*theta-b_function(theta))) + sum(log(eig_values1)) +
               sum(log(eig_values2)) + sum(log(eig_values3))
            return(Lu)
         }
         
         
         # initial obj value
         glm <- glm(as.logical(Y)~X%*%Gamma_init,family=binomial(link = "logit"))
         alpha <- glm$coefficients[1]
         eta <- glm$coefficients[2:(u+1)]
         obj_init <- obj_func(Gamma_init,alpha,eta)
         
         
         iter <- 0
         max_iter <- 100
         
         # TODO Preprocess data
         
         
         
         while (iter<max_iter) {
            iter <- iter + 1
            G <- Gamma
            GEidx <- gauss_elimination(G)
            GEidx1 <- GEidx[1:u]
            GEidx2 <- GEidx[(u+1):p]
            G1 <- G[GEidx1,]
            G2 <- G[GEidx2,]
            A <- G2%*%solve(G1) 
            # Find a G where G1 is symmetric and span(G) is the same.
            Iu <- diag(u) 
            eig1 <- eigen(Iu+t(A)%*%A)
            
            # Take G1 as the inverse of square root
            G1 <- eig1$vectors%*%diag(1/sqrt(eig1$values))%*%t(eig1$vectors)
            G2 <- A%*%G1
            G[GEidx1,] <- G1
            G[GEidx2,] <- G2
            
            
            
            
            # Give Gamma, optimize it over eta and alpha
            glm <- glm(as.logical(Y)~X%*%G,family=binomial(link = "logit"))
            alpha <- glm$coefficients[1]
            eta <- glm$coefficients[2:(u+1)]
            
            if(normv(eta)>1e7|abs(alpha)>1e7){
               break
            }
            
            
            # Give eta and alpha, optimize it over Gamma via reparameterization
            for (i in 1:(p-u)) {
               A1 <- A[setdiff(1:(p-u),i),]
               a_init <- A[i,]
               # row objective function
               obj_func_row <- function(a){
                  A <- matrix(0,p-u,u)
                  A[setdiff(1:(p-u),i),] <- A1
                  A[i,] <- a
                  
                  Iu <- diag(u)
                  eig1 <- eigen(Iu+t(A)%*%A)
                  # Take G1 as the inverse of square root
                  G1 <- eig1$vectors%*%diag(1/sqrt(eig1$values))%*%t(eig1$vectors)
                  G1_inv <- eig1$vectors%*%diag(sqrt(eig1$values))%*%t(eig1$vectors)
                  
                  G2 <- A%*%G1
                  G[GEidx1,] <- G1
                  G[GEidx2,] <- G2
                  CA <- G%*%G1_inv
                  CA1 <- CA[-GEidx2[i],]
                  
                  theta <- alpha+X%*%G%*%eta
                  
                  idxi <- GEidx2[i]
                  ridx <- setdiff(1:p,idxi)
                  
                  M <- SX
                  M12 <- as.matrix(M[idxi,ridx])
                  M11 <- M[ridx,ridx]
                  M22 <- M[idxi,idxi]
                  V <- SX_inv
                  V12 <- as.matrix(V[idxi,ridx])
                  V11 <- V[ridx,ridx]
                  V22 <- V[idxi,idxi]
                  
                  W1 <- t(CA1)%*%(M11-M12%*%t(M12)/M22)%*%CA1
                  W2 <- t(CA1)%*%(V11-V12%*%t(V12)/V22)%*%CA1
                  
                  
                  ######## DEBUG ########
                  if(min(abs(eigen(W1)$values))<1e-10){
                     stop("W1")
                  }
                  if(min(abs(eigen(W2)$values))<1e-10){
                     stop("W2")
                  }
                  if(min((eigen(W1)$values))<0){
                     stop("W1N")
                  }
                  if(min((eigen(W2)$values))<0){
                     stop("W2N")
                  }
                  
                  if(sum(a^2)>1e5){
                     # print(A1)
                     # print("2")
                  }
                  
                  
                  
                  
                  W1_inv <- chol2inv(chol(W1))
                  W2_inv <- chol2inv(chol(W2))
                  
                  a1 <- a + t(CA1)%*%M12/M22
                  a2 <- a + t(CA1)%*%V12/V22
                  
                  
                  
                  Lu <- -2/n*(sum(Y*theta-b_function(theta))) - 2*log(1+dot(a,chol2inv(chol(t(CA1)%*%CA1))%*%a)) +
                     log(1+M22*dot(a1,W1_inv%*%a1)) + log(1+V22*dot(a2,W2_inv%*%a2))
                  
                  return(Lu)
               }
               # gradient of 1d objective function
               obj_func_row_grad <- function(a){
                  Kuu <- matrixcalc::commutation.matrix(u,u)
                  D <- matrixcalc::D.matrix(u)
                  L <- matrixcalc::elimination.matrix(u) 
                  
                  A <- matrix(0,p-u,u)
                  A[setdiff(1:(p-u),i),] <- A1
                  A[i,] <- a
                  
                  Iu <- diag(u)
                  eig1 <- eigen(Iu+t(A)%*%A)
                  # Take G1 as the inverse of square root
                  G1 <-eig1$vectors%*%diag(1/sqrt(eig1$values))%*%t(eig1$vectors)
                  G1_inv <- eig1$vectors%*%diag(sqrt(eig1$values))%*%t(eig1$vectors)
                  
                  G2 <- A%*%G1
                  G[GEidx1,] <- G1
                  G[GEidx2,] <- G2
                  CA <- G%*%G1_inv
                  CA1 <- CA[-GEidx2[i],]
                  
                  idxi <- GEidx2[i]
                  ridx <- setdiff(1:p,idxi)
                  
                  # Gradient of theta w.r.t a
                  tmp1 <- Iu+t(A)%*%A
                  # tmp2 <- L%*%kronecker(tmp1,tmp1)%*%(kronecker(G1,Iu)+kronecker(Iu,G1)%*%Kuu)%*%D
                  
                  tmp2 <- G1%*%G1
                  tmp0 <- kronecker(G1,Iu)+kronecker(Iu,G1)
                  
                  
                  tmp3 <- -kronecker(t(eta),X%*%CA)%*%chol2inv(chol(tmp0))%*%kronecker(tmp2,tmp2)%*%
                     (kronecker(a,Iu)+kronecker(Iu,a))
                  C <- matrix(0,p*u,u)
                  C[cbind((1:u)*p-p+idxi,1:u)] <- 1
                  theta_grad <- tmp3+kronecker(t(G1%*%eta),X)%*%C
                  
                  theta <- alpha+X%*%G%*%eta
                  
                  # Gradient of second part
                  M <- SX
                  M12 <- as.matrix(M[idxi,ridx])
                  M11 <- M[ridx,ridx]
                  M22 <- M[idxi,idxi]
                  V <- SX_inv
                  V12 <- as.matrix(V[idxi,ridx])
                  V11 <- V[ridx,ridx]
                  V22 <- V[idxi,idxi]
                  
                  W1 <- t(CA1)%*%(M11-M12%*%t(M12)/M22)%*%CA1
                  W2 <- t(CA1)%*%(V11-V12%*%t(V12)/V22)%*%CA1
                  W1_inv <- chol2inv(chol(W1))
                  W2_inv <- chol2inv(chol(W2))
                  
                  ######## DEBUG ########
                  if(min(abs(eigen(W1)$values))<1e-10){
                     stop("W1")
                  }
                  if(min(abs(eigen(W2)$values))<1e-10){
                     stop("W2")
                  }
                  if(min((eigen(W1)$values))<0){
                     stop("W1N")
                  }
                  if(min((eigen(W2)$values))<0){
                     stop("W2N")
                  }
                  
                  if(sum(a^2)>1e8){
                     # print(A1)
                  }
                  
                  
                  ########################
                  
                  a1 <- a + t(CA1)%*%M12/M22
                  a2 <- a + t(CA1)%*%V12/V22
                  
                  
                  
                  Lu_grad <- -2/n*(t(theta_grad)%*%Y-t(theta_grad)%*%(1/(1+exp(-theta))))-4*chol2inv(chol(t(CA1)%*%CA1))%*%a/(1+dot(a,chol2inv(chol(t(CA1)%*%CA1))%*%a))+
                     1/(1+M22*dot(a1,W1_inv%*%a1))*2*M22*W1_inv%*%a1+1/(1+V22*dot(a2,W2_inv%*%a2))*2*V22*W2_inv%*%a2
                  
                  return(Lu_grad)
               }
               
               ######## DEBUG
               # od <- NULL
               obj <- function(a){
                  A <- matrix(0,p-u,u)
                  A[setdiff(1:(p-u),i),] <- A1
                  A[i,] <- a
                  Iu <- diag(u)
                  eig1 <- eigen(Iu+t(A)%*%A)
                  # Take G1 as the inverse of square root
                  G1 <- eig1$vectors%*%diag(1/sqrt(eig1$values))%*%t(eig1$vectors)
                  G2 <- A%*%G1
                  G[GEidx1,] <- G1
                  G[GEidx2,] <- G2
                  return(obj_func_G(G))
               }
               od <- NULL
               for (j in 1:10) {
                  a <- runif(u,-10,10)
                  if(is.null(od)){
                     od <- obj(a)-obj_func_row(a)
                  } else{
                     if(abs(od-(obj(a)-obj_func_row(a)))>1e-5*abs(od)){
                        # saveRDS(list(X=X,Y=Y,u=u),file = "SpecialCase.RData")
                        # stop("Wrong obj row function!!!!!")
                        print("1")
                     }
                  }
               }
               
               for (j in 1:10) {
                  a <- runif(u,-10,10)
               }
               
               # ######## DEBUG
               # numeric_grad <- function(f,x,k,dx=1e-8){
               #    d <- length(x)
               #    e <- rep(0,d)
               #    e[k] <- 1
               #    return((f(x+dx*e)-f(x))/dx)
               # }
               # numeric_gradv <- function(f,x,dx=1e-8){
               #    d <- length(x)
               #    grad <- rep(0,d)
               #    for (i in 1:d) {
               #       grad[i] <- numeric_grad(f,x,i,dx)
               #    }
               #    return(matrix(grad,ncol=1))
               # }
               # for (j in 1:10) {
               #    a <- runif(u,-5,5) 
               #    if(sum((numeric_gradv(obj_func_row,a,dx=1e-8)-obj_func_row_grad(a))^2)>1e-5*sqrt(sum(obj_func_row_grad(a)^2))){
               #       saveRDS(list(X=X,Y=Y,u=u),file = "SpecialCase.RData")
               #       stop("Gradient does not match the numerical gradient!!")
               #    }
               # } 
               
               
               
               
               
               res1 <- optim(a_init,obj_func_row,obj_func_row_grad, control=list(maxit=200),method = "L-BFGS")
               if(res1$convergence==1){
                  print(u)
                  print("Warning! Row optimization does not converge.")
               }
               if(sqrt(sum(res1$par^2))>1e4){
                  print("Much")
                  next
               }
               
               # if(obj_func_row(res1$par)>obj_func_row(a_init)){
               #    cat("Ascent",obj_func_row(res1$par)-obj_func_row(a_init),"\n")
               #    A[i,] <- a_init
               #    next
               # }
               
               if(obj(res1$par)-obj_init>1e-9){
                  # browser()
                  print(u)
                  print(eta)
                  print(alpha)
                  cat("Ascent",obj(res1$par)-obj_init,"\n")
                  A[i,] <- a_init
                  next
               }
               
               svd(t(G)%*%Gamma_init)
               A[i,] <- res1$par
            }
            Iu <- diag(u)
            eig1 <- eigen(Iu+t(A)%*%A)
            # Take G1 as the inverse of square root
            G1 <- eig1$vectors%*%diag(1/sqrt(eig1$values))%*%t(eig1$vectors)
            G2 <- A%*%G1
            
            G[GEidx1,] <- G1
            G[GEidx2,] <- G2
            
            # glm <- glm(as.logical(Y)~X%*%G,family=binomial(link='logit'))
            # alpha <- glm$coefficients[1]
            # eta <- glm$coefficients[2:(u+1)]
            
            # print(obj_func(G,alpha,eta)-obj_func(Gamma_t,alpha_t,eta_t))
            
            if(obj_func_G(G)-obj_init>1e-9){
               print(u)
               cat("1",obj_func_G(G)-obj_init,"\n")
            }
            if(abs(obj_func_G(G)-obj_func_G(Gamma))<1e-3*abs(obj_func_G(Gamma))){
               Gamma <- G
               break
            } else{
               Gamma <- G
            }
         }
         
         glm1 <- glm(as.logical(Y)~X%*%Gamma,family=binomial(link='logit'))
         alpha <- glm1$coefficients[1]
         eta <- glm1$coefficients[2:(u+1)]
         
         Gamma0 <- lmreg::compbasis(Gamma)
         beta <- Gamma%*%eta
         return(list(Gamma=Gamma,Gamma0=Gamma0,obj=obj_func_G(Gamma)*n/2,alpha=alpha,eta=eta,beta=beta))
      }
   }
}



# Select the dimension by BIC
u_env_glm <- function(X,Y,type="logit",asy=TRUE,init=NULL,tol_inc=3){
   
   # The asy and init is not available for now.
   
   # Check arguments
   X <- as.matrix(X)
   Y <- as.matrix(Y)
   n <- dim(X)[1]
   p <- dim(X)[2]
   
   bic_min <- Inf
   num_inc <- 0
   bic_seq <- NULL
   # tol_inc <- 3
   u <- NULL
   for (i in 0:p) {
      res <- env_glm(X,Y,i,asy = asy, init = init)
      # BIC
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





# Select dimension for gpls via cross validation
u_gpls <- function(X,Y,nfold=5){
   X <- as.matrix(X)
   Y <- as.matrix(Y)
   n <- dim(X)[1]
   p <- dim(X)[2]
   n_test <- n/nfold
   n_train <- n-n_test
   m <- min(n_train,p)
   
   mse <- rep(0,m)
   for (j in 1:nfold) {
      idx_test <- sample(1:n,n_test)
      idx_train <- setdiff(1:n,idx_test)
      X_train <- X[idx_train,]
      Y_train <- Y[idx_train,]
      X_test <- X[idx_test,]
      Y_test <- Y[idx_test,]
      for (i in 1:m) {
         gpls1 <- gpls::gpls(X_train,Y_train,i)
         # gpls(Y_train~X_train,K.prov = i)
         y_test_pred <- pred_logit(coef(gpls1)[-1],coef(gpls1)[1],X_test)
         mse[i] <- mse[i] + mean(Y_test!=y_test_pred)
      }
      
   }
   return(which.min(mse))
}

# Generalized function envelope linear model for categorical response
cfelmdir <- function(X_list,Y,ux,tx.list,x.knots.list,x.order.list,u.selected=FALSE,alpha=0.01,init=NULL,do.gpls=F){
   # Direct methods with spline basis
   # X_i and Y use the same spline bases
   # spbasis.x.list <- NULL
   p <- length(X_list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X_list, tx.list, and x.knots.
           list's length must be same.")
   }
   
   tmp <- get_coord_dir_sy_sp(X_list=X_list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
   
   X.cord <- tmp$X.cord
   OX <- tmp$OX
   basis.value.x.list <- tmp$basis.value.x.list
   
   # Coordinate of Y
   Y <- as.matrix(Y)
   Y.cord <- Y
   
   # return(list(X.cord,Y.cord))
   dx <- dim(X.cord)[2]
   dy <- dim(Y.cord)[2]
   
   if(u.selected){
      u.mat <- u_env_glm(X.cord,Y.cord)
      ux <- u.mat$u
      
      # u.mat1 <- u_env_glm1(X.cord,Y.cord)
      # ux1 <- u.mat1$u
      if(do.gpls){
         ncomp.gpls <- u_gpls(X.cord,Y.cord)
         mgpls <- gpls::gpls(X.cord,Y.cord,ncomp.gpls)
      }
      else{
         ncomp.gpls <- NULL
         mgpls <- NULL
      }
      
   }
   else{
      ncomp.gpls <- ux
      if(do.gpls){
         mgpls <- gpls::gpls(X.cord,Y.cord,ncomp.gpls)
      }
      else{
         mgpls <- NULL
      }
      ux1 <- ux
   }
   
   # return(list(dx,dy))
   menv <- env_glm(X.cord,Y.cord,ux,init=init)
   # menv1 <- env_glm1(X.cord,Y.cord,ux1,init=init)
   mfull <- env_glm(X.cord,Y.cord,dx,init=init)
   
   cv.reg <- cv.glmnet(x=X.cord,y=Y.cord,alpha=0,family="binomial",nfolds = 5)
   mglmnet_reg <- glmnet(x=X.cord,y=Y.cord,alpha=0,family="binomial",lambda = cv.reg$lambda.min)
   
   # return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$alpha,alphafull=mfull$alpha,
   #             basis.value.X=basis.value.x.list,ux=ux,OX=OX))
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$alpha,alphafull=mfull$alpha,
               basis.value.X=basis.value.x.list,ux=ux,OX=OX,mglmnet_reg=mglmnet_reg,mgpls=mgpls,ncomp.gpls=ncomp.gpls))
   # return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$alpha,alphafull=mfull$alpha,
   #             basis.value.X=basis.value.x.list,ux=ux,OX=OX,mglmnet_reg=mglmnet_reg,mgpls=mgpls,
   #             beta1=menv1$beta,alpha1=menv1$alpha,ux1=ux1))
}

# Multivariate Functional Predictor Envelope Linear Model with spline basis (KL Method)
mfelmKL <- function(X_list,Y,ux,uy,tx.list,ty,x.knots.list,y.knots,x.order.list,y.order,u.selected=FALSE,alpha=0.01){
   
   # KL expansion
   p <- length(X_list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X_list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Coordinate of X
   basis.value.x.list <- NULL
   GX12 <- NULL
   GX12inv <- NULL
   for (i in 1:p) {
      # Obtain an orthogonal basis list for x
      knots <- x.knots.list[[i]]
      tx <- tx.list[[i]]
      x.order <- x.order.list[[i]]
      
      # wrap the below as a function #
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
      if(is.null(GX12) | is.null(GX12inv)){
         GX12=Gx12
         GX12inv=Gx12inv
      }
      else{
         GX12 <- Matrix::bdiag(GX12,Gx12)
         GX12inv <- Matrix::bdiag(GX12inv,Gx12inv)
      }
      # spbasis.x.list <- append(spbasis.x.list,spbasis)
      basis.value.x <- orthogonalsplinebasis::evaluate(spbasis,tx)
      basis.value.x.list <- append(basis.value.x.list,list(basis.value.x))
   }
   
   X.cord <- bvfitl(X_list, basis.value.x.list)
   basis.value.X <- Matrix::bdiag(basis.value.x.list)
   
   # Coordinate of Y
   Y <- as.matrix(Y)
   
   if (min(y.knots) < 0) stop("Location of the knots cannot be negative.")
   if (is.null(y.knots)) {
      y.knots <- orthogonalsplinebasis::expand.knots(c(0, 0.25, 0.5, 0.75, 1), order = y.order)
   } else {
      y.knots <- y.knots / max(y.knots)
      y.knots <- orthogonalsplinebasis::expand.knots(y.knots, order = y.order)
   }
   
   spbasis.y <- orthogonalsplinebasis::SplineBasis(y.knots,order=y.order)
   Gy <- orthogonalsplinebasis::GramMatrix(spbasis.y)
   eigY <- eigen(Gy)
   GY12 <- eigY$vectors %*% diag(sqrt(eigY$values)) %*% t(eigY$vectors)
   GY12inv <- eigY$vectors %*% diag(1/sqrt(eigY$values)) %*% t(eigY$vectors)
   basis.value.Y <- orthogonalsplinebasis::evaluate(spbasis.y,ty)
   Y.cord <- bvfit(Y,basis.value.Y)
   
   
   tmpX <- eigen(GX12 %*% stats::cov(X.cord) %*% t(GX12))
   tmpY <- eigen(GY12 %*% stats::cov(Y.cord) %*% t(GY12))
   
   X.melm <- X.cord %*% GX12 %*% tmpX$vectors
   Y.melm <- Y.cord %*% GY12 %*% tmpY$vectors
   
   OX <- as.matrix(GX12 %*% tmpX$vectors)
   Oy <- as.matrix(GY12 %*% tmpY$vectors)
   X.melm <- as.matrix(X.melm)
   
   if(u.selected){
      
      u.mat <- Renvlp::u.stenv(X.melm,Y.melm,alpha = alpha)
      ux <- u.mat$u.bic[1]
      uy <- u.mat$u.bic[2]
      
      ncomp.pcr <- CV_pcr(X.melm,Y.melm)
      mpcr <- pls::pcr(Y.melm~X.melm,validation="CV",ncomp=ncomp.pcr)
      
      
      ncomp.pls <- CV_pls(X.melm,Y.melm)
      mpls <- pls::plsr(Y.melm~X.melm,method = "simpls",ncomp=ncomp.pls,validation="CV")
   }
   else{
      ncomp <- ux
      mpcr <- pls::pcr(Y.melm~X.melm,validation="CV",ncomp=ncomp)
      mpls <- pls::plsr(Y.melm~X.melm,method = "simpls",ncomp=ncomp,validation="CV")
   }
   
   phihat.cord <- GX12inv %*% tmpX$vectors
   psihat.cord <- GY12inv %*% tmpY$vectors
   phihat.X <- basis.value.X %*% phihat.cord
   psihat.Y <- basis.value.Y %*% psihat.cord
   
   
   dX <- dim(X.melm)[2]
   dY <- dim(Y.melm)[2]
   menv <- Renvlp::stenv(X.melm,Y.melm,ux,uy)
   mfull <- Renvlp::stenv(X.melm,Y.melm,dX,dY)
   # mpcr <- pls::pcr(Y.melm~X.melm)
   # mpls <- pls::plsr(Y.melm~X.melm)
   
   
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$mu,alphafull=mfull$mu,
               phihat.X=phihat.X,psihat.Y=psihat.Y,psihat.cord=psihat.cord
               ,mpcr=mpcr,mpls=mpls,Y.melm=Y.melm,X.melm=X.melm, ux=ux, uy=uy, OX=phihat.cord,
               Oy=psihat.cord))
}

# Multivariate Functional Predictor Envelope Linear Model with spline basis (KL Method)
mfpeKL <- function(X_list,Y,ux,tx.list,ty,x.knots.list,y.knots,x.order.list,y.order,u.selected=FALSE,alpha=0.01){
   
   # KL expansion
   p <- length(X_list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X_list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Coordinate of X
   basis.value.x.list <- NULL
   GX12 <- NULL
   GX12inv <- NULL
   for (i in 1:p) {
      # Obtain an orthogonal basis list for x
      knots <- x.knots.list[[i]]
      tx <- tx.list[[i]]
      x.order <- x.order.list[[i]]
      
      # wrap the below as a function #
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
      if(is.null(GX12) | is.null(GX12inv)){
         GX12=Gx12
         GX12inv=Gx12inv
      }
      else{
         GX12 <- Matrix::bdiag(GX12,Gx12)
         GX12inv <- Matrix::bdiag(GX12inv,Gx12inv)
      }
      # spbasis.x.list <- append(spbasis.x.list,spbasis)
      basis.value.x <- orthogonalsplinebasis::evaluate(spbasis,tx)
      basis.value.x.list <- append(basis.value.x.list,list(basis.value.x))
   }
   
   X.cord <- bvfitl(X_list, basis.value.x.list)
   basis.value.X <- Matrix::bdiag(basis.value.x.list)
   
   # Coordinate of Y
   Y <- as.matrix(Y)
   
   if (min(y.knots) < 0) stop("Location of the knots cannot be negative.")
   if (is.null(y.knots)) {
      y.knots <- orthogonalsplinebasis::expand.knots(c(0, 0.25, 0.5, 0.75, 1), order = y.order)
   } else {
      y.knots <- y.knots / max(y.knots)
      y.knots <- orthogonalsplinebasis::expand.knots(y.knots, order = y.order)
   }
   
   spbasis.y <- orthogonalsplinebasis::SplineBasis(y.knots,order=y.order)
   Gy <- orthogonalsplinebasis::GramMatrix(spbasis.y)
   eigY <- eigen(Gy)
   GY12 <- eigY$vectors %*% diag(sqrt(eigY$values)) %*% t(eigY$vectors)
   GY12inv <- eigY$vectors %*% diag(1/sqrt(eigY$values)) %*% t(eigY$vectors)
   basis.value.Y <- orthogonalsplinebasis::evaluate(spbasis.y,ty)
   Y.cord <- bvfit(Y,basis.value.Y)
   
   GX12 <- as.matrix(GX12)
   GX12inv <- as.matrix(GX12inv)
   
   tmpX <- eigen(GX12 %*% stats::cov(X.cord) %*% t(GX12))
   tmpY <- eigen(GY12 %*% stats::cov(Y.cord) %*% t(GY12))
   
   X.melm <- X.cord %*% GX12 %*% tmpX$vectors
   Y.melm <- Y.cord %*% GY12 %*% tmpY$vectors
   
   OX <- as.matrix(GX12 %*% tmpX$vectors)
   Oy <- as.matrix(GY12 %*% tmpY$vectors)
   X.melm <- as.matrix(X.melm)
   
   if(u.selected){
      
      u.mat <- Renvlp::u.xenv(X.melm,Y.melm,alpha = alpha)
      ux <- u.mat$u.bic
      
      ncomp.pcr <- CV_pcr(X.melm,Y.melm)
      mpcr <- pls::pcr(Y.melm~X.melm,validation="CV",ncomp=ncomp.pcr)
      
      
      ncomp.pls <- CV_pls(X.melm,Y.melm)
      mpls <- pls::plsr(Y.melm~X.melm,method = "simpls",ncomp=ncomp.pls,validation="CV")
   }
   else{
      ncomp <- ux
      mpcr <- pls::pcr(Y.melm~X.melm,validation="CV",ncomp=ncomp)
      mpls <- pls::plsr(Y.melm~X.melm,method = "simpls",ncomp=ncomp,validation="CV")
   }
   
   phihat.cord <- GX12inv %*% tmpX$vectors
   psihat.cord <- GY12inv %*% tmpY$vectors
   phihat.X <- basis.value.X %*% phihat.cord
   psihat.Y <- basis.value.Y %*% psihat.cord
   
   
   dX <- dim(X.melm)[2]
   dY <- dim(Y.melm)[2]
   menv <- Renvlp::xenv(X.melm,Y.melm,ux)
   mfull <- Renvlp::xenv(X.melm,Y.melm,dX)
   # mpcr <- pls::pcr(Y.melm~X.melm)
   # mpls <- pls::plsr(Y.melm~X.melm)
   
   
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$mu,alphafull=mfull$mu,
               phihat.X=phihat.X,psihat.Y=psihat.Y,psihat.cord=psihat.cord
               ,mpcr=mpcr,mpls=mpls,Y.melm=Y.melm,X.melm=X.melm, ux=ux, OX=phihat.cord,
               Oy=psihat.cord))
}


# Calculate the mse for new sample for PCR and PLS
MSEP <- function(model,ncomp,basis.value,X_new,Y_new){
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   if(is.null(ncomp))
      ncomp <- model$ncomp
   Y_hat <- predict(model,X_new,ncomp=ncomp)[,,1]%*%t(basis.value)
   return(sum(colMeans((Y_hat-Y_new)^2)))
}



# Calculate the mse for FFFR and FELM
MSEF <- function(beta,mu,basis.value,X_new,Y_new){
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   Y_hat <- (X_new %*% beta+matrix(1,n,1)%*%t(mu))%*%t(basis.value)
   return(sum(colMeans((Y_hat-Y_new)^2)))
}



# Calculate the mse for new sample for PCR and PLS for scalar
MSEPs <- function(model,ncomp,X_new,Y_new){
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   
   if(is.null(ncomp))
      ncomp <- model$ncomp
   Y_hat <- predict(model,X_new,ncomp=ncomp)[,,1]
   return(sum(colMeans((Y_hat-Y_new)^2)))
}

# Calculate the mse for FFFR and FELM for scalar
MSEFs <- function(beta,mu,X_new,Y_new){
   n <- dim(X_new)[1]
   if(n!=dim(Y_new)[1])
      print("X and must must have same number of row")
   Y_hat <- (X_new %*% beta+matrix(1,n,1)%*%t(mu))
   return(sum(colMeans((Y_hat-Y_new)^2)))
} 

pred <- function(beta,mu,X_new){
   # Return the coef of predicted Y w.r.t splines basis 
   n <- dim(X_new)[1]
   Y_hat.coef <- (X_new %*% beta+matrix(1,n,1)%*%t(mu))
   return(Y_hat.coef)
}


# The misclassifications accuracy
MMAC_F <- function(model,X_new,Y_new){
   y_pred <- predict(model,newdata=X_new)$class
   return(mean(Y_new!=y_pred))
}


# estimated prediction probability
pred_logit_p <- function(beta,mu,X_new){
   n <- dim(X_new)[1]
   if(is.null(beta)){
      eta <- rep(mu,n)
   } else{
      eta <- X_new%*%beta+mu
   }
   plot(1:n,1/(1+exp(-eta)))
   return(1/(1+exp(-eta)))
}

# Prediction for logitistics regression
pred_logit <- function(beta,mu,X_new){
   n <- dim(X_new)[1]
   if(is.null(beta)){
      eta <- rep(mu,n)
   } else{
      eta <- X_new%*%beta+mu
   }
   # plot(1:n,1/(1+exp(-eta)))
   return(1/(1+exp(-eta))>0.5)
}

generate_basis_value <- function(knots,norder,t,Oy){
   # wrap the below as a function #
   # TODO 
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


# Predicted Accuracy 
pred_acc <- function(beta,mu,X_new,Y_new){
   Y_pred <- pred_logit(beta,mu,X_new)
   Y_new <- as.logical(Y_new)
   return(mean(Y_new==Y_pred))
}
