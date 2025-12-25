################################################################################
# Functional Envelope-based Partial Least Squares (FEPLS) 
################################################################################
# This code provides an R implementation of Functional Envelope-based Partial 
# Least Squares (FEPLS), a dimension-reduction method leveraging envelope 
# methods to extract the essential variation in functional predictors relevant 
# to the response. FEPLS combines functional data analysis (FDA) with partial 
# least squares (PLS) regression, enhancing predictive accuracy by isolating 
# relevant functional information and minimizing noise.
#
# Key Features:
# - Applies envelope methods for dimension reduction, capturing essential 
#   variation in functional predictors
# - Supports three response types: binary, scalar, and functional
# - Uses natural cubic spline basis as the default for functional predictors
################################################################################
source("Basic_Functions.R")
################################################################################
# Generalized envelope linear model
################################################################################
env_glm <- function(X, Y, u, asy = asy, type="logit", init=NULL){
   #' Generalized envelope linear model for logistic regression
   #'
   #' Fits an envelope model for binary response using logistic regression with
   #' dimension reduction through envelope subspace estimation.
   #'
   #' @param X Predictor matrix (n x p)
   #' @param Y Binary response vector (n x 1)
   #' @param u Dimension of the envelope subspace
   #' @param asy Not implemented.
   #' @param type Model type (default: "logit")
   #' @param init Optional initial value for Gamma matrix
   #'
   #' @return List containing:
   #'   \item{Gamma}{Basis for envelope subspace (p x u)}
   #'   \item{Gamma0}{Basis for orthogonal complement}
   #'   \item{obj}{Objective function value}
   #'   \item{alpha}{Intercept parameter}
   #'   \item{eta}{Coefficients in envelope coordinates (u x 1)}
   #'   \item{beta}{Regression coefficients in original space (p x 1)}
   #'
   #' @examples
   #' X <- matrix(rnorm(100*5), 100, 5)
   #' Y <- rbinom(100, 1, 0.5)
   #' fit <- env_glm(X, Y, u=2)
   X <- as.matrix(X)
   Y <- as.matrix(Y)
   n <- dim(X)[1]
   p <- dim(X)[2]
   
   SX <- stats::cov(X)*(n-1)/n
   SX_inv <- chol2inv(chol(SX))
   
   normv <- function(z){
      return(sqrt(sum(z^2)))
   }
   
   
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
         # Full model case (u = p)
         
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
         # Null model case (u = 0)
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
         # One-dimensional envelope (u = 1)
         
         # Initialize using glmnet
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
         
         
         # Initialize parameters
         glm <- glm(as.logical(Y)~X%*%Gamma_init,family=binomial(link = "logit"))
         alpha <- glm$coefficients[1]
         eta <- glm$coefficients[2:(u+1)]
         obj_init <- obj_func(Gamma_init,alpha,eta)
         
         # Alternating optimization
         iter <- 0
         max_iter <- 30
         
         while (iter<max_iter) {
            iter <- iter + 1
            G <- Gamma
            
            # Update alpha and eta given Gamma
            glm <- glm(as.logical(Y)~X%*%G,family=binomial(link = "logit"))
            alpha <- glm$coefficients[1]
            eta <- glm$coefficients[2:(u+1)]
            
            if(normv(eta)>1e7|abs(alpha)>1e7){
               break
            }
            
            # Update Gamma given alpha and eta via reparameterization
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
            res1 <- optim(G,obj_func_col,obj_func_col_grad,control=list(maxit=300),method = "L-BFGS")
            
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
         # Multi-dimensional envelope (u > 1)
         
         # Initialize using glmnet
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
         
         
         # Initialize parameters
         glm <- glm(as.logical(Y)~X%*%Gamma_init,family=binomial(link = "logit"))
         alpha <- glm$coefficients[1]
         eta <- glm$coefficients[2:(u+1)]
         obj_init <- obj_func(Gamma_init,alpha,eta)
         
         # Alternating optimization with Grassmann manifold optimization
         iter <- 0
         max_iter <- 100
         
         
         while (iter<max_iter) {
            iter <- iter + 1
            G <- Gamma
            
            # Perform Gaussian elimination to identify pivot rows
            GEidx <- gauss_elimination(G)
            GEidx1 <- GEidx[1:u]
            GEidx2 <- GEidx[(u+1):p]
            G1 <- G[GEidx1,]
            G2 <- G[GEidx2,]
            A <- G2%*%solve(G1) 
            
            # Reparameterize G using A matrix
            Iu <- diag(u) 
            eig1 <- eigen(Iu+t(A)%*%A)
            # Take G1 as the inverse of square root
            G1 <- eig1$vectors%*%diag(1/sqrt(eig1$values))%*%t(eig1$vectors)
            G2 <- A%*%G1
            G[GEidx1,] <- G1
            G[GEidx2,] <- G2
            
            
            
            
            # Update alpha and eta given Gamma
            glm <- glm(as.logical(Y)~X%*%G,family=binomial(link = "logit"))
            alpha <- glm$coefficients[1]
            eta <- glm$coefficients[2:(u+1)]
            
            if(normv(eta)>1e7|abs(alpha)>1e7){
               break
            }
            
            
            # Update each row of A sequentially
            for (i in 1:(p-u)) {
               A1 <- A[setdiff(1:(p-u),i),]
               a_init <- A[i,]
               # Objective function for row i of A
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
                  

                  
                  W1_inv <- chol2inv(chol(W1))
                  W2_inv <- chol2inv(chol(W2))
                  
                  a1 <- a + t(CA1)%*%M12/M22
                  a2 <- a + t(CA1)%*%V12/V22

                  Lu <- -2/n*(sum(Y*theta-b_function(theta))) - 2*log(1+dot(a,chol2inv(chol(t(CA1)%*%CA1))%*%a)) +
                     log(1+M22*dot(a1,W1_inv%*%a1)) + log(1+V22*dot(a2,W2_inv%*%a2))
                  
                  return(Lu)
               }
               # Gradient of objective function for row i
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
                  
                  # Gradient of theta with respect to a
                  tmp1 <- Iu+t(A)%*%A
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
                  
                  
                  a1 <- a + t(CA1)%*%M12/M22
                  a2 <- a + t(CA1)%*%V12/V22
                  
                  
                  
                  Lu_grad <- -2/n*(t(theta_grad)%*%Y-t(theta_grad)%*%(1/(1+exp(-theta))))-4*chol2inv(chol(t(CA1)%*%CA1))%*%a/(1+dot(a,chol2inv(chol(t(CA1)%*%CA1))%*%a))+
                     1/(1+M22*dot(a1,W1_inv%*%a1))*2*M22*W1_inv%*%a1+1/(1+V22*dot(a2,W2_inv%*%a2))*2*V22*W2_inv%*%a2
                  
                  return(Lu_grad)
               }
               
               res1 <- optim(a_init,obj_func_row,obj_func_row_grad, control=list(maxit=200),method = "L-BFGS")
               
               svd(t(G)%*%Gamma_init)
               A[i,] <- res1$par
            }
            
            # Update G using optimized A
            Iu <- diag(u)
            eig1 <- eigen(Iu+t(A)%*%A)
            # Take G1 as the inverse of square root
            G1 <- eig1$vectors%*%diag(1/sqrt(eig1$values))%*%t(eig1$vectors)
            G2 <- A%*%G1
            
            G[GEidx1,] <- G1
            G[GEidx2,] <- G2
            
            
            # Check convergence
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
         
         # Gamma <- Gamma
         # eta <- eta
         
         beta <- Gamma%*%eta
         return(list(Gamma=Gamma,Gamma0=Gamma0,obj=obj_func_G(Gamma)*n/2,alpha=alpha,eta=eta,beta=beta))
      }
   }
}

################################################################################
# Functional Envelope-based PLS
################################################################################

mfpedir <- function(X.list,Y,ux,tx.list,ty,x.knots.list,y.knots,x.order.list,y.order,u.selected=FALSE,alpha=.01){
   #' Functional predictor envelope model with functional response (direct method)
   #'
   #' Fits an envelope model for functional response with functional predictors
   #' using spline basis representation and direct optimization.
   #'
   #' @param X.list List of functional predictor matrices
   #' @param Y Functional response matrix (n x m)
   #' @param ux Dimension of envelope subspace
   #' @param tx.list List of time points for each functional predictor
   #' @param ty Time points for functional response
   #' @param x.knots.list List of knot vectors for predictor splines
   #' @param y.knots Knot vector for response spline
   #' @param x.order.list List of spline orders for predictors
   #' @param y.order Spline order for response
   #' @param u.selected Logical; if TRUE, select u via BIC (default: FALSE)
   #' @param alpha Significance level for envelope dimension selection (default: 0.01)
   #'
   #' @return List containing fitted model components including beta coefficients,
   #'   basis functions, and comparison models (PCR, PLS, full model)

   # Check argument consistency
   p <- length(X.list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X.list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Estimate the coordinates for X
   tmp <- get_coord_dir_sy_sp(X.list=X.list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
   X.cord <- tmp$X.cord
   OX <- tmp$OX
   basis.value.x.list <- tmp$basis.value.x.list
   
   
   # Estimate the coordinates for Y
   Y <- as.matrix(Y)
   tmp <- get_coord_dir_sy_sp(list(Y),list(ty),list(y.knots),list(y.order))
   Y.cord <- tmp$X.cord
   Oy <- tmp$OX
   basis.value.Y <- tmp$basis.value.x.list[[1]]
   
   
   dx <- dim(X.cord)[2]
   dy <- dim(Y.cord)[2]
   
   if(u.selected){
      # Select envelope dimension using BIC
      u.mat <- Renvlp::u.xenv(X.cord,Y.cord,alpha = alpha)
      ux <- u.mat$u.bic
      
      # Fit comparison models with cross-validation
      ncomp.pcr <- CV_pcr(X.cord,Y.cord)
      mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
      
      
      ncomp.pls <- CV_pls(X.cord,Y.cord)
      mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp.pls,validation="CV")
   }
   else{
      ncomp <- ux
      ncomp.pcr <- ux
      ncomp.pls <- ux
      mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
      mpls <- pls::plsr(Y.cord~X.cord,method = "simpls",ncomp=ncomp.pls,validation="CV")
   }
   
   # Fit envelope and full models
   menv <- Renvlp::xenv(X.cord,Y.cord,ux)
   mfull <- Renvlp::xenv(X.cord,Y.cord,dx)
   
   
   
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$mu,alphafull=mfull$mu,
               mpcr=mpcr,mpls=mpls,basis.value.Y=basis.value.Y,
               basis.value.X=basis.value.x.list,ux=ux,OX=OX,Oy=Oy,ncomp.pcr=ncomp.pcr,ncomp.pls=ncomp.pls,menv=menv))
}

sfpedir <- function(X.list,Y,ux,tx.list,x.knots.list,x.order.list,u.selected=FALSE,alpha=0.01){
   #' Functional predictor envelope model with scalar response (direct method)
   #'
   #' Fits an envelope model for scalar (vector) response with functional predictors
   #' using spline basis representation and direct optimization.
   #'
   #' @param X.list List of functional predictor matrices
   #' @param Y Scalar response vector or matrix (n x r)
   #' @param ux Dimension of envelope subspace
   #' @param tx.list List of time points for each functional predictor
   #' @param x.knots.list List of knot vectors for predictor splines
   #' @param x.order.list List of spline orders for predictors
   #' @param u.selected Logical; if TRUE, select u via BIC (default: FALSE)
   #' @param alpha Significance level for envelope dimension selection (default: 0.01)
   #'
   #' @return List containing fitted model components including beta coefficients,
   #'   basis functions, and comparison models (PCR, PLS, full model)
   
   # Check argument consistency
   p <- length(X.list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X.list, tx.list, and x.knots.list's length must be same.")
   }
   
   # Estimate the coordinates for X
   tmp <- get_coord_dir_sy_sp(X.list=X.list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
   X.cord <- tmp$X.cord
   OX <- tmp$OX
   basis.value.x.list <- tmp$basis.value.x.list
   

   Y <- as.matrix(Y)
   
   Y.cord <- Y
   
   dx <- dim(X.cord)[2]
   dy <- dim(Y.cord)[2]
   
   if(u.selected){
      u.mat <- Renvlp::u.xenv(X.cord,Y.cord)
      ux <- u.mat$u.bic[1]
      
      
      if(dim(Y.cord)[2]==1){

         ncomp.pcr <- CV_pcr(X.cord,Y.cord)
         mpcr <- pls::pcr(Y.cord~X.cord,validation="CV",ncomp=ncomp.pcr)
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
   
   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$mu,alphafull=mfull$mu,
               mpcr=mpcr,mpls=mpls,
               basis.value.X=basis.value.x.list,ux=ux,OX=OX,ncomp.pcr=mpcr$ncomp,ncomp.pls=mpls$ncomp,X.cord=X.cord,Y.cord=Y.cord))
}

cfelmdir <- function(X.list,Y,ux,tx.list,x.knots.list,x.order.list,u.selected=FALSE,alpha=0.01,init=NULL,do.gpls=F){
   #' Functional predictor envelope model with categorical response (direct method)
   #'
   #' Fits an envelope model for binary response with functional predictors
   #' using spline basis representation and generalized envelope methods.
   #'
   #' @param X.list List of functional predictor matrices
   #' @param Y Binary response vector (n x 1)
   #' @param ux Dimension of envelope subspace
   #' @param tx.list List of time points for each functional predictor
   #' @param x.knots.list List of knot vectors for predictor splines
   #' @param x.order.list List of spline orders for predictors
   #' @param u.selected Logical; if TRUE, select u via custom criterion (default: FALSE)
   #' @param alpha Significance level (default: 0.01)
   #' @param init Optional initial value for Gamma matrix
   #' @param do.gpls Logical; if TRUE, fit generalized PLS model (default: FALSE)
   #'
   #' @return List containing fitted model components including beta coefficients,
   #'   basis functions, and comparison models (ridge regression, GPLS if requested)
   
   # Check argument consistency
   p <- length(X.list)
   if(p!=length(tx.list) | p!=length(x.knots.list)){
      stop("X.list, tx.list, and x.knots.
           list's length must be same.")
   }
   
   
   # Estimate the coordinates for X
   tmp <- get_coord_dir_sy_sp(X.list=X.list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
   X.cord <- tmp$X.cord
   OX <- tmp$OX
   basis.value.x.list <- tmp$basis.value.x.list
   
   # Coordinate of Y
   Y <- as.matrix(Y)
   
   Y.factor <- factor(Y)
   dx <- dim(X.cord)[2]
   
   if(u.selected){
      # Select envelope dimension using BIC
      u.mat <- u_env_glm(X.cord,Y)
      ux <- u.mat$u
      
      if(do.gpls){
         ncomp.gpls <- u_gpls(X.cord,Y)
         mgpls <- gpls::gpls(X.cord,Y,ncomp.gpls)
      }
      else{
         ncomp.gpls <- NULL
         mgpls <- NULL
      }
      
   }
   else{
      ncomp.gpls <- ux
      
      if(do.gpls){
         mgpls <- gpls::gpls(X.cord,Y.factor,ncomp.gpls)
      }
      else{
         mgpls <- NULL
      }
      ux1 <- ux
   }
   
   menv <- env_glm(X.cord,Y,ux,init=init)
   mfull <- env_glm(X.cord,Y,dx,init=init)
   
   cv.reg <- cv.glmnet(x=X.cord,y=Y.factor,alpha=0,family="binomial",nfolds = 5)
   mglmnet_reg <- glmnet(x=X.cord,y=Y.factor,alpha=0,family="binomial",lambda = cv.reg$lambda.min)
   

   return(list(beta=menv$beta,betafull=mfull$beta,alpha=menv$alpha,alphafull=mfull$alpha,
               basis.value.X=basis.value.x.list,ux=ux,OX=OX,mglmnet_reg=mglmnet_reg,mgpls=mgpls,ncomp.gpls=ncomp.gpls))
   

}






