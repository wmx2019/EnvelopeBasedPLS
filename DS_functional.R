source("FEPLS.R")
library(doFuture)
library(doRNG)
library(fda.usc)
Sys.setenv(LANG="en")

############################################################
# Descriptions: Data simulation for functional response
############################################################





set.seed(342)

registerDoFuture()
# Set the number of cores
plan(multisession,workers=23)

# Number of replications
n_sim <- 100




# Set parameters
ux <- 4
t <- seq(0,1.0,1/15)
tx.list <- list(t1=t,t2=t,t3=t)
# tx.list <- list(t1=t)
ty <- t
x.knots.list <- list(knots1=seq(0,1,1/4),knots2=seq(0,1,1/4),knots3=seq(0,1,1/4))
x.order.list <- list(order1=4,order2=4,order3=4)

# x.knots.list <- list(knots1=seq(0,1,1/4))
# x.order.list <- list(order1=4)

y.knots <- seq(0,1,1/5)
y.order <- 4

psi <- function(x){
   r <- 1
   for (i in 1:6) {
      r <- c(r,sqrt(2)*sin(2*i*pi*x),sqrt(2)*cos(2*i*pi*x))
   }
   return(r)
}


# The coordinate of operator B given the Fourier(TRUE) basis
B1 <- matrix(0,nrow = 13,ncol=13)
# B1[4,2] <- -1.2
# B1[4,3] <- -0.4
# B1[2,2] <- -0.1
# B1[2,3] <- 0.2
# B1[3,2] <- 0.3
# B1[3,3] <- -0.2
B1[2,2] <- -1.2
B1[4,3] <- 2.4

# B1[4,2] <- -1.2

# B_list <- list(B1)
B2 <- matrix(0,nrow = 13,ncol=13)
B2[2,3] <- -0.04
B2[1,3] <- -0.05

B3 <- matrix(0,nrow = 13,ncol=13)
B3[3,3] <- 0.03
B3[4,3] <- -0.01

B_list <- list(B1,B2,B3)
# tau1 <- c(6.4,1,5.6,2,4,3.2,2.4,1.6,0.8,0.5,0.3)^2
# tau1 <- c(7,4,2,2,5.6,2.2,3,1.6,1,0.8,0.5,0.4,0.2)^2
# tau1 <- c(4,4,1,5.6,2,2.2,3,1.6,1,0.8,0.5,0.4,0.2)^2


tau1 <- c(8,4,1.6,7.3,5.5,runif(8,0.2,0.3))
# tau_list <- list(tau1)
# tau2 <- c(6.4,1.4,5.6,1.1,4,1.2,4.4,1.6,0.8,0.4,0.4,0.4,0.2)^2
# tau3 <- c(6.4,1.2,5.6,1.8,4,1.2,4.4,1.6,0.9,0.6,0.3,0.4,0.2)^2
tau2 <- c(8,2,1,7.3,5.5,runif(8,0.2,0.3))
tau3 <- c(8,2,0.5,3,5.5,runif(8,0.2,0.3))

tau_list <- list(tau1,tau2,tau3)

# tau_list <- list(tau1)
# rho <- c(1,0.5,3,1,4,2,2.5,3,3.5,4,5)^2
# rho <- c(2,0.8,0.3,4,3.5,1,3.5,2.4,4,3.5,4.6,1,2)^2 

chi <- function(x){
   r <- 1
   for (i in 1:6) {
      r <- c(r,sqrt(2)*sin(2*i*pi*x),sqrt(2)*cos(2*i*pi*x))
   }
   return(r)
}

# rho <- c(3,1,1,8,5,4,2,2.5,3,3.5,4,5,4.6)^2
rho <- rep(1,13)
# rho <- sqrt(c(2,3,3,3,2,2,2,3,2,2,3,2,1))

# rho <- c(3,3,1,2,2,runif(8,0.1,0.2))^2
DG <- function(B_list,n,t){
   # Generating data for the scalar response
   # Args:
   #   n: sample size
   #   t: time series vector
   
   nt <- length(t)
   # tau_list <- list(tau1,tau2,tau3)
   # rho <- c(1,0.1,0.1,0.1,1,1,1,1,1,1,4)^2
   
   
   
   
   
   
   X_list <- list()
   xi_list <- list()
   p <- length(B_list)
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
   e <- matrix(0,n,nt)
   nu <- matrix(rnorm(13*n),nrow = n)
   for (i in 1:n) {
      for (j in 1:nt) {
         e[i,j] <- 0.2*sum(rho^{1/2}*nu[i,]*chi(t[j]))
      }
   }
   
   # Generate Y
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
   # Y <- Y+e
   return(list(X_list=X_list,Y=Y))
}


# x.knots.list <- list(knots1=seq(0,1,1/3),knots2=seq(0,1,1/3),knots3=seq(0,1,1/3))
# x.order.list <- list(order1=2,order2=2,order3=2)

mse <- NULL

for (i in 1:2) {
   # Set the sample sizes for sample set and test (validation) set
   n_samples <- 50*2^(i-1)
   n_validation <- 1000
   
   # Reproducibility #
   registerDoRNG(12)
   acc <- foreach(i=1:n_sim) %dopar% {
      tmp <- DG(B_list,n=n_samples,t)
      X_list <- tmp$X_list
      Y <- tmp$Y
      
      # Test Data
      tmp <- DG(B_list,n=n_validation,t)
      X_list_new <- tmp$X_list
      Y_new <- tmp$Y
      
      
      res_dir <- mfpedir(X_list,Y,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                         ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)
      # X_list = list(X1)
      # res_dir <- mfelm_dir(X_list,Y,ux,uy,t,t)
      
      
      res_KL <- mfpeKL(X_list,Y,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                       ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)
      
      
      # tmp <- get_coord_dir_sy_sp(X_list=X_list,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list)
      # X.cord <- tmp$X.cord
      # tmp <- get_coord_dir_sy_sp(list(Y),list(ty),list(y.knots),list(y.order))
      # Y.cord <- tmp$X.cord
      
      # X_fdata <- fdata(X_list[[1]],t)
      # Y_fdata <- fdata(Y,t)
      # X_fdata_new <- fdata(X_list_new[[1]],t)
      #
      # basisx <- create.bspline.basis(norder = 4,breaks=knots)
      # basisy <- create.bspline.basis(norder = 4,breaks=knots)
      # freg <- fregre.basis.fr(X_fdata,Y_fdata,basisx,basisy)
      # plot(Y_fdata,col=1)
      # lines(freg$fitted.values,col=2)
      
      # mean(colMeans((predict(freg,X_fdata_new)$data-Y_new)^2))
      
      
      # Coordinates of new data
      X_new_dir <- bvfitl(X_list_new,res_dir$basis.value.X)
      X_new_KL <- bvfit(cbindl(X_list_new),res_KL$phihat.X)
      
      
      
      # MSE
      # Direct estimation
      c(MSEF(res_dir$beta,res_dir$alpha,res_dir$basis.value.Y,X_new_dir,Y_new),MSEF(res_dir$betafull,res_dir$alphafull,res_dir$basis.value.Y,X_new_dir,Y_new),
        MSEP(res_dir$mpcr,NULL,res_dir$basis.value.Y,X_new_dir,Y_new),MSEP(res_dir$mpls,NULL,res_dir$basis.value.Y,X_new_dir,Y_new),res_dir$ux,
        MSEF(res_KL$beta,res_KL$alpha,res_KL$psihat.Y
             ,X_new_KL,Y_new),MSEF(res_KL$betafull,res_KL$alphafull,res_KL$psihat.Y,X_new_KL,Y_new),MSEP(res_KL$mpcr,NULL,res_KL$psihat.Y,X_new_KL,Y_new),
        MSEP(res_KL$mpls,NULL,res_KL$psihat.Y,X_new_KL,Y_new),res_KL$ux)
      
   }
   
   
   acc_M_matrix <- matrix(unlist(acc),byrow = TRUE,ncol = 10)
   colnames(acc_M_matrix) <- c("env_dir","full_dir","pcr_dir","pls_dir","ux_dir","env_KL","full_KL","pcr_KL","pls_KL","uy_KL")
   
   if(is.null(mse))
      mse <- signif(colMeans(acc_M_matrix),digits = 3)
   else
      mse <- cbind(mse,signif(colMeans(acc_M_matrix),digits = 3))
   
}
print(mse)



#### PLOT

tmp <- DG(B_list,n=n_samples,t)
X_list <- tmp$X_list
Y <- tmp$Y


res_dir <- mfpedir(X_list,Y,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                   ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)
# X_list = list(X1)
# res_dir <- mfelm_dir(X_list,Y,ux,uy,t,t)


res_KL <- mfpeKL(X_list,Y,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                 ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)

# Test Data
tmp <- DG(B_list,n=3,t)
X_list_new <- tmp$X_list
Y_new <- tmp$Y
X_new_dir <- bvfitl(X_list_new,res_dir$basis.value.X)

# basis.value.y <- generate_basis_value(y.knots,4,seq(0,1,1/50),res_dir$Oy)

Y_pred <- pred(res_dir$beta,res_dir$alpha,X_new_dir)%*%t(res_dir$basis.value.Y)
Y_fd_pred <- Data2fd(t,y=as.matrix(t(Y_pred)),basisobj = create.bspline.basis(breaks = y.knots,norder=4))
Y_fdata_pred <- fdata(Y_fd_pred)
Y_fdata_new <- fdata(Y_new,argvals = t)


plot(Y_fd_pred,col=2)
lines(Y_fdata_new,col=1)
####
