source("FEPLS.R")
library(doFuture)
library(doRNG)
library(fda.usc)
Sys.setenv(LANG="en")

############################################################
# Descriptions: Data simulation for scalar response
############################################################



# Set the sample sizes for sample set and test (validation) set
n_samples <- 25
n_validation <- 1000

set.seed(8)

registerDoFuture()
# Set the number of cores
plan(multisession,workers=23)


# Number of replications
n_sim <- 100


# Reproducibility #
registerDoRNG(144)
# registerDoRNG(12)

psi <- function(x){
   r <- 1
   for (i in 1:6) {
      r <- c(r,sqrt(2)*sin(2*i*pi*x),sqrt(2)*cos(2*i*pi*x))
   }
   return(r)
}

# Set the parameters
ux <- 2
t <- seq(0,1.0,1/12)

knots <- seq(0,1,1/4)
order <- 4

tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

# The coordinate of operator B given the Fourier(TRUE) basis
B1 <- matrix(0,nrow = 10,ncol = 13)
B1[1,2] <- -1.2
B1[1,3] <- -0.4
B1[2,2] <- 0.4
B1[2,3] <- 0.6
B1[3,2] <- 0.3
B1[3,3] <- -0.4
B1[4,2] <- -1.2
B1[4,3] <- 2.4

B_list <- list(B1=B1)
# tau1 <- c(2.4,1,5.6,2,4,3.2,2.4,1.6,.8,.5,.3)^2
# tau1 <- c(1,8,4,1,2,1,2,1,2,1)
tau1 <- c(8,4,1.6,7.3,5.5,runif(8,0.2,0.3))
# tau1 <- c(0.2,3,0.15,0.9,0.5,2,0.1,0.15,0.2,0.13,0.14,0.15,0.1)
# tau1 <- c(0.2,3,0.15,0.9,runif(7,0.1,0.2),runif(2,0.01,0.03))
tau_list <- list(tau1)


# X_fdata <- fdata(mdata = X_list[[1]],argvals = t)
# plot(X_fdata)
# tau1 <- c(5,3,1,1,runif(7,0.2,1))

DG_s <- function(B_list,n,t){
   # Generating data for the scalar response
   # Args:
   #   n: sample size
   #   t: time series vector
   # Returns:
   #   The list which contains the generated X (list) and Y 
   nt <- length(t)
   X_list <- list()
   xi_list <- list()
   p <- length(B_list)
   # TODO 
   # 
   r <- dim(B_list[[1]])[1]
   
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
   
   Y <- matrix(0,n,r)
   e <- matrix(rnorm(n*r,sd=1),n,r)
   # Generate Y
   for (l in 1:p) {
      xi <- xi_list[[l]]
      B <- B_list[[l]]
      for (i in 1:n) {
         Y[i,] <- Y[i,] + sum((tau_list[[l]]^{1/2}*xi[i,])%*%t(B))
      }
   }
   Y <- Y+1.1*e
   return(list(X_list=X_list,Y=Y))
}



acc <- foreach(i=1:n_sim) %dopar% {
   
   
   # Generating Data
   tmp <- DG_s(B_list,n=n_samples,t)
   X_list <- tmp$X_list
   X.list <- X_list
   Y <- tmp$Y
   # Test Data
   tmp <- DG_s(B_list,n=n_validation,t)
   X_list_new <- tmp$X_list
   Y_new <- tmp$Y
   
   res_dirs <- sfpedir(X_list,Y,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
   X_new_dir <- bvfitl(X_list_new,res_dirs$basis.value.X)
   X_dir <- bvfitl(X_list,res_dirs$basis.value.X)
   # bspl1 <- create.bspline.basis(breaks = knots,norder = order)
   # plot(bspl1)
   # 
   # X_fdata <- fdata(X_list[[1]],argvals = t)
   # X_fdata_new <- fdata(X_list_new[[1]],argvals = t)
   # ldat <- ldata("df"=as.data.frame(Y),"x"=x)
   # fda_lm <- fregre.lm(y~x,ldat,basis.b = bspl1)
   # fregre.lm(f,ldat,  basis.b=basis.b)
   # res <- fregre.pls(X_fdata,as.data.frame(Y),c(1:4))
   # fda_lm <- fregre.lm(y~x,ldata("df"=as.data.frame(Y),"x"=X_fdata),basis.b = bspl1)
   # 
   # fregre.pls.cv(X_fdata,Y,8)
   
   c(MSEFs(res_dirs$beta,res_dirs$alpha,X_new_dir,Y_new),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_new_dir,Y_new),
     MSEPs(res_dirs$mpcr,res_dirs$ncomp.pcr,X_new_dir,Y_new),MSEPs(res_dirs$mpls,res_dirs$ncomp.pls,X_new_dir,Y_new),
     res_dirs$ux,res_dirs$ncomp.pcr,res_dirs$ncomp.pls)
}




acc_M_matrix <- matrix(unlist(acc),byrow = TRUE,ncol = 7)


colnames(acc_M_matrix) <- c("env","full","pcr","pls","u env","u pcr","u pls")
signif(colMeans(acc_M_matrix),digits = 4)
















