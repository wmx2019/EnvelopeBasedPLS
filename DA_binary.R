# library(conflicted)
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
source("FEPLS.R")
Sys.setenv(LANG="en")
library(glmnet)
library(gpls)
library(dplyr)
# library(tidyverse)

# data("COVID19")

data(DTI)
Y <- DTI$case
X <- DTI$rcst
X.fd <- fdata(X)
plot(X.fd)

idx_row <- (1:dim(X)[1])[rowSums(is.na(X))==0]
X <- X[idx_row,]
Y <- Y[idx_row]
X_list <- list(X=X)
t <- seq(0,1,1/(dim(X)[2]-1))
knots <- seq(0,1,1/18)
order <- 4

tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)


# res <- cfelmdir(X_list,Y,NULL,tx.list,x.knots.list,x.order.list,u.selected = T)

pred_mis_cv <- function(X_list,Y,nfold=10){
   n <- length(Y)
   n_val <- as.integer(n/nfold)
   
   
   set.seed(4534)
   registerDoFuture()
   # Set the number of cores
   plan(multisession,workers=23)
   
   # Number of replications
   n_sim <- 100
   # Reproducibility #
   registerDoRNG(31235)
   
   pmr <- foreach(i=1:n_sim) %dopar% {
      idx_val <- sample(1:n,n_val)
      idx_train <- setdiff(1:n,idx_val)
      # Split the data into train and validation 
      Y_train <- Y[idx_train]
      Y_val <- as.matrix(Y[idx_val])
      X_list_train <- list()
      X_list_val <- list()
      for (X in X_list) {
         X_train <- X[idx_train,]
         X_val <- X[idx_val,]
         X_list_train <- append(X_list_train,list(X_train))
         X_list_val <- append(X_list_val,list(X_val))
      }
      
      
      res_dirc <- cfelmdir(X_list_train,Y_train,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
      
      # Coordinates of validation data
      X_val_dir <- bvfitl(X_list_val,res_dirc$basis.value.X)
      X.fd <- fdata(X_list_train[[1]],argvals = t)
      basis.x <- create.bspline.basis(breaks = knots,norder = order)
      a1 <- classif.glm(y ~ x,ldata("df"=data.frame(y=Y_train),"x"=X.fd)
                        ,family = binomial(link = "logit"),basis.x = basis.x,basis.b = basis.x)
      X_fd_val <- fdata(X_list_val[[1]],argvals = t)
      p1 <- predict(a1,ldata("df"=data.frame(y=Y_val),"x"=X_fd_val))
      
      
      
      c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
            pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
            sum(p1==Y_val)/n_val,
            pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_val_dir,Y_val)))
   }
   pmr_M_matrix <- matrix(unlist(pmr),byrow = TRUE,ncol = 5)
   colnames(pmr_M_matrix) <- c("env","fglm","fglmnet","classif.glm","fgpls")
   
   return(list(pmr=pmr_M_matrix))
}



nfold <- 15
n <- length(Y)
n_val <- as.integer(n/nfold)


set.seed(4534)
registerDoFuture()
# Set the number of cores
plan(multisession,workers=23)

# Number of replications
n_sim <- 100
# Reproducibility #
registerDoRNG(31235)

pmr <- foreach(i=1:n_sim) %dopar% {
   idx_val <- sample(1:n,n_val)
   idx_train <- setdiff(1:n,idx_val)
   # Split the data into train and validation 
   Y_train <- Y[idx_train]
   Y_val <- as.matrix(Y[idx_val])
   X_list_train <- list()
   X_list_val <- list()
   for (X in X_list) {
      X_train <- X[idx_train,]
      X_val <- X[idx_val,]
      X_list_train <- append(X_list_train,list(X_train))
      X_list_val <- append(X_list_val,list(X_val))
   }
   
   
   res_dirc <- cfelmdir(X_list_train,Y_train,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
   
   # Coordinates of validation data
   X_val_dir <- bvfitl(X_list_val,res_dirc$basis.value.X)
   X.fd <- fdata(X_list_train[[1]],argvals = t)
   basis.x <- create.bspline.basis(breaks = knots,norder = order)
   a1 <- classif.glm(y ~ x,ldata("df"=data.frame(y=Y_train),"x"=X.fd)
                     ,family = binomial(link = "logit"),basis.x = basis.x,basis.b = basis.x)
   X_fd_val <- fdata(X_list_val[[1]],argvals = t)
   p1 <- predict(a1,ldata("df"=data.frame(y=Y_val),"x"=X_fd_val))
   
   
   
   c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
         pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
         sum(p1==Y_val)/n_val,
         pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_val_dir,Y_val)))
}
pmr_M_matrix <- matrix(unlist(pmr),byrow = TRUE,ncol = 5)
colnames(pmr_M_matrix) <- c("env","fglm","fglmnet","classif.glm","fgpls")

colMeans(pmr_M_matrix)
