library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
source("FEPLS.R")
data("gasoline")


X <- gasoline$NIR
X.fd <- fdata(X)
plot(X.fd)
t <- seq(0,1,1/400)
X_list <- list(X=X)
Y <- as.matrix(gasoline$octane)


knots <- seq(0,1,1/20)
order <- 4
tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

res_dirs <- sfpedir(X_list,Y,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
pred_error_cv <- function(X_list,Y,nfold=10){
   n <- dim(Y)[1]
   n_val <- as.integer(n/nfold)
   
   
   set.seed(4534)
   registerDoFuture()
   # Set the number of cores
   plan(multisession,workers=23)
   
   # Number of replications
   n_sim <- 200
   # Reproducibility #
   registerDoRNG(31235)
   
   pse <- foreach(i=1:n_sim) %dopar% {
      idx_val <- sample(1:n,n_val)
      idx_train <- setdiff(1:n,idx_val)
      # Split the data into train and validation 
      Y_train <- Y[idx_train,]
      Y_val <- as.matrix(Y[idx_val,])
      X_list_train <- list()
      X_list_val <- list()
      for (X in X_list) {
         X_train <- X[idx_train,]
         X_val <- X[idx_val,]
         X_list_train <- append(X_list_train,list(X_train))
         X_list_val <- append(X_list_val,list(X_val))
      }
      
      
      res_dir <- sfpedir(X_list_train,Y_train,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
      
      # Coordinates of validation data
      X_val_dir <- bvfitl(X_list_val,res_dir$basis.value.X)
      
      c(MSEFs(res_dirs$beta,res_dirs$alpha,X_val_dir,Y_val),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_val_dir,Y_val),
        MSEPs(res_dirs$mpcr,res_dirs$ncomp.pcr,X_val_dir,Y_val),MSEPs(res_dirs$mpls,res_dirs$ncomp.pls,X_val_dir,Y_val),
        res_dirs$ux,res_dirs$ncomp.pcr,res_dirs$ncomp.pls)
   }
   pse_M_matrix <- matrix(unlist(pse),byrow = TRUE,ncol = 7)
   colnames(pse_M_matrix) <- c("env","full","pcr","pls","u env","u pcr","u pls")
   
   return(list(pse=pse_M_matrix))
}


pse <- pred_error_cv(X_list,Y,nfold = 10)
colMeans(pse$pse)
