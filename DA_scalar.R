############ 
#DA_scalar
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
library(fds)
source("FEPLS.R")
source("Data_Preprocess.R")

data(labp)
data(labc)
data(nirp)
data(nirc)

t_rescale <- function(t0){
   return((t0-min(t0))/(max(t0)-min(t0)))
}
tidx <- (1:700)

X <- t(nirc$y)[,tidx]
Y <- as.matrix(labc[2,])


# Centering data
X <- data_centered(X)
Y <- data_centered(Y)


X.fd <- fdata(X)
plot(X.fd)
t <- t_rescale(nirc$x[tidx])



X_list <- list(X=X)




knots <- seq(0,1,1/12)
order <- 4
tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

# Coordinates w.r.t the given basis functions
res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
X.cord <- res_cord$X.cord
summary(lm(Y~X.cord))



pred_error_cv <- function(X_list,Y,nfold=10){
   n <- dim(Y)[1]
   n_val <- as.integer(n/nfold)
   
   
   set.seed(4534)
   registerDoFuture()
   # Set the number of cores
   plan(multisession,workers=23)
   
   # Number of replications
   n_sim <- 100
   # Reproducibility #
   registerDoRNG(31235)
   
   pse <- foreach(i=1:n_sim) %dopar% {
      idx_val <- sample(1:n,n_val)
      idx_train <- setdiff(1:n,idx_val)
      # Split the data into train and validation 
      Y_train <- as.matrix(Y[idx_train,])
      Y_val <- as.matrix(Y[idx_val,])
      X_list_train <- list()
      X_list_val <- list()
      for (X in X_list) {
         X_train <- X[idx_train,]
         X_val <- X[idx_val,]
         X_list_train <- append(X_list_train,list(X_train))
         X_list_val <- append(X_list_val,list(X_val))
      }
      
      
      res_dirs <- sfpedir(X_list_train,Y_train,NULL,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
      
      # Coordinates of validation data
      X_val_dir <- bvfitl(X_list_val,res_dirs$basis.value.X)
      
      # ncomp.pls <- CV_pls(X_train,Y_train)
      # mpls <- pls::plsr(Y_train~X_train,method = "simpls",ncomp=ncomp.pls,validation="CV")
      # 
      # ncomp.pcr <- CV_pls(X_train,Y_train)
      # mpcr <- pls::pcr(Y_train~X_train,validation="CV",ncomp=ncomp.pcr)
      # 
      # 
      # c(MSEFs(res_dirs$beta,res_dirs$alpha,X_val_dir,Y_val),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_val_dir,Y_val),
      #   MSEPs(mpcr,ncomp.pcr,X_val,Y_val),MSEPs(mpls,ncomp.pls,X_val,Y_val),
      #   res_dirs$ux,ncomp.pcr,ncomp.pls)
      
      c(MSEFs(res_dirs$beta,res_dirs$alpha,X_val_dir,Y_val),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_val_dir,Y_val),
        MSEPs(res_dirs$mpcr,res_dirs$ncomp.pcr,X_val_dir,Y_val),MSEPs(res_dirs$mpls,res_dirs$ncomp.pls,X_val_dir,Y_val),
        res_dirs$ux,res_dirs$ncomp.pcr,res_dirs$ncomp.pls)
   }
   pse_M_matrix <- matrix(unlist(pse),byrow = TRUE,ncol = 7)
   colnames(pse_M_matrix) <- c("env","full","pcr","pls","u env","u pcr","u pls")
   
   return(list(pse=pse_M_matrix))
}


pse <- pred_error_cv(X_list,Y,nfold = 5)
colMeans(pse$pse)