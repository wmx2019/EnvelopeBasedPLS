library(doFuture)
library(doRNG)
library(fda.usc)
Sys.setenv(LANG="en")
data(CanadianWeather)
source("FEPLS.R")

X <- t(CanadianWeather$dailyAv[1:360,,2])
Y <- t(CanadianWeather$dailyAv[1:360,,1])



X <- vapply(1:24, function(i) rowMeans(X[,(15*(i-1)+1):(15*i)]), numeric(35))
Y <- vapply(1:24, function(i) rowMeans(Y[,(15*(i-1)+1):(15*i)]), numeric(35))

X <- log(X)

X_list <- list(X=X)





# X <- t(CanadianWeather$monthlyPrecip)
# X_list <- list(X=X)
# Y <- t(CanadianWeather$monthlyTemp)


x.knots.list <- list(knots1=seq(0,1,1/8))
x.order.list <- list(order1=4)
y.knots <- seq(0,1,1/4)
y.order <- 4
t <- seq(0,1.0,1/(dim(X)[2]-1))
tx.list <- list(t1=t)
ty <- seq(0,1.0,1/(dim(Y)[2]-1))


res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
X.cord <- res_cord$X.cord
res_cord <- get_coord_dir_sy_sp(list(Y=Y),list(ty=ty),list(y.knots=y.knots),list(y.order=y.order))
Y.cord <- res_cord$X.cord
summary(lm(Y.cord~X.cord))


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
      Y_val <- Y[idx_val,]
      X_list_train <- list()
      X_list_val <- list()
      for (X in X_list) {
         X_train <- X[idx_train,]
         X_val <- X[idx_val,]
         X_list_train <- append(X_list_train,list(X_train))
         X_list_val <- append(X_list_val,list(X_val))
      }
      
      
      res_dir <- mfpedir(X_list_train,Y_train,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                         ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)
      res_KL <- mfpeKL(X_list_train,Y_train,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                       ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)
      
      # Coordinates of validation data
      X_val_dir <- bvfitl(X_list_val,res_dir$basis.value.X)
      X_val_KL <- bvfit(cbindl(X_list_val),res_KL$phihat.X)
      
      # ncomp.pls <- CV_pls(X_train,Y_train)
      # mpls <- pls::plsr(Y_train~X_train,method = "simpls",ncomp=ncomp.pls,validation="CV")
      # 
      # ncomp.pcr <- CV_pls(X_train,Y_train)
      # mpcr <- pls::pcr(Y_train~X_train,validation="CV",ncomp=ncomp.pcr)
      # 
      # c(MSEF(res_dir$beta,res_dir$alpha,res_dir$basis.value.Y,X_val_dir,Y_val),MSEF(res_dir$betafull,res_dir$alphafull,res_dir$basis.value.Y,X_val_dir,Y_val)
      #   ,res_dir$ux,MSEF(res_KL$beta,res_KL$alpha,res_KL$psihat.Y
      #        ,X_val_KL,Y_val),MSEF(res_KL$betafull,res_KL$alphafull,res_KL$psihat.Y,X_val_KL,Y_val),MSEPs(mpcr,NULL,X_val,Y_val),
      #   MSEPs(mpls,NULL,X_val,Y_val),res_KL$ux,mpcr$ncomp,mpls$ncomp)
      
      c(MSEF(res_dir$beta,res_dir$alpha,res_dir$basis.value.Y,X_val_dir,Y_val),MSEF(res_dir$betafull,res_dir$alphafull,res_dir$basis.value.Y,X_val_dir,Y_val),
        MSEP(res_dir$mpcr,NULL,res_dir$basis.value.Y,X_val_dir,Y_val),MSEP(res_dir$mpls,NULL,res_dir$basis.value.Y,X_val_dir,Y_val),res_dir$ux,res_dir$mpcr$ncomp,res_dir$mpls$ncomp,
        MSEF(res_KL$beta,res_KL$alpha,res_KL$psihat.Y
             ,X_val_KL,Y_val),MSEF(res_KL$betafull,res_KL$alphafull,res_KL$psihat.Y,X_val_KL,Y_val),MSEP(res_KL$mpcr,NULL,res_KL$psihat.Y,X_val_KL,Y_val),
        MSEP(res_KL$mpls,NULL,res_KL$psihat.Y,X_val_KL,Y_val),res_KL$ux,res_KL$mpcr$ncomp,res_KL$mpls$ncomp)
   }
   
   pse_M_matrix <- matrix(unlist(pse),byrow = TRUE,ncol = 14)
   colnames(pse_M_matrix) <- c("env_dir","full_dir","pcr","pls","u_env_dir","u_pcr","u_pls","env_KL","full_KL","pcr","pls","u_env_KL","u_pcr","u_pls")
   
   return(list(pse=pse_M_matrix))
}


pse <- pred_error_cv(X_list,Y,nfold = 5)
colMeans(pse$pse)

nfold = 5
set.seed(36)
n <- dim(Y)[1]
n_val <- as.integer(n/nfold)
idx_val <- sample(1:n,n_val)
idx_train <- setdiff(1:n,idx_val)
# Split the data into train and validation 

idx_val <- idx_val[c(3,4,5)]
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

res_dir <- mfpedir(X_list_train,Y_train,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                   ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)
res_KL <- mfpeKL(X_list_train,Y_train,ux,tx.list=tx.list,ty=ty,x.knots.list=x.knots.list
                 ,y.knots=y.knots,x.order.list=x.order.list,y.order = y.order,u.selected=T)


X_val_dir <- bvfitl(X_list_val,res_dir$basis.value.X)
Y_pred <- pred(res_dir$beta,res_dir$alpha,X_val_dir)%*%t(res_dir$basis.value.Y)
Y_fd_pred <- Data2fd(t,y=as.matrix(t(Y_pred)),basisobj = create.bspline.basis(breaks = y.knots,norder=4))
Y_fdata_pred <- fdata(Y_fd_pred)
Y_fdata_val <- fdata(Y_val,argvals = t)



plot(Y_fdata_val,col=1)
lines(Y_fd_pred, col=2)

