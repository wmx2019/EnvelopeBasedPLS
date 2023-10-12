############ 
# Binary
Sys.setenv(LANG="en")
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
library(glmnet)
source("FEPLS.R")
data("gasoline")
source("Data_Preprocess.R")

X <- gasoline$NIR
Y <- as.matrix(gasoline$octane)



# Centering data
# X <- data_centered(X)

Y <- Y>88


t <- seq(0,1,1/400)
X_list <- list(X=X)

X.fd <- fdata(t(X))
plot(X.fd, col=factor(Y))

# Set the parameters for the spline basis
knots <- seq(0,1,1/12)
# knots <- seq(0,1,1/12)
# knots <- seq(0,1,1/16)
order <- 4
tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

# # Coordinates w.r.t the given basis functions
# res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
# X.cord <- res_cord$X.cord
# summary(glm(Y~X.cord,family =binomial()))


nfold <- 5
n <- length(Y)
n_val <- as.integer(n/nfold)


set.seed(4534)
registerDoFuture()
# Set the number of cores
plan(multisession,workers=23)

# Number of replications
n_sim <- 100
# Reproducibility #
# registerDoRNG(31235)
registerDoRNG(31236)

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
   
   
   res_dirc <- cfelmdir(X_list_train,Y_train,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T,do.gpls = T)
   
   # Coordinates of validation data
   X_val_dir <- bvfitl(X_list_val,res_dirc$basis.value.X)
   X.fd <- fdata(X_list_train[[1]],argvals = t)
   basis.x <- create.bspline.basis(breaks = knots,norder = order)
   a1 <- classif.glm(y ~ x,ldata("df"=data.frame(y=Y_train),"x"=X.fd)
                     ,family = binomial(link = "logit"),basis.x = basis.x,basis.b = basis.x)
   X_fd_val <- fdata(X_list_val[[1]],argvals = t)
   p1 <- predict(a1,ldata("df"=data.frame(y=Y_val),"x"=X_fd_val))
   
   # ncomp.gpls <- u_gpls(X_train,Y_train)
   # mgpls <- gpls::gpls(X_train,Y_train,ncomp.gpls)
   
   # c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
   #       pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
   #       sum(p1==Y_val)/n_val,
   #       pred_acc(coef(mgpls)[-1],coef(mgpls)[1],X_val,Y_val)),res_dirc$ux,ncomp.gpls)
   
   c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
         pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
         sum(p1==Y_val)/n_val,
         pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_val_dir,Y_val)),res_dirc$ux,res_dirc$ncomp.gpls)
}
pmr_M_matrix <- matrix(unlist(pmr),byrow = TRUE,ncol = 7)
colnames(pmr_M_matrix) <- c("env","fglm","fglmnet","classif.glm","gpls","u_env","u_gpls")

colMeans(pmr_M_matrix)
