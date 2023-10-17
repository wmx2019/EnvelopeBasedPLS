############ 
# dataset that is unable to preprocessed
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
source("FEPLS.R")
Sys.setenv(LANG="en")
library(dplyr)
library(tidyverse)
data("content")

infant_data <- content

t.days <- as.integer(names(table(infant_data$agedays))[table(infant_data$agedays)>50])
infant_data <- infant_data[infant_data$agedays %in% t.days,1:6] 






t.days <- c(28,56,84,98,112,126,140,154,168,196,224)
infant_data <- infant_data[infant_data$agedays %in% t.days,] 

t.days <- c(28,56,84,112,140,168,196,224)
infant_data <- infant_data[infant_data$agedays %in% t.days,] 

data.split <- split(infant_data[,1:6],infant_data$id)

df.len <- NULL
for (df in data.split) {
   df.len <- c(df.len,dim(df)[1])
}
table(df.len)


id_keep <- NULL
for (df in data.split) {
   if(dim(df)[1]==length(t.days)){
      id_keep <- c(id_keep,df$id)
   }
}
infant_data <- infant_data[infant_data$id %in% id_keep,] 

df <- bind_cols(data.split,.id="agedays")

df <- data.split %>% reduce(merge, by="agedays")

t.days <- idx.set[table(content$agedays)>100]
t.days <- c(46)
df <- content
select_common <- function(df,t.days){
   idx.set <- unique(df$id)
   id.keep <- NULL
   df.list <- split(df[,1:6],df$id)
   
   for (i in 1:length(idx.set)) {
      df.i <- df.list[[i]]
      if(sum(df.i$agedays %in% t.days)==length(t.days)){
         id.keep <- c(id.keep,i)
      }
   }
   return(id.keep)
}
select_common(df,t.days)






############ 
#DA_scalar
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
source("FEPLS.R")
source("Data_Preprocess.R")

data("tecator")
t_rescale <- function(t0){
   return((t0-min(t0))/(max(t0)-min(t0)))
}


X <- tecator$absorp.fdata$data
Y <- as.matrix(tecator$y$Water)


# Centering data
X <- data_centered(X)
Y <- data_centered(Y)


X.fd <- fdata(X)
plot(X.fd)
t <- t_rescale(tecator$absorp.fdata$argvals)



X_list <- list(X=X)




knots <- seq(0,1,1/7)
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
   n_sim <- 200
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
      
      ncomp.pls <- CV_pls(X_train,Y_train)
      mpls <- pls::plsr(Y_train~X_train,method = "simpls",ncomp=ncomp.pls,validation="CV")
      
      ncomp.pcr <- CV_pls(X_train,Y_train)
      mpcr <- pls::pcr(Y_train~X_train,validation="CV",ncomp=ncomp.pcr)
      
      
      c(MSEFs(res_dirs$beta,res_dirs$alpha,X_val_dir,Y_val),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_val_dir,Y_val),
        MSEPs(mpcr,ncomp.pcr,X_val,Y_val),MSEPs(mpls,ncomp.pls,X_val,Y_val),
        res_dirs$ux,ncomp.pcr,ncomp.pls)
      
      # c(MSEFs(res_dirs$beta,res_dirs$alpha,X_val_dir,Y_val),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_val_dir,Y_val),
      #   MSEPs(res_dirs$mpcr,res_dirs$ncomp.pcr,X_val_dir,Y_val),MSEPs(res_dirs$mpls,res_dirs$ncomp.pls,X_val_dir,Y_val),
      #   res_dirs$ux,res_dirs$ncomp.pcr,res_dirs$ncomp.pls)
   }
   pse_M_matrix <- matrix(unlist(pse),byrow = TRUE,ncol = 7)
   colnames(pse_M_matrix) <- c("env","full","pcr","pls","u env","u pcr","u pls")
   
   return(list(pse=pse_M_matrix))
}


pse <- pred_error_cv(X_list,Y,nfold = 5)
colMeans(pse$pse)


















############ 
# Binary
Sys.setenv(LANG="en")
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
library(glmnet)
library(fds)
source("FEPLS.R")
source("Data_Preprocess.R")

data("poblenou")
t_rescale <- function(t0){
   return((t0-min(t0))/(max(t0)-min(t0)))
}


X <- poblenou$nox$data

X.fd <- fdata(X)
plot(X.fd)


t <- t_rescale(poblenou$nox$argvals)
X_list <- list(X=X)

Y <- as.matrix(as.numeric(poblenou$df$day.festive))-1

# Remove NA
idx_row <- na_keep(X_list,Y)

for (i in 1:length(X_list)) {
   X <- X_list[[i]]
   X_list[[i]] <- X[idx_row,] 
}

Y <- Y[idx_row,]

knots <- seq(0,1,1/8)
order <- 4

tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

# Coordinates w.r.t the given basis functions
res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
X.cord <- res_cord$X.cord
summary(glm(Y~X.cord,family =binomial()))


res <- cfelmdir(X_list,Y,NULL,tx.list,x.knots.list,x.order.list,u.selected = T)

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



nfold <- 40
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










############ 
# Binary
Sys.setenv(LANG="en")
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
library(glmnet)
library(fds)
source("FEPLS.R")
source("Data_Preprocess.R")

data(Moisturespectrum)
data(Moisturevalues)
plot(Moisturespectrum)

idx <- sample(1:100,40)
X <- t(Moisturespectrum$y)[idx,]
Y <- as.matrix(Moisturevalues>14.5)[idx,]
Y <- as.matrix(Y)

X.fd <- fdata(X)
plot(X.fd)
t <- t_rescale(Moisturespectrum$x)
X_list <- list(X=X)

# Centering data
X <- data_centered(X)

knots <- seq(0,1,1/12)
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
   
   ncomp.gpls <- u_gpls(X_train,Y_train)
   mgpls <- gpls::gpls(X_train,Y_train,ncomp.gpls)
   
   c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
         pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
         sum(p1==Y_val)/n_val,
         pred_acc(coef(mgpls)[-1],coef(mgpls)[1],X_val,Y_val)),res_dirc$ux,ncomp.gpls)
   
   # c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
   #       pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
   #       sum(p1==Y_val)/n_val,
   #       pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_val_dir,Y_val)))
}
pmr_M_matrix <- matrix(unlist(pmr),byrow = TRUE,ncol = 7)
colnames(pmr_M_matrix) <- c("env","fglm","fglmnet","classif.glm","fgpls","u_env","u_gpls")

colMeans(pmr_M_matrix)


############ 
# scalar
Sys.setenv(LANG="en")
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
library(glmnet)
library(fds)
source("FEPLS.R")
source("Data_Preprocess.R")

data(Moisturespectrum)
data(Moisturevalues)
plot(Moisturespectrum)

set.seed(5646156)

idx <- sample(1:100,30)
X <- t(Moisturespectrum$y)[idx,]
Y <- Moisturevalues[idx]
Y <- as.matrix(Y)

X.fd <- fdata(X)
plot(X.fd)

t <- t_rescale(Moisturespectrum$x)


# Centering data
X <- data_centered(X)
Y <- data_centered(Y)


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

# res_cord <- get_coord_dir_sy_sp(X_list_train,tx.list,x.knots.list,x.order.list)
# X.cord <- res_cord$X.cord
# Y.cord <- as.matrix(Y[idx_train,])


res_dirs <- sfpedir(X_list,Y,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
pred_error_cv <- function(X_list,Y,nfold=4){
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
   
   # idx <- readRDS(file = "case1.RData")
   # idx_val <- idx$idx_val
   # idx_train <- idx$idx_train
   
   
   
   pse <- foreach(i=1:n_sim) %dopar% {
      
      idx_val <- sample(1:n,n_val)
      idx_train <- setdiff(1:n,idx_val)
      
      # if(i==2){
      #    saveRDS(list(idx_val=idx_val,idx_train=idx_train),file = "case1.RData")
      # }
      
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
      
      
      res_dir <- sfpedir(X_list_train,Y_train,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
      
      # Coordinates of validation data
      X_val_dir <- bvfitl(X_list_val,res_dir$basis.value.X)
      
      
      ncomp.pls <- CV_pls(X_train,Y_train)
      mpls <- pls::plsr(Y_train~X_train,method = "simpls",ncomp=ncomp.pls,validation="CV")
      
      ncomp.pcr <- CV_pls(X_train,Y_train)
      mpcr <- pls::pcr(Y_train~X_train,validation="CV",ncomp=ncomp.pcr)
      c(MSEFs(res_dirs$beta,res_dirs$alpha,X_val_dir,Y_val),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_val_dir,Y_val),
        MSEPs(mpcr,ncomp.pcr,X_val,Y_val),MSEPs(mpls,ncomp.pls,X_val,Y_val),
        res_dirs$ux,ncomp.pcr,ncomp.pls)
      # c(MSEFs(res_dirs$beta,res_dirs$alpha,X_val_dir,Y_val),MSEFs(res_dirs$betafull,res_dirs$alphafull,X_val_dir,Y_val),
      #   MSEPs(res_dirs$mpcr,res_dirs$ncomp.pcr,X_val_dir,Y_val),MSEPs(res_dirs$mpls,res_dirs$ncomp.pls,X_val_dir,Y_val),
      #   res_dirs$ux,res_dirs$ncomp.pcr,res_dirs$ncomp.pls)
   }
   pse_M_matrix <- matrix(unlist(pse),byrow = TRUE,ncol = 7)
   colnames(pse_M_matrix) <- c("env","full","pcr","pls","u env","u pcr","u pls")
   
   return(list(pse=pse_M_matrix))
}


pse <- pred_error_cv(X_list,Y,nfold = 10)
colMeans(pse$pse)





############ 
# Binary
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
source("Data_Preprocess.R")
# library(tidyverse)

# data("COVID19")

data(DTI)
Y <- DTI$case
X <- DTI$cca


idx_row <- (1:dim(X)[1])[rowSums(is.na(X))==0]
X <- X[idx_row,]
Y <- Y[idx_row]

X.fd <- fdata(t(X))
plot(X.fd)

idx <- sample(1:100,40)
X <- X[idx,]
Y <- Y[idx]
Y <- as.matrix(Y)

# Centering data
X <- data_centered(X)

X_list <- list(X=X)




t <- seq(0,1,1/(dim(X)[2]-1))
knots <- seq(0,1,1/20)
order <- 4

tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)

# Coordinates w.r.t the given basis functions
res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
X.cord <- res_cord$X.cord
summary(glm(Y~X.cord,family =binomial()))


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
registerDoRNG(31235)

pmr <- foreach(i=1:n_sim) %dopar% {
   idx_val <- sample(1:n,n_val)
   idx_train <- setdiff(1:n,idx_val)
   # Split the data into train and validation 
   Y_train <- as.matrix(Y[idx_train])
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
   
   ncomp.gpls <- u_gpls(X_train,Y_train)
   mgpls <- gpls::gpls(X_train,Y_train,ncomp.gpls)
   
   c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
         pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
         sum(p1==Y_val)/n_val,
         pred_acc(coef(mgpls)[-1],coef(mgpls)[1],X_val,Y_val)),res_dirc$ux,ncomp.gpls)
   
   # c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
   #       pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
   #       sum(p1==Y_val)/n_val,
   #       pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_val_dir,Y_val)),res_dirc$ux,res_dirc$ncomp.gpls)
}
pmr_M_matrix <- matrix(unlist(pmr),byrow = TRUE,ncol = 7)
colnames(pmr_M_matrix) <- c("env","fglm","fglmnet","classif.glm","fgpls","u_env","u_gpls")

colMeans(pmr_M_matrix)

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
tidx <- (1:700)[seq(5,700,5)]

X <- t(nirc$y)[,tidx]
Y <- as.matrix(labc[2,])


# Centering data
X <- data_centered(X)
Y <- data_centered(Y)


X.fd <- fdata(X)
plot(X.fd)
t <- t_rescale(nirc$x[tidx])



X_list <- list(X=X)




knots <- seq(0,1,1/11)
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
############
# scalar
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
source("FEPLS.R")
data("gasoline")
source("Data_Preprocess.R")

X <- gasoline$NIR
Y <- as.matrix(gasoline$octane)

X.fd <- fdata(t(X))
plot(X.fd)

# Centering data
X <- data_centered(X)
Y <- data_centered(Y)


t <- seq(0,1,1/400)
X_list <- list(X=X)


# Set the parameters for the spline basis
knots <- seq(0,1,1/12)
order <- 4
tx.list <- list(t1=t)
x.knots.list <- list(knots1=knots)
x.order.list <- list(order1=order)





# Coordinates w.r.t the given basis functions
res_cord <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)
X.cord <- res_cord$X.cord
summary(lm(Y~X.cord))


# res_dirs <- sfpedir(X_list,Y,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
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
      
      
      res_dirs <- sfpedir(X_list_train,Y_train,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
      
      # Coordinates of validation data
      X_val_dir <- bvfitl(X_list_val,res_dirs$basis.value.X)
      
      ncomp.pls <- CV_pls(X_train,Y_train)
      mpls <- pls::plsr(Y_train~X_train,method = "simpls",ncomp=ncomp.pls,validation="CV")
      
      ncomp.pcr <- CV_pls(X_train,Y_train)
      mpcr <- pls::pcr(Y_train~X_train,validation="CV",ncomp=ncomp.pcr)
      
      
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
   
   # ncomp.gpls <- u_gpls(X_train,Y_train)
   # mgpls <- gpls::gpls(X_train,Y_train,ncomp.gpls)
   # 
   # c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
   #       pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
   #       sum(p1==Y_val)/n_val,
   #       pred_acc(coef(mgpls)[-1],coef(mgpls)[1],X_val,Y_val)),res_dirc$ux,ncomp.gpls)
   # 
   c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_val_dir,Y_val),pred_acc(res_dirc$betafull,res_dirc$alphafull,X_val_dir,Y_val),
         pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_val_dir,Y_val),
         sum(p1==Y_val)/n_val,
         pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_val_dir,Y_val)))
}
pmr_M_matrix <- matrix(unlist(pmr),byrow = TRUE,ncol = 7)
colnames(pmr_M_matrix) <- c("env","fglm","fglmnet","classif.glm","gpls","u_env","u_gpls")

colMeans(pmr_M_matrix)

############ 
# Binary
Sys.setenv(LANG="en")
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
library(glmnet)
library(fds)
source("FEPLS.R")
source("Data_Preprocess.R")
data(aa)
data(ao)
data(dcl)
data(iy)
data(sh)

plot(iy)


set.seed(5646156)


# idx1 <- sample(1:400,100)
# idx2 <- sample(1:400,100)


X <- rbind(t(sh$y)[,],t(iy$y)[,])

Y <- as.matrix(as.logical(c(rep(1,400)[],rep(0,400)[])))

idx <- sample(1:800,45)
X <- X[idx,]
Y <- as.matrix(Y[idx,])

X.fd <- fdata(X)
plot(X.fd)

t <- t_rescale(aa$x)


# Centering data
# X <- data_centered(X)



X_list <- list(X=X)
knots <- seq(0,1,1/12)
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

############
# Functional
library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
source("FEPLS.R")
source("Data_Preprocess.R")

Sys.setenv(LANG="en")

dat1 <- read.csv("complete_dataset.csv")

date1 <- head(dat1$date,30)
dat1 <- dat1[dat1$date %in% date1,]

set.seed(50)

vol <- vapply(1:500, function(i) mean(dat1$volume[(30*(i-1)+1):(30*i)]), numeric(1))

name1 <- unique(dat1$Name)[sort(vol,decreasing = F,index.return=T)$ix[1:50]]



dat1 <- dat1[dat1$Name %in% name1,]




sl <- split(dat1,f = dat1$Name)

X1 <- NULL
for (i in 1:50) {
   X1 <- rbind(X1,sl[[i]]$volume)
}

# X2 <- NULL
# for (i in 1:50) {
#    X2 <- rbind(X2,sl[[i]]$open)
# }

Y <- NULL
for (i in 1:50) {
   Y <- rbind(Y,sl[[i]]$close)
}

# # Centering data
# X <- data_centered(X)
# Y <- data_centered(Y)



X_list <- list(X1=X1)

X1.fd <- fdata(X1)
plot(X1.fd)




x.knots.list <- list(knots1=seq(0,1,1/5))
x.order.list <- list(order1=4)

y.knots <- seq(0,1,1/8)
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
