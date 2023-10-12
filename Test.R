library(doFuture)
library(doRNG)
library(fda.usc)
Sys.setenv(LANG="en")
data(CanadianWeather)
source("FEPLS.R")
cimis <- read.csv("cimis.csv")
length(unique(cimis$Stn.Name))
dat1 <- data.frame(id=cimis$Stn.Id,temp=cimis$Avg.Air.Temp..F.,dew=cimis$Avg.Dew.Point..F.
                   ,wind=cimis$Avg.Wind.Speed..mph.,precip = cimis$Total.Precip..in.,time=cimis$Month.Year,hum=cimis$Avg.Rel.Hum....)

sl <- split(dat1,f=dat1$id)

idxRomve <- NULL
for (i in 1:length(sl)) {
   if(length(sl[[i]]$time)<12){
      idxRomve <- c(idxRomve,i)
   }
}

sl <- sl[-idxRomve]

ns <- length(sl)
X1 <- NULL
for (i in 1:ns) {
   X1 <- rbind(X1,sl[[i]]$temp)
}
row.names(X1) <- names(sl)

X2 <- NULL
for (i in 1:ns) {
   X2 <- rbind(X2,sl[[i]]$dew)
}

row.names(X2) <- names(sl)

X3 <- NULL
for (i in 1:ns) {
   X3 <- rbind(X3,sl[[i]]$wind)
}

row.names(X3) <- names(sl)

X4 <- NULL
for (i in 1:ns) {
   X4 <- rbind(X4,sl[[i]]$hum)
}

row.names(X4) <- names(sl)

Y <- NULL
for (i in 1:ns) {
   Y <- rbind(Y,sl[[i]]$precip)
}
row.names(Y) <- names(sl)

dl <- remove_na(list(X1=X1,X2=X2,X3=X3,X4=X4,Y=Y))

X1 <- dl$X1
X2 <- dl$X2
X3 <- dl$X3
X4 <- dl$X4
Y <- dl$Y


idxRomve <- (1:nrow(Y))[rowSums(Y<1e-8)>2]
X1 <- X1[-idxRomve,]
X2 <- X2[-idxRomve,]
X3 <- X3[-idxRomve,]
X4 <- X4[-idxRomve,]
Y <- Y[-idxRomve,]
Y <- log(Y+1e-2)


t <- seq(0,1.0,1/(dim(X1)[2]-1))

X_list <- list(X1=X1,X2=X2,X3=X3,X4=X4)
x.knots.list <- list(knots1=seq(0,1,1/6),knots2=seq(0,1,1/6),knots3=seq(0,1,1/6),knots4=seq(0,1,1/6))
x.order.list <- list(order1=4,order2=4,order3=4,order4=4)
tx.list <- list(t1=t,t2=t,t3=t,t4=t)

X_list <- list(X2=X2,X3=X3,X4=X4)
x.knots.list <- list(knots2=seq(0,1,1/5),knots3=seq(0,1,1/5),knots4=seq(0,1,1/5))
x.order.list <- list(order2=4,order3=4,order4=4)
tx.list <- list(t2=t,t3=t,t4=t)

y.knots <- seq(0,1,1/5)
y.order <- 4
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
set.seed(13123)
n <- dim(Y)[c(1,2)]
n_val <- as.integer(n/nfold)
idx_val <- sample(1:n,n_val)
idx_train <- setdiff(1:n,idx_val)
# Split the data into train and validation 

idx_val <- idx_val[1]
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

