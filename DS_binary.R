source("FEPLS.R")
library(doFuture)
library(doRNG)
library(glmnet)
library(fda.usc)

############################################################
# Descriptions: Data simulation for binary response
############################################################

Sys.setenv(LANG="en")

# Set the sample sizes for sample set and test (validation) set
# n_samples <- 40
n_validation <- 1000

set.seed(5)

# Set the number of cores #
registerDoFuture()
plan(multisession, workers = 23)



# Number of replications
n_sim <- 100
# tau1 <- c(0.2,3,0.15,0.9,runif(9,0.1,0.2))
# tau1 <- c(3,4,0.9,4,4,runif(8,0.1,0.2))
tau1 <- c(8,5,0.9,7.3,2,runif(8,0.1,0.2))
# tau1 <- c(1,3,2,1,runif(9,0.1,0.2))

# tau1 <- c(8,5,0.9,7.3,2,runif(8,0.1,0.2))

# tau1 <-c(2,4,7.3,1.6,3.5,runif(8,0.2,0.3))/2

# The coordinate of operator B given selected basis function
B1 <- matrix(0,nrow = 1,ncol=13)
# B1[1,2] <- -1.1
# B1[1,3] <- 3.6
B1[1,2] <- -0.8
B1[1,3] <- 4.4
# B1[1,2] <- -1.2
# B1[1,4] <- 2.4
B_list <- list(B1=B1)


ux <- 2
t <- seq(0,1.0,1/14)

# knots <- seq(0,1,1/5)
# order <- 4



tx.list <- list(t1=t)


psi <- function(x){
   r <- 1
   for (i in 1:6) {
      r <- c(r,sqrt(2)*sin(2*i*pi*x),sqrt(2)*cos(2*i*pi*x))
   }
   return(r)
}   


DG_s <- function(B_list,n,t){
   tau_list <- list(tau1)
   
   nt <- length(t)
   X_list <- list()
   xi_list <- list()
   p <- length(B_list)
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
   e <- matrix(rnorm(n),n,1)
   # Generate Y
   
   
   eta <- matrix(0,n,1)
   for (l in 1:p) {
      xi <- xi_list[[l]]
      B <- B_list[[l]]
      for (i in 1:n) {
         eta[i] <- eta[i] + sum((tau_list[[l]]^{1/2}*xi[i,])%*%t(B))
      }
   }
   # print(eta)
   lp <- exp(eta)/(1+exp(eta))
   
   Y <- rbinom(n,1,lp)
   
   return(list(X_list=X_list,Y=Y,lp=lp))
}


mis_rate <- NULL

for (i in 1:4) {
   
   n_samples <- 40*2^(i-1)
   
   knots <- seq(0,1,1/4)
   order <- 4
   x.knots.list <- list(knots1=knots)
   x.order.list <- list(order1=order)
   
   # selected seed
   registerDoRNG(5)
   
   acc <- foreach(i=1:n_sim) %dopar% {
      
      # psi <- function(x){
      #    r <- vapply(1:5, function(i) c(sqrt(2)*sin(2*i*pi*x),sqrt(2)*cos(2*i*pi*x)), numeric(2))
      #    return(c(1,r))
      # }
      
      
      
      tmp <- DG_s(B_list,n=n_samples,t)
      X_list <- tmp$X_list
      Y <- tmp$Y
      if(is.null(Y)){
         stop("NULL")
      }
      
      
      # Test Data
      tmp <- DG_s(B_list,n=n_validation,t)
      X_list_new <- tmp$X_list
      Y_new <- tmp$Y
      
      # plot(1:n_validation,tmp$lp)
      
      
      
      
      
      tmp1 <- vapply(t, psi, FUN.VALUE = numeric(13))[c(2,3),]
      
      Gamma <- t(get_coord_dir_sy_sp(list(tmp=tmp1),tx.list,x.knots.list,x.order.list)$X.cord)
      
      # get_coord_dir_sy(X_list,Y,tx.list,x.knots.list,x.order.list)$X.cord
      
      eig_g <- eigen(t(Gamma)%*%Gamma)
      Gamma <- Gamma%*%eig_g$vectors%*%diag(1/sqrt(eig_g$values))%*%t(eig_g$vectors)
      
      Gamma_t <- Gamma
      X <- get_coord_dir_sy_sp(X_list,tx.list,x.knots.list,x.order.list)$X.cord
      
      Y <- as.matrix(Y)
      
      
      
      # env_glm(X,Y,2,type="logit",asy=TRUE,init=Gamma)
      glm1 <- glm(as.logical(Y)~X%*%Gamma_t,family = binomial(link = "logit"))
      alpha <- coef(glm1)[1]
      eta <- coef(glm1)[-1]
      
      
      res_dirc <- cfelmdir(X_list,Y,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=T)
      
      # res_dirc1 <- cfelmdir(X_list,Y,ux,tx.list=tx.list,x.knots.list=x.knots.list,x.order.list=x.order.list,u.selected=F,init=Gamma)
      # X.cord <- res_dirc[[1]]
      # Y.cord <- res_dirc[[2]]
      
      
      # Coordinates of new data
      X_new_dir <- bvfitl(X_list_new,res_dirc$basis.value.X)
      
      # X_new_KL <- bvfit(cbindl(X_list_new),res_KL$phihat.X)
      
      # X_dir <- bvfitl(X_list,res_dir$basis.value.X)
      # X_KL <- bvfit(cbindl(X_list),res_KL$phihat.X)
      # reg <- res_dirc$mglmnet_reg 
      # glm1 <- res_dirc$mglm
      
      # roc_LE <- roc(Y_new,pred_logit_p(res_dirc$beta,res_dirc$alpha,X_new_dir)[,1])
      # roc_glm <- roc(Y_new,pred_logit_p(res_dirc$betafull,res_dirc$alphafull,X_new_dir)[,1])
      # roc_glmnet <- roc(Y_new,pred_logit_p(coef(reg)[-1],coef(reg)[1],X_new_dir)[,1])
      
      
      X.fdata <- fdata(X_list[[1]],argvals = t)
      # plot(X.fdata)
      basis.x <- create.bspline.basis(breaks = knots,norder = order)
      a1 <- classif.glm(y ~ x,ldata("df"=data.frame(y=Y),"x"=X.fdata)
                        ,family = binomial(link = "logit"),basis.x = basis.x,basis.b = basis.x)
      X.fdata.new <- fdata(X_list_new[[1]],argvals = t)
      p1 <- predict(a1,ldata("df"=data.frame(y=Y_new),"x"=X.fdata.new))
      
      
      
      
      
      
      # if(abs(sum(X_new_dir^2)-sum(X_new_KL^2))>1e-10)
      #    print(abs(mean(X_new_dir^2)-mean(X_new_KL^2)))
      
      # 1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_new_dir,Y_new),pred_acc(res_dirc1$beta,res_dirc1$alpha,X_new_dir,Y_new),
      #   pred_acc(res_dirc$betafull,res_dirc$alphafull,X_new_dir,Y_new),1-MMAC_F(res_dirc$mgpls,X_new_dir,Y_new),
      #   pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_new_dir,Y_new),pred_acc(Gamma%*%eta,alpha,X_new_dir,Y_new),sum(p1==Y_new)/n_validation
      # )
      
      # c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_new_dir,Y_new)),res_dirc$ux)
      c(1-c(pred_acc(res_dirc$beta,res_dirc$alpha,X_new_dir,Y_new),
            pred_acc(res_dirc$betafull,res_dirc$alphafull,X_new_dir,Y_new),
            pred_acc(coef(res_dirc$mglmnet_reg)[-1],coef(res_dirc$mglmnet_reg)[1],X_new_dir,Y_new),
            pred_acc(Gamma_t%*%eta,alpha,X_new_dir,Y_new),sum(p1==Y_new)/n_validation,
            pred_acc(coef(res_dirc$mgpls)[-1],coef(res_dirc$mgpls)[1],X_new_dir,Y_new)
            
      ),res_dirc$ux)
      
   }
   acc_M_matrix <- matrix(unlist(acc),byrow = TRUE,ncol = 7)
   
   colnames(acc_M_matrix) <- c("envlp","fglm","fglmnet","true","classif.glm","fgpls","ux")
   
   
   if(is.null(mis_rate)){
      mis_rate <- signif(colMeans(acc_M_matrix),3)
   }
   else{
      mis_rate <- cbind(mis_rate,signif(colMeans(acc_M_matrix),3))
   }
}

print(mis_rate)

# acc_M_matrix <- matrix(unlist(acc),byrow = TRUE,ncol = 7)
# 
# colnames(acc_M_matrix) <- c("envlp","fglm","fglmnet","true","classif.glm","fgpls","ux")
# signif(colMeans(acc_M_matrix),3)
# warnings()
