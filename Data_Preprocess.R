t_rescale <- function(t0){
   return((t0-min(t0))/(max(t0)-min(t0)))
}



na_keep <- function(X_list,Y){
   idx_row <- NULL
   for (X in X_list) {
      idx_row <- union(idx_row,(1:dim(X)[1])[rowSums(is.na(X))!=0])
   }
   Y <- as.matrix(Y)
   idx_row <- union(idx_row, (1:dim(Y)[1])[rowSums(is.na(Y))!=0])
   
   idx_row <- setdiff(1:dim(X)[1],idx_row)
   return(idx_row)
}



data_centered <- function(X){
   X.mean <- colMeans(X)
   p <- length(X.mean)
   n <- dim(X)[1]
   return(vapply(1:p,function(i) X[,i]-X.mean[i],numeric(n)))
}



remove_na <- function(data_list){
   nd <- nrow(data_list[[1]])
   for (x in data_list) {
      if(nd!=nrow(x)){
         raise("All the data should have the same number of rows!")
      }
   }

   
   idxRomve <- NULL
   for (x in data_list) {
      idxRomve <- c(idxRomve,(1:nd)[rowSums(is.na(x))>0])
   }
   idxRomve <- unique(idxRomve)
   
   data_list_new <- list()
   
   for (x in data_list) {
      data_list_new <- append(data_list_new,list(x[-idxRomve,]))
   }
   names(data_list_new) <- names(data_list)
   return(data_list_new)
}
