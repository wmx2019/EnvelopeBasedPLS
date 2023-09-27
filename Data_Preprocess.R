t_rescale <- function(t0){
   return((t0-min(t0))/(max(t0)-min(t0)))
}

na_remove <- function(X_list,Y){
   idx_row <- NULL
   for (X in X_list) {
      idx_row <- c(idx_row,(1:dim(X)[1])[rowSums(is.na(X))==0])
   }
   Y <- as.matrix(Y)
   idx_row <- c(idx_row, (1:dim(Y)[1])[rowSums(is.na(Y))==0])
   return(idx_row)
}
