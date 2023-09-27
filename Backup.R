library(refund)
library(fda.usc)
library(doFuture)
library(doRNG)
source("FEPLS.R")
data("content")
Sys.setenv(LANG="en")
library(dplyr)
library(tidyverse)


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
