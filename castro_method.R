#Testing quadrature stuff
library(dplyr)
#Load in data
setwd("/Users/Jake/Documents/SAMSI_Astro/")

mgCl <- log(t(read.table('marsfollowupstuff/MgCl_data.txt')))
mgCl_x <- as.numeric(as.matrix(read.table('marsfollowupstuff/MgCl_eV.txt')))

real <- read.table("marsfollowupstuff/Mg_May21_190030_2013_CCS.csv.dat")
exp_x <- real[which(real[,1]<max(mgCl_x) & real[,1]>min(mgCl_x)),1]
exp_amp <- real[which(real[,1]<max(mgCl_x) & real[,1]>min(mgCl_x)),2]

##First:omit some row of data
holdout <- 18
(subsample_rate = floor(length(mgCl_x)/length(exp_x)))
t_subset <- seq(1,length(mgCl_x),by = subsample_rate)

mgCl_trim <- mgCl[,t_subset]
mgCl_trim_train <- mgCl_trim[-holdout,]
mgCl_trim_test <- mgCl_trim[holdout,]

##Hmmm how do you apply the PCA rotation if the number of realizations 
##(i.e. points) are different?
###Some interpolation, as suggested in Castro_1986

trap_me <- function(x,y){
  n = length(y)
  idx = 2:(n-1)
  s <- c(y[1]*(x[2]-x[1])/2, y[idx]*(x[idx+1]-x[idx-1])/2, y[n]*(x[n]-x[n-1])/2)
  return(s)
}

#Find phi's for all paths using average W
W <- trap_me(mgCl_x[t_subset],apply(mgCl_trim_train,2,mean))

mhat <-apply(mgCl_trim_train,2,mean)
  
demean_train <- sweep(mgCl_trim_train,2,mhat)
C_hat <- cov(demean_train)

A = sweep(C_hat,2,sqrt(W),FUN="*") %>%
  sweep(1,sqrt(W),FUN="*")
e_A <- eigen(A)

psi <- e_A$vectors
phi <- sweep(psi,1,W^(-0.5),FUN="*")


alpha_list2 <- lapply(seq_along(1:dim(demean_train)[1]),function(sp){
  #W_sp <- trap_me(mgCl_x[t_subset],mgCl_trim_train[sp,])
  
  alphas <- t(phi)%*%diag(W)%*%matrix(demean_train[sp,])
  alphas
}) %>%
  do.call(cbind,.) %>%
  t()

alpha_mat <- t(alpha_list)

sp <- 15
W_sp <- trap_me(mgCl_x[t_subset],mgCl_trim_train[sp,])

alphas <- t(phi)%*%diag(W_sp)%*%matrix(demean_train[sp,])

#Recover sample path
r <- 2
test <- as.numeric(mhat + phi[,1:r]%*%matrix(alphas[1:r]))
plot(mgCl_x[t_subset],test,pch = 19,cex = .5,xlim = c(500,550))
lines(mgCl_x[t_subset],mgCl_trim_train[sp,],xlim = c(500,550))
#Get alpha's for each individual path by 
