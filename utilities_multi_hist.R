##Utility Functions for
##assess_sampler.R
##hist_calibration.R
##hier_GP_indv_updates.R

library(truncnorm)
library(mvtnorm)
library(dplyr)
library(caTools)
library(Matrix)
library(fields)
library(coda)

#GP between across the input space - covariance and cross-covariance matrices
GP_cov <- function(d,lambda,ell,type = 'sqr.exp', nu = 3/2, nugget = 0.){
  d_mat <- as.matrix(dist(d))
  if(type=="sqr.exp"){
    out_mat <- lambda^(-1)*exp(-d_mat^2/ell^2) + nugget*diag(length(d))
  }else if(type=="matern"){
    out_mat <- lambda^(-1)*Matern(d_mat,range=ell,smoothness = nu) + nugget*diag(length(d))
  }else{
    stop('Wrong covariance function type')
  } 
  
  
  return(out_mat)
}

GP_cross_cov <- function(d,d_star,lambda,ell,type = 'sqr.exp',nu = 3/2){
  inds <- 1:length(d)
  d_mat <- as.matrix(dist(c(d,d_star)))
  if(type=="sqr.exp"){
    out_mat <- lambda^(-1)*exp(-d_mat[inds,-inds]^2/ell^2)
  } else if(type=="matern"){
    out_mat <- lambda^(-1)*Matern(d_mat[inds,-inds],range = ell, smoothness = nu)
  }else{ 
    stop('Wrong covariance function type')
  }
  
  return(out_mat)
}


###Assessing the sampler

#Estimate the predicted bin probabilities given histogram components X_i
est_probs_i <- function(X_i){
  p_out <- sqrt(2)*C%*%sweep(t(X_i),1,r_vec,FUN = "*") + b/t_star
  return(p_out)
}

#Create conditioning matrices - mostly covariance matrices that predictions condition on
#(Sigma_22)^{-1} doesn't need to be computed more than once for each [t,n] pair
#Currently this doesn't save any time, but if we need to predict in an MCMC loop it'll be necessary
create_cond_mats <- function(X_mat_all,lam_mat,ell_mat){
  T_samp <- dim(X_mat_all)[1]
  N <- dim(X_mat_all)[2]
  I <- dim(X_mat_all)[3]
  
  #(Sigma_22)^{-1}
  cov_inv_mats <-  vector("list",N)
  
  #(Sigma_22)^{-1}(X_condition)
  inv_vec_mult <- array(0,dim=c(T_samp,N,I))#T x N x I
  
  #Create the matrices, store them in a big list of lists
  for(n in 1:N){
    cov_inv_mats[[n]] <- vector("list",T_samp)
    for(t in 1:T_samp){
      lam_tn <- lam_mat[t,n]
      ell_tn <- ell_mat[t,n]
      cov_inv_mats[[n]][[t]] <- solve(GP_cov(d_scaled,lambda = lam_tn,ell = ell_tn,nugget = 1E-6))
      inv_vec_mult[t,n,] <- cov_inv_mats[[n]][[t]]%*%X_thin[t,n,]
    }
  }
  
  return(list(cov_inv_mats = cov_inv_mats,inv_vec_mult = inv_vec_mult,
              ell_mat = ell_mat,lam_mat = lam_mat))
}

#Predict components for calculating predictive distribution of density
#Essentially integrate them away
pred_x <- function(cond_mats, d_prime, d_cond=d_scaled, verbose = FALSE,
                   mean_only=TRUE){
  #Loop through N components, collect values
  #500k x 9 matrix
  
  sig_22_inv_list <- cond_mats$cov_inv_mats
  sig_inv_x <- cond_mats$inv_vec_mult
  
  lam = cond_mats$lam_mat
  ell = cond_mats$ell_mat
  
  n_iters <- length(sig_22_inv_list[[1]])
  pred_components <- matrix(0,n_iters,Nx)
  
  
  for(n in 1:Nx){
    pred_mean <- pred_var <- n_iters
    time_start <- proc.time()
    for(t in 1:n_iters){
      sig_22_inv <- sig_22_inv_list[[n]][[t]]
      sig_12 <- t(GP_cross_cov(d_cond, d_prime,lambda = lam[t,n],ell = ell[t,n]))
      
      pred_mean[t] <- sig_12%*%sig_inv_x[t,n,]
      pred_var[t] <- 1/lam[t,n] - sig_12%*%sig_22_inv%*%t(sig_12)
      if(pred_var[t]<0)stop(paste("pred_var is",pred_var[t],"nt is",n,t)) 
      
    }
    
    
    #Swap between these two to perharps make run faster
    if(mean_only) pred_components[,n] <- pred_mean
    else pred_components[,n] <- rnorm(n_iters,pred_mean,sqrt(pred_var))
    
    
  }
  
  return(pred_components)
}


