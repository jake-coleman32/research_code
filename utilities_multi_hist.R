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

##Make the coefficient matrix
make_coefs <- function(type = 'sin'){
  #This is the "J" in the function for GP_hist_partial
  make_C_J <- length(alpha)
  if(type=='sin'){
    A <- matrix(0,N,make_C_J-1)
    B <- matrix(0,N,make_C_J-1)
    for(n in 1:N){
      for(j in 2:make_C_J){
        A[n,j-1] <- t_star*(sin(2*pi*n*alpha[j]/t_star) - sin(2*pi*n*alpha[j-1]/t_star))/(2*pi*n)
        B[n,j-1] <- t_star*(cos(2*pi*n*alpha[j-1]/t_star) - cos(2*pi*n*alpha[j]/t_star))/(2*pi*n)
      }
    }
    
    C <- cbind(t(A),t(B))
  }else if(type=='cos'){
    C <- matrix(0,N,(make_C_J-1))
    for(n in 1:N){
      for(j in 2:make_C_J){
        C[n,j-1] <- t_star*(sin(pi*n*alpha[j]/(t_star)) - sin(pi*n*alpha[j-1]/t_star))/(pi*n)
      }
    }
    C <- t(C)
  }else(stop('Need to choose sin or cos for type'))
  return(C)
}


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

#####################
##Running the Sampler
#####################

##Metropolis Hastings Sampler
hier_gp_mh_i <- function(iters = 1E4, burnin_prop = 0.1,
                         delta_lam = rep(0,Nx),
                         delta_ell = rep(0.3,Nx),
                         X_kap = replicate(Nx,rep(5E-1,I)), #column n is diagonal of proposal
                         #for GP n
                         verbose = FALSE
){
  burnin <- iters*burnin_prop
  
  X_array <- array(0,dim = c(iters,Nx,I)) #Columns are GPs, rows are histograms
  lam_mat <- ell_mat <- matrix(0,iters,Nx)
  
  #Current values
  X_cur <- matrix(rnorm(Nx*I),Nx,I) #currently Nx x I
  ell_cur <- rgamma(Nx,ell_a,ell_b)
  lam_cur <- rep(1,Nx)
  
  p_cur <- t(sqrt(2)*C%*%sweep(X_cur,1,r_vec,"*") + replicate(I,b/t_star))
  
  if(sum(p_cur<0)) stop("Unlucky, try again")
  l_cur <- sum(Y_trunc*log(p_cur))# + sum(Y_trunc[,j_trim]*log(1-apply(p_cur,1,sum)))
  
  ell_acc <- lam_acc  <- numeric(Nx)
  x_acc <- matrix(0,Nx,I)
  
  time_start <- proc.time()
  for(t in 1:(burnin + iters)){
    for(n in 1:Nx){
      if(!t%%100){
        if(verbose){
          flush.console()
          cat("\r t =", t, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
          #print(paste0("t = ", t, ", n = ", n))
        }
      }
      
      #Update X_n
      auto_reject = FALSE
      (x_cur_n <- X_cur[n,])
      
      #Need to find upper and lower bounds for each histogram
      #  l <- -Inf*(numeric(I) + 1)
      # u <- Inf*(numeric(I) + 1)
      for(i in 1:I){
        
        
        #All the constraints
        constr <- (-b/t_star - sqrt(2)*C[,-n]%*%(r_vec[-n]*X_cur[-n,i]))/(sqrt(2)*r_vec[n]*C[,n])
        
        #Lower constraints
        (a_constr <- constr[which(round(C[,n],10)>0)])
        
        #Upper constraints
        (b_constr <- constr[which(round(C[,n],10)<0)])
        
        #Bad constraints
        (bad_constr <- which(round(C[,n],10)==0))
        
        #If any bad constraints
        if(length(bad_constr)){
          for(j in bad_constr){
            
            #For each bad constraint, check if the rest of X[i,] take care of it
            #If they don't auto-reject the whole vector
            if(b[j]/t_star + sqrt(2)*C[j,-n]%*%(r_vec[-n]*X_cur[-n,i])<=0){
              print(paste("Bad times: i=",i,"n =",n,"j =",j))
              auto_reject = TRUE
              l_star = -Inf
            }
          }
        }
        
        #Update X_n[i]
        (l <- ifelse(length(a_constr),max(a_constr),-Inf))
        (u <- ifelse(length(b_constr),min(b_constr),Inf))
        x_cur_i <- x_cur_n[i]
        (x_st_i <- rtruncnorm(1,mean = x_cur_i,a = l, b = u,sd = X_kap[i,n]))
        
        (x_cond <- x_cur_n[-i])
        (sig_11 = 1/lam_cur[n])
        (sig_22 <- GP_cov(type = 'sqr.exp', d_scaled[-i],lam_cur[n],ell_cur[n],nugget = 1E-5))
        (sig_12 <- t(matrix(GP_cross_cov(type = 'sqr.exp', d_scaled[i],d_scaled[-i],lam_cur[n],ell_cur[n]))))
        
        (prior_mean = sig_12%*%(solve(sig_22)%*%x_cond))
        (prior_var = sig_11 - sig_12%*%solve(sig_22)%*%t(sig_12))
        
        X_star <- X_cur
        X_star[n,i] <- x_st_i
        
        
        p_star <- t(sqrt(2)*C%*%(sweep(X_star,1,r_vec,"*")) + replicate(I,b/t_star))
        if(sum(p_star<0)){
          print(paste("j = ",j,"i = ",i,"t = ",t,"n = ",n))
          stop("Dammit this shouldn't happen")
        }
        
        l_star <- sum(Y_trunc*log(p_star))
        
        adjust <- dtruncnorm(x=x_cur_i,a = l,b=u,mean = x_st_i,sd = X_kap[i,n]) -
          dtruncnorm(x=x_st_i,a = l,b=u,mean = x_cur_i,sd = X_kap[i,n])
        
        ratio <- l_star - l_cur +
          
          #Prior
          dnorm(x_st_i,mean = prior_mean, sd = sqrt(prior_var)) - 
          dnorm(x_cur_i,mean = prior_mean, sd = sqrt(prior_var)) + 
          
          #Adjustment ratio
          adjust
        
        if(length(ratio)>1){stop(paste("wtf n=",n))}
        if(runif(1)<exp(ratio)){
          x_acc[n,i] <- x_acc[n,i] + 1
          l_cur <- l_star
          X_cur[n,i] <- x_st_i
        }
        
        
      }# i loop
      
      #Calculate proposed likelihood based on proposed x_star_n
      
      #############
      #Update ell_n
      #############
      z_ell <- rnorm(1)
      ell_st <- ell_cur[n]*exp(delta_ell[n]*z_ell)
      adjust_ratio <- delta_ell[n]*z_ell
      
      #MH Ratio
      
      ratio <- 
        #"Likelihood" - X prior
        dmvnorm(X_cur[n,],sigma = GP_cov(type = 'sqr.exp', d_scaled,lam_cur[n],ell_st,nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[n,],sigma = GP_cov(type = 'sqr.exp', d_scaled,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
        #Priors
        dgamma(ell_st,ell_a[n],ell_b[n],log = TRUE) -
        dgamma(ell_cur[n],ell_a[n],ell_b[n],log = TRUE) +
        
        #Proposal Adjustment
        adjust_ratio
      
      if(ratio + rexp(1)>0){
        ell_acc[n] <- ell_acc[n] + 1
        ell_cur[n] <- ell_st
      }
      
      ############
      #Update lam_n
      ############
      z_lam <- rnorm(1)
      lam_st <- lam_cur[n]*exp(delta_lam[n]*z_lam)
      adjust_ratio <- delta_lam[n]*z_lam
      
      ratio <- 
        #"Likelihood" - X prior
        dmvnorm(X_cur[n,],sigma = GP_cov(type = 'sqr.exp', d_scaled,lam_st,ell_cur[n],nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[n,],sigma = GP_cov(type = 'sqr.exp', d_scaled,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
        #Priors
        dgamma(lam_st,lam_a[n],lam_b[n],log = TRUE)-
        dgamma(lam_cur[n],lam_a[n],lam_b[n],log = TRUE) +
        
        #Proposal Adjustment
        adjust_ratio
      
      
      if(ratio + rexp(1)>0){
        lam_acc[n] <- lam_acc[n] + 1
        lam_cur[n] <- lam_st
      }
      
      
    }# n loop
    
    if(t > burnin){
      X_array[t-burnin,,] <- X_cur
      lam_mat[t-burnin,] <- lam_cur
      ell_mat[t-burnin,] <- ell_cur
    }
  }# t loop
  
  return(list(X = X_array,lam = lam_mat,ell = ell_mat,
              x_acc = x_acc/(iters+burnin),
              lam_acc = lam_acc/(iters+burnin),
              ell_acc = ell_acc/(iters+burnin)))
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


