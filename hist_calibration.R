#Trying to do some calibration based on posterior GP estiamted values
#assume we're on the school computer

source('/home/grad/jrc71/Documents/Research/Computer_Emulation/research_code/utilities_multi_hist.R')

#File to show pictures/look at plots easily
setwd('/home/grad/jrc71/Documents/Research/Computer_Emulation/JETSCAPE/JETSCAPE-STAT/run_calibration')
load('data_list.Rdata')

#Change to what you need
setwd('N_9_cos_only/')
load('params_list.Rdata')
#######################
#Load necessary values
##########################

###Load parameters
list2env(params,envir = parent.frame())
# c - constant in r_vec: c*r^n
# Nx - total number of coefficient parameters per histogram
# r - ratio in r_vec: c*r^n
# Nz - index vector of coefficients to include for cosine basis values
# Nw - index vector of coefficients to include for sine basis value
# C - J x Nx vector of integrals of basis coefficients across bins
# basis_type - character of "sin" or "cos" determining which basis to use
# X_kap - I x Nx matrix of individual MH sigma tuning values


###Load data - shouldn't change between runs
list2env(data,envir = parent.frame())
# num_counts - number of total counts in histogram
# holdout - which q_hat was held out in training GP
# Y - all data with all q_hats
# Y_new - Y[holdout, ]
# Y_trunc - Y[-holdout,], truncated to t_star (0.5, because the data stops there)
# Y_new_trunc - Y_new truncated to t_star (0.5)
# t_star - truncation value on A_j - data only goes to 0.5 here
# d_scaled - q_hat design scaled to [0,1]
# d_new_s - holdout q_hat on [0,1] scale
# trunc_cols - index list of columns retained in Y_trunc
# I - number of histograms
# J - number of bins in Y_trunc
# alpha - bin edges (vector of size J+1)
# b - bin sizes (vector of size J)



##Sampler values for validated GP
load('sampler_vals.Rdata')
X <- hope_i$X
apply(hope_i$x_acc,1,mean)
dim(X)

#Extra burnin if necessary - currently none
extra_burn = 1
X_burn <- X[extra_burn:dim(X)[1],,]
lam_burn <- hope_i$lam[extra_burn:dim(X)[1],]
ell_burn <- hope_i$ell[extra_burn:dim(X)[1],]

#Thinning if necessary
thin = 200
iters <- 1:dim(X_burn)[1]
X_thin <- X_burn[!iters%%thin,,]
lam_thin <- lam_burn[!iters%%thin,]
ell_thin <- ell_burn[!iters%%thin,]


##Calibration function
cal_qhat <- function(kap, niters = 1E3){
  #Setup
  auto_reject = FALSE
  accept <- 0 
  q_cals <- numeric(niters)
  
  #Current values
  q_cur <- runif(1)
  
  #Only if we take the mean in pred_x, rather than rnorm
  X_pred_cur <- pred_x(run_cond_mats, d_prime = q_cur, d_cond = d_scaled)
  p_cur <- apply(est_probs_i(X_pred_cur),1,mean)
  
  
  time_start <- proc.time()
  #Go through loop
  for(i in 1:niters){
    
    #Propose new qhat
    q_st <- rnorm(1,q_cur,kap)
    if(q_st<0 | q_st > 1) auto_reject = TRUE
    
    #predict Y's with qhat
    X_pred_st <- pred_x(run_cond_mats, d_prime = q_st, d_cond = d_scaled)
    p_st <- apply(est_probs_i(X_pred_cur),1,mean)
    
    #Also get the same for the current q_hat
    #Only have this if we do rnorm in pred_x, rather than just take the mean
    #X_pred_cur <- pred_x(run_cond_mats, d_prime = q_cur, d_cond = d_scaled)
    #p_cur <- apply(est_probs_i(X_pred_cur),1,mean)
    
    #calculate likelihood of predicted Y's compared to actual
    log_ratio <- dmultinom(Y_new_trunc,prob = p_st,log = TRUE) - 
      dmultinom(Y_new_trunc,prob = p_cur,log = TRUE) 
    
    if(runif(1)<exp(log_ratio) & !auto_reject){
      p_cur <- p_st
      q_cur <- q_st
      accept <- accept + 1
    }
    
    q_cals[i] <- q_cur
    auto_reject <- FALSE
    
    if(!i%%(niters*0.1)){
      flush.console()
      cat("\r i = ", i, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
    }    
  }
  print(accept/niters)
  return(list(q_cals = q_cals, acc_ratio = accept/niters))
}

cal_try <- cal_qhat(kap = 0.5,niters=100)
save(cal_try,save_file = 'cal_try.Rdata')
stop('We\'re done here')


cal_try$acc_ratio

q_cal <- cal_try$q_cals
hist(q_cal, main = expression('Posterior Distribution of '~hat(q)),prob = TRUE,
     xlab = expression(hat(q)))
plot(q_cal,type = 'l')
