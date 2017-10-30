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
#Plot functions
##########################
list2env(params,envir = parent.frame())
list2env(data,envir = parent.frame())


#Load conditional matrices

#Calibration, given those draws from predictive posterior
get_cal_draws <- function(X_pred){
  pred_probs <- est_probs_i(X_pred)
  Y_pred <- apply(pred_probs,2,function(x){
    return(rmultinom(1,num_counts,x))
  })
  return(Y_pred)
}

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
hist(q_cal)
plot(q_cal,type = 'l')
