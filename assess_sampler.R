#File to show pictures/look at plots easily
setwd('/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/long_run/')
load('data_list.Rdata')

#Change to what you need
setwd('good_r_p5_c_1/')
load('params_list.Rdata')
#######################
#Plot functions
##########################
list2env(params,envir = parent.frame())
list2env(data,envir = parent.frame())

GP_cov <- function(d,lambda,ell,nugget = 0.){
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(d))^2/ell^2) + nugget*diag(length(d))
  
  return(out_mat)
}

GP_cross_cov <- function(d,d_star,lambda,ell){
  inds <- 1:length(d)
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(c(d,d_star)))[inds,-inds]^2/ell^2)
  
  return(out_mat)
}

#Trace plots for components of a given histogram
plot_traces <- function(x_hist,save_pics = FALSE){
  
  for(n in 1:N){
    if(n<=floor(N/2)){
      var_type = 'Z'
      index = n
    }
    else{
      var_type = 'W'
      index = n-ceiling(N/2) 
    }
    if(save_pics) pdf(paste0(path,var_type,n,suffix,'.pdf'))
    plot(x_hist[,n],type = 'l',ylab = bquote(.(var_type)[.(index)]),
         main = bquote('Trace Plot of '~.(var_type)[.(index)]~', Histogram'~.(i)),
         xlab = 'Iteration')
    if(save_pics) dev.off()
  }
}

#Ell and Lambda parameter trace plots
#Currently lambda set to 1, so only plotting ell
plot_thetas <- function(save_pics = FALSE){
  for(n in 1:N){
    #Lambda
    #if(save_pics) pdf(file = paste0(path,'lam',n,suffix,'.pdf'))
    #plot(lam_thin[,n],xlab = "Iteration",ylab = bquote(lambda[.(n)]),type = 'l',
    #     main = bquote('Trace Plot of '~lambda[.(n)]))
    #if(save_pics) dev.off()
    
    #Ell
    if(save_pics) pdf(file = paste0(path,'ell',n,suffix,'.pdf'))
    plot(ell_thin[,n],xlab = "Iteration",ylab = bquote('ell'[.(n)]),type = 'l',
         main = bquote('Trace Plot of ell'[.(n)]))
    if(save_pics) dev.off()
  }
}

#Estimate density with mean and credible intervals, given a histogram's parameters
#X_mat should be T x N matrix
est_dens_i <- function(x_mat){
  f_mat <-matrix(0,dim(x_mat)[1],length(T_out))
  Nx <- dim(x_mat)[2]
  for(t in 1:length(T_out)){
    cos_vec <- cos(2*pi*T_out[t]*Nz/t_star)
    sin_vec <- sin(2*pi*T_out[t]*Nw/t_star)
    cos_sin_vec <- c(cos_vec,sin_vec)*r_vec
    
    f_mat[,t] <- sqrt(2)*t(matrix(cos_sin_vec))%*%t(x_mat) + gam/t_star
  }
  return(f_mat)
}

add_hist <- function(Y = Y_new_trunc,add = TRUE){
  #Requires alpha
  
  Y_fake_dat <- c()
  for(j in 1:length(Y)){
    Y_fake_dat <- c(Y_fake_dat, runif(Y[j], alpha[j],alpha[j+1]))
  }
  hist(Y_fake_dat, breaks = alpha,probability = 1,add = add,
       col = rgb(1,0,0,0.1))
}

#Plot the predicted density, given the T x N matrix of components
plot_dens_i <- function(x_mat, save_pics = FALSE,legend_side = 'topright',...){
  f_est <- est_dens_i(x_mat)
  mean_est <- apply(f_est,2,mean)
  
  if(save_pics) pdf(paste0(meeting_parent,meeting_folder,'mh_dens',suffix,'.pdf'))
  plot(T_out,mean_est,type = 'l', main = 'GP Density Estimate',
       ylab = 'Density',xlab = 'y',lwd = 2,
       ylim = c(0,max(apply(f_est,2,quantile,0.975))),
       ...)
  lines(T_out,apply(f_est,2,quantile,0.025),col = 'blue',lty=2)
  lines(T_out,apply(f_est,2,quantile,0.975),col = 'blue',lty=2)
  # legend(legend_side,c('Post Mean','Post 95% Cred'),
  #       lwd = 2,lty = c(1,2),col = c('black','blue'))
  
  add_hist()
  
  colvec <- c(rgb(t(col2rgb('black'))),
              rgb(t(col2rgb('blue')),max = 255),
              rgb(t(col2rgb('red')),max = 255,alpha = 50))
  legend('topright',c('Post Mean', '95% Interval', 'Truth'),
         lty = c(1,2,0), lwd = c(1,1,0),bty = "n",
         col = colvec,
         pch = c(NA,NA, 15),
         pt.cex = 2)
  if(save_pics) dev.off()
}
plot_dens_i(X_pred,save_pics = FALSE) #Density

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
pred_x <- function(cond_mats, d_prime, d_cond=d_scaled, verbose = FALSE){
  #Loop through N components, collect values
  #500k x 9 matrix
  
  sig_22_inv_list <- cond_mats$cov_inv_mats
  sig_inv_x <- cond_mats$inv_vec_mult
  
  lam = cond_mats$lam_mat
  ell = cond_mats$ell_mat
  
  n_iters <- length(sig_22_inv_list[[1]])
  pred_components <- matrix(0,n_iters,N)
  
  
  for(n in 1:N){
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
    pred_components[,n] <- rnorm(n_iters,pred_mean,sqrt(pred_var))
    #pred_components[,n] <- pred_mean
    
  }
  
  return(pred_components)
}


##
##Bin probability comparison
##

#Counts with cred intervals compared to truth
bin_prob_comp <- function(X_pred,save_pics = FALSE){
  if(save_pics) pdf(paste0(path,'pred_comp',suffix,'.pdf'))
  plot(as.numeric(colnames(Y_trunc)),apply(est_probs_i(X_pred),1,mean),
       cex = 0.9, ylim = c(min(apply(est_probs_i(X_pred),1,quantile,probs = 0.025)),
                           max(apply(est_probs_i(X_pred),1,quantile,probs = 0.975))),main = 'Predicted Bin Probabilities - Full Model',
       xlab = 'Bin',ylab = 'Predicted Probability',col = rgb(1,0,0,0.3),pch = 19,
       las = 1)
  arrows(as.numeric(colnames(Y_trunc)), apply(est_probs_i(X_pred),1,quantile,probs = 0.025),
         as.numeric(colnames(Y_trunc)), apply(est_probs_i(X_pred),1,quantile,probs = 0.975), length=0.05, angle=90, code=3)
  points(as.numeric(colnames(Y_trunc)),Y_new_trunc/num_counts,
         col = rgb(0,0,1,0.3),cex = 0.8,pch = 19)
  legend('bottomright',c('Predicted','Truth'),col = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),pch = 19,
         pt.cex = c(0.9,0.8))
  if(save_pics) dev.off()
}

#Residuals with cred intervals compared to zero
resid_prob_comp <- function(X_pred,save_pics = FALSE){
  if(save_pics) pdf(paste0(path,'pred_resid',suffix,'.pdf'))
  plot(as.numeric(colnames(Y_trunc)),
       apply(est_probs_i(X_pred),1,mean) - Y_new_trunc/num_counts,
       cex = 0.9, 
       ylim = c(min(apply(est_probs_i(X_pred) - Y_new_trunc/num_counts,1,quantile,probs = 0.025)),
                max(apply(est_probs_i(X_pred) - Y_new_trunc/num_counts,1,quantile,probs = 0.975))),
       main = bquote(A[j]~'Bin Probability Residuals - Full Model'),
       xlab = bquote(A[j]~'Bin'),ylab = 'Residual',col = rgb(1,0,0,0.7),pch = 19,
       las = 1, yaxt = 'n')
  axis(2, at=c(-0.004,0,0.004),las = 1,labels = FALSE)
  text(y=c(-0.004,0,0.004),par("usr")[1],labels = c(-0.004,0,0.004),pos = 2,xpd = TRUE)
  
  arrows(as.numeric(colnames(Y_trunc)),apply(sweep(est_probs_i(X_pred),1,Y_new_trunc/num_counts,FUN = "-"),1,quantile,probs = 0.025),
         as.numeric(colnames(Y_trunc)), apply(sweep(est_probs_i(X_pred),1,Y_new_trunc/num_counts,FUN = "-"),1,quantile,probs = 0.975), length=0.05, angle=90, code=3)
  abline(h = 0, col = rgb(0,0,1,0.7))
  if(save_pics) dev.off()
}

##
##Connected histograms, truth vs sampler prediction
##

#Plotting the data points, connecting bins of same histogram
#Marking the holdout sample
plot_true_bins <- function(save_pics = FALSE){
  cols = rainbow(I)
  if(save_pics) pdf(file = paste0(path,'all_data',suffix,'.pdf'))
  plot(as.numeric(colnames(Y_trunc)),
       Y_trunc[1,]/num_counts,col = cols[1],type = 'l',
       xlab = 'Aj',
       ylab = 'Fraction of Bin',
       main = 'Data Across Input')
  for(i in 2:I){
    lines(as.numeric(colnames(Y_trunc)),
          Y_trunc[i,]/num_counts,col = cols[i],type = 'l')
  }
  
  lines(as.numeric(colnames(Y_trunc)),Y_new_trunc/num_counts,lwd = 2)
  legend('bottomright','Out of Sample Histgram',lwd = 2)
  if(save_pics) dev.off()
}

#Plotting predicted bins, connecting for same hist
#Marking out-of-sample
plot_pred_bins <- function(X_mats_all, X_pred,save_pics = FALSE){
  cols = rainbow(I)
  
  #Now lets plot our estiamtes
  #Doesn't vary more than the 95% confidence intervals
  if(save_pics) pdf(file = paste0(path,'all_emulations',suffix,'.pdf'))
  plot(as.numeric(colnames(Y_trunc)),
       apply(est_probs_i(X_mats_all[,,1]),1,mean),col = cols[1],type = 'l',
       xlab = 'Aj',
       ylab = 'Probability of Bin',
       main = 'Emulated Values Across Input')
  for(i in 2:I){
    lines(as.numeric(colnames(Y_trunc)),
          apply(est_probs_i(X_mats_all[,,i]),1,mean),col = cols[i],type = 'l')
  }
  
  lines(as.numeric(colnames(Y_trunc)),apply(est_probs_i(X_pred),1,mean),
        lwd = 2)
  legend('bottomright','Out of Sample Prediction',lwd = 2)
  if(save_pics) dev.off() 
}

#Look at the coefficients, see if the 95% credible interval
# overlaps with zero
plot_coefs <- function(xmat){
  means <- apply(xmat,2,mean)
  plot(means,pch = 19, col = 'blue', xlab = "Component",ylab = 'Value',
       main = 'Components 95% Credible Intervals',
       ylim = c(min(apply(xmat,2,quantile,probs = 0.025)),
                max(apply(xmat,2,quantile,probs = 0.975))))
  abline(h = 0,col = 'red')
  
  arrows(1:dim(xmat)[2], apply(xmat,2,quantile,probs = 0.025),
         1:dim(xmat)[2], apply(xmat,2,quantile,probs = 0.975), 
         length=0.05, angle=90, code=3)
  
}

##########
#Load what you need
#########

load('sampler_vals.Rdata')
X <- hope_i$X
apply(hope_i$x_acc,1,mean)

#Extra burnin if necessary
#extra_burn = min(which(X[,9,]<(-1)))
extra_burn = 1
X_burn <- X[extra_burn:dim(X)[1],,]
lam_burn <- hope_i$lam[extra_burn:dim(X)[1],]
ell_burn <- hope_i$ell[extra_burn:dim(X)[1],]

#Thinning if necessary
thin = 500
#thin = 200
iters <- 1:dim(X_burn)[1]
X_thin <- X_burn[!iters%%thin,,]
lam_thin <- lam_burn[!iters%%thin,]
ell_thin <- ell_burn[!iters%%thin,]

#Which in-sample histogram to look at
i <- 4
X_i <- X_thin[,,i]


meeting_parent <- '/Users/Jake/Dropbox/Research/Computer_Emulation/meetings/2017/'
meeting_folder <- 'meeting_4_13/'
path <- paste0(meeting_parent,meeting_folder)
save_pics = FALSE
suffix = '_old_good'


#Look for evidence of non-convergence
plot_traces(x_hist = X_i,save_pics = FALSE)
plot_thetas(save_pics = FALSE)

#Get predictive Xs
run_cond_mats <- create_cond_mats(X_mat_all = X_thin,lam_mat = lam_thin,ell_mat = ell_thin)
system.time(X_pred <- pred_x(run_cond_mats, d_prime = d_new_s, d_cond = d_scaled))

#Make some predictive plots
gam=1
T_out=seq(0,t_star,length=101)

plot_dens_i(X_pred,save_pics = FALSE) #Density
bin_prob_comp(X_pred = X_pred,save_pics = FALSE) #Bin probs
resid_prob_comp(X_pred = X_pred,save_pics = FALSE) #Residuals

plot_true_bins(save_pics = FALSE)#True bins
plot_pred_bins(X_mats_all = X_thin,X_pred = X_pred,save_pics = FALSE) #predicted bins


##############################
##Scratch work

#Calibration, given those draws from predictive posterior
get_cal_draws <- function(X_pred){
  pred_probs <- est_probs_i(X_pred)
  Y_pred <- apply(pred_probs,2,function(x){
    return(rmultinom(1,num_counts,x))
  })
  return(Y_pred)
}

cal_qhat <- function(kap, num_cals = 1E3){
  #Setup
  auto_reject = FALSE
  accept <- 0 
  q_cals <- numeric(num_cals)
  
  #Current values
  q_cur <- runif(1)
  X_pred_cur <- pred_x(run_cond_mats, d_prime = q_cur, d_cond = d_scaled)
  p_cur <- apply(est_probs_i(X_pred_cur),1,mean)

  time_start <- proc.time()
  #Go through loop
  for(i in 1:num_cals){
    
    #Propose new qhat
    q_st <- rnorm(1,q_cur,kap)
    if(q_st<0 | q_st > 1) auto_reject = TRUE
    
    #predict Y's with qhat
    X_pred_st <- pred_x(run_cond_mats, d_prime = q_st, d_cond = d_scaled)
    p_st <- apply(est_probs_i(X_pred_cur),1,mean)
    
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
    
    if(!i%%100){
      flush.console()
      cat("\r i = ", i, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
      #print(paste0("t = ", t, ", n = ", n))
    }    
  }
  print(accept/num_cals)
  return(list(q_cals = q_cals, acc_ratio = accept/num_cals))
}

cal_try <- cal_qhat(kap = 0.1)
cal_try$acc_ratio

q_cal <- cal_try$q_cals
hist(q_cal)
plot(q_cal,type = 'l')



#Checking out singularity in C
for(n in 1:dim(C)[2]){
  View(round(sweep(C,1,C[,n] + 1E-5,"/"),3))
}

View(round(sweep(C,1,C[,4],"/"),3))

#Plotting predictive probabilities by the predicted density
new_t <- vector("list",J)
pred_dens <- est_dens_i(X_pred)
pred_ps <- matrix(0,dim(pred_dens)[1],J)
for(a in 1:(length(alpha)-1)){
  new_t <- which(T_out>=alpha[a]&T_out<alpha[a+1])
  if(a ==J){new_t = c(new_t,which(T_out==alpha[a+1]))}
  print(new_t)
  pred_ps[,a] <- apply(pred_dens[,new_t],1,function(x){trapz(T_out[new_t],x)})
}
test <- apply(pred_dens,1,function(x){trapz(T_out,x)})

summary(test)
summary(apply(pred_ps,1,sum))#why does this not give 1


sum(pred_ps)


