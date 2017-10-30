library(truncnorm)
library(mvtnorm)
library(dplyr)
library(caTools)
library(Matrix)
library(fields)
library(coda)

#Load utility functions
source('/Users/Jake/Dropbox/Research/Computer_Emulation/R_code/utilities_multi_hist.R')

#File to show pictures/look at plots easily
setwd('/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/long_run/')
load('data_list.Rdata')

#Change to what you need
setwd('N_9_cos_only/')
load('params_list.Rdata')
#######################
#Plot functions
##########################
list2env(params,envir = parent.frame())
list2env(data,envir = parent.frame())


#Trace plots for components of a given histogram
plot_traces <- function(x_hist,save_pics = FALSE){
  
  for(n in 1:Nx){
    if(n<=length(Nz)){
      var_type = 'Z'
      index = n
    }
    else{
      var_type = 'W'
      index = n-length(Nz)
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
  for(n in 1:Nx){
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
est_dens_i <- function(x_mat,basis = 'sin'){
  f_mat <-matrix(0,dim(x_mat)[1],length(T_out))
  Nx <- dim(x_mat)[2]
  for(t in 1:length(T_out)){
    if(basis=='sin'){
      cos_vec <- cos(2*pi*T_out[t]*Nz/t_star)
      sin_vec <- sin(2*pi*T_out[t]*Nw/t_star)
      cos_sin_vec <- c(cos_vec,sin_vec)*r_vec
    }else if(basis=='cos'){
      cos_sin_vec <- cos(pi*Nz*(T_out[t])/(t_star))*r_vec
    }
    
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
plot_dens_i <- function(x_mat,basis = 'sin', save_pics = FALSE,legend_side = 'topright',...){
  f_est <- est_dens_i(x_mat,basis = basis)
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
#plot_dens_i(X_pred,basis = basis_type, save_pics = FALSE) #Density


#est_probs_i() is a utility function
#create_cond_mats() is a utility function
#pred_x() is a utility function


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
  cols = heat.colors(I)
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
dim(X)

#Extra burnin if necessary
#extra_burn = min(which(X[,9,]<(-1)))
extra_burn = 1
X_burn <- X[extra_burn:dim(X)[1],,]
lam_burn <- hope_i$lam[extra_burn:dim(X)[1],]
ell_burn <- hope_i$ell[extra_burn:dim(X)[1],]

#Thinning if necessary
#thin = 500
thin = 200
iters <- 1:dim(X_burn)[1]
X_thin <- X_burn[!iters%%thin,,]
lam_thin <- lam_burn[!iters%%thin,]
ell_thin <- ell_burn[!iters%%thin,]

#Which in-sample histogram to look at
i <- 4
X_i <- X_thin[,,i]
dim(X_i)
effectiveSize(X_i)


meeting_parent <- '/Users/Jake/Dropbox/Research/Computer_Emulation/meetings/2017/'
meeting_folder <- 'meeting_10_26/'
path <- paste0(meeting_parent,meeting_folder)
save_pics = FALSE
suffix = '_cos_only'


#Look for evidence of non-convergence
plot_traces(x_hist = X_i,save_pics = save_pics)
plot_thetas(save_pics = save_pics)

#Get predictive Xs
run_cond_mats <- create_cond_mats(X_mat_all = X_thin,lam_mat = lam_thin,ell_mat = ell_thin)
system.time(X_pred <- pred_x(run_cond_mats, d_prime = d_new_s, d_cond = d_scaled))

#Make some predictive plots
gam=1
T_out=seq(0,t_star,length=101)

basis_type = 'cos'
plot_dens_i(X_pred,basis=basis_type, save_pics = save_pics) #Density
bin_prob_comp(X_pred = X_pred,save_pics = save_pics) #Bin probs
resid_prob_comp(X_pred = X_pred,save_pics = save_pics) #Residuals

plot_true_bins(save_pics = save_pics)#True bins
plot_pred_bins(X_mats_all = X_thin,X_pred = X_pred,save_pics = save_pics) #predicted bins


##############################
##Scratch work





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


