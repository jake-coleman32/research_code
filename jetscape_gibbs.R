library(MASS)
library(cluster)
library(emdbook)
library(dplyr)
library(reshape2)
library(RobustGaSP)
library(mvtnorm)
library(Matrix)
library(stringr)

#Covariance Function - Currently Squared Exponential
cov_mat <- function(pred_pts,
                    design_pts,
                    ell_vec,
                    lambda,
                    alpha_vec = rep(1.9,dim(design_pts)[2]), nugget=0.){
  
  M = dim(design_pts)[1]
  N = dim(pred_pts)[1]
  
  if(dim(design_pts)[2]!=dim(design_pts)[2]) stop('You got dimension problems, big fella')
  
  K = dim(pred_pts)[2]
  
  sig_12 = exp(-as.matrix(1/ell_vec[1]*abs(outer(pred_pts[,1],design_pts[,1],FUN = "-")))^alpha_vec[1])
  if(K>1){
    for(k in 2:K){
      sig_12 = sig_12*exp(-as.matrix(1/ell_vec[k]*abs(outer(pred_pts[,k],design_pts[,k],FUN = "-")))^alpha_vec[k])
    }
  }
  
  if(dim(pred_pts)[1] == dim(design_pts)[1]){
    sig_12 = sig_12 + nugget*diag(M)
  }
  sig_12 = (1/lambda_eta)*sig_12
  return(sig_12)
  
}



setwd("/Users/Jake/Dropbox/Research/JETSCAPE/second_project_discrepancy/")
folder = "forSTAT/"
save_path <- "/Users/Jake/Dropbox/Research/JETSCAPE/second_project_discrepancy/output_1_15"
save_pics = FALSE


dataset <- "PbPb5020-cen-30-50"


(design <- read.table("latin_hc_design.txt",header = TRUE))
load('ranges.Rdata')

scaled_d <- matrix(0,dim(design)[1],dim(design)[2]) %>%
  as.data.frame()
for(j in 1:dim(design)[2]){
  scaled_d[,j] <- (design[,j] - ranges[[j]][1])/(ranges[[j]][2] - ranges[[j]][1])
}
colnames(scaled_d) <- c('A','B','C','D')

n_design <- dim(design)[1]

current_dset = dataset
(current_experiment = strsplit(dataset,'-')[[1]][1])

all_info <- read.table(paste0(folder,current_dset,".dat"))
(names(all_info) <- c("pT","RAA_exp","Stat_err","Sys_err",paste0("RAA_",as.character(1:n_design))))

#Separate output, change from wide to long
mod_dat <- dplyr::select(all_info, -c(RAA_exp,Stat_err,Sys_err)) %>%
  melt(id = "pT") %>%
  arrange(by = pT)%>%
  cbind(do.call("rbind", rep(list(scaled_d), 3)))%>%
  mutate(pT = (pT - min(pT))/(max(pT) - min(pT)))



###############
##Start labeling variables
################
n = length(all_info$pT)
x_E = with(all_info,(pT - min(pT))/(max(pT) - min(pT)))

m = dim(mod_dat)[1]
D_c = dplyr::select(mod_dat,pT,A,B,C,D)

y_E = all_info$RAA_exp
y_c = mod_dat$value

comp_mod <- rgasp(D_c,mod_dat$value,kernel_type = 'pow_exp', nugget.est = T)

theta_st_test = c(0.1, 0.2, 0.3, 0.4)
pred_mat = cbind(x_E,t(replicate(12,theta_st_test)))


ell_eta = 1/comp_mod@beta_hat
lambda_eta = 1/comp_mod@sigma2_hat

L_inv =  backsolve(r = comp_mod@L, x = diag(ncol(comp_mod@L)),
                   upper.tri = FALSE)

eta_sig_22_inv = lambda_eta*t(L_inv)%*%L_inv

test2 = cov_mat(D_c,
                D_c,
                ell_eta,
                lambda_eta,
                nugget = comp_mod@nugget)

View(test1%*%eta_sig_22_inv)

sig_12_test = cov_mat(pred_pts = pred_mat,
                      design_pts = D_c,
                      ell_vec = ell_eta,
                      lambda = lambda_eta)

calc_mu_st <- function(theta_st,
                       exp_control_inputs = x_E,
                       design_pts = D_c,
                       design_output = y_c,
                       ell_vec = ell_eta,
                       lambda = lambda_eta,
                       sigma_22_inv = eta_sig_22_inv,
                       intercept_hat = comp_mod@theta_hat){
  
  pred_pts = cbind(exp_control_inputs,t(replicate(12,theta_st)))
  
  sig_12 = cov_mat(pred_pts = pred_pts,
                   design_pts = design_pts,
                   ell_vec = ell_vec,
                   lambda = lambda)
  
  return(intercept_hat + sig_12%*%sigma_22_inv%*%(design_output - intercept_hat))
  
}

calc_sig_st <- function(theta_st,
                        exp_control_inputs = x_E,
                        design_pts = D_c,
                        design_output = y_c,
                        ell_vec = ell_eta,
                        lambda = lambda_eta,
                        sigma_22_inv = eta_sig_22_inv,
                        intercept_hat = comp_mod@theta_hat){
  
  pred_pts = cbind(exp_control_inputs,t(replicate(length(x_E),theta_st)))
  
  sig_11 = cov_mat(pred_pts = pred_pts,
                   design_pts = pred_pts,
                   ell_vec = ell_vec,
                   lambda = lambda,
                   nugget = comp_mod@nugget)
  
  sig_12 = cov_mat(pred_pts = pred_pts,
                   design_pts = design_pts,
                   ell_vec = ell_vec,
                   lambda = lambda)
  sig_21 = t(sig_12)
  
  return(sig_11 - sig_12%*%sigma_22_inv%*%sig_21)
  
}

comp_pred = RobustGaSP::predict(comp_mod,pred_mat)

comp_pred$mean - calc_mu_st(theta_st_test)

comp_pred$sd^2 - diag(calc_sig_st(theta_st_test))

#Gamma prior for lambda_E
a <- 2; b <- 5

#Gamma prior for lambda_delta
c <- 2; d <- 5;

jet_gibbs <- function(niters, burnin = 0.3, ell_kap = 0.5, theta_kap = 1E-1){
  
  #Set up output contains  
  burn_iters = floor(niters*burnin)
  total_iters = niters + burn_iters
  
  eta_mat <- matrix(0,total_iters, length(x_E))
  colnames(eta_mat) = paste0('eta',1:length(x_E))
  
  delta_mat <- matrix(0,total_iters, length(x_E))
  colnames(delta_mat) = paste0('delta',1:length(x_E))
  
  lambda_E_mat <- matrix(0,total_iters,1)
  colnames(lambda_E_mat) = paste0('lambda_E')
  
  lambda_delta_mat <- matrix(0,total_iters,1)
  colnames(lambda_delta_mat) <- paste0('lambda_delta')
  
  ell_delta_mat <- matrix(0,total_iters,1)
  colnames(ell_delta_mat) <- paste0('ell_delta')
  
  theta_mat <- matrix(0,total_iters, dim(scaled_d)[2])
  colnames(theta_mat) <- colnames(scaled_d)
  
  
  #More MCMC setup
  accept_ratio_theta = 0
  accept_ratio_ell = 0
  ob_ratio_theta = 0
  time_start <- proc.time()
  
  
  #Draw from priors
  theta_mat[1,] <- runif(dim(theta_mat)[2])
  lambda_E_mat[1] <- rgamma(1,shape = a, rate = b)
  lambda_delta_mat[1] <- rgamma(1,shape = a, rate = b)
  ell_delta_mat[1] <- runif(1)
  
  mu_st <- calc_mu_st(theta_mat[1,])
  sig_st <- calc_sig_st(theta_mat[1,])
  sig_st <- sig_st%*%t(sig_st)/2 #non-symmetric rounding issue
  eta_mat[1,] <- rmvnorm(1,mean = calc_mu_st(theta_mat[1,]),
                         sigma = sig_st) 
  
  
  
  delta_mat[1,] <- rmvnorm(1,mean = rep(0,length(x_E)),
                           sigma = cov_mat(pred_pts = matrix(x_E),
                                           design_pts = matrix(x_E),
                                           ell_vec = ell_delta_mat[1],
                                           lambda = lambda_delta_mat[1],
                                           nugget = 1E-12))
  
  for(i in 2:total_iters){
    
    #Print elapsed time every ten percent of iterations
    if(!i%%(niters*0.1)){
      flush.console()
      cat("\n i = ", i, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
      cat("\n accept_theta: ", accept_ratio_theta/i)
      cat("\n accept_ell_delta: ", accept_ratio_ell/i)
    }
    
    #################
    ###Eta Update###
    ##################
    sig_st <- calc_sig_st(theta_mat[i-1,])
    mu_st <- calc_mu_st(theta_mat[i-1,])
    V_eta <- solve(lambda_E_mat[i-1]*diag(n) + solve(sig_st))
    mean_eta <- V_eta%*%(lambda_E_mat[i-1]*(matrix(y_E) - delta_mat[i-1,]) +
                           solve(sig_st)%*%mu_st)
    
    V_eta <- V_eta%*%t(V_eta)/2 #symmetric rounding issue
    eta_mat[i,] <- rmvnorm(1,mean_eta, V_eta) 
    
    
    ##################
    ###Delta Update###
    ###################
    V_delta <- solve(lambda_E_mat[i-1]*diag(n) + solve(cov_mat(pred_pts = matrix(x_E),
                                                      design_pts = matrix(x_E),
                                                      ell_vec = ell_delta_mat[i-1],
                                                      lambda = lambda_delta_mat[i-1],
                                                      nugget = 1E-12)))
    mean_delta <- V_delta%*%(lambda_E_mat[i-1]*(matrix(y_E) - eta_mat[i,]))
    delta_mat[i,] <- rmvnorm(1,mean_delta,V_delta)#This doesn't need rounding?
    
    #####################
    ###lambda_E Update######
    ########################
    lambda_E_mat[i] <- rgamma(1,a + n/2,
                              b + 0.5*t(matrix(y_E) - eta_mat[i,] - delta_mat[i,])%*%
                                (matrix(y_E) - eta_mat[i,] - delta_mat[i,]))
    
    #####################
    ###lambda_delta Update######
    ########################    
    sig_delta <- cov_mat(pred_pts = matrix(x_E),
                         design_pts = matrix(x_E),
                         ell_vec = ell_delta_mat[i-1],
                         lambda = 1,
                         nugget = 1E-12)
    lambda_delta_mat[i] <- rgamma(1,c + n/2,
                                  b + 0.5*t(matrix(delta_mat[i,]))%*%solve(sig_delta)%*%
                                    (matrix(delta_mat[i,])))
   
    
    ######################
    ####Theta update#####
    #####################
    #Currently proposing together, could do individual updates
    
    theta_prop <- rnorm(dim(theta_mat)[2],theta_mat[i-1,],theta_kap)
    auto_reject_theta = FALSE
    if(sum(theta_prop<0) + sum(theta_prop>1)){
      #print('oops')
      auto_reject_theta = TRUE
      ob_ratio_theta = ob_ratio_theta + 1
    }
    
    mu_st_prop <- calc_mu_st(theta_prop)
    sig_st_prop <- calc_sig_st(theta_prop)
    
    mu_st_cur <- calc_mu_st(theta_mat[i-1,])
    sig_st_cur <- calc_sig_st(theta_mat[i-1,])
    
    post_prop_theta <- -0.5*determinant(sig_st_prop,logarithm = TRUE)$modulus +
      -0.5*t(matrix(eta_mat[i,]) - mu_st_prop)%*%solve(sig_st_prop)%*%
               (matrix(eta_mat[i,]) - mu_st_prop)
    
    post_cur_theta <- -0.5*determinant(sig_st_cur,logarithm = TRUE)$modulus +
      -0.5*t(matrix(eta_mat[i,]) - mu_st_cur)%*%solve(sig_st_cur)%*%
      (matrix(eta_mat[i,]) - mu_st_cur)
    
    ratio_theta <- post_prop_theta - post_cur_theta
    
    if(rexp(1) > -ratio_theta & !auto_reject_theta){
      theta_mat[i,] <- theta_prop
      accept_ratio_theta = accept_ratio_theta + 1
    }else{
      theta_mat[i,] <- theta_mat[i-1,]
    }
    
    ############################
    #####Ell_delta############
    #####################
    ell_delta_Z <- rnorm(1)
    ell_delta_prop <- ell_delta_mat[i-1]*exp(ell_kap*ell_delta_Z)
    
    sig_delta_prop <- cov_mat(pred_pts = matrix(x_E),
                              design_pts = matrix(x_E),
                              ell_vec = ell_delta_prop,
                              lambda = 1,
                              nugget = 1E-12)
    sig_delta_cur <- cov_mat(pred_pts = matrix(x_E),
                             design_pts = matrix(x_E),
                             ell_vec = ell_delta_mat[i-1],
                             lambda = 1,
                             nugget = 1E-12)
    
    prop_post_ell_delta <- -0.5*determinant(sig_delta_prop,logarithm = TRUE)$modulus +
      -(lambda_delta_mat[i]/2)*t(matrix(delta_mat[i,]))%*%solve(sig_delta_prop)%*%matrix(delta_mat[i,])
    
    cur_post_ell_delta <- 0.5*determinant(sig_delta_cur,logarithm = TRUE)$modulus + 
      -(-lambda_delta_mat[i]/2)*t(matrix(delta_mat[i,]))%*%solve(sig_delta_prop)%*%matrix(delta_mat[i,])
    

    ratio <- prop_post_ell_delta - cur_post_ell_delta
    
    if(rexp(1)> (-ratio - ell_delta_Z*ell_kap)){
      ell_delta_mat[i] <- ell_delta_prop
      accept_ratio_ell <- accept_ratio_ell + 1 
    }else{
      ell_delta_mat[i] <- ell_delta_mat[i-1]
    }
    
  }
  
  
  print(paste0('Final accept_ratio_theta: ', round(accept_ratio_theta/total_iters,3)))
  print(paste0('Final theta OB rate: ', round(ob_ratio_theta/total_iters,3)))
  print(paste0('Final accept_ratio_ell: ', round(accept_ratio_ell/total_iters,3)))
  print(paste0('Total minutes: ', round(((proc.time() - time_start)[3])/60,3)))
  out_mat <- cbind(eta_mat, delta_mat, lambda_E_mat, lambda_delta_mat, ell_delta_mat, theta_mat)
  return(out_mat)
}

Rprof('halp.txt')
test_gibbs <- jet_gibbs(100,theta_kap = 5E-2)
Rprof()

summaryRprof('halp.txt')


