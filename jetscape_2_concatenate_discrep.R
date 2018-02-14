library(MASS)
library(cluster)
library(emdbook)
library(dplyr)
library(reshape2)
library(RobustGaSP)
library(mvtnorm)
library(Matrix)
library(stringr)

###NOTE
#CURRENTLY USING WHITENED PCA, 

setwd("/Users/Jake/Dropbox/Research/JETSCAPE/second_project_discrepancy/")
folder = "all_data/"
save_path <- "/Users/Jake/Dropbox/Research/JETSCAPE/second_project_discrepancy/output_2_4/jonah_method/"
save_pics = FALSE


all_dsets <- c(
  "AuAu200-cen-00-10"
  ,"AuAu200-cen-40-50"
  ,"PbPb2760-cen-00-05"
  ,"PbPb2760-cen-30-40"
  # ,"PbPb5020-cen-00-10"
  ,"PbPb5020-cen-30-50"
)

design_file = paste0(folder,'all_design')
range_file = paste0(folder,'ranges.Rdata')

systems_to_calibrate = c(#"AuAu200"
  #, "PbPb2760"
  "PbPb5020"
)
if(length(systems_to_calibrate)==3){
  (hist_main = "All Datasets Simultaneous Calibration")
}else{
  (hist_main = paste(paste0(systems_to_calibrate,collapse = " & "),"Calibration"))
}

(design <- read.table(design_file,header = TRUE))
load(range_file)
scaled_d <- matrix(0,dim(design)[1],dim(design)[2]) %>%
  as.data.frame()
for(j in 1:dim(design)[2]){
  scaled_d[,j] <- (design[,j] - ranges[[j]][1])/(ranges[[j]][2] - ranges[[j]][1])
}
colnames(scaled_d) <- c('Lambda','alpha_s')

# scaled_d <- lapply(design,function(x){
#   (x-min(x))/(max(x) - min(x))
# }) %>% 
#   do.call(cbind,.) %>%
#   as.data.frame()

n_design <- dim(design)[1]


#Covariance Function - Currently Squared Exponential
cov_mat <- function(x,ell,lambda,alpha = 2, nugget=0.){
  #inds <- 1:length(x)
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(x)/ell)^alpha) + 
    nugget*diag(length(x))
  
  return(out_mat)
}
(dsets_to_use = which(str_detect(all_dsets, 
                                 paste(systems_to_calibrate,collapse = "|"))
))
covs <- au <- vector('list',length(dsets_to_use))

first_concat = TRUE
write_csvs = FALSE
for(i in 1:length(dsets_to_use)){
  dset = dsets_to_use[i]
  (current_dset = all_dsets[dset])
  (current_experiment = strsplit(current_dset,'-')[[1]][1])
  
  au[[i]] <- read.table(paste0(folder,current_dset,".dat"),header = TRUE)
  (names(au[[i]]) <- c("pT","RAA_exp","Stat_err","Sys_err",paste0("RAA_",as.character(1:n_design))))
  
  if(write_csvs){
    #Write csv of experimental data for easily manipulation later
    write.csv(au[[i]][,1:4],file=paste0(folder,current_dset,"_exp.csv"),row.names = FALSE)
  }
  #Separate output, change from wide to long
  mod_dat <- dplyr::select(au[[i]], -c(RAA_exp,Stat_err,Sys_err)) %>%
    melt(id = "pT") %>%
    arrange(by = pT)
  
  (exp_dat <- au[[i]]$RAA_exp)
  stat_cov <- diag(au[[i]]$Stat_err^2)
  
  #Not sure this is right - scale or don't scale pT?
  pT_scaled <- (au[[i]]$pT - min(au[[i]]$pT))/(max(au[[i]]$pT) - min(au[[i]]$pT))
  err_mod <- rgasp(pT_scaled,au[[i]]$RAA_exp, kernel_type = "pow_exp",alpha = 1.9)
  1/err_mod@beta_hat
  
  sys_cov <- outer(au[[i]]$Sys_err,au[[i]]$Sys_err)*
    cov_mat(x = pT_scaled,ell = 1/err_mod@beta_hat, lambda = 1,alpha = 1.9)
  
  covs[[i]] <- stat_cov + sys_cov
  
  #Take the experimental data, train a GP to it
  #Find the MLE length parameter when the standard deviations are fixed
  
  cur_concat <- dcast(mod_dat,variable~pT)
  colnames(cur_concat) = c('design',paste0(current_experiment,'_',colnames(cur_concat)[-1]))
  
  if(write_csvs){
    #Write csv of each dataset for easily manipulation later
    write.csv(cur_concat[,-1],file=paste0(folder,current_dset,".csv"),row.names = FALSE)
  }    
  if(first_concat){
    first_concat = FALSE
    all_mod_dat = cur_concat
    all_exp_dat = exp_dat
  }else{
    all_mod_dat = cbind(all_mod_dat,cur_concat[,-1])
    all_exp_dat = c(all_exp_dat,exp_dat)
  }
}
length(all_exp_dat)

block_covs <- as.matrix(bdiag(covs))
#Check normality?

#Get experimental vector - all_exp_dat
#Build block covariance matrices for experimental error - block_covs

#Rotate with PCA
Y <- all_mod_dat[,-1]
sd_vec = apply(Y,2,sd)
#sd_vec = rep(1,dim(Y)[2])
y_means = apply(Y,2,mean)
Y_final <- as.matrix(sweep(Y,2,y_means)) %*% diag(1/sd_vec)
Y_svd <- svd(Y_final)

#Test for Normality
#None are <0.05 even without accounting for multiple testing
for(j in 1:dim(Y_final)[2]){
  print(shapiro.test(Y_final[,j])$p.value)
}

#How many PCs?
eigs <- Y_svd$d^2
(V_q <- cumsum(eigs)/sum(eigs))

if(save_pics) pdf(paste0(save_path,'var_explained.pdf'))
plot(V_q[1:6],type = 'o',
     xlab = "Number of PC r",
     main = "Fraction of Variance Explained",
     ylab = expression(F[r]),
     cex.lab = 2,
     cex.axis = 1.3,
     cex.main = 1.4,
     mgp = c(2.3,1,0),
     pch = 19)
if(save_pics) dev.off()

q <- 5
V = Y_svd$v[,1:q]
S = diag(Y_svd$d[1:q])

m = dim(Y_final)[1]
Z <-as.matrix(Y_final)%*%V %>%
  as.data.frame()

##Extra variation
V_b = Y_svd$v[,(q+1):dim(Y_final)[2]]
S_b = diag(Y_svd$d[(q+1):dim(Y_final)[2]])



#############
##Prediction
############

final_mod <- lapply(Z,rgasp,design = scaled_d, nugget.est=TRUE)#, nugget = 1E-5)
#final_mod <- lapply(Z,km,design=scaled_d,formula=~1,nugget.estim=TRUE)

#(sigma2hats <- lapply(final_mod,function(x)x@sigma2_hat) %>% unlist())




#Calibration


#Rotate the experimental data with same PCA values
Y_exp_final <- t(all_exp_dat - y_means)%*%diag(1/sd_vec)
Y_cov_final <- diag(1/sd_vec)%*%block_covs %*%diag(1/sd_vec)

(Z_exp <- Y_exp_final%*%V)

(Z_cov <-t(V)%*%Y_cov_final%*%V)

#Extra variation
cov_extra_phys_cal = 1/(m)*diag(sd_vec)%*%V_b%*%S_b^2%*%t(V_b)%*%diag(sd_vec)

#cov_extra_pca_cal = t(V_b)%*%Y_cov_final%*%V_b

lmu_ell = log(0.3)
lsd_ell = 0.5

a_lam <- 1
b_lam <- 1
mh_cal <- function(niters = 1E4,burnin_prop = 0.3,
                   t_kap = 0.1,
                   proposal_cor_mat =diag(dim(design)[2]),
                   delta_kap = 1E-1,
                   kap_lam = 0.3,
                   kap_ell = 0.3){
  
  #Do various sampler setups
  time_start <- proc.time()
  
  burnin <- burnin_prop*niters
  
  ##Set up containers
  t_out <- matrix(0,niters,dim(design)[2])
  delta_out <- matrix(0,niters,length(all_exp_dat))
  lambda_out <- ell_out <- numeric(niters)
  t_ratio <- ob_ratio <- delta_ratio <- lambda_ratio <- ell_ratio <- 0
  
  
  lambda_cur <- b_lam/a_lam
  ell_cur <- exp(lmu_ell + lsd_ell/2)
  delta_cur <- rmvnorm(1,mean=rep(0,length(all_exp_dat)),
                       sigma = cov_mat(x=pT_scaled,
                                       ell = ell_cur,
                                       lambda = lambda_cur,
                                       alpha = 1.9))
  #draw from priors
  t_cur <- runif(dim(design)[2])
  
  #Get current values
  (X_cur <- t(as.matrix(t_cur)))
  colnames(X_cur) <- colnames(scaled_d)
  cur_mod <- lapply(final_mod, RobustGaSP::predict, testing_input = X_cur)

  track_ell <- numeric(niters+burnin)
  for(i in 1:(niters + burnin)){
    #Get the time every ten percent of iterations
    if(!i%%(niters*0.1)){
      flush.console()
      cat("\r i = ", i, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
    }
    track_ell[i] <- ell_cur
    
    #####################
    ### Theta Update ####
    #####################
    #Propose new t_st
    (t_st <- mvtnorm::rmvnorm(1,mean=t_cur,sigma = t_kap*proposal_cor_mat))
    
    #Check to make sure parameter is within [0,1]
    auto_reject = FALSE
    if(sum(t_st<0) + sum(t_st>1)){
      auto_reject = TRUE
      ob_ratio = ob_ratio + 1
    }
    
    #(X_st_sqrt <- t(as.matrix(t_st)))
    X_st_sqrt <- t_st
    (X_st <- X_st_sqrt)#^2
    colnames(X_st) <- colnames(scaled_d)
    #Get the predictive mean/sd for t_st
    st_mod <- lapply(final_mod, RobustGaSP::predict, testing_input = X_st)
    #st_mod <- lapply(final_mod,DiceKriging::predict, newdata=data.frame(X_st),type ='UK')
    
    round(z_st_means <- do.call(c,lapply(st_mod,function(x){x$mean})),3)
    round(z_st_sd <- do.call(c,lapply(st_mod,function(x){x$sd})),3)
    
    round(z_st_var <- z_st_sd^2,3)
    
    round(z_cur_means <- do.call(c,lapply(cur_mod,function(x){x$mean})), 3)
    round(z_cur_sd <- do.call(c,lapply(cur_mod,function(x){x$sd})), 3)
    
    round(z_cur_var <- z_cur_sd^2,3)
    
    (y_st_means <- t(z_st_means)%*%t(V)%*%diag(sd_vec) + y_means)
    (y_st_cov <- diag(sd_vec)%*%V%*%diag(z_st_var)%*%t(V)%*%diag(sd_vec))
    
    (y_cur_means <- t(z_cur_means)%*%t(V)%*%diag(sd_vec) + y_means)
    (y_cur_cov <- diag(sd_vec)%*%V%*%diag(z_cur_var)%*%t(V)%*%diag(sd_vec))
    
    
    #Get likelihoods
    
    (like_cur <- mvtnorm::dmvnorm(x = as.numeric(all_exp_dat),
                                  mean = as.numeric(y_cur_means) + delta_cur,
                                  sigma = block_covs + y_st_cov + cov_extra_phys_cal,
                                  log = TRUE))
    (like_st <- mvtnorm::dmvnorm(x = as.numeric(all_exp_dat),
                                 mean = as.numeric(y_st_means) + delta_cur,
                                 sigma = block_covs + y_cur_cov + cov_extra_phys_cal,
                                 log = TRUE))
    #Includes constant prior 
    (ratio <- sum(like_st) - sum(like_cur)) #+ sum(log(2*X_st)) - sum(log(2*X_cur))
    
    
    if(rexp(1) > -ratio & !auto_reject){
      cur_mod <- st_mod
      t_cur <- t_st
      t_ratio <- t_ratio + 1
      y_cur_means <- y_st_means
    }
    
    
    ###################
    ### Delta Update ###
    #####################
    
    #Joint update?
    (delta_prop <- rnorm(length(all_exp_dat),delta_cur,delta_kap))
    
    (like_cur <- mvtnorm::dmvnorm(x = as.numeric(all_exp_dat),
                                  mean = as.numeric(y_cur_means) + delta_cur,
                                  sigma = block_covs + y_st_cov + cov_extra_phys_cal,
                                  log = TRUE))
    
    (like_prop <- mvtnorm::dmvnorm(x = as.numeric(all_exp_dat),
                                   mean = as.numeric(y_cur_means) + delta_prop,
                                   sigma = block_covs + y_st_cov + cov_extra_phys_cal,
                                   log = TRUE))
    
    (prior_cur <- mvtnorm::dmvnorm(x = delta_cur,
                                  mean = rep(0,length(all_exp_dat)),
                                  ######### MAKE x_E ###########
                                  sigma = cov_mat(x = pT_scaled,
                                                  ell = ell_cur,
                                                  lambda = lambda_cur,
                                                  alpha = 1.9),
                                  log = TRUE))
    
    (prior_prop <- mvtnorm::dmvnorm(x = delta_prop,
                                   mean = rep(0,length(all_exp_dat)),
                                   ######### MAKE x_E ###########
                                   sigma = cov_mat(x = pT_scaled,
                                                   ell = ell_cur,
                                                   lambda = lambda_cur,
                                                   alpha = 1.9),
                                   log = TRUE))
    
    (r_delta <- like_prop + prior_prop - like_cur - prior_cur)
    
    if(rexp(1)>-r_delta){
      delta_cur <- delta_prop
      delta_ratio <- delta_ratio + 1
    }
    
    ###########################
    #### Delta lambda update #####
    ##############################
    
    Z_lam <- rnorm(1)
    lambda_prop <- lambda_cur*exp(kap_lam*Z_lam)
    
    (like_cur <- mvtnorm::dmvnorm(x = delta_cur,
                                 mean = rep(0,length(all_exp_dat)),
                                 ######### MAKE x_E ###########
                                 sigma = cov_mat(x = pT_scaled,
                                                 ell = ell_cur,
                                                 lambda = lambda_cur,
                                                 alpha = 1.9),
                                 log = TRUE))
    (like_prop <- mvtnorm::dmvnorm(x = delta_prop,
                                 mean = rep(0,length(all_exp_dat)),
                                 ######### MAKE x_E ###########
                                 sigma = cov_mat(x = pT_scaled,
                                                 ell = ell_cur,
                                                 lambda = lambda_prop,
                                                 alpha = 1.9),
                                 log = TRUE))
    
    (prior_cur <- dgamma(lambda_cur, a_lam, b_lam,log = TRUE))
    (prior_prop <- dgamma(lambda_prop, a_lam, b_lam, log = TRUE))
    
    (r_lam <- like_prop + prior_prop - like_cur - prior_cur + kap_lam*Z_lam)
    
    
    if(rexp(1)>-r_lam){
      lambda_cur <- lambda_prop
      lambda_ratio <- lambda_ratio + 1
    }
    
    
    ##########################
    #### Ell Update ##########
    ##########################
    
    Z_ell <- rnorm(1)
    (ell_prop <- ell_cur*exp(kap_ell*Z_ell))
    
    (like_cur <- mvtnorm::dmvnorm(x = delta_cur,
                                 mean = rep(0,length(all_exp_dat)),
                                 ######### MAKE x_E ###########
                                 sigma = cov_mat(x = pT_scaled,
                                                 ell = ell_cur,
                                                 lambda = lambda_cur,
                                                 alpha = 1.9),
                                 log = TRUE))
    (like_prop <- mvtnorm::dmvnorm(x = delta_cur,
                                  mean = rep(0,length(all_exp_dat)),
                                  ######### MAKE x_E ###########
                                  sigma = cov_mat(x = pT_scaled,
                                                  ell = ell_prop,
                                                  lambda = lambda_cur,
                                                  alpha = 1.9),
                                  log = TRUE))
    
    (prior_cur <- dlnorm(ell_cur,meanlog = lmu_ell,sdlog = lsd_ell,log= TRUE))
    (prior_prop <- dlnorm(ell_prop,meanlog = lmu_ell,sdlog = lsd_ell,log = TRUE))
    
    
    (r_ell <- like_prop + prior_prop - like_cur - prior_cur + kap_ell*Z_ell)
    
    
    if(rexp(1)>-r_ell){
      ell_cur <- ell_prop
      ell_ratio <- ell_ratio + 1
    }
    
    
    if(i > burnin){
      t_out[i-burnin,] <- t_cur
      delta_out[i-burnin,] <- delta_cur
      lambda_out[i-burnin] <- lambda_cur
      ell_out[i-burnin] <- ell_cur
    }
    


  }#i loop
  
  print("t_ratio is")
  print(t_ratio/(niters + burnin))
  print("ob ratio is")
  print(ob_ratio/(niters + burnin))
  
  print('delta_ratio is')
  print(delta_ratio/(niters + burnin))
 
  print('lambda_ratio is')
  print(lambda_ratio/(niters + burnin))
  
  print('ell_ratio is')
  print(ell_ratio/(niters + burnin))
  
  
  return(list(params = t_out,delta=delta_out,lambda=lambda_out,ell=ell_out))
  
}


res <- mh_cal(niters = 1E4,t_kap = 1E-4,delta_kap = 5E-3)#,
#proposal_cor_mat = matrix(c(1,-.8,-.8,1),ncol = 2))
#save(res,file = paste0(save_path,'res_old.Rdata'))

param_plot <- matrix(0,dim(res$params)[1],dim(res$params)[2])
#param_plot <- matrix(0,dim(res$params)[1],dim(res$params)[2])
for(j in 1:dim(param_plot)[2]){
  param_plot[,j] <- res$params[,j]*(ranges[[j]][2] - ranges[[j]][1]) + ranges[[j]][1]
  plot(param_plot[,j],type = 'l',main = colnames(design)[j])
}
#param_plot_sqrt <- sqrt(param_plot)

plot(res$lambda,type = 'l')
plot(res$ell,type = 'l')
for(j in 1:length(all_exp_dat)){
  plot(res$delta[,j],type = 'l',main = paste('delta',j))
}

param_names = colnames(design)
for(k in 1:(dim(design)[2]-1)){
  for(l in (k+1):dim(design)[2]){
    f1 <- kde2d(param_plot[,k], (param_plot[,l]), n = 100,
                lims = c(ranges[[k]],ranges[[l]]))
    if(save_pics) pdf(paste0(save_path,'heatmap_',param_names[k],'_',param_names[l],'.pdf'))
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    image(f1,
          xlab = param_names[k],
          ylab = "",
          main = hist_main,
          cex.lab = 2,
          cex.axis = 1.3,
          cex.main = 1.4,
          mgp = c(3,1,0)#,
          # xlim = ranges[[1]],
          # ylim = ranges[[2]]
    )
    title(ylab = param_names[l], 
          mgp = c(2.1,1,0),
          cex.lab = 2,
          cex.axis = 1.3)
    
    perc_lvl = c(.6,.75,.9)
    HPDregionplot(param_plot[,c(k,l)], prob = perc_lvl,
                  col=c("black"), lty = c(1,5,3), add=TRUE)
    
    
    # legend('topright',paste0(perc_lvl*100,"%"),title = "Highest Density Kernel Estimate",
    #        lty = c(1,5,3),
    #        bty="n")
    if(save_pics) dev.off()
  }
}




plot(param_plot[,1],type = 'l', ylab = expression(alpha), xlab = 'Iteration',
     main = expression(alpha~'Trace Plot'))
plot(param_plot[,2],type = 'l', ylab = expression(beta), xlab = 'Iteration',
     main = expression(~beta~'Trace Plot'))
# plot(param_plot[,3],type = 'l',ylab = 'Gamma')
# plot(param_plot[,4],type = 'l',ylab = 'Delta')
# 
# plot(res$params[,1],type = 'l', ylab = expression(sqrt(Lambda^jet)), xlab = 'Iteration',
#      main = expression(sqrt(Lambda^jet)~'Trace Plot'))
# plot(res$params[,2],type = 'l', ylab = expression(sqrt(alpha[s]^med)), xlab = 'Iteration',
#      main = expression(sqrt(alpha[s]^med)~'Trace Plot'))

#write.table(param_plot,file = paste0(save_path,'post_draws.txt'),
#           row.names = FALSE)

#Now there are 12 bivariate normal likelihoods, rather than 12 univariate normal likelihoods
#Everything else is the same


panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y,xlim = c(0,0.35),  ...)
}

panel.image <- function(x,y,...){
  f1 <- kde2d(x, y, n = 100)
  image(f1,add = TRUE)
  #points(x,y, col = rgb(1,0,0,.1))
}

if(save_pics){ pdf(paste0(save_path,
                          paste0(systems_to_calibrate,collapse = "_"),
                          ".pdf"))}
pairs(param_plot, 
      panel = panel.image,
      diag.panel = panel.hist,
      #pch = 19,
      #cex = .3,
      col = rgb(1,0,0,.1),
      labels = c('A','B','C','D'),
      #expression(alpha),
      #expression(beta),
      #expression(gamma),
      #expression(delta)),
      upper.panel = NULL,
      cex.lab = 1.5,
      cex.axis = 1.3,
      cex.main = 1.4,
      mgp = c(2.3,1,0),
      las = 2,
      main = hist_main)
if(save_pics) dev.off()

effectiveSize(param_plot)






#####
#Extra visualization/analysis
#####

plot_lines <- function(dat,pT_col = 'pT',exp_col = 'RAA_exp',
                       err_col_stat = 'exp_err_stat', 
                       err_col_sys = 'exp_err_sys',
                       col_start = 5,
                       alpha_val = 0.01,
                       ...){
  
  plot(dat[,pT_col],dat[,col_start],
       cex = 0.5,
       col = rgb(1,0,0,alpha_val),type = 'l',
       # ylim = c(min(dat[,exp_col] - 2*dat[,err_col_sys],
       #              dat[,exp_col] - 2*dat[,err_col_stat]),
       #          max(dat[,exp_col] + 2*dat[,err_col_sys],
       #              dat[,exp_col] + 2*dat[,err_col_stat])),
       ylim = c(0,1),
       xlab = expression(p[T]),
       ylab = expression(R[AA]),
       cex.main = 1.3,
       cex.lab = 2,
       mgp = c(2.3,1,0),
       ...)
  
  x_range = max(dat[,pT_col]) - min(dat[,pT_col])
  arrow_offset = x_range/130
  
  for(i in (col_start+1):dim(dat)[2]){
    lines(dat[,pT_col],dat[,i],cex = 0.5,
          col = rgb(1,0,0,alpha_val))
  }
  
  points(dat[,pT_col],dat[,exp_col],
         pch = 19)
  #ylim = c(0,1))
  arrows(dat[,pT_col]-arrow_offset, dat[,exp_col] - 1.96*dat[,err_col_stat],
         dat[,pT_col]-arrow_offset, dat[,exp_col] + 1.96*dat[,err_col_stat],
         length=0.05, angle=90, code=3,col = 'blue')
  arrows(dat[,pT_col]+arrow_offset, dat[,exp_col] - 1.96*dat[,err_col_sys],
         dat[,pT_col]+arrow_offset, dat[,exp_col] + 1.96*dat[,err_col_sys],
         length=0.05, angle=90, code=3, col = 'black')
  legend('bottomright',c('Systematic Errors',
                         'Statistical Errors'),
         col = c('blue','black'), lwd =1)
  
}

##First, go from concatenated to many datasets, with p_T in rows

#1. Cut into 6 original datasets

plot_all_dsets <- function(all_dset_mat,
                           orig_dset_list = au,
                           dset_strings = all_dsets[dsets_to_use],
                           title_end = "",
                           alpha_val = 0.01,
                           save_pics = FALSE){
  cur_col <- 1
  for(j in 1:length(orig_dset_list)){
    (num_pT <- dim(orig_dset_list[[j]])[1])
    cur_dset <- t(all_dset_mat[,cur_col:(cur_col + num_pT - 1)]) %>%
      t() %>%
      cbind('pT' = orig_dset_list[[j]]$pT,
            'RAA_exp' = orig_dset_list[[j]]$RAA_exp,
            'exp_err_stat' = orig_dset_list[[j]]$Stat_err,
            'exp_err_sys' =  orig_dset_list[[j]]$Sys_err,
            .)
    (cur_col = cur_col + num_pT)
    
    
    #Plot 'em
    if(save_pics) pdf(paste0(save_path,dset_strings[j],'_',title_end,'.pdf'))
    plot_lines(cur_dset,
               main = paste(dset_strings[j],title_end),
               alpha_val = alpha_val)
    if(save_pics) dev.off()
  }
}


#Predict PCA vals
plot_draws <- function(draws,
                       train_mod = final_mod,
                       rot_mat = V,
                       scale_mat = Y,
                       num_samples = 1E3,
                       title_end = "",
                       alpha_val = 0.01,
                       save_pics = FALSE){
  
  colnames(draws) <- colnames(scaled_d)
  post_pred_mod <-lapply(train_mod, RobustGaSP::predict,
                         testing_input = draws)
  
  # post_pred_mod <-lapply(train_mod, DiceKriging::predict,
  #                        newdata = draws,type = 'UK')
  
  pred_means <- lapply(post_pred_mod,function(x)x$mean) %>%
    do.call(cbind,.)
  
  #Rotate
  Y_pred_scaled <- pred_means %*%t(rot_mat)
  
  #Scale
  Y_pred <- sweep(Y_pred_scaled,2,sd_vec,FUN = "*") %>%
    sweep(2,apply(scale_mat,2,mean),FUN = '+') 
  #Y_pred <- sweep(Y_pred_scaled,2,apply(scale_mat,2,mean),FUN = '+') 
  
  if(dim(Y_pred)[1]<num_samples){
    stop('Num samples bigger than number of draws')
  }
  
  Y_pred <- Y_pred[sample(1:dim(Y_pred)[1],num_samples),]
  
  plot_all_dsets(Y_pred,title_end = title_end,
                 alpha_val = alpha_val,
                 save_pics = save_pics)
}

plot_draws(res$params,title_end = 'Posterior',
           alpha_val = 0.03,
           rot_mat = V,
           save_pics = FALSE)
plot_draws(scaled_d[-80,],title_end = "Design",
           alpha_val = 1,
           rot_mat = V,
           num_samples = dim(design)[1]-1,
           save_pics = FALSE)


###########
##Validation Plot
##########

holdout = 35
train_d = scaled_d[-holdout,]
test_d = scaled_d[holdout,]

train_Y = Y_final[-holdout,]
test_Y = Y[holdout,]

V_train = svd(train_Y)$v[,1:q]

train_Z = as.matrix(train_Y)%*%V_train
test_Z = as.matrix(test_Y)%*%V_train

train_mod <- lapply(as.data.frame(train_Z),rgasp,design = train_d,nugget = 1E-5)
test_mod <- lapply(train_mod,RobustGaSP::predict,testing_input = test_d)

pred_Z <- lapply(test_mod,function(x)x$mean) %>%
  do.call(cbind,.) 
pred_err_Z <- lapply(test_mod,function(x)x$sd^2) %>%
  do.call(cbind,.) %>%
  as.numeric()%>%
  diag() 

pred_Y <- pred_Z %*%t(V_train) %>%
  sweep(2,sd_vec,FUN = '*') %>%
  sweep(2,apply(Y,2,mean),FUN = '+')

#This gives negative values
pred_err_Y <- V_train%*%pred_err_Z %*%t(V_train) %>%
  sweep(2,sd_vec^2,FUN = '*') %>%
  diag()

#This is what I argue for in Meeting 6/29 IT'S THE SAME!!!!!
# pred_err_Y2 = diag(pred_err_Z)%*%t(V^2) %>%
#   sweep(2,apply(Y,2,sd)^2,FUN = '*')

plot(as.numeric(test_Y),as.numeric(pred_Y),pch = 19,cex = .5,
     xlab = 'Holdout Values',
     ylab = 'Predicted Values',
     main = 'Emulator Prediction',
     cex.lab = 2,
     cex.axis = 1.3,
     cex.main = 1.4,
     mgp = c(2.4,1,0),
     ylim = c(min(as.numeric(pred_Y) - 2*sqrt(pred_err_Y)),
              max(as.numeric(pred_Y) + 2*sqrt(pred_err_Y))))
arrows(as.numeric(test_Y), as.numeric(pred_Y) - 2*sqrt(pred_err_Y),
       as.numeric(test_Y),as.numeric(pred_Y) + 2*sqrt(pred_err_Y),
       length=0, angle=90, code=3)


abline(a = 0,b = 1)








samples <- param_plot

##Not needed
# perc_lvl = c(.6,.75,.9)
# 
# fit <- cov.mve(samples, quantile.used = nrow(samples) * perc_lvl,
#                nsamp = 1E4)
# points_in_ellipse <- samples[fit$best, ]
# ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
# lines(ellipse_boundary, col="black", lwd=3)
# legend('topright',paste0(perc_lvl*100,"% Highest Density Ellipse"), lwd = 3,
#        bty = "n")



arrows(median(param_plot[,1]), hdi(log(param_plot),.9)[1,2],
       median(param_plot[,1]),hdi(log(param_plot),.9)[2,2],
       length=0.05, angle=90, code=3)
arrows(hdi(param_plot,.9)[1,1],median(log(param_plot[,2])),
       hdi(param_plot,.9)[2,1],median(log(param_plot[,2])),
       length=0.05, angle=90, code=3)
legend('topright',c('Posterior Medians','90% HPD'),
       pch = c(15,NA),
       lwd = c(NA, 1))

unimode_x = which.max(density(param_plot[,1])$y) %>%
  density(param_plot[,1])$x[.]
unimode_y = which.max(density(param_plot[,2])$y) %>%
  density(param_plot[,2])$x[.]
points(unimode_x,unimode_y)

hdi(param_plot,.6)[,1]
