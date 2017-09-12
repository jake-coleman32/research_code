library(MASS)
library(cluster)
library(emdbook)
library(dplyr)
library(reshape2)
library(RobustGaSP)
library(mvtnorm)
library(Matrix)

setwd("/Users/Jake/Dropbox/Research/JETSCAPE/second_project/")
folder = "forSTAT 2/"
save_path <- "/Users/Jake/Dropbox/Research/JETSCAPE/second_project/forSTAT_old_long_run/"
save_pics = FALSE


all_dsets <- c(
  "AuAu200-cen-00-10"
  ,"AuAu200-cen-40-50"
  ,"PbPb2760-cen-00-05"
  ,"PbPb2760-cen-30-40"
  ,"PbPb5020-cen-00-10"
  ,"PbPb5020-cen-30-50"
)

systems_to_calibrate = c("AuAu200"
                        ,  "PbPb2760"
                         , "PbPb5020"
                         )
if(length(systems_to_calibrate)==3){
  hist_main = "All Datasets Simultaneous Calibration"
}else{
  hist_main = paste(paste0(systems_to_calibrate,collapse = " & "),"Calibration")
}


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
for(i in 1:length(dsets_to_use)){
  dset = dsets_to_use[i]
  current_dset = all_dsets[dset]
  (current_experiment = strsplit(current_dset,'-')[[1]][1])
  
  au[[i]] <- read.table(paste0(folder,current_dset,".dat"))
  names(au[[i]]) <- c("pT","RAA_exp","Stat_err","Sys_err",paste0("RAA_",as.character(1:24)))
  
  #Separate output, change from wide to long
  mod_dat <- dplyr::select(au[[i]], -c(RAA_exp,Stat_err,Sys_err)) %>%
    melt(id = "pT") %>%
    arrange(by = pT)
  
  exp_dat <- au[[i]]$RAA_exp
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
Y_final <- sweep(Y,2,apply(Y,2,mean)) %>%
  sweep(2,apply(Y,2,sd),FUN = '/') 
Y_svd <- svd(Y)

#Test for Normality
#Not egregiously off
# for(j in 1:dim(Y_final)[2]){
#   qqnorm(Y_final[,j],main = j)
#   qqline(Y_final[,j])
# }

#None are <0.05 even without accounting for multiple testing
for(j in 1:dim(Y_final)[2]){
  print(shapiro.test(Y_final[,j])$p.value)
}

#How many PCs?
eigs <- Y_svd$d^2
V_q <- cumsum(eigs)/sum(eigs)

if(save_pics) pdf(paste0(save_path,'au_var_explained.pdf'))
plot(V_q[1:6],type = 'o',
     xlab = "Number of PC r",
     main = "Fraction of Variance Explained",
     ylab = expression(F[r]),
     cex.lab = 2,
     cex.axis = 1.3,
     cex.main = 1.4,
     pch = 19)
if(save_pics) dev.off()

q <- 3

V = Y_svd$v[,1:q]
S = diag(Y_svd$d[1:q])


Z <- as.matrix(Y_final)%*%V %>%
  as.data.frame()

(design <- read.table("latin_hc_design.txt",header = TRUE))
scaled_d <- lapply(design,function(x){
  (x-min(x))/(max(x) - min(x))
}) %>% 
  do.call(cbind,.) %>%
  as.data.frame()

final_mod <- lapply(Z,rgasp,design = scaled_d)

#############
##Prediction
############
holdout = 15
train_d = scaled_d[-holdout,]
test_d = scaled_d[holdout,]

train_Y = Y_final[-holdout,]
test_Y = Y[holdout,]

train_Z = as.matrix(train_Y)%*%V
test_Z = as.matrix(test_Y)%*%V

train_mod <- lapply(as.data.frame(train_Z),rgasp,design = train_d)
test_mod <- lapply(train_mod,RobustGaSP::predict,testing_input = test_d)

pred_Z <- lapply(test_mod,function(x)x$mean) %>%
  do.call(cbind,.) 
pred_err_Z <- lapply(test_mod,function(x)x$sd) %>%
  do.call(cbind,.) %>%
  as.numeric()%>%
  diag() 

pred_Y <- pred_Z %*%t(V) %>%
  sweep(2,apply(Y,2,sd),FUN = '*') %>%
  sweep(2,apply(Y,2,mean),FUN = '+')

#This gives negative values
pred_err_Y <- V%*%pred_err_Z %*%t(V) %>%
  sweep(2,apply(Y,2,sd)^2,FUN = '*') %>%
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
     cex.main = 1.4)
arrows(as.numeric(test_Y), as.numeric(pred_Y) - 2*sqrt(pred_err_Y),
       as.numeric(test_Y),as.numeric(pred_Y) + 2*sqrt(pred_err_Y),
       length=0, angle=90, code=3)

abline(a = 0,b = 1)

#Calibration
#Rotate the experimental data with same PCA values
Y_exp_final <- (all_exp_dat -apply(Y,2,mean))/apply(Y,2,sd)

(Z_exp <- Y_exp_final%*%V)

(Z_cov <- t(V)%*%diag(1/apply(Y,2,sd))%*%block_covs %*%diag(1/apply(Y,2,sd)) %*%V)

mh_cal <- function(niters = 1E4,burnin_prop = 0.3,
                   t_kap = c(0.1,0.1)
                   ){
  
  #Do various sampler setups
  time_start <- proc.time()
  
  burnin <- burnin_prop*niters
  
  t_out <- matrix(0,niters,2)
  t_ratio <- ob_ratio <- 0
  
  #draw from priors
  t_cur <- runif(2)

  #Get current values
  (X_cur <- (t(as.matrix(t_cur))))
  cur_mod <- lapply(final_mod, RobustGaSP::predict, testing_input = X_cur)
  
  #(descrep_cov_cur <- cov_mat(pT_scaled,ell_cur,lam_cur))
  
  for(i in 1:(niters + burnin)){
    #Get the time every ten percent of iterations
    if(!i%%(niters*0.1)){
      flush.console()
      cat("\r i = ", i, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
    }
    
    #Propose new t_st
    t_st <- rnorm(2,t_cur,sd = t_kap)
    
    #Check to make sure parameter is within [0,1]
    auto_reject = FALSE
    if(sum(t_st<0) + sum(t_st>1)){
      auto_reject = TRUE
      ob_ratio = ob_ratio + 1
    }
    
    (X_st_sqrt <- t(as.matrix(t_st)))
    X_st <- X_st_sqrt^2
    #Get the predictive mean/sd for t_st
    st_mod <- lapply(final_mod, RobustGaSP::predict, testing_input = X_st)
    
    (z_st_means <- do.call(c,lapply(st_mod,function(x){x$mean})))
    z_st_sd <- do.call(c,lapply(st_mod,function(x){x$sd}))
    
    z_cur_means <- do.call(c,lapply(cur_mod,function(x){x$mean}))
    z_cur_sd <- do.call(c,lapply(cur_mod,function(x){x$sd}))
    
    (zmod_st <- rnorm(length(z_st_means),z_st_means,z_st_sd))
    zmod_cur <- rnorm(length(z_cur_means),z_cur_means,z_cur_sd)
    
    #Draw y^M for current and st predictive distributions
    #(ymod_st <- rnorm(length(st_mod$mean),st_mod$mean,st_mod$sd))
    #(ymod_cur <- rnorm(length(cur_mod$mean),cur_mod$mean,cur_mod$sd))
    
    #Get likelihoods
 
      like_cur <- mvtnorm::dmvnorm(x = Z_exp,
                          mean = zmod_cur,
                          sigma = Z_cov,
                          log = TRUE)
      like_st <- mvtnorm::dmvnorm(x = Z_exp,
                         mean = zmod_st,
                         sigma = Z_cov,
                         log = TRUE)

    #Includes constant prior 
    ratio <- sum(like_st) - sum(like_cur) + sum(log(2*X_st)) - sum(log(2*X_cur))
    
    #Put exponential prior on first input parameter?
    # if(j==1) ratio <- ratio + dexp(tj_st,log = TRUE) - dexp(t_cur[j],log = TRUE)
    
    if(runif(1) < exp(ratio) & !auto_reject){
      cur_mod <- st_mod
      t_cur <- t_st
      t_ratio <- t_ratio + 1
    }
    
    if(i > burnin){
      t_out[i-burnin,] <- t_cur
    }
  }#i loop
  
  print("t_ratio is")
  print(t_ratio/(niters + burnin))
  print("ob ratio is")
  print(ob_ratio/(niters + burnin))
  
  # print(paste("lam ratio is", lam_ratio/(niters + burnin)))
  #print(paste("ell ratio is", ell_ratio/(niters + burnin)))
  
  return(list(params = t_out))
  
}


res <- mh_cal(niters = 1E5,t_kap = c(8E-2,8E-2))
#save(res,file = paste0(save_path,'res_old.Rdata'))

param_plot <- matrix(0,dim(res$params)[1],dim(res$params)[2])
#param_plot <- matrix(0,dim(res$params)[1],dim(res$params)[2])
for(j in 1:2){
  param_plot[,j] <- res$params[,j]^2*(max((design[,j])) - min((design[j]))) + min((design[j]))
}
param_plot_sqrt <- sqrt(param_plot)

plot(param_plot[,1],type = 'l', ylab = expression(Lambda^jet), xlab = 'Iteration',
     main = expression(Lambda^jet~'Trace Plot'))
plot(param_plot[,2],type = 'l', ylab = expression(~alpha[s]^med), xlab = 'Iteration',
     main = expression(~alpha[s]^med~'Trace Plot'))

plot(res$params[,1],type = 'l', ylab = expression(sqrt(Lambda^jet)), xlab = 'Iteration',
     main = expression(sqrt(Lambda^jet)~'Trace Plot'))
plot(res$params[,2],type = 'l', ylab = expression(sqrt(alpha[s]^med)), xlab = 'Iteration',
     main = expression(sqrt(alpha[s]^med)~'Trace Plot'))

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
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

panel.image <- function(x,y,...){
  #f1 <- kde2d(x, y, n = 100)
  #image(f1,add = TRUE)
  points(x,y, col = rgb(1,0,0,.1))
}

if(save_pics){ pdf(paste0(save_path,
                         paste0(systems_to_calibrate,collapse = "_"),
                         ".pdf"))}
pairs(res$params, 
      panel = panel.image,
      diag.panel = panel.hist,
      #pch = 19,
      #cex = .3,
      col = rgb(1,0,0,.1),
      labels = c(expression(sqrt(Lambda^jet)),
                 expression(sqrt(alpha[s]^med))),
      upper.panel = NULL,
      cex.lab = 1.5,
      cex.axis = 1.3,
      cex.main = 1.4,
      mgp = c(2.3,1,0),
      las = 2,
      main = hist_main)
if(save_pics) dev.off()

effectiveSize(param_plot)


f1 <- kde2d(param_plot[,1], (param_plot[,2]), n = 100)
par(mar=c(5.1, 5.1, 4.1, 2.1))
image(f1,
      xlab = expression(Lambda^jet),
      ylab = "",
      main = "Heat Map of Posterior Draws",
      cex.lab = 2,
      cex.axis = 1.3,
      cex.main = 1.4,
      mgp = c(3,1,0)
)
title(ylab = expression(alpha[s]^med), 
      mgp = c(2.1,1,0),
      cex.lab = 2,
      cex.axis = 1.3)

perc_lvl = c(.6,.75,.9)
HPDregionplot(param_plot, prob = perc_lvl,
              col=c("black"), lty = c(1,5,3), add=TRUE)


legend('topright',paste0(perc_lvl*100,"%"),title = "Highest Density Kernel Estimate",
       lty = c(1,5,3),
       bty="n")


#####
#Extra visualization/analysis
#####

plot_lines <- function(dat,pT_col = 'pT',exp_col = 'RAA_exp',
                       err_col_stat = 'exp_err_stat', 
                       err_col_sys = 'exp_err_sys',
                       col_start = 3,
                       alpha_val = 0.01,
                       ...){
  
  plot(dat[,pT_col],dat[,i],
       cex = 0.5,
       col = rgb(1,0,0,alpha_val),type = 'l',
       ylim = c(min(dat[,exp_col] - 2*dat[,err_col_sys],
                    dat[,exp_col] - 2*dat[,err_col_stat]),
                max(dat[,exp_col] + 2*dat[,err_col_sys],
                    dat[,exp_col] + 2*dat[,err_col_stat])),
       xlab = expression(p[T]),
       ylab = expression(R[AA]),
       cex.main = 1.3,
       cex.lab = 1.5,
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
  arrows(dat[,pT_col]-arrow_offset, dat[,exp_col] - 2*dat[,err_col_stat],
         dat[,pT_col]-arrow_offset, dat[,exp_col] + 2*dat[,err_col_stat],
         length=0.05, angle=90, code=3,col = 'blue')
  arrows(dat[,pT_col]+arrow_offset, dat[,exp_col] - 2*dat[,err_col_sys],
         dat[,pT_col]+arrow_offset, dat[,exp_col] + 2*dat[,err_col_sys],
         length=0.05, angle=90, code=3, col = 'black')
  legend('topleft',c('Systematic Errors',
                     'Statistical Errors'),
         col = c('blue','black'), lwd =1)
  
}

##First, go from concatenated to many datasets, with p_T in rows

#1. Cut into 6 original datasets

plot_all_dsets <- function(all_dset_mat,
                           orig_dset_list = au,
                           dset_strings = all_dsets,
                           title_end = "",
                           alpha_val = 0.01,
                           save_pics = FALSE){
  cur_col <- 1
  for(j in 1:length(orig_dset_list)){
    num_pT <- dim(orig_dset_list[[j]])[1]
    cur_dset <- all_dset_mat[,cur_col:(cur_col + num_pT - 1)] %>%
      t() %>%
      cbind('pT' = orig_dset_list[[j]]$pT,
            'RAA_exp' = orig_dset_list[[j]]$RAA_exp,
            'exp_err_stat' = orig_dset_list[[j]]$Stat_err,
            'exp_err_sys' =  orig_dset_list[[j]]$Sys_err,
            .)
    cur_col = cur_col + num_pT
    
    
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
  
  post_pred_mod <-lapply(train_mod, RobustGaSP::predict,
                         testing_input = draws)
  
  pred_means <- lapply(post_pred_mod,function(x)x$mean) %>%
    do.call(cbind,.)
  
  #Rotate
  Y_pred_scaled <- pred_means %*%t(rot_mat)
  
  #Scale
  Y_pred <- sweep(Y_pred_scaled,2,apply(scale_mat,2,sd),FUN = "*") %>%
    sweep(2,apply(scale_mat,2,mean),FUN = '+') 
  
  if(dim(Y_pred)[1]<num_samples){
    stop('Num samples bigger than number of draws')
  }
  
  Y_pred <- Y_pred[sample(1:dim(Y_pred)[1],num_samples),]
  
  plot_all_dsets(Y_pred,title_end = title_end,
                 alpha_val = alpha_val,
                 save_pics = save_pics)
}

plot_draws(res$params,title_end = "Posterior",
           alpha_val = 0.03,
           save_pics = TRUE)
plot_draws(matrix(runif(2E3),ncol = 2),title_end = "Prior",
           alpha_val = 0.03,
           save_pics = TRUE)






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
