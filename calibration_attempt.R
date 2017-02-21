library(coda)
#draw x_st from multivariate normal
#find z_st via GP emulation
#plug into likelihood, do MH (prior is uniform)
#let's just try it with one well-behaved value

draws <- 1000
param_vec <- c(
  'norm' = .2,
  'alp' = .2,
  'tau' = .2,
  'etas' = .5,
  'kpi' = .2
)

#Figure 3, model outputs
#ALICE points
alice <- read.csv("Dropbox/Research/JETSCAPE/mtd-paper/data/exp/alice_means.csv")[,-1]
alice$cent <- as.factor(cent_names_m)
(y_exp <- with(alice,matrix(c(sqrt(mult),v2,v3),nrow = 1)))
colnames(y_exp) <- paste(rep(c("sqrt_mult","v2","v3"),each = 6),rep(as.character(alice$cent),3),sep = "_")

#regular standardization
y_exp_stand = (y_exp - unlist(lapply(Y[-c(1,4),],mean)))/unlist(lapply(Y[-c(1,4),],sd))

#jonah's standardization
y_exp_stand = (y_exp/alice_m)*weights_m - unlist(lapply(Y2,mean))


(z_exp = sqrt(m)*(y_exp_stand)%*%V%*%solve(S))
sig_z <- 0.06
cov_z <- sig_z^2*abs(diag(as.numeric(z_exp)))
#Ok how do we calculate z_exp?
#z_st comes from data divided by alice, multiplied by weights, mean subtracted, PCs

get_z <- function(pred_mod_list,param_vec){
  pred_z_list <- lapply(pred_mod_list,RobustGaSP::predict,#from cleaner_jonah.R
                        testing_input = matrix(param_vec,nrow = 1))
  means = do.call(cbind,lapply(pred_z_list,function(x){x$mean}))
  covs = do.call(cbind,lapply(pred_z_list,function(x){x$sd}))
  out_z <- as.numeric(rmvnorm(1,mean = means,sigma = diag(as.numeric(covs))))
  return(out_z)
}

(old_z <- get_z(run_robust$rob_mod,param_vec))

metro_cal <- function(n_iter,k_alp){
  param_mat <- matrix(numeric(n_iter*length(param_vec)),ncol = length(param_vec))
  param_mat[1,] <- param_vec
  accept <- ob <- rep(0,n_iter)
  for(i in 2:n_iter){
    #alp_star <- rnorm(1,alp[i-1],sqrt(k_alp))
    
    params_star <- rmvnorm(1,param_mat[i-1,],diag(rep(k_alp),length(param_vec)))
    #params_star['alp'] <- alp_star
    
    z_star <- get_z(run_robust$rob_mod,params_star)
    
    (log_star <- dmvnorm(as.numeric(z_star),mean = as.numeric(z_exp),sigma = cov_z,log = TRUE))
    (log_old <- dmvnorm(as.numeric(old_z),mean = as.numeric(z_exp),sigma = cov_z,log = TRUE))
    
    test_log <- log_star - log_old
    
    if(log(runif(1)) < test_log & !sum(params_star<0 | params_star >1)){
      param_mat[i,] <- params_star
      accept[i] <- 1
    }else{
      if(!!sum(params_star<0 | params_star >1)){
        ob[i] = 1
      }
      param_mat[i,] <- param_mat[i-1,]
    }
    old_z = get_z(run_robust$rob_mod,param_mat[i,])
    if(!i%%(n_iter/100)){
      perc = n_iter/100
      print(paste("Simulation is",i/perc,"percent done"))
    }
  }
  print(paste("Out of bounds is",sum(ob)/n_iter))
  print(paste("Acceptance fraction is",sum(accept)/n_iter))
  return(list(param_mat = param_mat,ob_rate = sum(ob)/n_iter, accept_rate = sum(accept)/n_iter))
}

start_t <- proc.time()
param_cal_s <- metro_cal(1000000,0.0005)
proc.time() - start_t

param_burnin <- param_cal_s$param_mat[2000:dim(param_cal_s$param_mat)[1],]
indices <- 1:dim(param_burnin)[1]
param_1000 <- param_burnin[!indices%%1000,]


design_mins <- do.call(cbind,lapply(design_m,function(x){min(x)}))
design_maxes <- do.call(cbind,lapply(design_m,function(x){max(x)}))
params_unscaled <- matrix(numeric(prod(dim(param_1000))),ncol = dim(param_1000)[2])
for(j in 1:dim(param_1000)[2]){
  params_unscaled[,j] <- param_1000[,j]*(design_maxes[j] - design_mins[j]) + design_mins[j]
}

acf(params_unscaled[,1],main = "ACF Norm")
acf(params_unscaled[,2],main = expression("ACF "~ alpha))
acf(params_unscaled[,3],main = expression("ACF "~ tau [0]))
acf(params_unscaled[,4],main = expression("ACF "~ eta~"/ s"))
acf(params_unscaled[,5],main = expression("ACF "~ "k" [pi]))

xtable(t(as.matrix(coda::effectiveSize(param_1000))))

par(mar = c(5.1,4.1,1.1,2.1))

hist(params_unscaled[,1],xlab = "Norm",main = "",breaks = 20,
     xlim = c(design_mins[1],design_maxes[1]))
hist(params_unscaled[,2],xlab = expression(alpha),main = "",breaks = 20,
     xlim = c(design_mins[2],design_maxes[2]))
hist(params_unscaled[,3],xlab = expression(tau [0]),main = "",breaks = 20,
     xlim = c(design_mins[3],design_maxes[3]))
hist(params_unscaled[,4],xlab = expression(eta~"/ s"),main = "",breaks = 20,
     xlim = c(design_mins[4],design_maxes[4]))
hist(params_unscaled[,5],xlab = expression("k" [pi]),main = "",breaks = 20,
     xlim = c(design_mins[5],design_maxes[5]))

plot(params_unscaled[,1],ylab = "Norm",xlab = "Iteration",main = "Trace Plot",type= 'l')
plot(params_unscaled[,2],ylab = expression(alpha),xlab = "Iteration",main = "Trace Plot",type= 'l')
plot(params_unscaled[,3],ylab = expression(tau [0]),xlab = "Iteration",main = "Trace Plot",type= 'l')
plot(params_unscaled[,4],ylab = expression(eta~"/ s"),xlab = "Iteration",main = "Trace Plot",type= 'l')
plot(params_unscaled[,5],ylab = expression("k" [pi]),xlab = "Iteration",main = "Trace Plot",type= 'l')

