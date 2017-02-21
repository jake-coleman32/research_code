setwd("/Users/Jake/Dropbox/Research/Computer_Emulation/meeting_8_3")

load(file = "jonah_stand_param_cal_xs.Rdata")
all_runs <- list("jonah_xs" = param_cal_s)
rm(param_cal_s)

load(file = "j_stand_param_cal_s.Rdata")
all_runs[[length(all_runs) + 1]] = param_cal_s
names(all_runs)[length(all_runs)] = "jonah_s"
rm(param_cal_s)

load(file = "j_stand_xs_alp_2.Rdata")
all_runs[[length(all_runs) + 1]] = param_cal_s
names(all_runs)[length(all_runs)] = "jonah_xs_alp2"
rm(param_cal_s)

setwd("/Users/Jake/Dropbox/Research/Computer_Emulation/meeting_7_27")

load(file = "param_cal_xs.Rdata")
all_runs[[length(all_runs) + 1]] = param_cal_s
names(all_runs)[length(all_runs)] = "regular_xs"
rm(param_cal_s)

load(file = "param_cal_s.Rdata")
all_runs[[length(all_runs) + 1]] = param_cal_s
names(all_runs)[length(all_runs)] = "regular_s"
rm(param_cal_s)

# load(file = "param_cal.Rdata")
# all_runs[[length(all_runs) + 1]] = param_cal
# names(all_runs)[length(all_runs)] = "reg"
# rm(param_cal)


accept_rates <- unlist(lapply(all_runs, function(x){
    return(x$accept_rate)
}))

ob_rates <- unlist(lapply(all_runs, function(x){
  x$ob_rate
}))
sigma_sq_z = c(0.0036,0.2,0.0036,0.0036,0.2)
rate_mat <- cbind(accept_rates,ob_rates,sigma_sq_z)
xtable(rate_mat,digits = 4)

design_mins <- do.call(cbind,lapply(design_m,function(x){min(x)}))
design_maxes <- do.call(cbind,lapply(design_m,function(x){max(x)}))
params_unscaled <- matrix(numeric(prod(dim(param_1000))),ncol = dim(param_1000)[2])
for(j in 1:dim(param_1000)[2]){
  params_unscaled[,j] <- param_1000[,j]*(design_maxes[j] - design_mins[j]) + design_mins[j]
}

param_mats <- lapply(seq_along(1:dim(all_runs[[1]]$param_mat)[2]),function(i){
  temp <- do.call(cbind,lapply(all_runs,function(x){
    x$param_mat[,i]
  }))
  temp_burnin <- temp[2000:dim(temp)[1],]
  indices <- 1:dim(temp_burnin)[1]
  scaled <- temp_burnin[!indices%%500,]
  unscaled <- scaled
  for(j in 1:dim(scaled)[2]){
    unscaled[,j] <- scaled[,j]*(design_maxes[j] - design_mins[j]) + design_mins[j]
  }
  return(scaled)
})

par(mai=c(1.92,1.12,0.82,0.42))

param_names = c('Norm',expression(alpha),expression(tau [0]),expression(eta~"/ s"),
                expression("k" [pi]))
boxplot(param_mats[[1]],main = "Boxplot for Norm",las = 2,
        ylab = "Parameter_value")
boxplot(param_mats[[2]],main = expression("Boxplots for"~ alpha),las = 2,
        ylab = "Parameter_value")
boxplot(param_mats[[3]],main = expression("Boxplots for"~tau [0]),las = 2,
        ylab = "Parameter_value")
boxplot(param_mats[[4]],main = expression("Boxplots for"~eta~"/ s"),las = 2,
        ylab = "Parameter_value")
boxplot(param_mats[[5]],main = expression("Boxplots for"~"k" [pi]),las = 2,
        ylab = "Parameter_value")

e_sample_sizes <- do.call(cbind, lapply(param_mats,function(x){
  effectiveSize(x)
}))
xtable(e_sample_sizes,digits = 0)

par(mai=c(1.02,0.82,0.82,0.42))