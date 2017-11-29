##Simulation Study

#Exploring betas
x <- rbeta(1E4,3,7)

source('/Users/Jake/Dropbox/Research/Computer_Emulation/R_code/utilities_multi_hist.R')



for(i in 1:12){
 # x <- rbeta(1E4,i+1,4)
  #hist(x,main = i)
  x <- seq(0.01,0.99,by =0.01)
  plot(x,dbeta(x,6,i+0.5),type = 'l',main = i)
}

#Make the histogram data
#My own function makes it easier to control bins
my_hist <- function(dat,bounds){
  bins <- numeric(length(alpha)-1)
  for(i in 1:length(bins)){
    bins[i] <- length(which(dat>bounds[i] & dat <=bounds[i+1]))
  }
  return(bins)
}

n_q = 12
n_bins <- 10

alpha <- seq(0,1,length = n_bins+1)

Y_beta <- matrix(0,n_q,n_bins)
num_counts = 1E4
beta_a = 6
for(i in 1:n_q){
  x_dat <- rbeta(num_counts,beta_a,i+0.5)
  Y_beta[i,] <- my_hist(x_dat,bounds=alpha)
  plot(Y_beta[i,]/num_counts)
}
(pretend_q_vals <- 1:n_q + 0.5)#basically just don't want beta(x,1)


(colnames(Y_beta) <- alpha[-length(alpha)])

holdout <- 3

Y_new <- Y_beta[holdout,]
Y <- Y_beta[-holdout,]
d <- pretend_q_vals[-holdout]
d_new <-pretend_q_vals[holdout]
(d_scaled <- (d-min(d))/(max(d)-min(d)))
d_new_s <- (d_new-min(d))/(max(d)-min(d))
t_minus = 0
t_star <- 1
(trunc_cols <- which(as.numeric(colnames(Y))<t_star))
Y_trunc <- Y[,trunc_cols]
Y_new_trunc <- Y_new[trunc_cols]

#data and parameters
data <- list()
data$num_counts <- num_counts
data$holdout <- holdout
data$beta_a <- beta_a
data$Y <- Y
data$Y_new <- Y_new
data$Y_trunc <- Y_trunc
data$Y_new_trunc <- Y_new_trunc
data$t_star <- t_star
data$d_scaled <- d_scaled
data$d_new_s <- d_new_s
data$trunc_cols <-trunc_cols


(I <- dim(Y)[1])
(J <- dim(Y_trunc)[2])

jitter = FALSE
if(jitter){
  alpha[-c(1,J+1)] <- alpha[-c(1,J+1)] + rnorm(J-1,0,1E-3)
}
(b <- (alpha - lag(alpha))[-1])

data$I <- I
data$J <- J
data$alpha <- alpha
data$b <- b

ave_data = FALSE #Change this if you change something
if(save_data|jitter){
  save(data, file = "data_list.Rdata")
}

#Parameter things
params <- list()

N <- 10 #One fewer than number of bins: the max number for N
#N is the top of the summation - Nx will be the number of parameters

#N <- 20 #One fewer than number of bins: the max number for N
#More in B matrix than A matrix



# if(N%%2){
#   n = ceiling(N/2)
#   for(j in 2:(J+1)){
#     B[n,j-1] <- t_star*(cos(2*pi*n*alpha[j-1]/t_star) - cos(2*pi*n*alpha[j]/t_star))/(2*pi*n)
#   }
#}

basis_type = 'cos'

#make_coefs is in utilities_multi_hist.R
C <- make_coefs(type= basis_type)
apply(C,2,sum)

#Get the columns that have all zeros
(bad_cols <- which(apply(C,2,function(x){sum(abs(x)<1E-5)})==dim(C)[1]))
if(length(bad_cols)) C <- C[,-bad_cols]

#Nx is the number of parameters
(Nx <- dim(C)[2])
(bin_sizes  <- (alpha-lag(alpha))[-1])

r = 0.95
(sig <- needed_sig(1E-5,t_st = t_star, t_min = t_minus))
(c <- sig_to_c(sig = sig,r=r))

#This won't cause problems for the cos-only vector, because we'll never use Nw
##and bad_cols will only be in 1:N
Nz <- Nw <-  1:N
if(length(bad_cols)){
  bad_z <- bad_cols[which(bad_cols<=N)]
  bad_w <- (bad_cols-N)[which(bad_cols>N)]
  if(length(bad_z)) Nz <- Nz[-bad_z]
  if(length(bad_w)) Nw <- Nw[-bad_w]
}

if(basis_type == 'sin'){
  (r_vec <- c*c(r^Nz,r^Nw))
}else if(basis_type=='cos'){
  (r_vec <- c*r^Nz)
}else stop('You got basis problems, big fella')



params$c <- c
params$Nx <- Nx
params$r <- r
params$Nz <- Nz
params$Nw <- Nw
params$r_vec <- r_vec
params$C <- C
params$basis_type = basis_type


(run_description <- paste0("r_vec is c*r^n, with r = ",r,", c=",round(c,3),
                           ", basis type = ",basis_type, ", Nx = ",Nx))
write(run_description,file="model_description.txt")


#Check the ranks of our coefficient matrices
print(paste("Nx =",Nx))
print(paste("r =",r))
print(paste("c =",round(c,3)))

#Starting place for X
#Use fact that the sum of probabilities add to one
#Could be issue if we trim before observations end
#I.e. if gam != 1

#j_trim <- 1
#C_trim <- C[-j_trim,]
#P_hat <- t(Y_trunc[,-j_trim]/num_counts) #P_hat is (J-1)xI

#Different starting points
#X_mle <- sqrt(0.5)*solve(R)%*%solve(C_trim)%*%(P_hat - replicate(I,b[-j_trim]/t_star))
#X_0 <- replicate(I,rep(0,Nx))



#Priors
lam_a <- rep(1,Nx)
lam_b <- rep(1,Nx)

ell_a <- rep(1,Nx)
ell_b <- rep(1,Nx)


###hier_gp_mh_i() is in utilities_multi_hist.R

print('pre-X_kap')
X_kap_run = matrix(nrow = I, byrow = FALSE,data = c(
  rep(5E-2,I)#1
  ,rep(1E-1,I)#2
  ,rep(1E-1,I) #3
  ,rep(1E-1,I) #4
  ,rep(1E-1,I)#5
  ,rep(1E-1,I)#6
  ,rep(1E-1,I)#7
  ,rep(1E-1,I)#8
  ,rep(1E-1,I)##9
  # ,rep(5E-1,I)#10
  # ,rep(5E-1,I)#11
  # ,rep(5E-1,I)#12
  # ,rep(5E-1,I)#13
  # ,rep(1,I)#14
  # ,rep(1,I)#15
  # ,rep(3,I)#16
  # ,rep(3,I)#17
  # ,rep(5,I)))#18
))

params$X_kap <- X_kap_run

#Save the parameters to use later
save(params,file = "params_list.Rdata")
print('starting run')
start_t <- proc.time()
hope_i <- hier_gp_mh_i(iters = 5E5,verbose = TRUE, burnin_prop = 0.4,X_kap = X_kap_run)
write(paste((proc.time() - start_t)[3]/60,'minutes'),file = 'model_time.txt')

save(hope_i, file = "sampler_vals.Rdata")
